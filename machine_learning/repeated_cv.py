from dataclasses import dataclass
from math import sqrt
from random import shuffle

import numpy as np
from pandas import DataFrame, Series
import pandas as pd
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.model_selection import train_test_split
from tqdm.auto import tqdm

from helpers.p_values import gpd_p_value
from helpers.r import p_adjust
from .preprocessing import (
    StaticFilter, DynamicFilter, RSideNormalizer
)
from .roc_auc import roc_auc_plot_data, compare_roc_curves


def is_instance_by_name(child, parent):
    """Like isinstance, but works with autoreloader.
    
    Unfotunately, is sensitive to location (will yield false positives).
    """
    return parent.__name__ in {x.__name__ for x in child.__class__.mro()}


class ExtendedPipeline(Pipeline):
    
    def apply_steps(self, kind, data):
        for step, transformer in self.steps:
            if is_instance_by_name(transformer, kind):
                data = transformer.fit_transform(data)
        return data

    def exclude(self, to_exclude):
        return make_extended_pipeline(
            *[
                transformer
                for step, transformer in self.steps
                # again, a workaround for autoreloader
                # if type(transformer) not in to_exclude
                if transformer.__class__.__name__ not in {
                    t.__name__
                    for t in to_exclude
                }
            ]
        )


def make_extended_pipeline(*args, **kwargs):
    pipeline = make_pipeline(*args, **kwargs)
    pipeline.__class__ = ExtendedPipeline
    return pipeline


@dataclass
class RepeatedCVResult:
    coefficients: DataFrame
    test_y: Series
    cv_accuracy: Series
    predicted_probabilities: Series
    class_imbalance: float
    omic: DataFrame
    patient_to_class: Series
    case_class: Series
    #pipeline: ExtendedPipeline

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def roc_auc_data(self):
        return roc_auc_plot_data(self.predicted_probabilities, self.test_y)


def subset_and_scale(omic, pipeline, patients_filter, early_normalization, scale='warn'):
    omic_filtered = omic.iloc[omic.index.str.contains(patients_filter)]

    if scale not in ['warn', None]:
        omic_scaled = scale.fit_transform(omic_filtered)
        omic_filtered = DataFrame(
            omic_scaled,
            index=omic_filtered.index,
            columns=omic_filtered.columns
        )
    elif scale ==' warn':
        print('Warning: not scaling!')

    if omic_filtered.isnull().any().any():
        print(
            'Rejecting n =',
            len(omic_filtered.isnull().any(axis=1)),
            'genes with nulls'
        )
        omic_filtered = omic_filtered.loc[:, ~omic_filtered.isnull().any(axis=0)]

    omic_filtered = pipeline.apply_steps(StaticFilter, omic_filtered)
    pipeline = pipeline.exclude({StaticFilter})
    
    if early_normalization:
        omic_filtered = pipeline.apply_steps(RSideNormalizer, omic_filtered)
        pipeline = pipeline.exclude({RSideNormalizer})
        
    return omic_filtered, pipeline


def summarize_coefficients(
    coeffs_across_cv: RepeatedCVResult, as_contributions: bool,
    null_distributions=None, p_value_kwargs={},
    abundance=None
):
    
    coeffs = coeffs_across_cv.mean(axis=1).to_frame(name='mean')
    coeffs['selected_in'] = (coeffs_across_cv != 0).mean(axis=1)
    coeffs['positive_in'] = (coeffs_across_cv > 0).mean(axis=1)
    coeffs['negative_in'] = (coeffs_across_cv < 0).mean(axis=1)

    coeffs['volatile'] = (
        2 * coeffs[['positive_in', 'negative_in']].min(axis=1)
        / coeffs['selected_in']
    )

    if as_contributions:
        coeffs_across_cv = coeffs_across_cv.abs() / coeffs_across_cv.abs().sum()

    coeffs['mean'] = coeffs_across_cv.mean(axis=1)
    coeffs['stdev'] = coeffs_across_cv.std(axis=1)
    coeffs['ci'] = 1.96 * coeffs_across_cv.std(axis=1) / sqrt(len(coeffs_across_cv))

    coeffs['quantile_0.25'] = coeffs_across_cv.quantile(q=0.25, axis=1)
    coeffs['quantile_0.50'] = coeffs_across_cv.quantile(q=0.50, axis=1)
    coeffs['quantile_0.75'] = coeffs_across_cv.quantile(q=0.75, axis=1)

    if null_distributions is not None:

        if as_contributions:
            null_distribution = null_distributions['mean_contributions']
        else:
            null_distribution = null_distributions['mean_coefficients']

        # reorder
        null_distribution = null_distribution.loc[coeffs.index]

        positive = coeffs[coeffs['mean'] >= 0]
        negative = coeffs[coeffs['mean'] < 0]
        
        p_values = pd.concat([
            gpd_p_value(
                coeffs['mean'][positive.index],
                null_distribution.loc[positive.index],
                **p_value_kwargs
            ),
            gpd_p_value(
                -coeffs['mean'][negative.index],
                -null_distribution.loc[negative.index],
                **p_value_kwargs
            )
        ]).loc[coeffs.index]

        assert (p_values.index == coeffs.index).all()

        coeffs[p_values.columns] = p_values
        coeffs['fdr'] = p_adjust(coeffs['p_value'], 'BH')

        # originally the idea was to look at the differences,
        # but it does not seem worth it
        # abs_p_values = gpd_p_value(
        #     abs(coeffs['mean']),
        #     abs(null_distribution),
        #     **p_value_kwargs
        # )
        # coeffs['absolute_p_value'] = abs_p_values['p_value']
        # coeffs['absolute_fdr'] = p_adjust(coeffs['absolute_p_value'], 'BH')

    if abundance is not None:
        coeffs['mean_abundance'] = abundance.sum().loc[coeffs.index].values

    coeffs['gene'] = pd.Categorical(
        coeffs.index,
        categories=coeffs['mean'].sort_values().index,
        ordered=True
    )
        
    return coeffs


def select_coefficients(
    coeffs, kind, non_zero_ratio_required=0.05,
):    
    selected_coeffs = coeffs[coeffs.selected_in >= non_zero_ratio_required]
    
    print('Selected', len(selected_coeffs), 'coefficients diferent from zero in at least', non_zero_ratio_required, 'repeats')

    return selected_coeffs


def repeated_cross_validation(
    pipeline: ExtendedPipeline, omic, patients_filter,
    response, data_comes_from_normalizer=True,
    case_class='Tuberculosis',
    n=1000, use_seed=False, stratify=False,
    permute=False, progress_bar=True, _subseted_and_scaled=False,
    only_coefficients=False,
    # this will produce over-optimistic prediction assessment
    # when compared to real world performance, due to moderate
    # data leakage in normalization and filtering steps. However,
    # the early normalization improves significantly on the 
    # speed of execution and generally makes the permutation test
    # feasible. Moreover, it closely corresponds to differential
    # expression analysis pipeline where normalization and filtering
    # are often not only using the entire dataset, but are sometimes
    # also not blind to group assignment.
    # This does NOT moves the scaling step before the cross-validation!
    # Scaling may have larger data leakage (depending on the scaler)
    # and thus may completely trash the performance assessment!
    early_normalization=False,
    verbose=True
) -> RepeatedCVResult:
    cv_scores = []
    cv_coefs = []
    cv_proba = []
    cv_test_y = []

    if not _subseted_and_scaled:
        omic_filtered, pipeline = subset_and_scale(
            omic, pipeline, patients_filter, early_normalization,
            scale=None
        )
    else:
        omic_filtered = omic
        # usually subsetting and scaling is perfromed here, but for the performance reasons
        # is is done beforehand for permutation test and for it only.
        assert permute
        
    classes = omic_filtered.index.map(response)
    assert len(set(classes)) == 2
    # to binary vector
    patients_classes = classes == case_class

    # shuffle split
    seed = 0
    
    iterable = tqdm(range(n), total=n) if progress_bar else range(n)
    
    for i in iterable:
        
        if permute:
            shuffle(patients_classes)

        ok = False
        
        while not ok:
            (
                train_X, test_X,
                train_y, test_y,
                train_patients, test_patients
            ) = train_test_split(
                omic_filtered,
                patients_classes,
                omic_filtered.index,
                stratify=patients_classes if stratify else None,
                random_state=seed if use_seed else None
            )
            seed += 1
            
            if len(set(train_y)) == 2 and len(set(test_y)) == 2:
                # we have so few patients that we need to check if split
                # has at least one in each class...
                ok = True
            elif verbose:
                print('Skipping a split with too few response samples')

        pipeline.fit(train_X, train_y)

        if data_comes_from_normalizer:
            train_X = pipeline.apply_steps(
                DynamicFilter,
                omic_filtered[train_X.columns]
            )

        estimator = pipeline._final_estimator
        cv_coefs.append(Series(estimator.coef_[0], index=train_X.columns))
        
        if not only_coefficients:
            cv_scores.append(pipeline.score(test_X, test_y))
            proba = pipeline.predict_proba(test_X)
            cv_proba.append(Series(proba[:, 1], test_patients))
            cv_test_y.append(test_y)

    return RepeatedCVResult(**{
        'coefficients': pd.concat(cv_coefs, axis=1),
        'test_y': cv_test_y,
        'cv_accuracy': cv_scores,
        'predicted_probabilities': cv_proba,
        'class_imbalance': np.mean(patients_classes),
        'omic': omic,
        'patient_to_class': Series(classes, index=omic_filtered.index),
        'case_class': case_class
    })


def null_distributions_over_cv(
    pipeline, omic, patients_filter, response, early_normalization,
    permutations=100, **kwargs
):
    """Null distributions for an alternative hypothesis that the:
    a) mean coefficient value over cross validations is not random
    b) mean contribution value over cross validations is not random
    
    Arguments:
        pipeline: the extended pipeline
        *args: passed to cross_validated_lasso
        permutations: number of permutations
        **kwargs: passed to cross_validated_lasso
    """
    omic_filtered, pipeline = subset_and_scale(
        omic, pipeline, patients_filter, early_normalization,
        scale=None
    )
    
    mean_coefficients = []
    mean_contributions = []
    
    for i in tqdm(range(permutations)):
        coefficients = repeated_cross_validation(
            pipeline, omic_filtered, patients_filter,
            data_comes_from_normalizer=True, response=response,
            permute=True, progress_bar=False, _subseted_and_scaled=True,
            only_coefficients=True,
            **kwargs
        )['coefficients']
        
        mean_coefficients.append(coefficients.mean(axis=1))
        mean_contributions.append((coefficients / coefficients.abs().sum()).mean(axis=1))
        
    return {
        'mean_coefficients': pd.concat(mean_coefficients, axis=1),
        'mean_contributions': pd.concat(mean_contributions, axis=1)
    }


def response_and_mean_prediction(result: RepeatedCVResult):
    mean_prediction = (
        pd.concat(result['predicted_probabilities'])
        .rename_axis('patient')
        .groupby('patient')
        .mean()
    )
    response = (
        result['patient_to_class'].loc[mean_prediction.index]
        ==
        result['case_class']
    )
    return response, mean_prediction


def performance_comparison(
    pipeline, omic, datasets, patients_filter,
    datasets_subset=None,
    **kwargs
):
    assert omic in {'rna', 'protein'}
    comparison_auc = []
    comparison_coef = []
    models = {}
    p_values = {}
    powers = []
    
    for name, (norm_rna, norm_protein) in tqdm(datasets.items()):
        if datasets_subset and name not in datasets_subset:
            continue

        r = repeated_cross_validation(
            pipeline,
            norm_rna if omic == 'RNA' else norm_protein,
            patients_filter, use_seed=True, **kwargs
        )

        auc = r.roc_auc_data.assign(group=name)
        
        for other, r2 in models.items():

            a_response, a_prediction = response_and_mean_prediction(r)
            b_response, b_prediction = response_and_mean_prediction(r2)

            paired = (
                len(r['patient_to_class']) == len(r2['patient_to_class'])
                and all(r['patient_to_class'] == r2['patient_to_class'])
            )
            
            p, power_report = compare_roc_curves(
                a_response, a_prediction,
                b_response, b_prediction,
                paired=paired,
                compute_power=paired
            )

            if paired:
                powers.append({'roc1': name, 'roc2': other, **power_report})
            else:
                # generally, this should not happen for most cases
                print('Skipped power computation for non-paired ROC curves')

            p_values[(name, other)] = p
        
        comparison_auc.append(auc)
        comparison_coef.append(r['coefficients'])
        models[name] = r

    return {
        'auc': pd.concat(comparison_auc, axis=0),
        'coefficients': pd.concat(comparison_coef),
        'power': DataFrame(powers),
        'auc_differ_p_value': DataFrame({
            'p_value': p_values,
        }).T
    }