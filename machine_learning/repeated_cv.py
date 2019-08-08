from random import shuffle

from pandas import DataFrame, Series
import pandas as pd
from sklearn.model_selection import train_test_split
from tqdm.auto import tqdm

from .coefficients import Coefficients, Contributions
from .flexible_pipeline import FlexiblePipeline

from .data_classes import DataSet, Result, CrossValidationResult
from .preprocessing import (
    StaticFilter, DynamicFilter, RSideNormalizer
)
from .roc_auc import compare_roc_curves


def subset_and_scale(omic: DataFrame, pipeline: FlexiblePipeline, patients_filter, early_normalization):
    omic_filtered = omic.iloc[omic.index.str.contains(patients_filter)]

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


def repeated_cross_validation(
    pipeline: FlexiblePipeline, omic, patients_filter,
    response, data_comes_from_normalizer=True,
    case_class='Tuberculosis',
    n=1000, use_seed=False, stratify=False,
    permute=False, progress_bar=True, _subseted_and_scaled=False,
    only_coefficients=False, stripped=False,
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
    verbose=False,
    validation_data=None
) -> CrossValidationResult:
    cv_coeffs = []
    full_pipeline = pipeline

    if not _subseted_and_scaled:
        omic_filtered, pipeline = subset_and_scale(
            omic, pipeline, patients_filter, early_normalization
        )
    else:
        omic_filtered = omic
        # usually subsetting and scaling is performed here, but for the performance reasons
        # it is done beforehand for permutation test and for it only.
        assert permute

    train_data = DataSet(data=omic_filtered, case_class=case_class, _response=response)
    validation_data = DataSet(data=validation_data, case_class=case_class, _response=response)

    cv_results = Result(
        train_data=train_data,
        test_data=train_data,
        predicted_probabilities=[],
        binary_true_responses=[]
    )

    sub_sampling_validation_result = Result(
        train_data=train_data,
        test_data=validation_data,
        predicted_probabilities=[],
        binary_true_responses=[]
    )

    patients_classes = train_data.binary_response

    # shuffle split
    seed = 0

    if verbose:
        print('Fitting cross-validation models...')

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
            train_X = pipeline.apply_steps(DynamicFilter, train_X)

        cv_coeffs.append(
            pipeline.get_coefficients(train_X.columns)
        )

        if not only_coefficients:
            proba = pipeline.predict_proba(test_X)

            cv_results.predicted_probabilities.append(
                Series(proba[:, 1], index=test_patients)
            )
            cv_results.binary_true_responses.append(test_y)

    if verbose:
        print('Re-fit on the entire dataset...')

    pipeline.fit(omic_filtered, patients_classes)

    validation_result = Result.from_validation_set(
        pipeline=pipeline,
        validation_set=validation_data,
        train_set=train_data
    ) if validation_data.data else None

    cv_coeffs = pd.concat(cv_coeffs, axis=1)

    if data_comes_from_normalizer:
        omic_filtered = pipeline.apply_steps(DynamicFilter, omic_filtered)

    coeffs = DataFrame(pipeline.get_coefficients(index=omic_filtered.columns))

    return CrossValidationResult(
        # coefficients and contributions
        coefficients=Coefficients(coeffs, abundance=omic_filtered, skip_compute=stripped),
        contributions=Contributions(coeffs, abundance=omic_filtered, skip_compute=stripped),
        cross_validation_coefficients=Coefficients(cv_coeffs, abundance=omic_filtered, skip_compute=stripped),
        cross_validation_contributions=Contributions(cv_coeffs, abundance=omic_filtered, skip_compute=stripped),

        # pipelines
        stripped_pipeline=pipeline, full_pipeline=full_pipeline,

        # results
        cross_validation_results=cv_results,
        validation_result=validation_result,
        sub_sampling_validation_results=sub_sampling_validation_result
    )


def null_distributions_over_cv(
    pipeline, omic, patients_filter, response, early_normalization,
    permutations=100,
    **kwargs
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

    kinds = {
        'coefficients', 'contributions',
        'cross_validation_coefficients',
        'cross_validation_contributions',
    }

    nulls = {kind: [] for kind in kinds}

    omic_filtered, pipeline = subset_and_scale(
        omic, pipeline, patients_filter, early_normalization
    )

    for i in tqdm(range(permutations)):
        result = repeated_cross_validation(
            pipeline, omic_filtered, patients_filter,
            data_comes_from_normalizer=True, response=response,
            permute=True, progress_bar=False, _subseted_and_scaled=True,
            only_coefficients=True, stripped=True,
            **kwargs
        )

        for kind in kinds:
            nulls[kind].append(getattr(result, kind))

    return {
        f'mind_{kind}': pd.concat(values, axis=1)
        for kind, values in nulls.items()
    }


def performance_comparison(
    pipeline, omic, datasets, patients_filter,
    datasets_subset=None,
    **kwargs
):
    # TODO: this is defunct
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

            a_response, a_prediction = r.cv_response_and_mean_prediction()
            b_response, b_prediction = r2.response_and_mean_prediction()

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
