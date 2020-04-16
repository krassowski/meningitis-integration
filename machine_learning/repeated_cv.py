from collections import defaultdict
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import List, Dict, Union, Tuple
from warnings import warn

from pandas import DataFrame, Series
import pandas as pd
from sklearn import clone
from sklearn.model_selection import BaseCrossValidator, StratifiedShuffleSplit
from sklearn.model_selection._search import BaseSearchCV
from sklearn.model_selection._split import BaseShuffleSplit, ShuffleSplit
from tqdm.auto import tqdm

from machine_learning.model_selection.split import (
    IndependentMultiSizeShuffleSplits, ValidatingSplitter,
    classification_split_validator,
)
from .data_classes import AttributesStore, BlocksWrapper, RandomizeCallback
from .data_classes import MultiBlockDataSet
from .coefficients import CoefficientsManager
from .coefficients import Coefficients, Contributions
from .multi_block_pipeline import MultiBlockPipeline, LeakyPipeline

from .predictions import Result
from .roc_auc import compare_roc_curves

SklearnSplitter = Union[BaseCrossValidator, BaseShuffleSplit]


@dataclass
class CrossValidation:
    """Perform Cross validation, cross-validation on permuted samples,
    and validation of holdout set on the models trained for cross-validation.

    Note: CV does not care about pre-processing, this is the responsibility of the pipeline.
    """

    splits: List[Dict[str, Union[MultiBlockPipeline, float]]] = field(init=False)
    train_dataset: MultiBlockDataSet
    coefficients: AttributesStore[Coefficients] = field(init=False)
    contributions: AttributesStore[Contributions] = field(init=False)

    def __post_init__(self):
        self.splits = []

    @staticmethod
    def are_blocks_indexes_aligned(train_data):
        return BlocksWrapper(blocks=train_data).verify_index_integrity(require_ordered=True)

    def cross_validate(
        self, pipeline: MultiBlockPipeline,
        cv: Union[ValidatingSplitter, SklearnSplitter],
        only_coefficients: bool = False, stripped: bool = False,
        verbose: bool = True, coefficients: Dict[str, str] = {'x': 'coef_'},
    ) -> Result:

        if not isinstance(cv, ValidatingSplitter):
            cv = ValidatingSplitter(cv)

        train_dataset = self.train_dataset

        cv_results = Result(
            train_data=train_dataset,
            test_data=train_dataset,
            predicted_probabilities=[],
            binary_true_responses=[]
        )

        if set(train_dataset.data.keys()) - set(pipeline.block_pipelines.keys()) and verbose:
            print('Note: the pipeline does not utilize all of the data blocks')

        train_data = {
            block: train_dataset.data[block]
            for block in pipeline.block_pipelines.keys()
        }

        assert self.are_blocks_indexes_aligned(train_data) is True

        # note: one should compare the abundances from early and late normalization;
        # if they differ too much, this might explain low performance (if there is such)
        # > abundance_from_early_normalization = pipeline.transformed_blocks
        # later normalization:
        abundance_transformed = clone(pipeline).fit_transform_blocks(train_data).transformed_blocks

        coefficients_manager = CoefficientsManager(coefficients, abundance_matrices=abundance_transformed)

        patients_classes = train_dataset.binary_response

        if verbose:
            print('Fitting cross-validation models...')

        for train_indices, test_indices in cv.split(x=train_dataset.observations, y=patients_classes):
            train: Dict[str, DataFrame] = {}
            test: Dict[str, DataFrame] = {}
            for block, data in train_data.items():
                # TODO: use .values to speed up?
                train[block] = train_data[block].iloc[train_indices]
                test[block] = train_data[block].iloc[test_indices]

            split_pipeline: MultiBlockPipeline = clone(pipeline)
            split_pipeline.fit(train)

            test_ratio = len(test_indices) / (len(train_indices) + len(test_indices))

            self.splits.append({
                'test_ratio': test_ratio,
                'train_ratio': 1 - test_ratio,
                'pipeline': split_pipeline
            })

            coefficients_manager.add(split_matrices=split_pipeline.transformed_blocks, pipeline=split_pipeline)

            if only_coefficients:
                continue

            cv_dataset = MultiBlockDataSet(test, case_class=train_dataset.case_class, response=train_dataset.response)

            prediction = split_pipeline.predict(split_pipeline, cv_dataset)
            cv_results.predicted_probabilities.append(prediction)
            assert len(prediction) == len(cv_dataset.binary_response)

            cv_results.binary_true_responses.append(cv_dataset.binary_response)

        coefficients_manager.concatenate()

        self.coefficients = coefficients_manager.to_store(Coefficients, skip_compute=stripped)
        self.contributions = coefficients_manager.to_store(Contributions, skip_compute=stripped)

        return cv_results

    def validate(self, dataset: MultiBlockDataSet) -> Result:

        result = Result(
            train_data=self.train_dataset,
            test_data=dataset,
            predicted_probabilities=[],
            binary_true_responses=[]
        )

        for split in tqdm(self.splits):
            pipeline = split['pipeline']
            prediction = pipeline.predict(pipeline, dataset)
            result.predicted_probabilities.append(prediction)
            result.binary_true_responses.append(dataset.binary_response)

        return result


@dataclass
class CrossValidationResult:
    # CrossValidation stores all estimators, coefficients and contributions obtained across cross-validation splits
    cross_validation: CrossValidation

    # Coefficients obtained by the model trained on a full training set
    coefficients: AttributesStore[Coefficients]
    contributions: AttributesStore[Contributions]

    # results for multiple models trained on subsets of training set and tested on subsets of the same training set
    cross_validation_results: Result

    # results for multiple models trained on subsets of training set and tested on the holdout set
    sub_sampling_test_results: Result

    # result for the model trained on the ALL training set and tested on the holdout set
    test_result: Result

    # NOTE: the pipeline is re-trained on the full training set using chosen_params (if any);
    # as such they should not be used to obtain predictions of the training set for the
    # performance assessment, as this is better served by the cross validated results.
    pipeline: MultiBlockPipeline

    # original pipeline, not fitted
    original_pipeline: MultiBlockPipeline

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
    early_normalization: bool

    chosen_params: dict
    good_params: DataFrame

    @property
    def validation_dataset(self) -> MultiBlockDataSet:
        return self.test_result.test_data

    @property
    def training_dataset(self) -> MultiBlockDataSet:
        return self.cross_validation_results.train_data

    @property
    def case_class(self) -> str:
        # TODO: ideally, this should go over the results to verify integrity
        return self.cross_validation_results.train_data.case_class

    @property
    def _response(self) -> Series:
        # TODO: ideally, this should go over the results to verify integrity
        return self.cross_validation_results.train_data.response

    def validate(self, dataset: MultiBlockDataSet) -> Result:

        result = Result(
            train_data=self.training_dataset,
            test_data=dataset,
            predicted_probabilities=[],
            binary_true_responses=[]
        )

        prediction = self.pipeline.predict(self.pipeline, dataset)
        result.predicted_probabilities.append(prediction)
        result.binary_true_responses.append(dataset.binary_response)

        return result


def choose_params(cross_validation, pipeline, verbose) -> Tuple[dict, DataFrame]:

    best_params = None
    good_params = None

    # a little bit of duck-typing to check if any parameters were estimated:
    splits = list(cross_validation.splits)
    if isinstance(pipeline.model, BaseSearchCV) or (splits and hasattr(splits[0]['pipeline'].model, 'best_params_')):
        if verbose:
            print('Choosing best parameters out of the ', type(pipeline.model))

        # the best params set is chosen in voting:
        good_params_list = []
        good_scores = defaultdict(list)

        for split in cross_validation.splits:
            internal_cv: BaseSearchCV = split['pipeline'].model
            params = internal_cv.best_params_
            good_params_list.append(params)
            good_scores[frozenset(params.items())].append(internal_cv.best_score_)

        good_params = DataFrame(good_params_list)
        vote_result = (
            good_params
            .groupby(list(good_params.columns))
            .size()
            .sort_values(ascending=False)
            .rename('vote_')
        )

        # TODO: sort by vote and then by score
        good_params_sorted_by_vote = vote_result.reset_index().drop('vote_', axis=1)

        best_params = good_params_sorted_by_vote.iloc[0].to_dict()
        second_best_params = good_params_sorted_by_vote.iloc[1].to_dict()

        scores_of_winner = Series(good_scores[frozenset(best_params.items())])
        scores_of_second = Series(good_scores[frozenset(second_best_params.items())])
        print(
            f'Params chosen with CV vote: {best_params},'
            f' based on voting: {vote_result.head().to_dict()} '
            f' with average {pipeline.model.refit}={scores_of_winner.mean()}±{scores_of_winner.std()},'
            f' compared to {scores_of_second.mean()}±{scores_of_second.std()} of the next best choice'
            f' ({second_best_params}).'
        )

    return best_params, good_params


def cross_validate_and_test(
    pipeline: MultiBlockPipeline, train_data: Dict[str, DataFrame],
    response, case_class, n=1000,
    use_seed=False, stratify=True,
    randomize: Union[str, RandomizeCallback] = False, progress_bar=True,
    only_coefficients=False, stripped=False,
    early_normalization=False,
    verbose=False,
    test_data: Dict[str, DataFrame] = None,
    coefficients={'x': 'coef_'},
    multi_scale=True, test_size_min=0.3, test_size_max=0.3,
    min_class_members=2, min_classes_number=2, max_class_number=2
) -> CrossValidationResult:
    assert case_class in set(response)

    train_dataset = MultiBlockDataSet(data=train_data, case_class=case_class, response=response)
    test_dataset = MultiBlockDataSet(data=test_data or [], case_class=case_class, response=response)

    if randomize:
        if verbose:
            print('Shuffling the response vector of the train dataset...')
        train_dataset = train_dataset.randomized(method=randomize)

    if early_normalization:
        if not isinstance(pipeline, LeakyPipeline):
            warn('Using early normalization without leaky pipeline - is this intended?')
        pipeline.fit_transform_blocks(train_dataset.data)

    cross_validation = CrossValidation(train_dataset=train_dataset)

    if multi_scale:
        assert test_size_max != test_size_min
    else:
        assert test_size_max == test_size_min

    validate_split = classification_split_validator(
        min_class_members_n=min_class_members,
        min_classes_n=min_classes_number,
        max_classes_n=max_class_number
    )

    cv = IndependentMultiSizeShuffleSplits(
        shuffle_cv=StratifiedShuffleSplit if stratify else ShuffleSplit,
        n_sizes=n, test_sizes=(test_size_min, test_size_max),
        validator=validate_split, progress_bar=progress_bar,
        # ShuffleSplit arguments
        random_state=0 if use_seed else None, n_splits=1  # splits for each split size
    )

    cv_results = cross_validation.cross_validate(
        pipeline, cv=cv, only_coefficients=only_coefficients,
        stripped=stripped, verbose=verbose,
        coefficients=coefficients
    )

    sub_sampling_test_result = cross_validation.validate(test_dataset) if test_data else None

    best_params, good_params = choose_params(cross_validation, pipeline, verbose)

    if verbose:
        if best_params:
            print(f'Re-fitting on the entire dataset using parameters from CV: {best_params}')
        else:
            print('Re-fitting on the entire dataset (no parameters estimated)...')

    fitted_pipeline = clone(pipeline)
    if best_params:
        fitted_pipeline.model = fitted_pipeline.model.estimator.set_params(**best_params)

    # NOTE: this pipeline should not be used blindly; if any parameter estimation was performed,
    # the best estimator should be chosen using the cross_validation
    fitted_pipeline.fit(train_dataset.data)

    test_result = Result.from_test_set(
        pipeline=fitted_pipeline,
        test_set=test_dataset,
        train_set=train_dataset
    ) if test_data else None

    coefficients_manager = CoefficientsManager(coefficients, abundance_matrices=fitted_pipeline.transformed_blocks)
    coefficients_manager.add(split_matrices=fitted_pipeline.transformed_blocks, pipeline=fitted_pipeline)
    coefficients_manager.concatenate()

    return CrossValidationResult(
        # TODO: move those into results or move CV result into CV?
        # coefficients and contributions
        coefficients=coefficients_manager.to_store(subclass=Coefficients, skip_compute=stripped),
        contributions=coefficients_manager.to_store(subclass=Contributions, skip_compute=stripped),

        # pipelines
        original_pipeline=pipeline,
        pipeline=fitted_pipeline,
        cross_validation=cross_validation,

        # results
        cross_validation_results=cv_results,
        test_result=test_result,
        sub_sampling_test_results=sub_sampling_test_result,

        # parameters
        early_normalization=early_normalization,
        chosen_params=best_params,
        good_params=good_params
    )


def null_distributions_over_cv(
    pipeline: MultiBlockPipeline, train_data, response, block, method='permute', permutations=100,
    groups=('test', 'cross_validation'),
    kinds=('coefficients', 'contributions'),
    **kwargs
):
    """Null distributions for an alternative hypothesis that the:
    a) mean coefficient value over cross validations is not random
    b) mean contribution value over cross validations is not random

    Arguments:
        pipeline: a multi-block pipeline
        method: randomization method (permute, bootstrap or custom callback)
        block: x or y
        permutations: number of permutations
        **kwargs: passed to cross_validate_and_test method
    """

    nulls = {
        group: {
            kind: []
            for kind in kinds
        }
        for group in groups
    }

    for i in tqdm(range(permutations)):
        result = cross_validate_and_test(
            pipeline, train_data, response=response,
            randomize=method, progress_bar=False,
            only_coefficients=True, stripped=True,
            **kwargs
        )

        for group in groups:
            parent = result if group == 'test' else result.cross_validation
            for kind in kinds:
                store: AttributesStore = getattr(parent, kind)
                nulls[group][kind].append(
                    getattr(store, block).data
                )

    return SimpleNamespace(**{
        group: SimpleNamespace(**{
            f'{kind}': pd.concat(values, axis=1)
            for kind, values in kinds.items()
        })
        for group, kinds in nulls.items()
    })


# TODO move out of here
def performance_comparison(
    pipeline, omic, datasets, patients_filter, to_compare,
    datasets_subset=None,
    **kwargs
):
    assert omic in {'rna', 'protein'}
    assert to_compare in {'cross_validation_results', 'test_result', 'sub_sampling_test_results'}

    comparison_auc = []
    # comparison_coef = []
    models = {}
    p_values = {}
    powers = []

    for name, (norm_rna, norm_protein) in tqdm(datasets.items()):
        if datasets_subset and name not in datasets_subset:
            continue

        results = cross_validate_and_test(
            pipeline,
            norm_rna if omic == 'RNA' else norm_protein,
            patients_filter, use_seed=True, **kwargs
        )
        a = getattr(results, to_compare)

        auc = a.roc_auc.assign(group=name)

        for other_name, other in models.items():

            b = getattr(results, to_compare)

            a_response, a_prediction = a.consensus_result()
            b_response, b_prediction = b.consensus_result()

            paired = (
                len(a.test_data.observations) == len(b.test_data.observations)
                and all(a.test_data.observations == b.test_data.observations)
                and all(a.test_data.binary_response == b.test_data.binary_response)
            )

            p, power_report = compare_roc_curves(
                a_response, a_prediction,
                b_response, b_prediction,
                paired=paired,
                compute_power=paired
            )

            if paired:
                powers.append({'roc1': name, 'roc2': other_name, **power_report})
            else:
                # generally, this should not happen for most cases
                print('Skipped power computation for non-paired ROC curves')

            p_values[(name, other_name)] = p

        comparison_auc.append(auc)
        models[name] = results

    return {
        'auc': pd.concat(comparison_auc, axis=0),
        'power': DataFrame(powers),
        'results': models,
        'auc_differ_p_value': DataFrame({
            'p_value': p_values,
        }).T
    }
