from dataclasses import dataclass, field
from random import shuffle
from types import SimpleNamespace
from typing import List, Dict

from pandas import DataFrame, Series
import pandas as pd
from sklearn import clone
from sklearn.model_selection import train_test_split
from tqdm.auto import tqdm

from .data_classes import AttributesStore
from .data_classes import MultiBlockDataSet
from .coefficients import CoefficientsManager
from .coefficients import Coefficients, Contributions
from .multi_block_pipeline import TwoBlockPipeline

from .predictions import Result
from .roc_auc import compare_roc_curves


@dataclass
class CrossValidation:
    """Perform Cross validation, cross-validation on permuted samples,
    and validation of holdout set on the models trained for cross-validation.

    Note: CV does not care about pre-processing, this is the responsibility of the pipeline.
    """

    estimators: List[TwoBlockPipeline] = field(init=False)
    train_dataset: MultiBlockDataSet
    coefficients: AttributesStore[Coefficients] = field(init=False)
    contributions: AttributesStore[Contributions] = field(init=False)

    def __post_init__(self):
        self.estimators = []

    def cross_validate(
        self, pipeline: TwoBlockPipeline, n=100,
        use_seed: bool = True, stratify: bool = True,
        only_coefficients: bool = False, stripped: bool = False,
        permute: bool = False, progress_bar: bool = True, verbose: bool = True,
        coefficients: Dict[str, str] = {'x': 'coef_'}, early_normalization=False
    ) -> Result:

        cv_results = Result(
            train_data=self.train_dataset,
            test_data=self.train_dataset,
            predicted_probabilities=[],
            binary_true_responses=[]
        )
        patients_classes = self.train_dataset.binary_response
        train_data = {
            block: getattr(self.train_dataset, block)
            for block in pipeline.block_pipelines.keys()
        }

        some_block = next(iter(train_data.values()))
        assert all(
            all(train_block.index == some_block.index)
            for train_block in train_data.values()
        )

        coefficients_manager = CoefficientsManager(coefficients, abundance_matrices=train_data)

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
                    train_binary, test_binary,
                    train_patients, test_patients
                ) = train_test_split(
                    patients_classes,
                    self.train_dataset.observations,
                    stratify=patients_classes if stratify else None,
                    random_state=seed if use_seed else None
                )
                seed += 1

                if len(set(train_binary)) == 2 and len(set(test_binary)) == 2:
                    # we have so few patients that we need to check if split
                    # has at least one in each class...
                    ok = True
                elif verbose:
                    print('Skipping a split with too few response samples')

            train: Dict[str, DataFrame] = {}
            test: Dict[str, DataFrame] = {}
            for block, data in train_data.items():
                train[block] = train_data[block].loc[train_patients]
                test[block] = train_data[block].loc[test_patients]

            split_pipeline: TwoBlockPipeline = clone(pipeline)

            if early_normalization:
                split_pipeline.transformed_blocks = {
                    name: block.loc[train_patients]
                    for name, block in pipeline.transformed_blocks.items()
                }
                for name, p in pipeline.block_pipelines.items():
                    setattr(split_pipeline, name, p)
                split_pipeline.combine.fit(*split_pipeline.transformed_blocks.values())
            else:
                split_pipeline.fit(*train.values())

            self.estimators.append(split_pipeline)

            if verbose > 1:
                print('Applying dynamic filters to train get filtered columns')

            coefficients_manager.add(split_matrices=split_pipeline.transformed_blocks, pipeline=split_pipeline)

            if only_coefficients:
                continue

            cv_dataset = MultiBlockDataSet(
                [*test.values()],
                case_class=self.train_dataset.case_class,
                response=self.train_dataset.response
            )

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

        for estimator in tqdm(self.estimators):
            prediction = estimator.predict(estimator, dataset)
            result.predicted_probabilities.append(prediction)
            result.binary_true_responses.append(dataset.binary_response)

        return result


@dataclass
class CrossValidationResult:
    # CrossValidation stores all estimators, coefficients and contributions obtained across cross-validation splits
    cross_validation: CrossValidation

    # Coefficients obtained by the model trained on full training set
    coefficients: AttributesStore[Coefficients]
    contributions: AttributesStore[Contributions]

    # results for multiple models trained on subsets of training set and tested on subsets of the same training set
    cross_validation_results: Result

    # results for multiple models trained on subsets of training set and tested on the holdout set
    sub_sampling_test_results: Result

    # result for the model trained on the ALL training set and tested on the holdout set
    test_result: Result

    # NOTE: the pipeline is re-trained on the full training set;
    # as such they should not be used to obtain predictions of the training set for the
    # performance assessment, as this is better served by the cross validated results.
    pipeline: TwoBlockPipeline

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


def repeated_cross_validation(
    pipeline: TwoBlockPipeline, train_data: List[DataFrame],
    response, case_class='Tuberculosis',
    n=1000, use_seed=False, stratify=True,
    permute=False, progress_bar=True,
    only_coefficients=False, stripped=False,
    early_normalization=False,
    verbose=False,
    test_data=None,
    coefficients={'x': 'coef_'}
) -> CrossValidationResult:
    assert case_class in set(response)

    train_dataset = MultiBlockDataSet(data=train_data, case_class=case_class, response=response)
    test_dataset = MultiBlockDataSet(data=test_data if test_data is not None else [], case_class=case_class, response=response)

    if early_normalization:
        pipeline.fit_transform_blocks(train_dataset.x, train_dataset.y)

    cross_validation = CrossValidation(train_dataset=train_dataset)

    cv_results = cross_validation.cross_validate(
        pipeline,
        n=n, use_seed=use_seed, stratify=stratify, only_coefficients=only_coefficients,
        stripped=stripped, permute=permute, progress_bar=progress_bar, verbose=verbose,
        coefficients=coefficients, early_normalization=early_normalization
    )

    sub_sampling_test_result = cross_validation.validate(test_dataset) if test_data is not None else None

    if verbose:
        print('Re-fitting on the entire dataset...')

    # TODO: make it work with 1D
    if early_normalization:
        pipeline.combine.fit(*pipeline.transformed_blocks.values())
    else:
        pipeline.fit(train_dataset.x, train_dataset.y)

    test_result = Result.from_test_set(
        pipeline=pipeline,
        test_set=test_dataset,
        train_set=train_dataset
    ) if test_data is not None else None

    coefficients_manager = CoefficientsManager(
        # TODO: the abundance matrices should be taken straight from the dataset...
        #  the internal data storage could be refactored to Dict[str, DataFrame]
        coefficients, abundance_matrices={'x': train_dataset.x, 'y': train_dataset.y}
    )
    coefficients_manager.add(split_matrices=pipeline.transformed_blocks, pipeline=pipeline)
    coefficients_manager.concatenate()

    return CrossValidationResult(
        # TODO: move those into results or move CV result into CV?
        # coefficients and contributions
        coefficients=coefficients_manager.to_store(subclass=Coefficients, skip_compute=stripped),
        contributions=coefficients_manager.to_store(subclass=Contributions, skip_compute=stripped),

        # pipelines
        pipeline=pipeline,
        cross_validation=cross_validation,

        # results
        cross_validation_results=cv_results,
        test_result=test_result,
        sub_sampling_test_results=sub_sampling_test_result,

        # parameters
        early_normalization=early_normalization
    )


def null_distributions_over_cv(pipeline, omic, response, block, permutations=100, **kwargs):
    """Null distributions for an alternative hypothesis that the:
    a) mean coefficient value over cross validations is not random
    b) mean contribution value over cross validations is not random
    
    Arguments:
        pipeline: the extended pipeline
        omic: passed to cross_validated_lasso
        response: passed to cross_validated_lasso
        block: x or y
        permutations: number of permutations
        **kwargs: passed to cross_validated_lasso
    """

    groups = {'test', 'cross_validation'}
    kinds = {'coefficients', 'contributions'}

    nulls = {
        group: {
            kind: []
            for kind in kinds
        }
        for group in groups
    }

    for i in tqdm(range(permutations)):
        result = repeated_cross_validation(
            pipeline, omic, response=response,
            permute=True, progress_bar=False,
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

        results = repeated_cross_validation(
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
