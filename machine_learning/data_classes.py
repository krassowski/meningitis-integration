from dataclasses import dataclass, field
from functools import partial
from types import SimpleNamespace, FunctionType
from typing import List

import pandas as pd
from pandas import DataFrame, Series

from .coefficients import Coefficients, Contributions
from .flexible_pipeline import FlexiblePipeline
from .roc_auc import roc_auc_plot_data, compute_cv_statistics


@dataclass
class DataSet:
    data: DataFrame
    case_class: str

    # this may contains responses from patients outside of this dataset,
    # therefore is private and should not be accessed directly
    _response: Series

    @property
    def binary_response(self) -> Series:
        return self.binary_response_for(observations=self.observations)

    @property
    def observations(self) -> Series:
        """For example, patients in DA problems"""
        return self.data.index

    def binary_response_for(self, observations):
        classes = observations.map(self._response)
        assert len(set(classes)) == 2
        return classes == self.case_class

    @property
    def class_imbalance(self) -> float:
        return self.binary_response.mean()


@dataclass
class Result:
    """Can hold a single or multiple results of predictions"""
    # predictions based on test set
    predicted_probabilities: List[Series]
    binary_true_responses: List[Series]

    # the model can be trained on one dataset and tested on another;
    # however, the train and test set may be the same (for sake of
    # assessment of performance on the trained data)
    test_data: DataSet
    train_data: DataSet  # do we really need that here?

    metrics: SimpleNamespace = field(init=False)

    def __post_init__(self):
        from sklearn import metrics
        self.metrics = SimpleNamespace(**{
            name: partial(self.compute_metric, metric=func)
            for name in dir(metrics)
            for func in [getattr(metrics, name)]
            if isinstance(func, FunctionType)
        })

    def compute_metric(self, metric, round=False):
        """Compute an sklearn metric.

        Args:
            round: may be used as a cheap proxy of classifier.predict()
                   for functions which require class labels
        """
        return Series(
            metric(true, (predicted.round() if round else predicted))
            for true, predicted in zip(self.binary_true_responses, self.predicted_probabilities)
        )

    @classmethod
    def from_validation_set(cls, pipeline, validation_set: DataSet, train_set=None):
        return cls(
            predicted_probabilities=[Series(
                pipeline.predict_proba(validation_set.data)[:, 1],
                index=validation_set.observations
            )],
            binary_true_responses=[
                validation_set.binary_response
            ],
            test_data=validation_set,
            train_data=train_set
        )

    def consensus_result(self):
        # WARNING: this - when used for performance assessment
        # gives too optimistic a result (as it averages the
        # probabilities from all repeats). It is very easy to
        # get an AUC of 1 here - but it is kind of incorrect.
        mean_probability = (
            pd.concat(self.predicted_probabilities)
            .rename_axis('observation')
            .groupby('observation')
            .mean()
        )

        # reorder to match the test data:
        # mean_probability = mean_probability.loc[self.test_data.observations]
        # NOPE - this could cause an issue when only a subset of full test data
        # was evaluated (producing nans)

        if set(mean_probability.index) != set(self.test_data.observations):
            missing = set(self.test_data.observations) - set(mean_probability.index)
            print(
                f'Warning: not all of the test observations were predicted'
                f' prior to consensus ranking: {len(missing)} missing out of'
                f' {len(mean_probability)}'
            )

        binary_response = self.test_data.binary_response_for(
            mean_probability.index
        )

        return Result(
            predicted_probabilities=[mean_probability],
            binary_true_responses=[binary_response],
            test_data=self.test_data,
            train_data=None
        )

    @property
    def roc_auc(self):
        return roc_auc_plot_data(
            self.predicted_probabilities,
            self.binary_true_responses
        )

    @property
    def cv_auc(self):
        return compute_cv_statistics(
            self.predicted_probabilities,
            self.binary_true_responses
        )[0]


@dataclass
class CrossValidationResult:
    # Coefficients obtained across cross-validations
    cross_validation_coefficients: Coefficients
    cross_validation_contributions: Contributions

    # Coefficients obtained by the model trained on full training set
    coefficients: Coefficients
    contributions: Contributions

    # results for multiple models trained on subsets of training set and tested on subsets of the same training set
    cross_validation_results: Result

    # results for multiple models trained on subsets of training set and tested on the holdout set
    sub_sampling_validation_results: Result

    # results for the model trained on the ALL training set and tested on the holdout set
    validation_result: Result

    # NOTE: the pipelines are re-trained on the full training set;
    # as such they should not be used to obtain predictions of the training set for the
    # performance assessment, as this is better served by the cross validated results.

    # pipeline without the pre-filtering steps
    stripped_pipeline: FlexiblePipeline
    # full pipeline
    full_pipeline: FlexiblePipeline

    @property
    def validation_dataset(self) -> DataSet:
        return self.validation_result.test_data

    @property
    def training_dataset(self) -> DataSet:
        return self.cross_validation_results.train_data

    @property
    def case_class(self) -> Series:
        # this should go over the Results and after integrity check
        # return the case_class extracted from those
        raise NotImplementedError

    @property
    def _response(self) -> Series:
        raise NotImplementedError
