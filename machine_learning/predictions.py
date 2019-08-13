from dataclasses import field
from functools import partial
from types import SimpleNamespace, FunctionType
from typing import List

import pandas as pd
from pandas import Series

from .data_classes import MultiBlockDataSet, dataclass
from .multi_block_pipeline import TwoBlockPipeline
from .roc_auc import roc_auc_plot_data, compute_cv_statistics

# TODO: something is amiss with pydantic here?
from dataclasses import dataclass


@dataclass
class Result:
    """Can hold a single or multiple results of predictions"""
    # predictions based on test set
    predicted_probabilities: List[Series]
    binary_true_responses: List[Series]

    # the model can be trained on one dataset and tested on another;
    # however, the train and test set may be the same (for sake of
    # assessment of performance on the trained data)
    test_data: MultiBlockDataSet
    train_data: MultiBlockDataSet  # do we really need that here?

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
    def from_test_set(cls, pipeline: TwoBlockPipeline, test_set: MultiBlockDataSet, train_set=None):
        return cls(
            predicted_probabilities=[pipeline.predict(pipeline, test_set)],
            binary_true_responses=[test_set.binary_response],
            test_data=test_set,
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
