from abc import ABC, abstractmethod
from dataclasses import field
from functools import partial
from types import SimpleNamespace, FunctionType
from typing import List, Callable, Union
from warnings import warn

import pandas as pd
from pandas import Series

from .data_classes import MultiBlockDataSet, Block
from .multi_block_pipeline import MultiBlockPipeline
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
    def from_test_set(cls, pipeline: MultiBlockPipeline, test_set: MultiBlockDataSet, train_set=None):
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

        # NOTE: no reordering here
        if set(mean_probability.index) != set(self.test_data.observations):
            missing = set(self.test_data.observations) - set(mean_probability.index)
            warn(
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


PredictionFunction = Callable[[MultiBlockPipeline, Block], Series]


class PipelinePredictor(ABC):

    @abstractmethod
    def fit(self, pipeline: MultiBlockPipeline):
        """Fit to given pipeline."""

    @abstractmethod
    def predict(self, dataset: MultiBlockDataSet):
        """Predict (probabilities or classes) for given dataset."""

    def fit_predict(self, pipeline: MultiBlockPipeline, dataset: MultiBlockDataSet):
        """Fit to the pipeline and predict for the dataset."""
        self.fit(pipeline)
        return self.predict(dataset)

    def __call__(self, pipeline: MultiBlockPipeline, dataset: MultiBlockDataSet):
        return self.fit_predict(pipeline, dataset)


class RobustTwoBlockPredictor(PipelinePredictor):
    """A utility wrapper around prediction function for use in real-world pipelines,

    which use filtering and may receive data with nans, which needs to be accounted for.
    """
    fitted_block: Block = None
    pipeline: MultiBlockPipeline = None

    def __init__(self, predict: Union[PredictionFunction, str], verbose=False):
        """
        Args:
            predict: a function callback returning predictions given a MultiBlockPipeline and a dataset,
                or name of a method of the pipeline, which is a prediction function (e.g. 'predict_proba' for sklearn)
            verbose: whether diagnostic messages should be displayed
        """
        self._predict = predict
        self.verbose = verbose

    def remove_nans(self, data: Block, response: Block):
        if data.isnull().any().any():
            if self.verbose:
                print('Dropping nans: ', data.isnull().any(axis=1).sum())
            not_dropped = data.index[~data.isnull().any(axis=1)]
            data = data.dropna(axis=0)
            response = response.loc[not_dropped]
        return data, response

    def select_variables_by_fitted_filter(self, data: Block):
        if len(self.fitted_block.columns) != len(data.columns) or (self.fitted_block.columns != data.columns).all():
            if self.verbose:
                print('Taking only the variables as selected in the previously fitted model')
            difference = self.fitted_block.columns.difference(data.columns)
            assert not difference.any()  # this can happen if outliers are excluded from one block only
            data = data.loc[:, self.fitted_block.columns].fillna(0)
        return data

    def fit(self, pipeline: MultiBlockPipeline):
        full_blocks = pipeline.combine.transform(pipeline.transformed_blocks)
        full_blocks = pipeline.adapter.transform(full_blocks)
        self.fitted_block = full_blocks['X']
        self.pipeline = pipeline

    def predict(self, dataset: MultiBlockDataSet):
        blocks = self.pipeline.partial_transform(blocks=dataset.data)
        if len(blocks) == 2:
            data = blocks['X']
            response = blocks['y']

            data, response = self.remove_nans(data, response)
            data = self.select_variables_by_fitted_filter(data)

            if callable(self._predict):
                p = self._predict(self.pipeline, data)
            else:
                p = self.pipeline.call(self._predict, data)
                p = p[:, 1]
            return Series(p, index=response)
        else:
            raise ValueError(
                'predict_proba can only handle two blocks after combination'
                '(with one being the outcome);'
                f' {len(blocks)} provided'
            )


def predict_proba(pipeline: MultiBlockPipeline, dataset: MultiBlockDataSet):
    """This is intended to work with standard sklearn API
    (usually for a single data block + single response block).
    """
    return RobustTwoBlockPredictor(predict='predict_proba').fit_predict(pipeline, dataset)
