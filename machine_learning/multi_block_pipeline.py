from abc import abstractmethod, ABC
from dataclasses import field

from pydantic import BaseConfig
from pydantic.dataclasses import dataclass
from typing import Callable, Union, Dict, Type

from pandas import Series, DataFrame
from sklearn.base import BaseEstimator
from sklearn.pipeline import Pipeline

from .data_classes import MultiBlockDataSet


class ValidationConfig(BaseConfig):
    arbitrary_types_allowed = True


def find_attribute(estimator, attribute='coef_', to_traverse=('best_estimator_', '_final_estimator'), self=True):
    """It is assumed that estimator has no more than one of the 'to_traverse' attributes"""
    if self and hasattr(estimator, attribute):
        return getattr(estimator, attribute)
    for field in to_traverse:
        if hasattr(estimator, field):
            next_estimator = getattr(estimator, field)
            return find_attribute(next_estimator, attribute)
    return None


BlockId = Type[str]


@dataclass(config=ValidationConfig)
class GeneralizedBlockPipeline(ABC, BaseEstimator):

    # preprocess: Pipeline
    combine: Union[Pipeline, BaseEstimator]
    predict: Callable[['GeneralizedBlockPipeline', MultiBlockDataSet], Series]
    transformed_blocks: Dict[str, DataFrame] = field(init=False)

    def __post_init__(self):
        self.transformed_blocks: Dict[str, DataFrame] = {}

    @property
    @abstractmethod
    def block_pipelines(self) -> Dict[BlockId, Pipeline]:
        pass

    def fit_transform_blocks(self, *args):
        for (block_id, block_pipeline), block_data in zip(self.block_pipelines.items(), args):
            self.transformed_blocks[block_id] = block_pipeline.fit_transform(block_data)
        return self

    def fit(self, *args):
        self.fit_transform_blocks(*args)
        self.combine.fit(*self.transformed_blocks.values())
        return self

    def attribute(self, attribute):
        """Locate an attribute/property in the combined pipeline"""
        return find_attribute(self.combine, attribute)

    def call(self, attribute, *args, **kwargs):
        """Locate a function"""
        method = self.attribute(attribute)
        return method(*args, **kwargs)

    def get_coefficients(self, index: Series, coefficient='coef_'):
        coef = find_attribute(self.combine, attribute=coefficient)
        # TODO: this is a hack, revise API instead!
        if coefficient.endswith('_'):
            coef = coef[0]
        else:
            coef = coef[:, 0]
        return Series(coef, index=index)


@dataclass(config=ValidationConfig)
class MultiBlockPipeline(GeneralizedBlockPipeline):
    """Intended for N-block integration"""
    block_pipelines: Dict[BlockId, Pipeline]


@dataclass(config=ValidationConfig)
class TwoBlockPipeline(GeneralizedBlockPipeline):
    """Intended for regression, classification, discriminant analysis etc."""
    x: Pipeline
    y: Pipeline

    # workaround for https://github.com/samuelcolvin/pydantic/issues/739
    def __post_init__(self):
        # workaround for a re-loader issue
        super(self.__class__, self).__post_init__()

    @property
    def block_pipelines(self):
        return {'x': self.x, 'y': self.y}


@dataclass(config=ValidationConfig)
class OneBlockPipeline(GeneralizedBlockPipeline):
    """Intended for PCA, unsupervised clustering etc."""
    x: Pipeline

    # workaround for https://github.com/samuelcolvin/pydantic/issues/739
    def __post_init__(self):
        # workaround for a re-loader issue
        super(self.__class__, self).__post_init__()

    @property
    def block_pipelines(self):
        return {'x': self.x}

    @property
    def transformed(self):
        return self.transformed_blocks['x']


def predict_proba(pipeline: GeneralizedBlockPipeline, dataset: MultiBlockDataSet):
    """This is intended to work with standard sklearn API (usually for a single block)"""
    p = pipeline.call('predict_proba', pipeline.x.transform(dataset.x))
    return Series(p[:, 1], index=pipeline.y.transform(dataset.y))
