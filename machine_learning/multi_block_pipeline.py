from warnings import warn

from pydantic import BaseModel
from typing import Callable, Union, Dict

from pandas import Series
from sklearn.base import BaseEstimator
from sklearn.pipeline import Pipeline

from machine_learning.adapter import BlocksAdapter, SklearnAdapter
from machine_learning.combine import BlocksCombiner
from .data_classes import MultiBlockDataSet, BlockId, Blocks, Block


class ValidatedModel(BaseModel):
    class Config:
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


class BlocksRelay(BlocksCombiner):

    def fit(self, *args):
        return self

    def transform(self, blocks: Blocks) -> Blocks:
        return blocks


class MultiBlockPipeline(BaseEstimator, ValidatedModel):
    """Intended for N-block integration"""

    block_pipelines: Dict[BlockId, Pipeline]
    model: Union[Pipeline, BaseEstimator]
    predict: Callable[['MultiBlockPipeline', MultiBlockDataSet], Series]
    combine: BlocksCombiner = BlocksRelay()
    adapter: BlocksAdapter = BlocksRelay()
    verbose: bool = False
    transformed_blocks: Dict[BlockId, Block] = None

    def __init__(self, **kwargs):
        ValidatedModel.__init__(self, **kwargs)
        # workaround for pydantic rewriting the dict and for sklearn checking for identity
        # to verify if still needed: assert kwargs['block_pipelines'] is not self.block_pipelines
        assert kwargs['block_pipelines'] == self.block_pipelines
        self.block_pipelines = kwargs['block_pipelines']

    def __post_init__(self):
        self.transformed_blocks = {}

    __custom_params__ = set()
    __ignored_params__ = {'transformed_blocks'}

    @property
    def _params(self):
        params = set()
        for cls in self.__class__.mro():
            if hasattr(cls, '__annotations__'):
                params.update(cls.__annotations__.keys())
        return (params - self.__ignored_params__) | self.__custom_params__

    def get_params(self, deep=True):
        return {
            key: getattr(self, key)
            for key in self._params
        }

    def set_params(self, **params):
        for parameter, value in params.items():
            setattr(self, parameter, value)
        return self

    def fit_transform_blocks(self, blocks: Blocks):
        self.transformed_blocks = self.transform_blocks(blocks, fit=True)
        return self

    def partial_transform(self, blocks: Blocks, fit=False):
        """Transform given blocks preparing them for the use in the final model."""
        blocks = self.transform_blocks(blocks, fit=fit)
        blocks = self.combine.transform(blocks)
        blocks = self.adapter.transform(blocks)
        return blocks

    def fit_from_transformed(self, transformed_blocks: Blocks):
        """Fit model using pre-fitted blocks.
        Please be very-careful when using this, as it might lead to an unintended leakage.
        """
        blocks = self.combine.transform(transformed_blocks)
        blocks = self.adapter.transform(blocks)
        return self.model.fit(**blocks)

    def fit(self, blocks: Blocks):
        self.fit_transform_blocks(blocks)
        combined_blocks = self.combine.fit_transform(self.transformed_blocks)
        adapted_blocks = self.adapter.fit_transform(combined_blocks)
        self.model.fit(**adapted_blocks)
        return self

    def transform_blocks(self, blocks: Blocks, fit=False):
        transformed_blocks = {}
        if set(blocks.keys()) != set(self.block_pipelines):
            warn(f'Blocks ({blocks.keys()}) differ from the set of available block pipelines ({self.block_pipelines})')
        for block_id, block_data in blocks.items():
            block_pipeline = self.block_pipelines[block_id]
            if hasattr(block_data, 'is_transformed') and block_data.is_transformed:
                if self.verbose:
                    print('skipping', block_id, 'marked as already transformed')
                transformed = block_data
            else:
                if fit:
                    transformed = block_pipeline.fit_transform(block_data)
                else:
                    transformed = block_pipeline.transform(block_data)
            transformed_blocks[block_id] = transformed
        return transformed_blocks

    def attribute(self, attribute):
        """Locate an attribute/property in the combined pipeline"""
        return find_attribute(self.model, attribute)

    def call(self, attribute, *args, **kwargs):
        """Locate a function"""
        method = self.attribute(attribute)
        return method(*args, **kwargs)

    def get_coefficients(self, index: Series, coefficient='coef_'):
        coef = find_attribute(self.model, attribute=coefficient)
        # TODO: this is a hack, revise API instead!
        if coefficient.endswith('_'):
            coef = coef[0]
        else:
            coef = coef[:, 0]
        return Series(coef, index=index)


class TwoBlockPipeline(MultiBlockPipeline):
    """Intended for regression, classification, discriminant analysis etc."""
    __custom_params__ = {'x', 'y'}
    __ignored_params__ = {'transformed_blocks', 'block_pipelines'}

    def __init__(self, x: Pipeline, y: Pipeline, adapter=SklearnAdapter(x='X', y='y'), **kwargs):
        kwargs['block_pipelines'] = {'x': x, 'y': y}
        kwargs['adapter'] = adapter
        super().__init__(**kwargs)

    @property
    def x(self) -> Pipeline:
        return self.block_pipelines['x']

    @property
    def y(self) -> Pipeline:
        return self.block_pipelines['y']


class OneBlockPipeline(MultiBlockPipeline):
    """Intended for PCA, unsupervised clustering etc."""
    __custom_params__ = {'x'}
    __ignored_params__ = {'transformed_blocks', 'block_pipelines'}

    def __init__(self, x: Pipeline, adapter=SklearnAdapter(x='X'), **kwargs):
        kwargs['block_pipelines'] = {'x': x}
        kwargs['adapter'] = adapter
        super().__init__(**kwargs)

    @property
    def x(self) -> Pipeline:
        return self.block_pipelines['x']


def predict_proba(pipeline: MultiBlockPipeline, dataset: MultiBlockDataSet):
    """This is intended to work with standard sklearn API (usually for a single block)"""
    blocks = pipeline.partial_transform(blocks=dataset.data)
    if len(blocks) == 2:
        data = blocks['X']
        response = blocks['y']
        if data.isnull().any().any():
            print('Dropping nans: ', data.isnull().any(axis=1).sum())
            not_dropped = data.index[~data.isnull().any(axis=1)]
            data = data.dropna(axis=0)
            response = response.loc[not_dropped]

        full_blocks = pipeline.combine.transform(pipeline.transformed_blocks)
        full_blocks = pipeline.adapter.transform(full_blocks)
        fitted_to_block = full_blocks['X']
        if len(fitted_to_block.columns) != len(data.columns) or (fitted_to_block.columns != data.columns).all():
            # print('Taking only the variables as selected in the previously fitted model')
            difference = fitted_to_block.columns.difference(data.columns)
            assert not difference.any()  # this can happen if outliers are excluded from one block only
            data = data.loc[:, fitted_to_block.columns].fillna(0)

        p = pipeline.call('predict_proba', data)
        return Series(p[:, 1], index=response)
    else:
        raise ValueError(
            'predict_proba can only handle two blocks after combination'
            '(with one being the outcome);'
            f' {len(blocks)} provided'
        )
