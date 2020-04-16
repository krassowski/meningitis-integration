import random
from typing import TypeVar, Generic, Union, Dict, Callable

from pandas import DataFrame, Series, Index
from pydantic import BaseConfig

from pydantic.dataclasses import dataclass
from sklearn.pipeline import Pipeline


class ValidationConfig(BaseConfig):
    arbitrary_types_allowed = True


dataclass = dataclass(config=ValidationConfig)

BlockId = str
Block = Union[DataFrame, Series]
Blocks = Dict[BlockId, DataFrame]


@dataclass
class BlocksWrapper:
    """A thin wrapper around Dict[str, Union[DataFrame, Series]]
    facilitating easy manipulation of blocks.
    Use it for block manipulation if needed and then just pass .blocks to the pipeline.
    """
    blocks: Dict[BlockId, Block]

    @property
    def aligned_index(self) -> Index:
        self.verify_index_integrity(require_ordered=True)
        blocks_iterator = iter(self.blocks.values())
        first_block = next(blocks_iterator)
        return first_block.index

    def consensus_index(self, operation: str):
        assert operation in {'intersection', 'union', 'difference'}
        # TODO use reduce instead for better readability?
        blocks_iterator = iter(self.blocks.values())
        first_block = next(blocks_iterator)
        index = first_block.index
        for block in blocks_iterator:
            index = getattr(index, operation)(block.index)
        return list(index)

    def take_by_consensus_index(self, operation: str) -> 'BlocksWrapper':
        index = self.consensus_index(operation)
        return BlocksWrapper({
            block_id: block.loc[index]
            for block_id, block in self.blocks.items()
        })

    def take_union(self) -> 'BlocksWrapper':
        return self.take_by_consensus_index('union')

    def take_intersection(self) -> 'BlocksWrapper':
        return self.take_by_consensus_index('intersection')

    def align_indices(self) -> 'BlocksWrapper':
        """Assuming that all blocks have identical indices (up to a permutation),

        align them so that the indices have exact same ordering."""
        self.verify_index_integrity(require_ordered=False)
        blocks_iterator = iter(self.blocks.values())
        first_block = next(blocks_iterator)
        index = first_block.index
        new = BlocksWrapper({
            block_id: block.loc[index]
            for block_id, block in self.blocks.items()
        })
        new.verify_index_integrity(require_ordered=True)
        return new

    def verify_index_integrity(self, require_ordered: bool):
        blocks_iterator = iter(self.blocks.values())
        first_block = next(blocks_iterator)
        for block in blocks_iterator:
            if require_ordered:
                ok = (block.index == first_block.index).all()
            else:
                ok = set(block.index) == set(first_block.index)
            if not ok:
                print(block.index, first_block.index)
                raise ValueError()
        return True

    def without(self, block_id: BlockId) -> 'BlocksWrapper':
        new_blocks = self.blocks.copy()
        del new_blocks[block_id]

        return BlocksWrapper(new_blocks)

    def nullify_exclusively(
        self, to_nullify: BlockId,
        index='intersection'
    ) -> 'BlocksWrapper':

        nullified = DataFrame(
            index=self.without(to_nullify).consensus_index(index),
            columns=self.blocks[to_nullify].columns
        ).fillna(0)

        new_blocks = self.blocks.copy()
        new_blocks[to_nullify] = nullified

        return BlocksWrapper(new_blocks)

    def transform_single_block(self, block_id: BlockId, pipeline: Pipeline, mark=True) -> 'BlocksWrapper':

        new_blocks = self.blocks.copy()
        new_blocks[block_id] = pipeline.fit_transform(new_blocks[block_id])

        if mark:
            new_blocks[block_id].is_transformed = True

        return BlocksWrapper(new_blocks)

    def add_supervision_block(self, conditions_vector: Series, name: str) -> 'BlocksWrapper':
        self.verify_index_integrity(require_ordered=True)
        first_block = next(iter(self.blocks.values()))
        index = first_block.index

        conditions_list = list(conditions_vector)
        conditions = Series(conditions_list, index=conditions_list)
        conditions_block = conditions.loc[index]

        return BlocksWrapper({
            name: conditions_block,
            **self.blocks
        })

    def to_dataset(self, **kwargs) -> 'MultiBlockDataSet':
        return MultiBlockDataSet(self.blocks, **kwargs)


RandomizeCallback = Callable[[Series], Series]
randomization_methods = {
    # random.sample samples without replacement
    'permute': lambda x: random.sample(x, len(x)),
    # random.sample samples with replacement
    'bootstrap': lambda x: random.choices(x, k=len(x))
}


@dataclass
class MultiBlockDataSet:
    data: Dict[str, Block]
    case_class: str

    # this may contains responses from patients outside of this dataset,
    # therefore is private and should not be accessed directly
    response: Series

    @property
    def binary_response(self) -> Series:
        return Series(self.binary_response_for(observations=self.observations))

    def binary_response_for(self, observations):
        classes = observations.map(self.response)
        assert len(set(classes)) == 2
        return classes == self.case_class

    def randomized(self, method: Union[str, RandomizeCallback]) -> 'MultiBlockDataSet':
        """Randomize the data using provided method ('permute, 'bootstrap', or custom callback)."""

        if isinstance(method, str):
            method = randomization_methods[method]

        # when randomizing, we have to make sure we do not introduce labels from outside of the dataset
        # (having an outcome for those does not mean that we have the corresponding data)
        response_subset = list(set(self.observations) & set(self.response.index))

        return MultiBlockDataSet(
            data=self.data,
            case_class=self.case_class,
            response=Series(data=method(list(self.response.loc[response_subset])), index=response_subset)
        )

    @property
    def class_imbalance(self) -> float:
        return self.binary_response.mean()

    def __getattr__(self, item):
        try:
            return self.data[item]
        except KeyError:
            raise AttributeError

    def apply(self, func):
        for block_id, block in self.data.items():
            self.data[block_id] = func(block)
        return self

    @property
    def observations(self) -> Series:
        """For example, patients in DA problems"""
        first_block = next(iter(self.data.values()))
        BlocksWrapper(blocks=self.data).verify_index_integrity(require_ordered=True)
        return first_block.index


ModelAttribute = TypeVar('ModelAttribute')


@dataclass
class AttributesStore(Generic[ModelAttribute]):
    """Stores model attributes such as coefficients per-block of the model.

    Probably a better design would build upon the blocks...
    """
    blocks: Dict[BlockId, ModelAttribute]

    def __getattr__(self, item):
        return self.blocks[item]

    @classmethod
    def from_dicts(cls, coefficients_values, abundance_matrices, subclass, skip_compute=False):
        return cls(blocks={
            matrix: subclass(values, abundance=abundance_matrices[matrix], skip_compute=skip_compute)
            for matrix, values in coefficients_values.items()
        })
