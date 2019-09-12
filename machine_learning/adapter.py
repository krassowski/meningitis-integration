from abc import ABC, abstractmethod

from sklearn.base import TransformerMixin

from machine_learning.data_classes import Blocks

ArgumentName = str


class BlocksAdapter(ABC, TransformerMixin):

    @abstractmethod
    def transform(self, blocks: Blocks) -> Blocks:
        pass


class SklearnAdapter(BlocksAdapter):
    """Use this to pass the combined block to sklearn.

    Sklearn can usually handle up to two blocks; for more blocks,
    consider concatenating your blocks or using a different method.

    Args:
        blocks_map - maps the block to the argument name
    """

    def __init__(self, **blocks_map: ArgumentName):
        self.blocks_argument_map = blocks_map

    def fit(self, *args):
        return self

    def transform(self, blocks: Blocks):
        assert set(blocks.keys()) == set(self.blocks_argument_map.keys())
        return {
            self.blocks_argument_map[block_id]: block
            for block_id, block in blocks.items()
        }
