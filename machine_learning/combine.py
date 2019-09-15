from abc import ABC, abstractmethod
from typing import List

from pandas import concat
from sklearn.base import TransformerMixin

from machine_learning.data_classes import Blocks


class BlocksCombiner(ABC, TransformerMixin):

    @abstractmethod
    def transform(self, blocks: Blocks) -> Blocks:
        pass


class BlocksConcatenation(BlocksCombiner):
    """Similar to sklearn.compose.FeatureUnion"""

    def __init__(self, blocks: List[str], axis=1, verbose=False):
        self.blocks = blocks
        self.verbose = verbose
        self.axis = axis

    def fit(self, *args):
        return self

    def transform(self, blocks: Blocks):
        blocks_to_concatenate = []
        remaining_blocks = {}

        for name, block in blocks.items():
            if name in self.blocks:
                blocks_to_concatenate.append(block)
            else:
                remaining_blocks[name] = block

        return {
            'Combined': concat(blocks_to_concatenate, axis=self.axis),
            **remaining_blocks
        }
