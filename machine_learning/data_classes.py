from dataclasses import field
from typing import List, TypeVar, Generic, Union

from pandas import DataFrame, Series
from pydantic import BaseConfig

from pydantic.dataclasses import dataclass


class ValidationConfig(BaseConfig):
    arbitrary_types_allowed = True


dataclass = dataclass(config=ValidationConfig)


@dataclass
class MultiBlockDataSet:
    data: List[Union[DataFrame, Series]]
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

    @property
    def class_imbalance(self) -> float:
        return self.binary_response.mean()

    @property
    def x(self):
        return self.data[0]

    @property
    def y(self):
        return self.data[1]

    def apply(self, func):
        for i in range(len(self.data)):
            self.data[i] = func(self.data[i])
        return self

    @property
    def observations(self) -> Series:
        """For example, patients in DA problems"""
        return self.data[0].index


ModelAttribute = TypeVar('ModelAttribute')


@dataclass
class AttributesStore(Generic[ModelAttribute]):
    """Stores model attributes such as coefficients per-block of the model.

    Probably a better design would build upon the blocks...
    """
    x: ModelAttribute
    y: ModelAttribute = field(default_factory=lambda: None)

    @classmethod
    def from_dicts(cls, coefficients_values, abundance_matrices, subclass, skip_compute=False):
        return cls(**{
            matrix: subclass(values, abundance=abundance_matrices[matrix], skip_compute=skip_compute)
            for matrix, values in coefficients_values.items()
        })
