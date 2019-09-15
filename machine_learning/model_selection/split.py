import random
from abc import ABC, abstractmethod
from typing import Callable, Iterable, TypeVar, Generic, Union, Tuple
from warnings import warn

import numpy
from numpy import linspace
from pandas import DataFrame, Series
from sklearn.model_selection import KFold, ShuffleSplit
from tqdm.auto import tqdm

SplitRatioSequenceGenerator = Callable[[float, float, int], Iterable[float]]
SplitRatioGenerator = Callable[[float, float], float]
ArrayLike = TypeVar('ArrayLike', DataFrame, numpy.array, Series)
Split = Tuple[numpy.ndarray, numpy.ndarray]
Validate = Callable[[Split, ArrayLike, ArrayLike], bool]


class SplitterInterface(ABC, Generic[ArrayLike]):

    @abstractmethod
    def __init__(self, split_specification: Union[float, int]):
        pass

    @abstractmethod
    def split(self, x: ArrayLike, y: ArrayLike = None) -> Iterable[Tuple[numpy.ndarray, numpy.ndarray]]:
        pass


def always_ok(split, x, y):
    return True


def at_least_n_in_each_class(values, n):
    return values.min() >= n


def class_number_between(values, min_n, max_n):
    return min_n <= len(values) <= max_n


def classification_split_validator(min_class_members_n: int, min_classes_n: int, max_classes_n: int) -> Validate:
    """
    Creates a validator for classification problems, checking the quantity of classes in the response vector (y)
    and the number of observations in each of the classes, for both train and test data.

    Useful if you very few patients and are in danger of having a split with no patients of given class.
    """
    def validate(split, x, y):
        train, test = split
        train_binary, test_binary = y.loc[train], y.loc[test]
        train_values = train_binary.value_counts(sort=False)
        test_values = test_binary.value_counts(sort=False)
        return all([
            class_number_between(train_values, min_classes_n, max_classes_n),
            class_number_between(test_values, min_classes_n, max_classes_n),
            at_least_n_in_each_class(train_values, min_class_members_n),
            at_least_n_in_each_class(test_values, min_class_members_n)
        ])
    return validate


class ValidatingSplitter:
    """Splitter with validation (allows to avoid splits with too few observations of any class)
    """

    is_valid: Validate

    def __init__(self, cv=None, validator: Validate = None):
        self.cv_instance = cv
        self.is_valid = validator or always_ok

    def split(self, x: ArrayLike, y: ArrayLike) -> Iterable[Split]:
        cv = self.cv_instance
        for split in cv.split(x, y):
            if not self.is_valid(split, x, y):
                warn(f'Invalid split - skipping')
            else:
                yield split


class _MultiSizeSplits(ValidatingSplitter):
    """Multi-scale repeated splits for an arbitrary randomized CV splitter."""

    def __init__(
        self, min_max: Tuple[float, float], n_values: int, progress_bar: bool = True,
        single_split_generator: SplitRatioGenerator = None,
        sequence_generator: SplitRatioSequenceGenerator = None,
        validator: Validate = None, n_retrials=100, **kwargs
    ):
        """
        Args:
            min_max: Minimum and maximum value for the generator
            progress_bar: whether to display a progress bar (requires tqdm)
            **kwargs: additional arguments to be passed to the CV splitter.

        Attributes:
            values: the generated split parameter values
        """
        super().__init__(validator=validator)
        assert len(min_max) == 2
        min_size, max_size = min_max
        assert max_size >= min_size
        assert single_split_generator or sequence_generator and not (single_split_generator and sequence_generator)
        self.progress_bar = progress_bar
        self.kwargs = kwargs
        self.min_max = min_max
        self.n_values = n_values
        if sequence_generator:
            self.values = sequence_generator(min_size, max_size, self.n_values)
        else:
            self.values = [
                single_split_generator(min_size, max_size)
                for _ in range(self.n_values)
            ]
        self.single_split_generator = single_split_generator
        self.n_retrials = n_retrials

    @abstractmethod
    def _cv(self, size, **kwargs) -> SplitterInterface:
        """Create the splitter"""

    def split(self, x: ArrayLike, y: ArrayLike) -> Iterable[Tuple[numpy.ndarray, numpy.ndarray]]:
        split_sizes = self.values
        if self.progress_bar:
            split_sizes = tqdm(split_sizes, total=self.n_values)
        for size in split_sizes:
            for split in self._split_and_check(x, y, size, self.n_retrials):
                yield split

    def _split_and_check(self, x: ArrayLike, y: ArrayLike, size: float, n_retrials: int, history=None):
        if not history:
            history = []
        if n_retrials == 0:
            raise ValueError(
                f'Could not find a split of size {history}'
                f'satisfying the validation ({n_retrials} attempts);'
            )
        cv = self._cv(size, **self.kwargs)
        splits = []
        for split in cv.split(x, y):
            if self.is_valid(split, x, y):
                splits.append(split)
            else:
                if not self.single_split_generator:
                    warn(f'Invalid split for size {size}, and no single split generator - skipping')
                else:
                    min_size, max_size = self.min_max
                    new_size = self.single_split_generator(min_size, max_size)
                    history.append(size)
                    return self._split_and_check(new_size, x, y, n_retrials - 1, history)
        return splits


class MultiSizeFolds(_MultiSizeSplits):

    def __init__(
        self, k_number=5, k_range=(3, 5), generator: SplitRatioSequenceGenerator = range, k_fold_cv=KFold,
        **kwargs
    ):
        """
        Args:
            k_number: the number of k values to be generated
            k_range: a tuple of `(min_k, max_k)` where min_k is the minimal and max_k is the maximal k,
                max_k > min_k, and min_k >= 2
            generator: a function that given `(min_k, max_k, k_number)` returns an iterable with chosen k-values, e.g.:
                - `range` for an exhaustive search
                - `lambda start, stop, n: random.choices(range(start, stop), k=n)` for random sampling of the fold sizes
            k_fold_cv: a CV splitter accepting n_splits argument representing the number of folds
            **kwargs: additional arguments passed to the k_fold_cv and to the parent class
        """
        min_size, max_size = k_range
        assert min_size >= 2
        super().__init__(min_max=k_range, n_values=k_number, sequence_generator=generator, **kwargs)
        self.k_fold_cv = k_fold_cv

    def _cv(self, k: int, **kwargs):
        return self.k_fold_cv(n_splits=k, **kwargs)


class IndependentMultiSizeFolds(_MultiSizeSplits):

    def __init__(
        self, k_number=5, k_range=(3, 5), generator: SplitRatioGenerator = random.randint, k_fold_cv=KFold,
        **kwargs
    ):
        """
        Args:
            k_number: the number of k values to be generated
            k_range: a tuple of `(min_k, max_k)` where min_k is the minimal and max_k is the maximal k,
                max_k > min_k, and min_k >= 2
            generator: a function that given `(min_k, max_k)` returns a random k-value
            k_fold_cv: a CV splitter accepting n_splits argument representing the number of folds
            **kwargs: additional arguments passed to the k_fold_cv and to the parent class
        """
        min_size, max_size = k_range
        assert min_size >= 2
        super().__init__(min_max=k_range, n_values=k_number, single_split_generator=generator, **kwargs)
        self.k_fold_cv = k_fold_cv

    def _cv(self, k: int, **kwargs):
        return self.k_fold_cv(n_splits=k, **kwargs)


class MultiSizeShuffleSplits(_MultiSizeSplits):
    def __init__(
        self, n_sizes=10, test_sizes=(0.1, 0.3), generator: SplitRatioSequenceGenerator = linspace,
        shuffle_cv=ShuffleSplit, **kwargs
    ):
        """
        Args:
            n_sizes: the number of different split ratios to be generated
            test_sizes: the range of test sizes to be used for split ratios generation;
                a tuple of `(min_test_size, max_test_size)` with both values in (0, 1) range
            generator: a function that given `(min_test_size, max_test_size, n_sizes)` returns an iterable with chosen
                test sizes, for example:
                    - `np.linspace` - exhaustive split space traversal with linear sampling
                    - `np.logspace` - exhaustive split space traversal with log sampling
                    - `lambda start, stop, n: np.random.standard_normal(n) * stop - start`
                      useful if any specific split ratio is more desired than others
                    - `lambda start, stop, n: random.choices(np.linspace(start, stop, num=x), n)`
                      allows to constraint sampling resolution by changing num, for example using:
                        `start=0.25, end=0.5, num=5`
                      the possible values are:
                        `{0.25, 0.3125, 0.375, 0.4375, 0.5}`
            shuffle_cv: a CV splitter accepting test_size argument
            **kwargs: additional arguments passed to the shuffle_cv and to the parent class
        """
        min_size, max_size = test_sizes
        assert 1 > min_size > 0 and 1 > max_size > 0
        super().__init__(min_max=test_sizes, n_values=n_sizes, sequence_generator=generator, **kwargs)
        self.shuffle_cv = shuffle_cv

    def _cv(self, test_size: float, **kwargs):
        return self.shuffle_cv(test_size=test_size, **kwargs)


class IndependentMultiSizeShuffleSplits(_MultiSizeSplits):
    def __init__(
        self, n_sizes=10, test_sizes=(0.1, 0.3), generator: SplitRatioGenerator = random.uniform,
        shuffle_cv=ShuffleSplit, **kwargs
    ):
        """
        Args:
            n_sizes: the number of different split ratios to be generated
            test_sizes: the range of test sizes to be used for split ratios generation;
                a tuple of `(min_test_size, max_test_size)` with both values in (0, 1) range
            generator: a function that given (min_test_size, max_test_size) returns a test sizes, for example:
                - `random.uniform` for continuous, uniform sampling
                - `lambda start, stop: random.choice(np.linspace(start, stop, num=x))`
            shuffle_cv: a CV splitter accepting test_size argument
            **kwargs: additional arguments passed to the shuffle_cv and to the parent class
        """
        min_size, max_size = test_sizes
        assert 1 > min_size > 0 and 1 > max_size > 0
        super().__init__(min_max=test_sizes, n_values=n_sizes, single_split_generator=generator, **kwargs)
        self.shuffle_cv = shuffle_cv

    def _cv(self, test_size: float, **kwargs):
        return self.shuffle_cv(test_size=test_size, **kwargs)
