from pandas import Series
from sklearn.base import TransformerMixin
from sklearn.base import BaseEstimator
from types import FunctionType


class ConditionsVectorizer(BaseEstimator):

    def __init__(self, mapping, case_class):
        self.mapping = mapping
        self.case_class = case_class

    def fit(self, x, y=None):
        return self

    def transform(self, x: Series):
        classes = x.map(self.mapping)
        assert len(set(classes)) == 2
        return classes == self.case_class


class Filter(TransformerMixin):
    
    def __init__(self, verbose=True):
        self.verbose = verbose
        self.to_filter_out = None

    def transform(self, x):
        to_filter_out_columns = x.columns.isin(self.to_filter_out)
        if self.verbose:
            name = self.__class__.__name__
            print(
                f'{name}: filtering out {sum(to_filter_out_columns)}/'
                f'{len(self.to_filter_out)} variables'
            )
        if self.verbose > 1 and sum(to_filter_out_columns):
            print(list(x.columns[to_filter_out_columns]))

        return x[x.columns[~to_filter_out_columns]]


class StaticFilter(Filter):
    """Will always filter out the same rows"""
    pass


class OutliersFilter(StaticFilter):
    
    def __init__(self, outlier_patients, verbose=True):
        self.outlier_patients = outlier_patients or []
        super().__init__(verbose=verbose)

    def fit(self, x, y=None):
        return self

    def transform(self, x):
        outlier_rows = x.index.isin(self.outlier_patients)
        
        if self.verbose:
            name = self.__class__.__name__
            print(f'{name}: filtering out {sum(outlier_rows)} outliers')

        return x.loc[~outlier_rows]

    
class DynamicFilter(Filter):
    """filter out different rows, depending on the provided data"""
    pass


class KeepFilter(StaticFilter):

    def __init__(self, to_keep, verbose=True):
        self.to_keep = to_keep
        super().__init__(verbose=verbose)

    def fit(self, x, y=None):
        self.to_filter_out = x.columns[
            ~x.columns.isin(self.to_keep)
        ]
        return self


class PreFilterLowestExpresion(KeepFilter):

    def __init__(self, func, data, smallest_group_n, verbose=True):
        self.n = smallest_group_n
        self.original_data = data
        to_keep = data.index[
            Series(func(data, self.n, return_frame=False)) == 1
            ]
        super().__init__(
            to_keep=to_keep,
            verbose=verbose
        )


class LowCountsFilter(DynamicFilter):

    def __init__(self, ratio=1/3, verbose=True):
        self.ratio = ratio
        super().__init__(verbose=verbose)
    
    def fit(self, x, y=None):
        self.to_filter_out = x.columns[
            (x != x.median()).sum() <= self.ratio * len(x)
        ]
        return self

    
class LowVarianceFilter(DynamicFilter):

    def __init__(self, quantile=0.001, verbose=True):
        self.quantile = quantile
        super().__init__(verbose=verbose)
    
    def fit(self, x, y=None):
        self.to_filter_out = x.columns[
            x.var() < x.var().quantile(self.quantile)
        ]
        return self


class RSideNormalizer(TransformerMixin):
    """Normalizer with the data living in R.
    
    Each copy procedure is costly, thus we do not copy the data into R;
    Instead the data are stored in R and subsetted in R,
    and only the subsetted and normalized data frame is copied to Python.

    As this might be prone to errors - caution is recommended.
    """

    def __init__(self, func: FunctionType, omic: str, **kwargs):
        self.omic = omic
        self.kwargs = kwargs
        self.func = func
    
    def fit(self, data, y=None):
        self.kwargs['subset'] = data.index
        self.kwargs['subset_rows'] = data.columns
        return self
    
    def transform(self, data):
        assert len(data.columns)
        self.kwargs['subset'] = data.index
        # NOTE: subset_rows was previously commented out for some reason!
        self.kwargs['subset_rows'] = data.columns
        # here we go to R and back
        normalized = self.func(self.omic, **self.kwargs).T
        return normalized
