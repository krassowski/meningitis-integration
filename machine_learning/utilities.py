from pandas import DataFrame
from sklearn.base import TransformerMixin


def df_keeping(estimator):
    class Wrapped(estimator):
        def transform(self, x):
            new = super().transform(x)
            return DataFrame(new, index=x.index, columns=x.columns)
    Wrapped.__name__ = 'Pandas' + estimator.__name__
    return Wrapped


class Suffixer(TransformerMixin):
    """Adding unique suffixes to the variable names from different blocks is a good way of preventing conflicts

    when concatenating blocks or analyzing joint coefficients."""
    def __init__(self, suffix=''):
        self.suffix = suffix

    def fit(self, *args):
        return self

    def transform(self, block):
        return block.add_suffix(self.suffix)
