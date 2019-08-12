from pandas import DataFrame


def df_keeping(estimator):
    class Wrapped(estimator):
        def transform(self, x):
            new = super().transform(x)
            return DataFrame(new, index=x.index, columns=x.columns)
    Wrapped.__name__ = 'Pandas' + estimator.__name__
    return Wrapped
