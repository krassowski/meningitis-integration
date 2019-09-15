from pandas import DataFrame
from sklearn.base import TransformerMixin
from sklearn.model_selection import KFold


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


def n_splits(cv, X):
    if isinstance(cv, int):
        cv = KFold(cv)
    return cv.get_n_splits(X)


def format_grid_results(reg, X=None, strip_last_step_prefix=True):
    statistics = {
        'split' + str(split) + '_test_' + statistic: statistic
        for statistic in reg.scoring.keys()
        for split in range(n_splits(reg.cv, X))
    }

    params = [k for k in reg.cv_results_.keys() if k.startswith('param_')]

    grid_results = DataFrame(reg.cv_results_).dropna()[
        [*params, *statistics.keys()]
    ].astype(float)
    grid_results = grid_results.melt(id_vars=params)

    if hasattr(reg.estimator, 'steps'):
        prefix = reg.estimator.steps[-1][0]
        grid_results = grid_results.rename(
            columns=lambda c: 'param_' + c[6 + len(prefix) + 2:] if c.startswith('param_' + prefix) else c
        )

    grid_results['scoring'] = grid_results.variable.map(statistics)
    grid_results['split'] = grid_results.variable.str.split('_').str[0].str.replace('split', '')
    return grid_results
