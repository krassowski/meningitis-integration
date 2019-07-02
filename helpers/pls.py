from pandas import DataFrame
import pandas as pd
from sklearn.cross_decomposition import PLSRegression, PLSCanonical


def pls_wrapper(pls):
    
    class PLSPandasMixin(pls):

        def fit(self, x, y):
            self.x = x
            self.y = y
            return super().fit(x, y)

        def transform(self, x, y):
            #assert all(x.index == self.x.index) and all(y.index == self.y.index)
            T, U = super().transform(x, y)
            return DataFrame(T, index=x.index), DataFrame(U, index=y.index)

        @property
        def W(self):
            return DataFrame(self.x_weights_, index=self.x.columns)

        @property
        def C(self):
            return DataFrame(self.y_weights_, index=self.x.columns)

        @property
        def P(self):
            return DataFrame(self.x_loadings_, index=self.x.columns)

        @property
        def Q(self):
            return DataFrame(self.y_loadings_, index=self.y.columns)
    
    PLSPandasMixin.__name__ = 'Pandas' + pls.__name__
    return PLSPandasMixin


PandasPLSRegression = pls_wrapper(PLSRegression)
PandasPLSCanonical = pls_wrapper(PLSCanonical)


def format_grid_results(reg):
    statistics = {
        'split' + str(split) + '_test_' + statistic: statistic
        for statistic in reg.scoring.keys()
        for split in range(reg.cv.get_n_splits())
    }
    grid_results = DataFrame(reg.cv_results_).dropna()[
        ['param_n_components', *statistics.keys()]
    ].astype(float)
    grid_results = grid_results.melt(id_vars='param_n_components')
    grid_results['scoring'] = grid_results.variable.map(statistics)
    grid_results['split'] = grid_results.variable.str.split('_').str[0].str.replace('split', '')
    return grid_results


def custom_scorer(func, greater_is_better=False):
    if not greater_is_better:
        return lambda e, x, y: -func(e, x, y)
    return func


def scores_by_component(scores_map):
    return pd.concat([
        pd.concat([
            DataFrame(scores[component]).rename(columns={component: kind})
            for kind, scores in scores_map.items()
        ], axis=1).rename_axis('gene').reset_index().assign(component=component)
        for component in range(3)
    ])

def format_test_train_scores(scores_train_X, scores_train_Y, scores_test_X, scores_test_Y):
    return pd.concat([
        scores_by_component({'T': scores_train_X, 'U': scores_train_Y}).assign(split='train'),
        scores_by_component({'T': scores_test_X, 'U': scores_test_Y}).assign(split='test')
    ])