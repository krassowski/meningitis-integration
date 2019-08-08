from dataclasses import dataclass
from functools import partial
from math import sqrt

import pandas as pd
from pandas import DataFrame, Series
from rpy2.robjects import r

from helpers.p_values import gpd_p_value
from helpers.r import p_adjust, r_function, r_function_numpy


@dataclass
class Coefficients:
    data: DataFrame
    non_zero_coeffs_count_by_cv: Series
    coeffs_across_cv: DataFrame

    def transform(self, df):
        return df

    def extract_null(self, null_distributions):
        return null_distributions['mean_contributions']

    def __init__(
        self, coeffs_across_cv: DataFrame, abundance: DataFrame = None,
        null_distributions: dict = None, p_value_kwargs: dict = None,
        skip_compute=False
    ):
        coeffs = coeffs_across_cv.mean(axis=1).to_frame(name='mean')
        if skip_compute:
            self.data = coeffs
            return

        self.non_zero_coeffs_count_by_cv = (coeffs_across_cv != 0).sum(axis=0)

        coeffs['selected_in'] = (coeffs_across_cv != 0).mean(axis=1)
        coeffs['positive_in'] = (coeffs_across_cv > 0).mean(axis=1)
        coeffs['negative_in'] = (coeffs_across_cv < 0).mean(axis=1)

        coeffs['volatile'] = (
            (2 * coeffs[['positive_in', 'negative_in']].min(axis=1) / coeffs['selected_in'])
        ).fillna(0)

        coeffs_across_cv = self.transform(coeffs_across_cv)
        self.coeffs_across_cv = coeffs_across_cv

        coeffs['mean'] = coeffs_across_cv.mean(axis=1)
        coeffs['stdev'] = coeffs_across_cv.std(axis=1)
        coeffs['ci'] = 1.96 * coeffs_across_cv.std(axis=1) / sqrt(len(coeffs_across_cv))

        coeffs['quantile_0.25'] = coeffs_across_cv.quantile(q=0.25, axis=1)
        coeffs['quantile_0.50'] = coeffs_across_cv.quantile(q=0.50, axis=1)
        coeffs['quantile_0.75'] = coeffs_across_cv.quantile(q=0.75, axis=1)

        self.data = coeffs

        if null_distributions is not None:
            null_distribution = self.extract_null(null_distributions)
            self.add_significance(null_distribution, **p_value_kwargs)
        else:
            assert not p_value_kwargs

        if abundance is not None:
            self.add_abundance(abundance)

        self.order_genes()

    def add_abundance(self, abundance):
        self.data['mean_abundance'] = abundance.sum().loc[self.data.index].values

    def order_genes(self, by='mean'):
        self.data['gene'] = pd.Categorical(
            self.data.index,
            categories=self.data[by].sort_values().index,
            ordered=True
        )

    def add_permutation_significance(self, null_distribution: DataFrame = None, **kwargs):

        coeffs = self.data
        # reorder
        null_distribution = null_distribution.loc[coeffs.index]

        positive = coeffs[coeffs['mean'] >= 0]
        negative = coeffs[coeffs['mean'] < 0]

        p_values = pd.concat([
            gpd_p_value(
                coeffs['mean'][positive.index],
                null_distribution.loc[positive.index],
                **kwargs
            ),
            gpd_p_value(
                -coeffs['mean'][negative.index],
                -null_distribution.loc[negative.index],
                **kwargs
            )
        ]).loc[coeffs.index]

        assert (p_values.index == coeffs.index).all()

        coeffs[p_values.columns] = p_values
        coeffs['fdr'] = p_adjust(coeffs['p_value'], 'BH')

    def add_hdi_significance(self, x, y, family='binomial', cores=4):
        r('library("hdi")')
        lasso_proj = partial(r_function_numpy, 'lasso.proj')
        as_matrix = partial(r_function, 'as.matrix')
        from rpy2 import robjects
        robjects.globalenv['y'] = y.tolist()
        r('y = as.numeric(y)')
        result = lasso_proj(
            x=as_matrix(x.values), y=r['y'],
            family=family, parallel=True, ncores=cores
        )
        self.data[f'{family}_p'] = Series(result.rx2('pval'), index=x.columns)
        self.data[f'{family}_FDR'] = p_adjust(self.data[f'{family}_p'], 'BH')

    def add_significance(self, method='hdi', **kwargs):
        assert method in {'hdi', 'permutation'}

        if method == 'permutation':
            self.add_permutation_significance(**kwargs)
        else:
            self.add_hdi_significance(**kwargs)

    def add_weighted_auc(self, cv_auc: DataFrame):
        coeffs = self.data
        contribution = self.coeffs_across_cv.abs() / self.coeffs_across_cv.abs().sum()
        coeffs['weighted_auc'] = (
            contribution.values * pd.concat([Series(cv_auc) for _ in range(len(coeffs))], axis=1).T.values
        ).mean(axis=1)

    def select(self, non_zero_ratio_required=0.05):
        coeffs = self.data
        selected_coeffs = coeffs[coeffs.selected_in >= non_zero_ratio_required]
        print('Selected', len(selected_coeffs), 'coefficients different from zero in at least', non_zero_ratio_required, 'repeats')
        return selected_coeffs


class Contributions(Coefficients):

    def transform(self, df):
        return df.abs() / df.abs().sum()

    def extract_null(self, null_distributions):
        return null_distributions['mean_coefficients']
