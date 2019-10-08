from collections import defaultdict
from functools import partial
from math import sqrt, log2
from typing import Type, Dict, Union, Callable

import pandas as pd
from numpy import sign
from pandas import DataFrame, Series
from rpy2.robjects import r
from rpy2 import robjects
from tqdm.auto import tqdm

from helpers.p_values import gpd_p_value
from helpers.r import p_adjust, r_function, r_function_numpy

from .data_classes import AttributesStore, dataclass


CoefficientsGetter = Union[str, Callable]


@dataclass
class Coefficients:
    data: DataFrame
    non_zero_coeffs_count_by_cv: Series
    coeffs_across_cv: DataFrame

    def transform(self, df):
        return df

    def extract_null(self, null_distributions):
        return null_distributions['contributions']

    def add_quantiles(self):
        coeffs = self.data
        coeffs_across_cv = self.coeffs_across_cv
        coeffs['quantile_0.25'] = coeffs_across_cv.quantile(q=0.25, axis=1)
        coeffs['quantile_0.50'] = coeffs_across_cv.quantile(q=0.50, axis=1)
        coeffs['quantile_0.75'] = coeffs_across_cv.quantile(q=0.75, axis=1)

    def add_above_abs_quantile(self, q=0.5):
        coeffs = self.data
        coeffs_across_cv = self.coeffs_across_cv

        abs_coeffs_cv = coeffs_across_cv.abs()
        coeffs[f'above_abs_iteration_quantile_{q}'] = (
            abs_coeffs_cv.ge(
                abs_coeffs_cv.quantile(q=q, axis=0),
                axis='columns'
            )
        ).mean(axis=1)

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

        self.data = coeffs

        if null_distributions is not None:
            null_distribution = self.extract_null(null_distributions)
            self.add_significance(null_distribution, **p_value_kwargs or {})
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

    def add_permutation_significance(self, null_distribution: DataFrame, **kwargs):

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

    def add_hdi_significance(self, x, y, family='binomial', cores=4, suppress_group_testing=False, **kwargs):
        r('library("hdi")')
        lasso_proj = partial(r_function_numpy, 'lasso.proj')
        as_matrix = partial(r_function, 'as.matrix')
        robjects.globalenv['y'] = y.tolist()
        r('y = as.numeric(y)')
        result = lasso_proj(
            x=as_matrix(x.values), y=r['y'],
            family=family, parallel=True, ncores=cores,
            **{'suppress.grouptesting': suppress_group_testing, **kwargs}
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
        """Pass result.sub_sampling_test_results.cv_auc here"""
        coeffs = self.data
        coeffs_across_cv = self.coeffs_across_cv
        null_models = self.coeffs_across_cv.sum() == 0
        null_genes = self.coeffs_across_cv.sum(axis=1) == 0
        if null_models.any():
            print(null_models.sum(), 'null models detected')
            coeffs_across_cv = self.coeffs_across_cv.fillna(0)
        if null_genes.any():
            print(null_genes.sum(), 'genes not selected in any model detected')
        contribution = coeffs_across_cv.abs() / coeffs_across_cv.abs().sum()
        if null_models.any():
            contribution = contribution.fillna(0)
        # otherwise this is useless
        assert len(set(cv_auc)) > 1
        weighted_auc = (
            contribution.values * pd.concat([Series(cv_auc) for _ in range(len(coeffs))], axis=1).T.values
        ).sum(axis=1) / (coeffs_across_cv != 0).sum(axis=1)
        if null_genes.any():
            weighted_auc = weighted_auc.fillna(0)
        coeffs['weighted_auc'] = weighted_auc

    def select(self, non_zero_ratio_required=None):
        coeffs = self.data
        if non_zero_ratio_required is None:
            print('Adjusting non-zero ratio required to keep all significant coeffs')
            non_zero_ratio_required = coeffs[coeffs.fdr < 0.05].selected_in.min()
        selected_coeffs = coeffs[coeffs.selected_in >= non_zero_ratio_required]
        print('Selected', len(selected_coeffs), 'coefficients different from zero in at least', non_zero_ratio_required, 'repeats')
        return selected_coeffs

    def co_selection_network(self, contribution_threshold=0, selection_threshold=0):
        # frequent interactions
        counts_co_selected_graph = defaultdict(int)
        # strong interactions
        weights_graph = defaultdict(float)
        sign_agreement_graph = defaultdict(int)

        coeffs = self.coeffs_across_cv

        for cv_split_id in tqdm(coeffs.columns):
            cv_split = coeffs[cv_split_id]
            for i, protein in enumerate(cv_split.index):
                a = cv_split.loc[protein]

                if a == 0 or abs(a) < contribution_threshold:
                    continue

                if selection_threshold and self.data['selected_in'].loc[protein] < selection_threshold:
                    continue

                a_sign = sign(a)
                for other in cv_split.index[i + 1:]:

                    b = cv_split.loc[other]

                    pair = protein, other

                    if selection_threshold and self.data['selected_in'].loc[other] < selection_threshold:
                        continue

                    weights_graph[pair] += abs(a) + abs(b)
                    counts_co_selected_graph[pair] += 1
                    sign_agreement_graph[pair] += (a_sign == sign(b))

        df = DataFrame([
            {
                'a': other,
                'b': protein,
                'count': count,
                'frequency': count / len(coeffs.columns),
                'sign_agreement': sign_agreement_graph[protein, other],
                # average contribution of the pair when selected together
                'weight': weights_graph[protein, other] / count,
            }
            for (protein, other), count in counts_co_selected_graph.items()
        ])

        average_pair_frequency = (df.b.map(self.data['selected_in']) * df.a.map(self.data['selected_in']))
        df['log2_frequency_ratio'] = (
            df['frequency'] / average_pair_frequency
        ).apply(log2)

        average_pair_weight = (df.b.map(self.data['mean']) + df.a.map(self.data['mean'])) / 2
        df['relative_weight'] = df['weight'] / average_pair_weight

        return df


class Contributions(Coefficients):

    def transform(self, df):
        return df.abs() / df.abs().sum() * df.apply(sign)

    def extract_null(self, null_distributions):
        return null_distributions['coefficients']


class CoefficientsManager:

    def __init__(self, coefficients: Dict[str, CoefficientsGetter], abundance_matrices: Dict):
        self.abundance_matrices = abundance_matrices.copy()
        self.matrix_to_attribute = coefficients
        assert not (coefficients.keys() - {'x', 'y', 'combined'})
        self.values = {matrix: [] for matrix in coefficients}
        self.concatenated = False
        self.abundance_combined_added = False

    def add(self, split_matrices, pipeline):
        for matrix, attribute in self.matrix_to_attribute.items():
            # requested for coefficients of a single matrix
            if matrix in split_matrices:
                matrix_data = split_matrices[matrix]
            # requested for coefficients of combined matrices
            elif matrix == 'combined':
                matrix_data = pipeline.combine.transform(split_matrices)['Combined']
                if not self.abundance_combined_added:
                    combined = pipeline.combine.transform(self.abundance_matrices)
                    assert len(combined) <= 2   # maximum of one combined block + one outcome block allowed
                    self.abundance_matrices['combined'] = combined['Combined']
                    self.abundance_combined_added = True
            else:
                raise ValueError('Unknown matrix for coefficients extraction: ', matrix)
            self.values[matrix].append(
                pipeline.get_coefficients(matrix_data.columns, coefficient=attribute)
            )

    def concatenate(self):
        for matrix in self.matrix_to_attribute:
            self.values[matrix] = pd.concat(self.values[matrix], axis=1)
        self.concatenated = True

    def to_store(self, subclass: Type, skip_compute: bool):
        assert self.concatenated
        return AttributesStore.from_dicts(
            self.values, self.abundance_matrices, subclass=subclass, skip_compute=skip_compute
        )
