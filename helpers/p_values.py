from itertools import starmap
from warnings import warn

import numpy as np
from pandas import Series, DataFrame
from skgof import ks_test, ad_test
from scipy.stats import genpareto


def ecdf_p_value(test_statistics: Series, null_distribution: DataFrame, ecdf_pseudocount=1):
    # note: using biased p-value estimator as in https://www.ncbi.nlm.nih.gov/pubmed/21044043
    does_test_exceed_null = null_distribution.ge(test_statistics, axis='rows')
    N = len(null_distribution.columns)
    return (ecdf_pseudocount + does_test_exceed_null.sum(axis=1)) / (N + 1)


def gpd_p_value(test_statistics: Series, null_distribution: DataFrame, n=250, m=10, decrease_n_by=10, ecdf_pseudocount=1):
    """Partially-vectorized calculation of p-values using Empirical Cumulative Distribution Function (ECDF)
    with Generalized Pareto Distribution (GPD) correction for extreme values,
    usign algorithm proposed by (Knijnenburg, et. al, 2009), see:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687965/
    
    Maximum Likelihood Estimate is used to fit GPD.
    
    Args:
        test_statistics: vector of values of the test statistics for observations 
        null_distribution: null distribution obtained from permutation testing,
            containing iterations in columns and observations in rows
        n: initial number of permutation exceedances ($N_{\text{exc}}$) used to determine exceedances threshold ($t$),
            used to select data for fitting of GDP 
        m: the threshold for test statistic exceedances $M$ below which GDP is used instead of ECDF
        decrease_n_by: how much should n be decreased in each iteration if the GDP fit is not significant
        ecdf_pseudocount: small value added to ECDF p-value estimate to prevent p_values equal to zero
        
    Details:
        a formula using (N+1) in denominator of biased p-value (ECDF-based) estimator
        is used (as in https://www.ncbi.nlm.nih.gov/pubmed/21044043) rather than N as
        in the (Knijnenburg, et. al, 2009).
    """
    # all tested observations have to be in the null
    assert not (set(test_statistics.index) - set(null_distribution.index))
    
    if (
        len(null_distribution.index) != len(test_statistics.index)
        or
        any(null_distribution.index != test_statistics.index)
    ):
        # reorder and limit null so that it matches tested observations
        null_distribution = null_distribution.loc[test_statistics.index]
    
    # TODO:
    #assert not any(null_distribution.isnull())

    # number of permutations
    N = len(null_distribution.columns)
    
    assert n < N

    # indicator vector
    does_test_exceed_null = null_distribution.ge(test_statistics, axis='rows')

    # count of test statistic exceedances
    M = does_test_exceed_null.sum(axis=1)

    if len(M) == 0:
        print('No observations for p-value estimation!')
        return DataFrame(
            columns={
                'p_value': np.nan,
                'is_gpd': np.nan,
                'gof': np.nan,
                'ecdf_p': np.nan
            }.keys()
        )

    if min(M) >= m:
        print(f'Got enough permutations, no observations with fewer exceedances than {m}, falling back to ECDF')
        return DataFrame({
            'p_value': ecdf_p_value(test_statistics, null_distribution, ecdf_pseudocount),
            'is_gpd': False,
            'gof': np.nan,
            'ecdf_p': np.nan
        })
    
    def _p(test_i, null_i, M_i, d_i):
        gpd_fit = None
        gpd_fit_p_value = None

        n_i = n
        
        # TODO: no need to sort as much as N numbers, do partial sort:
        #  but this requires some tests (both performance and unit)
        # null_i_partitioned = np.partition(null_i, n_i+1)
        # null_i_first_n_sorted = sorted(null_i_partitioned[:-n_i+1])
        null_i = sorted(null_i)
        t = None
        
        if all(np.isnan(null_i)):
            return np.nan, False, np.nan, np.nan
        
        # compute ecdf based, biased estimate of p-value
        raw_ecdf_estimate = (ecdf_pseudocount + d_i.sum()) / (N + 1)
        
        if M_i < m:
            # fit GDP, reducing $n$ until convergance
            while n_i > 0:
                
                # -1 because Python has 0-based indexing
                t = (null_i[-n_i-1] + null_i[-n_i-2]) / 2
                
                y_untill_n = null_i[-n_i:]
                exceedences = y_untill_n - t

                assert all(y_untill_n >= t)
                assert len(exceedences) == n_i
                
                fit = genpareto.fit(exceedences)
                fitted = genpareto(*fit)
                gpd_fit = fitted
                
                gpd_fit_p_value = ad_test(exceedences, fitted).pvalue

                if gpd_fit_p_value <= 0.05:
                    break
                else:
                    n_i -= decrease_n_by

        if gpd_fit and gpd_fit_p_value < 0.05:
            return n_i / N * (1 - gpd_fit.cdf(test_i - t)), True, gpd_fit_p_value, raw_ecdf_estimate
        else:
            if gpd_fit:
                # TODO: get index and highlight which observation could not be fitted!
                warn(f'A good GPD fit could not be reached, using ECDF estimate instead')
            
            return raw_ecdf_estimate, False, np.nan, raw_ecdf_estimate

    p_values = DataFrame(starmap(
        _p,
        list(zip(
            test_statistics.values, null_distribution.values,
            M.values, does_test_exceed_null.values
        ))
    ), columns=['p_value', 'is_gpd', 'gof', 'ecdf_p'], index=test_statistics.index)
    return p_values
