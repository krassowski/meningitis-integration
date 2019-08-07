from contextlib import contextmanager, redirect_stdout
import re
from functools import partial, reduce
from inspect import signature
from operator import mul
from statistics import mean
from types import SimpleNamespace
from typing import Optional
from io import StringIO
from pandas import DataFrame, Series

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import numpy as np

from helpers.r import r_function
from dataclasses import dataclass


# load the overrides
r_function('source', 'thirdparty/OmicsPLS_overrides.R')


@dataclass(frozen=True)
class Statistic:
    name: str
    symbol: str
    func: callable
    equation: Optional[str] = ''


def as_matrix(df):
    return (
        df.values
        if hasattr(df, 'values') else
        df
    )


@dataclass(frozen=True)
class O2PLSStatistic(Statistic):
    model: str = None  # X, Y or XY
    
    statistics = {}
    
    def calc(self, estimator, X, Y):
        if hasattr(estimator, '_final_estimator'):
            estimator = estimator._final_estimator
        if len(signature(self.func).parameters) > 1:
            return self.func(estimator, as_matrix(X), as_matrix(Y))
        else:
            return self.func(estimator)


def o2pls_statistic(*args, **kwargs):
    
    def decorator(func):
        statistic = O2PLSStatistic(*args, func=func, **kwargs)
        O2PLSStatistic.statistics[func.__name__] = statistic
        return statistic
    
    return decorator


x_o2pls_statistic = partial(o2pls_statistic, model='X')
y_o2pls_statistic = partial(o2pls_statistic, model='Y')


def ratio_of_sums_of_squares(nominator, denominator):
    return (
        sum(sum(pow(np.atleast_2d(as_matrix(nominator)), 2)))
        /
        sum(sum(pow(np.atleast_2d(as_matrix(denominator)), 2)))
    )


def modeled_variation(residuals, matrix):
    return 1 - ratio_of_sums_of_squares(residuals, matrix)


def combine_statistic(func, a, b, how='average', m1=True):
    if isinstance(a, DataFrame):
        a = a.values
    if isinstance(b, DataFrame):
        b = b.values
    a = np.atleast_2d(a)
    b = np.atleast_2d(b)
    try:
        if how == 'average':
            return mean([func(a[i, :], b[i, :]) for i in range(a.shape[0])])
        else:
            return 1 - product([1 - func(a[i, :], b[i, :]) for i in range(a.shape[0])])
    except:
        print(a)
        print(b)
        print(func)
        raise


def cum_modeled_variation(residuals, matrix):
    return cumulative_statistic(modeled_variation, residuals, matrix)


@x_o2pls_statistic('Modeled variation of $X$', 'R^2X', equation=r'1 - \sum(E)^2 / \sum X^2')
def modeled_variation_of_X(estimator):
    return modeled_variation(residuals=estimator.E, matrix=estimator.X)


@y_o2pls_statistic('Modeled variation of $Y$', 'R^2Y', equation=r'1 - \sum(F)^2 / \sum Y^2')
def modeled_variation_of_Y(estimator):
    return modeled_variation(residuals=estimator.F, matrix=estimator.Y)


def r2(true, pred, check_shape=True, detract_mean=False):
    if check_shape:
        assert true.shape == pred.shape

    return (
        1 - ratio_of_sums_of_squares(
            nominator=pred - true,
            # note: the original paper (Trygg, Wold 2002) did not detract mean
            # read more here: https://stats.stackexchange.com/a/26205/233900
            # on such a type of R2
            denominator=true - (true.mean(axis=0) if detract_mean else 0)
        )
    )


def product(seq):
    return reduce(mul, seq, 1)


def cum_r2(true, pred):
    # TODO: use # cumulative metric as in
    #  https://www.bioconductor.org/packages/devel/bioc/vignettes/ropls/inst/doc/ropls-vignette.html ?
    # return combine_statistic(r2, true, pred, 'cumulative')
    return r2(true, pred)


@y_o2pls_statistic('Predicted variation of $X$', 'R^2X_{hat}')
def predicted_variation_of_X(estimator):
    return cum_r2(
        true=estimator.X,
        pred=estimator.predict_x(estimator.Y)
    )


@x_o2pls_statistic('Predicted variation of $Y$', 'R^2Y_{hat}')
def predicted_variation_of_Y(estimator):
    return cum_r2(
        true=estimator.Y,
        pred=estimator.predict_y(estimator.X)
    )


@y_o2pls_statistic('Cross-validated predictions of $X$', 'Q^2X_{hat}')
def cross_validated_predictions_of_X(estimator, X, Y):
    return cum_r2(
        true=X,
        pred=estimator.predict_x(Y)
    )


@x_o2pls_statistic('Cross-validated predictions of $Y$', 'Q^2Y_{hat}')
def cross_validated_predictions_of_Y(estimator, X, Y):
    return cum_r2(
        true=Y,
        pred=estimator.predict_y(X)
    )


@x_o2pls_statistic('Modeled structured noise variation of $X$', 'R^2X_{YO}')
def modeled_structured_noise_variation_of_X(estimator):
    return ratio_of_sums_of_squares(
        nominator=estimator.T_Yosc @ estimator.P_Yosc.T,
        denominator=estimator.Y
    )


@y_o2pls_statistic('Modeled structured noise variation of $Y$', 'R^2Y_{XO}')
def modeled_structured_noise_variation_of_Y(estimator):
    return ratio_of_sums_of_squares(
        nominator=estimator.U_Xosc @ estimator.P_Xosc.T,
        denominator=estimator.Y
    )


@x_o2pls_statistic('Modeled joint X-Y covariation of $X$', 'R^2X_{corr}')
def modeled_joint_XY_covariation_of_X(estimator):
    return ratio_of_sums_of_squares(
        nominator=estimator.T @ estimator.W.T,
        denominator=estimator.X
    )


@y_o2pls_statistic('Modeled joint X-Y covariation of $Y$', 'R^2Y_{corr}')
def modeled_joint_XY_covariation_of_Y(estimator):
    return ratio_of_sums_of_squares(
        nominator=estimator.U @ estimator.C.T,
        denominator=estimator.Y
    )


# meta-statistics

@o2pls_statistic(
    'Average of cross-validated predictions of $X$ and $Y$',
    '\\frac{Q^2X_{hat} + Q^2Y_{hat}}{2}',
    model='XY'
)
def average_cv_predictions(estimator, X, Y):
    return (
        cross_validated_predictions_of_X.calc(estimator, X, Y)
        + cross_validated_predictions_of_Y.calc(estimator, X, Y)
    ) / 2


# o2m2 = partial(r_function, 'o2m')
o2m2 = partial(r_function, 'fixed_o2m2')
svd = partial(r_function, 'fast_svd')
o2m = partial(r_function, 'fixed_o2m')
# predict_o2m = partial(r_function, 'predict')
predict_o2m = partial(r_function, 'predict_o2m')
r_matrix = partial(r_function, 'as.matrix')
        

def named_vector_to_dict(x):
    assert len(set(x.names)) == len(x.names)
    with localconverter(ro.default_converter + pandas2ri.converter):
        return {
            k: ro.conversion.rpy2py(v)
            for k, v in x.items()
        }


@contextmanager
def silent_output(patterns):
    f = StringIO()
    silenced_lines = []
    try:
        with redirect_stdout(f):
            yield silenced_lines
    except Exception:
        raise
    finally:
        for line in f.getvalue().split('\n'):
            if not line:
                continue
            if not re.match('|'.join(patterns), line):
                print(line)


def maybe_add_columns(values, source):
    return DataFrame(
        values,
        columns=(
            source.columns
            if isinstance(source, DataFrame) else
            None
        )
    )


def enforce_2d(result):
    result = np.atleast_2d(result)
    if result.shape[0] == 1:
        result = result.T
    return result


class O2PLS:
    
    attributes_map = {
        'T': 'Tt',   # joint X-scores
        'U': 'U',    # joint Y-scores
        # not implemented for o2m2, see https://github.com/selbouhaddani/OmicsPLS/issues/11
        # 'E': 'E',    # joint X-residuals
        # 'F': 'Ff',    # joint Y-residuals
        'W': 'W.',   # joint X-loadings
        'C': 'C.',   # joint Y-loadings
        'P_Yosc': 'P_Yosc.',   # orthogonal X-loadings
        'P_Xosc': 'P_Xosc.',   # orthogonal Y-loadings
        'T_Yosc': 'T_Yosc',   # orthogonal X-scores
        'U_Xosc': 'U_Xosc',   # orthogonal Y-scores
        'B_U': 'B_U',  # weights/coefficients for X prediction
        'B_T': 'B_T.',  # weights/coefficients for Y prediction
        'C_Xosc': 'C_Xosc',
        'W_Yosc': 'W_Yosc'
    }
    
    def __init__(self, max_iterations=200, algorithm='nipals', verbose=False, save_memory=False, **params):
        """
        Args:
            max_iterations
        """
        self.params = {
            'max_iterations': max_iterations,
            'algorithm': algorithm,
            'verbose': verbose,
            **params
        }
        self.fitted = False
        self.metric = SimpleNamespace(**{
            name: partial(statistic.func, self)
            for name, statistic in O2PLSStatistic.statistics.items()
        })

    def calc_metrics(self, key='symbol', **kwargs):
        metrics = {}
        for name, statistic in O2PLSStatistic.statistics.items():
            try:
                k = getattr(statistic, key)
                if key == 'symbol':
                    k = '$' + k + '$'
                metrics[k] = statistic.func(self, **kwargs)
                # missing arguments
            except TypeError:
                pass
        return Series(metrics).to_frame().T
        
    def fit(self, X, Y):
        if self.params['verbose']:
            print('Fitting started...')
        Y = np.atleast_2d(Y)
        self.X = X
        self.Y = Y

        assert X.shape[0] == Y.shape[0]

        messages_to_silence = [
            r'Power Method \(comp \d+\) stopped after (\d+) iterations\.',
            'Data is not centered, proceeding...'
        ]
        algorithm = self.params['algorithm']
        assert algorithm in ['nipals', 'svd', 'default_svd']
        r_implementation = o2m2 if algorithm == 'nipals' else svd
        if algorithm == 'default_svd':
            r_implementation = o2m
        # TODO: test if \1 > max_iterations, warn about lack of convergence

        with silent_output(messages_to_silence) as silenced_lines:
            self.r_model = r_implementation(
                as_matrix(X), r_matrix(as_matrix(Y)),
                #stripped=True,
                n=self.params['joint_components'],
                nx=self.params['x_ortho_components'],
                ny=self.params['y_ortho_components'],
                max_iterations=self.params['max_iterations']
            )
        
        self._logs = silenced_lines
        self.model = named_vector_to_dict(self.r_model)
        
        for key, value in self.attributes_map.items():
            setattr(self, key, self.model[value])

        self.E = X - self.predict_x(Y)
        self.F = Y - self.predict_y(X)
        self.fitted = True
        return self
    
    @property
    def x_coefficients(self) -> DataFrame:
        return maybe_add_columns(self._x_coefficients, self.Y)            
    
    @property
    def y_coefficients(self) -> DataFrame:
        return maybe_add_columns(self._y_coefficients, self.X)

    @property
    def _x_coefficients(self):
        return self.W @ self.B_T @ self.C.T
    
    @property
    def _y_coefficients(self):
        return self.C @ self.B_U @ self.W.T
    
    @property
    def coef_(self):
        return {
            'x': self.x_coefficients,
            'y': self.y_coefficients
        }
    
    def predict_x_in_R(self, Y):
        Y = r_matrix(as_matrix(Y))
        return predict_o2m(self.r_model, Y, XorY='Y')
    
    def predict_y_in_R(self, X):
        X = r_matrix(as_matrix(X))
        return predict_o2m(self.r_model, X, XorY='X')

    def predict_y(self, X):
        X = as_matrix(X)

        # see predict_y for explanations
        if self.B_T.shape == (1, 1):
            return enforce_2d(self.predict_y_in_R(X))

        X_corrected = X - X @ self.W_Yosc @ self.W_Yosc.T
        return X_corrected @ self.W @ self.B_T @ self.C.T

    def predict_x(self, Y):
        Y = as_matrix(Y)

        # when there is only one LV, compute in R as:
        # - in Python it takes two A4 pages to write scalar-compatible multiplication (np.dot etc)
        # - it is not as expensive to transfer "small" matrices and matrices for single LV are the smallest possilbe
        if self.B_U.shape == (1, 1):
            return enforce_2d(self.predict_x_in_R(Y))

        Y_corrected = Y - Y @ self.C_Xosc @ self.C_Xosc.T
        return Y_corrected @ self.C @ self.B_U @ self.W.T

    def predict(self, X, Y):
        return self.predict_x(Y), self.predict_y(X)
    
    def get_params(self, deep=False):
        return self.params.copy()

    def set_params(self, **params):
        self.params.update(params)
        return self

    @property
    def joint_X(self):
        return self.T @ self.W.T

    @property
    def joint_Y(self):
        return self.U @ self.C.T
    
    def calculate_Y_scores_U(self, Y):
        return (Y - Y @ self.C_Xosc @ self.C_Xosc.T) @ self.C
    
    def calculate_X_scores_T(self, X):
        return (X - X @ self.W_Yosc @ self.W_Yosc.T) @ self.W
    
    def calculate_Y_scores_U_without_osc(self, Y):
        return Y @ self.C
    
    def calculate_X_scores_T_without_osc(self, X):
        return X @ self.W
    
    @property
    def r_metrics(self):
        return Series({
            k: v[0] for k, v in self.model.items()
            if len(v) == 1 and isinstance(v[0], float)
        }).to_frame().T

    def __repr__(self):
        n = self.params['joint_components']
        nx = self.params['x_ortho_components']
        ny = self.params['y_ortho_components']
        return (
            f'<{"fitted " if self.fitted else ""}O2PLS model with: '
            f'{n} joint components, '
            f'{nx} orthogonal X components, '
            f'and {ny} orthogonal Y components>'
        )


def add_grid_metadata(grid_results):
    grid_results['statistic'] = grid_results['scoring'].apply(
        lambda name: O2PLSStatistic.statistics[name]
    )
    grid_results['scoring_label'] = grid_results['statistic'].apply(
        lambda stat: f'${stat.symbol}$: {stat.name}'
    )
    grid_results['model'] = grid_results['statistic'].apply(lambda stat: stat.model)
    grid_results['symbol'] = grid_results['statistic'].apply(lambda stat: f'${stat.symbol}$')
    grid_results = grid_results.drop(columns=['statistic'])
    grid_results['joint_mean'] = grid_results.groupby(
        ['param_joint_components', 'scoring']
    ).value.transform('mean')

    by_components = grid_results.groupby(
        ['param_joint_components', 'param_x_ortho_components', 'param_y_ortho_components', 'scoring']
    )
    grid_results['mean'] = by_components.value.transform('mean')
    grid_results['std'] = by_components.value.transform('std')

    grid_results['ortho_gain_percent'] = (
        100 *
        (grid_results['mean'] - grid_results['joint_mean'])
        /
        abs(grid_results['joint_mean'])
    )
    return grid_results


def summary_table_for(grid_results, model, joint_components, x_ortho_components, y_ortho_components):
    grid_summary = grid_results.drop_duplicates([
        'param_joint_components', 'param_x_ortho_components', 'param_y_ortho_components',
        'scoring'
    ])

    if model == 'Y':
        orthogonal_query = """
        param_x_ortho_components == @x_ortho_components
        and
        param_y_ortho_components <= @y_ortho_components
        """

    if model == 'X':
        orthogonal_query = """
        param_x_ortho_components <= @x_ortho_components
        and
        param_y_ortho_components == @y_ortho_components
        """

    df = DataFrame(grid_summary.query(f"""(
    param_joint_components <= @joint_components
    and
    {orthogonal_query}
    and
    model == @model
    )""".replace('\n', '')))
    m = model.lower()

    df['value_and_sd'] = (
            df['value'].apply(lambda x: f'{x:.2f}').astype(str)
            + ' ± '
            + df['std'].apply(lambda x: f'{x:.2f}').astype(str)
    )

    df['component'] = (
            'n='
            + df['param_joint_components'].astype(int).astype(str)
            + f', {m}_orth='
            + df[f'param_{m}_ortho_components'].astype(int).astype(str)
    )
    constant = 'param_x_ortho_components' if model == 'Y' else 'param_y_ortho_components'

    df = df.drop(columns=[constant])

    df = df.pivot('component', 'symbol', 'value_and_sd')
    df.columns.name = f'${model}$ model'
    return df
