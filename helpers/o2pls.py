from contextlib import redirect_stdout
import re
from functools import partial
from typing import NamedTuple, Optional
from io import StringIO

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# from helpers.pls import custom_scorer
from helpers.r import r_function


class Statistic(NamedTuple):
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


class O2PLSStatistic(Statistic):
    statistics = {}
    
    def calc(self, estimator, X, Y):
        return self.func(estimator, as_matrix(X), as_matrix(Y))

    
def o2pls_statistic(name, *args, **kwargs):
    
    def decorator(func):
        statistic = O2PLSStatistic(name, *args, func=func, **kwargs)
        O2PLSStatistic.statistics[name] = statistic
        return statistic
    
    return decorator


def modeled_variation(residuals, matrix):
    return 1 - sum(sum(pow(as_matrix(residuals), 2))) / sum(sum(pow(as_matrix(matrix), 2)))


@o2pls_statistic('Modeled variation of $X$', 'R^2X', equation='1 - \sum(E)^2 / \sum X^2')
def modeled_variation_of_X(estimator, x, y):
    return modeled_variation(residuals=estimator.E, matrix=estimator.X)


@o2pls_statistic('Modeled variation of $Y$', 'R^2Y', equation='1 - \sum(F)^2 / \sum Y^2')
def modeled_variation_of_Y(estimator, x, y):
    return modeled_variation(residuals=estimator.F, matrix=estimator.Y)


def r2(true, pred):
    return 1 - sum(sum(as_matrix(pred))) / sum(sum(as_matrix(true)))


@o2pls_statistic('Predicted variation of $X$', 'R^2X_{hat}')
def predicted_variation_of_X(estimator, x, y):
    return r2(
        true=estimator.X,
        pred=estimator.predict_x(estimator.Y)
    )


@o2pls_statistic('Predicted variation of $Y$', 'R^2Y_{hat}')
def predicted_variation_of_Y(estimator, x, y):
    return r2(
        true=estimator.Y,
        pred=estimator.predict_y(estimator.X)
    )


@o2pls_statistic('Cross-validated predictions of $X$', 'Q^2X_{hat}')
def cross_validated_predictions_of_X(estimator, X, Y):
    return r2(
        true=X,
        pred=estimator.predict_x(Y)
    )


@o2pls_statistic('Cross-validated predictions of $Y$', 'Q^2Y_{hat}')
def cross_validated_predictions_of_Y(estimator, X, Y):
    return r2(
        true=Y,
        pred=estimator.predict_y(X)
    )


@o2pls_statistic('Modeled structured noise variation of $X$', 'R^2X_{YO}')
def modeled_structured_noise_variation_of_X(estimator, x, y):
    return modeled_variation(
        residuals=estimator.T_Yosc @ estimator.P_Yosc.T,
        matrix=estimator.Y
    )


@o2pls_statistic('Modeled structured noise variation of $Y$', 'R^2Y_{XO}')
def modeled_structured_noise_variation_of_Y(estimator, x, y):
    return modeled_variation(
        residuals=estimator.U_Xosc @ estimator.P_Xosc.T,
        matrix=estimator.Y
    )


@o2pls_statistic('Modeled joint X-Y covariation of $X$', 'R^2X_{corr}')
def modeled_joint_XY_covariation_of_X(estimator, x, y):
    return 1 - modeled_variation(residuals=estimator.T @ estimator.W.T, matrix=estimator.X)


@o2pls_statistic('Modeled joint X-Y covariation of $Y$', 'R^2Y_{corr}')
def modeled_joint_XY_covariation_of_Y(estimator, x, y):
    return 1 - modeled_variation(residuals=estimator.U @ estimator.C.T, matrix=estimator.Y)


# meta-statistics

@o2pls_statistic(
    'Average of cross-validated predictions of $X$ and $Y$',
    '(Q^2X_{hat} + Q^2Y_{hat}) / 2'
)
def average_cv_predictions(estimator, X, Y):
    return (
        cross_validated_predictions_of_X.calc(estimator, X, Y)
        + cross_validated_predictions_of_Y.calc(estimator, X, Y)
    ) / 2


o2m2 = partial(r_function, 'o2m2')
predict_o2m = partial(r_function, 'predict') 


def named_vector_to_dict(x):
    assert len(set(x.names)) == len(x.names)
    with localconverter(ro.default_converter + pandas2ri.converter):
        return {
            k: ro.conversion.rpy2py(v)
            for k, v in x.items()
        }

    
class O2PLS:
    
    attributes_map = {
        'T': 'Tt',   # joint X-scores
        'U': 'U',    # joint Y-scores
        # not implemented for o2m2, see https://github.com/selbouhaddani/OmicsPLS/issues/11
        # 'E': 'E',    # joint X-residuals
        # 'F': 'Ff',    # joint Y-residals
        'W': 'W.',   # joint X-loadings
        'C': 'C.',   # joint Y-loadings
        'P_Yosc': 'P_Yosc.',   # orthogonal X-loadings
        'P_Xosc': 'P_Xosc.',   # orthogonal Y-loadingss
        'T_Yosc': 'T_Yosc.',   # orthogonal X-scores
        'U_Xosc': 'U_Xosc.',   # orthogonal Y-scores
        'B_U': 'B_U',  # weights/coefficients for X predction
        'B_T': 'B_T.',  # weights/coefficients for Y predction
        'C_Xosc': 'C_Xosc',
        'W_Yosc': 'W_Yosc'
    }
    
    def __init__(self, max_iterations=100, **params):
        """
        Args:
            max_iterations
        """
        self.params = {
            'max_iterations': max_iterations,
            **params
        }

    def fit(self, X, Y):
        self.X = X
        self.Y = Y
        
        f = StringIO()
        
        with redirect_stdout(f):
            self.model_r = o2m2(
                X.values, Y.values,
                stripped=True,
                n=self.params['joint_components'],
                nx=self.params['x_ortho_components'],
                ny=self.params['y_ortho_components'],
                max_iterations=self.params['max_iterations']
            )

        self._logs = f.getvalue()
        self._check_logs()
        self.model = named_vector_to_dict(self.model_r)
        
        for key, value in self.attributes_map.items():
            setattr(self, key, self.model[value])

        self.E = X - self.predict_x(Y)
        self.F = Y - self.predict_y(X)
        return self
    
    def _check_logs(self):
        for line in self._logs.split('\n'):
            if not line:
                continue
            assert re.match(r'Power Method \(comp \d+\) stopped after \d+ iterations\.', line)

    def predict_x_in_R(self, Y):
        Y = as_matrix(Y)
        return predict_o2m(self.model_r, Y, XorY='X')
    
    def predict_y_in_R(self, X):
        X = as_matrix(X)
        return predict_o2m(self.model_r, X, XorY='Y')
            
    def predict_y(self, X):
        X = as_matrix(X)
        X_corrected = X - X @ self.W_Yosc @ self.W_Yosc.T
        return X_corrected @ self.W @ self.B_T @ self.C.T

    def predict_x(self, Y):
        Y = as_matrix(Y)
        Y_corrected = Y - Y @ self.C_Xosc @ self.C_Xosc.T
        return Y_corrected @ self.C @ self.B_U @ self.W.T

    def predict(self, X, Y):
        return self.predict_x(Y), self.predict_y(X)
    
    def get_params(self, deep=False):
        return self.params.copy()

    def set_params(self, **params):
        self.params.update(params)
        return self
