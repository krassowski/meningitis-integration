from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


def p_adjust(x, *args, **kwargs):
    if not isinstance(x, list):
        x = list(x)
    return r['p.adjust'](x, *args, **kwargs)


def r_function(__name__, *args, **kwargs):
    with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):
        return r[__name__](*args, **kwargs)

