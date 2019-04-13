from rpy2.robjects import r


def p_adjust(x, *args, **kwargs):
    if not isinstance(x, list):
        x = list(x)
    return r['p.adjust'](x, *args, **kwargs)
