from scipy.stats import spearmanr


def results_correlation(l, r, n=None, p_value_column='adj.P.Val'):
    subset = [p for p in l.index if p in r.index]
    if n: subset = subset[:n]
    return spearmanr(
        l[p_value_column][subset].rank(),
        r[p_value_column][subset].rank()
    )