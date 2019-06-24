from functools import partial
from . import _venn


def venn_n(data, n, names=None, fill='number'):
    labels = _venn.get_labels([
        subset
        for subset in data.values()
    ], fill=[fill])
    return getattr(_venn, f'venn{n}')(labels, names=names if names else list(data.keys()))


venn2 = partial(venn_n, n=2)
venn3 = partial(venn_n, n=3)
venn4 = partial(venn_n, n=4)

