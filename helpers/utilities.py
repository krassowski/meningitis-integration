from copy import deepcopy, copy
from itertools import chain
from functools import reduce, partial
from statistics import mean

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# from helpers.r import *

import pandas as pd

from jupyter_helpers.table import display_table
from jupyter_helpers.source import embed_source_styling, show_source
from jupyter_helpers.namespace import Namespace, NeatNamespace, HorizontalNamespace
import jupyter_helpers.rpy2_autocompletion

from IPython.display import HTML
from pandas import read_table, read_csv, read_excel, concat, Series, DataFrame
import numpy as np
from tqdm.auto import tqdm


pd.options.display.max_rows = 10
pd.options.display.max_columns = 10

# enable stable sort by default to minimize diffs on notebooks
DataFrame.sort_values.__defaults__ = tuple(
    'mergesort' if default == 'quicksort' else default
    for default in DataFrame.sort_values.__defaults__
)


show_table = display_table

full_table = partial(display_table, n_rows=None)


def keys(obj):
    return list(obj.keys())


# embed_source_styling()

T = True
F = False
