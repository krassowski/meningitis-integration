from copy import deepcopy, copy
from itertools import chain
from functools import reduce, partial
from os import chdir, getcwd
from pathlib import Path
from statistics import mean

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd

from jupyter_helpers.table import display_table
from jupyter_helpers.source import embed_source_styling, show_source
from jupyter_helpers.namespace import Namespace, NeatNamespace, HorizontalNamespace
from jupyter_helpers.selective_import import skip_on_import
import jupyter_helpers.rpy2_autocompletion

from IPython.display import HTML
from pandas import read_table, read_csv, read_excel, concat, Series, DataFrame
import numpy as np
from tqdm.auto import tqdm
from .io import create_paths, save_outputs, load_inputs

pd.options.display.max_rows = 10
pd.options.display.max_columns = 10

# enable stable sort by default to minimize diffs on notebooks
DataFrame.sort_values.__defaults__ = tuple(
    'mergesort' if default == 'quicksort' else default
    for default in DataFrame.sort_values.__defaults__
)

from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
rpy2_logger.addFilter(lambda record: 'notch went outside hinges' not in record.msg)

import rpy2.rinterface

def print_to_notebook(x):
    print(x, end='')

# rpy2.rinterface_lib.callbacks.consolewrite_print = print_to_notebook
rpy2.rinterface_lib.callbacks.consolewrite_warnerror = print_to_notebook


show_table = display_table
full_table = partial(display_table, n_rows=None)


def keys(obj):
    return list(obj.keys())

# embed_source_styling()


T = True
F = False


local_dir = getcwd()
top_level = Path(__file__).parent.parent

# always use the same, absolute paths - which makes
# moving the notebooks around easier in the future
chdir(top_level)


class Dummy:
    def __getattr__(self, key):
        return


dummy = Dummy()


def get_or_dummy(callback, *args, **kwargs):
    try:
        return callback(*args, **kwargs)
    except ValueError:
        return dummy

    
load_inputs = partial(
    load_inputs,
    main_loader=partial(read_csv, index_col=0)
)
