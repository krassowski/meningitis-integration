from typing import Dict

from pandas import DataFrame, Series

from machine_learning.preprocessing import OutliersFilter
from machine_learning.data_classes import BlocksWrapper


def subset(omic, subset, outliers=None):
    df = omic[omic.columns.intersection(subset)].T
    if outliers is not None:
        of = OutliersFilter(outlier_patients=outliers, verbose=True)
        df = of.fit_transform(df)
    return df


def add_supervision_block(
    blocks: Dict[str, DataFrame], conditions_vector: Series,
    name='y'
):
    return BlocksWrapper(blocks).add_supervision_block(
        conditions_vector=conditions_vector,
        name=name
    ).blocks