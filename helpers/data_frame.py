from inflection import underscore

from typing import List, Union
from pandas import DataFrame, Series, Index


def to_nested_series(df: DataFrame) -> Series:
    return Series(dict(zip(
        df.columns,
        df.T.values.tolist()
    )))


def explode_rows_with_lists(df: DataFrame) -> DataFrame:
    data = []
    c = 0
    for i, row in enumerate(df.itertuples(index=False)):
        row_dict = row._asdict()
        first_value = next(iter(row_dict.values()))
        for j in range(len(first_value)):
            c += 1
            sub_row_j = {
                key: row_dict[key][j]
                for key in row_dict.keys()
            }
            data.append({'index': c, 'group': i, **sub_row_j})

    return DataFrame(data).set_index(['group', 'index'])


def extract_duplicates(data: DataFrame, duplicate_columns: List[str], index_columns: List[str], fill_na='NA') -> DataFrame:
    # without filling nulls we would get false negatives as nan != nan in Python
    if fill_na:
        data = data.fillna(fill_na)
    # pre-filter to operate only on duplicates
    data = data[duplicate_columns + index_columns]
    data = data[data.duplicated(keep=False, subset=duplicate_columns)]
    by_values = (
        data.reset_index()
        .groupby(duplicate_columns)
        [index_columns]
        .apply(to_nested_series)
    )
    if by_values.empty:
        return DataFrame()
    some_index_column = by_values[index_columns[0]]
    df = by_values[some_index_column.apply(len) > 1].reset_index(drop=True)
    return explode_rows_with_lists(df)


def set_duplicates_group(duplicates, group, other):
    duplicates[group] = duplicates.merge(other)[group].values
    duplicates = (
        duplicates
        .reset_index('group', drop=True)
        .reset_index()
        .set_index([group, 'index'])
    )
    return duplicates


def to_underscore(data: Union[Index, List[str], Series], limit_to=None):
    return [
        underscore(name)
        if not limit_to or name in limit_to else
        name
        for name in data
    ]


def to_lowercase(data: Union[Index, List[str], Series], limit_to=None):
    return [
        name.lower()
        if not limit_to or name in limit_to else
        name
        for name in data
    ]


def select_columns(data, match=None, exclude=None):
    assert (match or exclude) and not (match and exclude)
    regex = match or exclude
    selection_vector = data.columns.str.contains(regex)
    if exclude:
        selection_vector = ~selection_vector
    return data[data.columns[selection_vector]]


def juxtapose_tables(a, b, suffixes, keep_index=False, only_keep=None):
    if not keep_index:
        a = a.reset_index(drop=True)
        b = b.reset_index(drop=True)
    if only_keep:
        a = a[only_keep]
        b = b[only_keep]
    df = a.join(b, lsuffix='|' + suffixes[0], rsuffix='|' + suffixes[1])
    df.columns = df.columns.map(lambda x: tuple(reversed(x.split('|'))))
    return df
