import pathlib
import numpy as np


def resolve_relative_path(file: str, path: str, parent_levels=0) -> pathlib.Path:
    """
    returns absolute path to a file given its relative :param path from :param file
    :param file: path from root to a file
    :param path: relative path to file2 from the last common point in the two paths
    :param parent_levels: number of levels in the path to go up to get to the last common point
    :return: path from root to file2
    """
    return pathlib.Path(file).parents[parent_levels].joinpath(path)


def calculate_distance(a: tuple, b: tuple) -> float:
    x1 = a[0]
    x2 = b[0]
    y1 = a[1]
    y2 = b[1]
    d = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    return d
