import logging

from typing import Callable, Dict, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

matplotlib_logger = logging.getLogger("matplotlib")
matplotlib_logger.setLevel(logging.INFO)
PIL_logger = logging.getLogger("PIL")
PIL_logger.setLevel(logging.INFO)


def scatter_separation_as_function_of_dv(
    ax: Axes, apsis_of_separation: str, dv_separations_map: Dict[float, float]
):

    ax.scatter(list(dv_separations_map.keys()), list(dv_separations_map.values()))
    ax.set_title(
        f"Distance due to impulse after 1 rev - {apsis_of_separation} to {apsis_of_separation}"
    )
    ax.set_xlabel("Impulse [m/s]")
    ax.set_ylabel("Distance [m]")

    return ax


def plot_fitted_curve(ax: Axes, xdata: np.array, func, popt):

    ax.plot(
        xdata,
        func(xdata, *popt),
        "g--",
        label="fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f" % tuple(popt),
    )

    return ax


def plot_root(ax: Axes, root: float, image: float):

    ax.scatter(
        [root],
        [image],
        c=["r"],
        marker="*",
    )
    return ax
