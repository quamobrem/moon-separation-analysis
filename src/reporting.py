import logging

from typing import Callable, Dict, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from org.orekit.propagation import SpacecraftState


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


def report_results(
    solution: float,
    desired_separation: float,
    apsis_of_separation: str,
    state_at_separation: SpacecraftState,
    final_state: SpacecraftState,
):
    print(
        f"Final Solution Reporting\n"
        f"========================\n"
        f"A prograde impulsive dV of {solution} m/sec, applied at {apsis_of_separation}, achieves the desired separation of {desired_separation} m "
        f"after a one orbital revolution.\n"
        f"========================\n"
        f"The states at separation and the final state follow, expressed in a ICRF-oriented, Moon Centered Frame:\n"
        f"The state at separation is {state_at_separation.getPVCoordinates()}\n"
        f"The final state is {final_state.getPVCoordinates()}"
    )
