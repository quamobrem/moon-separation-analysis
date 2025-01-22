import logging

from typing import Dict

import matplotlib.pyplot as plt

matplotlib_logger = logging.getLogger("matplotlib")
matplotlib_logger.setLevel(logging.INFO)
PIL_logger = logging.getLogger("PIL")
PIL_logger.setLevel(logging.INFO)


def plot_separation_as_function_of_dv(
    apsis_of_separation: str, dv_separations_map: Dict[float, float]
):

    fig, ax = plt.subplots()

    ax.scatter(list(dv_separations_map.keys()), list(dv_separations_map.values()))
    ax.set_title(
        f"Distance due to impulse after 1 rev - {apsis_of_separation} to {apsis_of_separation}"
    )
    ax.set_xlabel("Impulse [m/s]")
    ax.set_ylabel("Distance [m]")

    plt.show()
