import math
import numpy as np

from org.orekit.propagation import SpacecraftState


def wrap_to_2_pi(angle: float) -> float:
    return ((angle % (2 * math.pi)) + 2 * math.pi) % (2 * math.pi)


def get_true_anomaly_from_equinoctial(state: SpacecraftState) -> float:

    ex = state.getEquinoctialEx()
    ey = state.getEquinoctialEy()
    hx = state.getHx()
    hy = state.getHy()

    raan = math.atan2(hy, hx)
    omega = math.atan2(ey, ex) - raan
    ta = state.getLv() - raan - omega

    return ta
