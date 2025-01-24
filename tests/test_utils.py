import math

import pytest

import orekit

from main import *
from src.utils import *


def test_wrap2_pi():

    angle = 540.0
    result = wrap_to_2_pi(math.radians(angle))
    assert math.pi == pytest.approx(result)


def test_true_anomaly_calculation():

    epoch = datetime_to_absolutedate(dt.fromisoformat("2025-01-01T12:00:00"))
    frame = FramesFactory.getEME2000()
    basic_orbit = KeplerianOrbit(
        7555.0,
        0.00115,
        math.radians(28.5),
        math.radians(0.0),
        math.radians(0.0),
        math.radians(135.0),
        PositionAngleType.TRUE,
        frame,
        epoch,
        Constants.WGS84_EARTH_MU,
    )
    state1 = SpacecraftState(basic_orbit)

    assert math.radians(135.0) == pytest.approx(
        get_true_anomaly_from_equinoctial(state1)
    )
