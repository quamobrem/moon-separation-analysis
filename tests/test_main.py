import pytest

import orekit

from main import *


def test_root_finder():
    """
    Am I using the root finder correctly?
    Testing it against a simple linear function
    """

    def func(x, m, q):
        return m * x + q

    root = root_finder(func, [5.0, 12.0], 55.0, 62.0)

    assert root == pytest.approx(10.0)


def test_calc_separation():
    """Tests the calculation of the RIC.
    The test states have been manufactured to be at an exact 10km offset.
    """

    epoch1 = datetime_to_absolutedate(dt.fromisoformat("2025-01-18T13:28:13.539"))
    frame = FramesFactory.getEME2000()
    from org.orekit.utils import PVCoordinates

    coordinates1 = PVCoordinates(
        Vector3D(6608.710713 * 1000, -853.075402 * 1000, -433.605184 * 1000),
        Vector3D(1.108565 * 1000, 6.717702 * 1000, 3.652032 * 1000),
    )
    orbit1 = KeplerianOrbit(coordinates1, frame, epoch1, Constants.WGS84_EARTH_MU)
    state1 = SpacecraftState(orbit1)

    epoch2 = datetime_to_absolutedate(dt.fromisoformat("2025-01-18T13:27:45.000"))
    coordinates2 = PVCoordinates(
        Vector3D(6571.712998 * 1000, -1052.939504 * 1000, -542.277051 * 1000),
        Vector3D(1.360760 * 1000, 6.681405 * 1000, 3.633399 * 1000),
    )
    orbit2 = KeplerianOrbit(coordinates2, frame, epoch2, Constants.WGS84_EARTH_MU)
    state2 = SpacecraftState(orbit2)

    separation = calc_separation(state1, state2)

    assert pytest.approx(separation, rel=10.0) == 10.0 * 1000
