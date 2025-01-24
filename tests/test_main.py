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
