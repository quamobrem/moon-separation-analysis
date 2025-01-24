import pytest

import orekit

from main import *


@pytest.fixture(autouse=True, scope="session")
def init_orekit():
    setup_orekit()
