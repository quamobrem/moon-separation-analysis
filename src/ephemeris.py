from datetime import datetime as dt
import inspect
from typing import List

from orekit.pyhelpers import absolutedate_to_datetime


def create_ephemeris_header(object_name: str, initial_date: dt, final_date: dt) -> str:
    start_time = initial_date.strftime("%Y-%m-%dT%H:%M:%S.%f")
    end_time = final_date.strftime("%Y-%m-%dT%H:%M:%S.%f")
    header = f"""
    CCSDS_OEM_VERS = 1.0

    CREATION_DATE  = {dt.now().isoformat()}
    ORIGINATOR     = Myself

    META_START
    OBJECT_NAME          = {object_name}
    OBJECT_ID            = 2000-000A
    CENTER_NAME          = Moon
    REF_FRAME            = ICRF
    TIME_SYSTEM          = UTC
    START_TIME           = {start_time}
    USEABLE_START_TIME   = {start_time}
    USEABLE_STOP_TIME    = {end_time}
    STOP_TIME            = {end_time}
    INTERPOLATION        = Lagrange
    INTERPOLATION_DEGREE = 5
    META_STOP

    """

    return inspect.cleandoc(header)


def create_oem_data(states: List) -> str:
    output = ""
    line = ""
    for state in states:
        line = "{date_} {x_} {y_} {z_} {v_x} {v_y} {v_z}\n".format(
            date_=absolutedate_to_datetime(state.getDate()).strftime(
                "%Y-%m-%dT%H:%M:%S.%f"
            ),
            x_=state.getPosition().getX() / 1000,
            y_=state.getPosition().getY() / 1000,
            z_=state.getPosition().getZ() / 1000,
            v_x=state.getVelocity().getX() / 1000,
            v_y=state.getVelocity().getY() / 1000,
            v_z=state.getVelocity().getZ() / 1000,
        )
        output += line

    return output


def createOem(object_name: str, initial_date: dt, final_date: dt, states: List) -> str:

    header = create_ephemeris_header(object_name, initial_date, final_date)
    data = create_oem_data(states)

    return header + data


def save_oem_to_file(object_name: str, oem: str):

    with open(object_name + ".oem", "w+") as fp:
        fp.write(oem)
