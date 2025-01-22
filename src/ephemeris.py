from datetime import datetime as dt
import inspect
import os
from pathlib import Path
from typing import List

import numpy as np

from orekit.pyhelpers import absolutedate_to_datetime, datetime_to_absolutedate

from org.orekit.propagation import SpacecraftState, BoundedPropagator
from org.orekit.utils import TimeStampedPVCoordinates


EPHEMERIS_RECORDS_TIMESTEP = float(os.getenv("EPHEMERIS_RECORDS_TIMESTEP", 60))


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

    return inspect.cleandoc(header) + "\n"


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


def create_oem(object_name: str, initial_date: dt, final_date: dt, states: List) -> str:

    header = create_ephemeris_header(object_name, initial_date, final_date)
    data = create_oem_data(states)

    return header + data


def save_oem_to_file(prefix: Path | None, object_name: str, oem: str):

    filepath = Path(object_name).with_suffix(".oem")
    if prefix:
        filepath = prefix.joinpath(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w+") as fp:
        fp.write(oem)


def save_final_states(states_to_save: list[SpacecraftState]):

    with open("final_states.txt", "w+") as fp:
        for state in states_to_save:
            fp.write(
                f"{state.getDate().toString()} {state.getPVCoordinates().toString()}\n"
            )


def unravel_bounded_propagator(
    bounded_propagator: BoundedPropagator,
) -> List[TimeStampedPVCoordinates]:

    min_date = absolutedate_to_datetime(bounded_propagator.getMinDate())
    max_date = absolutedate_to_datetime(bounded_propagator.getMaxDate())

    step = EPHEMERIS_RECORDS_TIMESTEP
    propagation_time = (max_date - min_date).total_seconds()
    samples = int(propagation_time / step)
    linspace = np.linspace(
        min_date,
        max_date,
        samples + 1,
    )

    states = []
    for date_ in linspace:
        state = bounded_propagator.propagate(datetime_to_absolutedate(date_))
        states.append(state.getPVCoordinates())

    return states


def save_shooter_ephemeris(
    ephemeris_array_of_arrays: List[List[TimeStampedPVCoordinates]],
    prefix: str | None = None,
):

    for idx, timestamped_pv_coordinates_array in enumerate(ephemeris_array_of_arrays):

        min_date = absolutedate_to_datetime(
            timestamped_pv_coordinates_array[0].getDate()
        )
        max_date = absolutedate_to_datetime(
            timestamped_pv_coordinates_array[-1].getDate()
        )

        if idx == 0:
            object_name = "ReferenceSat"
        else:
            object_name = f"SeparatedSat{idx}"

        oem = create_oem(
            object_name,
            min_date,
            max_date,
            timestamped_pv_coordinates_array,
        )
        if prefix:
            prefix = Path(prefix)
        save_oem_to_file(prefix, object_name, oem)
