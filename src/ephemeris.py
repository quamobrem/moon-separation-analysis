from datetime import datetime as dt
import inspect


def create_ephemeris_header(initial_date: dt, final_date: dt):
    start_time = initial_date.strftime("%Y-%m-%dT%H:%M:%S.%f")
    end_time = final_date.strftime("%Y-%m-%dT%H:%M:%S.%f")
    header = f"""
    CCSDS_OEM_VERS = 1.0

    CREATION_DATE  = {dt.now().isoformat()}
    ORIGINATOR     = Myself

    META_START
    OBJECT_NAME          = Sat1
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
