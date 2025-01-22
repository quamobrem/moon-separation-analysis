import argparse
import logging
from math import degrees, radians
from datetime import datetime as dt
import numpy as np

import orekit
from orekit.pyhelpers import (
    setup_orekit_curdir,
    absolutedate_to_datetime,
    datetime_to_absolutedate,
)

from org.orekit.bodies import CelestialBodyFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants
from org.orekit.frames import FramesFactory, Transform, Frame

from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import SpacecraftState
from org.orekit.orbits import OrbitType
from org.orekit.forces.gravity import NewtonianAttraction, ThirdBodyAttraction
from org.orekit.propagation.events import ApsideDetector, EventsLogger
from org.orekit.propagation.events.handlers import StopOnEvent
from org.orekit.forces.maneuvers import SmallManeuverAnalyticalModel


from orekit import JArray_double


minStep = 0.001
maxstep = 1000.0
initStep = 60.0
positionTolerance = 1.0
step = 60


logging.basicConfig(level=logging.INFO)
module_logger = logging.getLogger(__name__)


def setup_orekit():

    try:
        vm = orekit.initVM()
        print("Java version:", vm.java_version)
        print("Orekit version:", orekit.VERSION)

        setup_orekit_curdir()
    except Exception as exc:
        print(
            f"Exception when setting up the Orekit data folder: {exc}! Cannot continue"
        )
        SystemExit(1)


def get_input_args():

    parser = argparse.ArgumentParser(
        prog="Given an initial condition for a Lunar Orbiter, calculate the minimum DV in m/s to achieve a given separation after 1 revolution"
    )

    parser.add_argument("--initial_epoch", default="2025-01-18T12:00:00")
    parser.add_argument("-apoalt", "--altitude_of_apoapsis", default=10_000)
    parser.add_argument("-perialt", "--altitude_of_periapsis", default=100)
    parser.add_argument("-i", "--inclination", default=28.5)
    parser.add_argument("-raan", "--right_ascension_of_ascending_node", default=0.0)
    parser.add_argument("-w", "--argument_of_periapsis", default=270.0)
    parser.add_argument("-ta", "--true_anomaly", default=0.0)
    parser.add_argument("--spacecraft_masses", default=100)
    parser.add_argument("--separation_to_achieve", default=10_000)
    parser.add_argument(
        "--impulse_at_which_apsis",
        choices=["apoapsis", "periapsis"],
        default="apoapsis",
    )

    args = parser.parse_args()
    return args


def create_moon_centered_ICRF_frame(
    initialDate: AbsoluteDate,
):

    moon = CelestialBodyFactory.getMoon()

    ## Inertial frame where the satellite is defined
    inertialMoonFrame = moon.getInertiallyOrientedFrame()

    ## ICRF Frame, for output
    ICRF = FramesFactory.getICRF()
    rotation_to_ICRF = inertialMoonFrame.getStaticTransformTo(
        ICRF, initialDate
    ).getRotation()
    transform_to_ICRF = Transform(initialDate, rotation_to_ICRF)
    moon_centered_ICRF = Frame(
        inertialMoonFrame, transform_to_ICRF, "moon_centered_ICRF", True
    )

    return moon_centered_ICRF


def create_initial_orbit_and_state(
    initialDate, frame, rp, ra, i, omega, raan, lv, satellite_mass
):

    a = (rp + ra + 2 * Constants.MOON_EQUATORIAL_RADIUS) / 2.0
    e = 1.0 - (rp + Constants.MOON_EQUATORIAL_RADIUS) / a

    ## Orbit construction as Keplerian
    initialOrbit = KeplerianOrbit(
        a,
        e,
        i,
        omega,
        raan,
        lv,
        PositionAngleType.TRUE,
        frame,
        initialDate,
        Constants.JPL_SSD_MOON_GM,
    )

    satellite_mass = 100.0  # The models need a spacecraft mass, unit kg.
    initialState = SpacecraftState(initialOrbit, satellite_mass)

    return initialOrbit, initialState


def setup_propagation_objects(initialOrbit, initialState):

    tolerances = NumericalPropagator.tolerances(
        positionTolerance, initialOrbit, initialOrbit.getType()
    )

    integrator = DormandPrince853Integrator(
        minStep,
        maxstep,
        JArray_double.cast_(
            tolerances[0]
        ),  # Double array of doubles needs to be casted in Python
        JArray_double.cast_(tolerances[1]),
    )
    integrator.setInitialStepSize(initStep)

    propagator_num = NumericalPropagator(integrator)
    propagator_num.setOrbitType(OrbitType.KEPLERIAN)
    propagator_num.setInitialState(initialState)

    lunar_gravity = NewtonianAttraction(Constants.JPL_SSD_MOON_GM)
    earth_gravity = ThirdBodyAttraction(CelestialBodyFactory.getEarth())
    sun_gravity = ThirdBodyAttraction(CelestialBodyFactory.getSun())

    propagator_num.addForceModel(lunar_gravity)
    propagator_num.addForceModel(earth_gravity)
    propagator_num.addForceModel(sun_gravity)

    apsis_detector = ApsideDetector(initialOrbit).withHandler(StopOnEvent())

    logger = EventsLogger()
    logged_detector = logger.monitorDetector(apsis_detector)

    propagator_num.addEventDetector(logged_detector)
    generator = propagator_num.getEphemerisGenerator()

    return propagator_num, generator, logger


def setup_propagation_context(initialDate: AbsoluteDate, initialOrbit: KeplerianOrbit):
    propagation_time = initialOrbit.getKeplerianPeriod()
    samples = int(propagation_time / step)
    linspace = np.linspace(
        absolutedate_to_datetime(initialDate),
        absolutedate_to_datetime(initialDate.shiftedBy(float(propagation_time))),
        samples + 1,
    )
    return linspace


def propagate_to_separation_apsis(
    propagator: NumericalPropagator,
    initialDate: AbsoluteDate,
    initialOrbit: KeplerianOrbit,
    apsis: str,
):

    state = propagator.propagate(
        initialDate.shiftedBy(initialOrbit.getKeplerianPeriod())
    )

    true_anomaly = degrees(
        round(
            KeplerianOrbit(
                OrbitType.KEPLERIAN.convertType(state.getOrbit())
            ).getTrueAnomaly()
        )
    )
    need_to_propagate_more = True
    if apsis == "periapsis":
        if true_anomaly == 0.0:
            module_logger.info("We are at periapsis, which is where we want to be")
            need_to_propagate_more = False
        else:
            module_logger.info(
                "We are at periapsis, but we want to start at apoapsis. Propagating until next stopping condition."
            )

    else:
        if true_anomaly == 180.0:
            module_logger.info("We are at apoapsis, which is where we want to be")
            need_to_propagate_more = False
        else:
            module_logger.info(
                "We are at periapsis, but we want to start at apoapsis. Propagating until next stopping condition."
            )

    if need_to_propagate_more:
        state = propagator.propagate(
            initialDate.shiftedBy(initialOrbit.getKeplerianPeriod())
        )

    module_logger.info(
        f"Starting state at separation is {state.getDate()}, {state.getOrbit()}"
    )
    return state


def main():
    try:
        setup_orekit()
        args = get_input_args()

        initial_datetime = dt.fromisoformat(args.initial_epoch)
        altitude_of_apoapsis = args.altitude_of_apoapsis * 1000
        altitude_of_periapsis = args.altitude_of_periapsis * 1000
        inclination = radians(args.inclination)
        right_ascension_of_ascending_node = radians(
            args.right_ascension_of_ascending_node
        )
        argument_of_periapsis = radians(args.argument_of_periapsis)
        true_anomaly = radians(args.true_anomaly)
        spacecraft_masses = args.spacecraft_masses
        separation_to_achieve = args.separation_to_achieve
        impulse_at_which_apsis = args.impulse_at_which_apsis

        initialDate = datetime_to_absolutedate(initial_datetime)

        moon_centered_ICRF_frame = create_moon_centered_ICRF_frame(initialDate)

        initialOrbit, initialState = create_initial_orbit_and_state(
            initialDate,
            moon_centered_ICRF_frame,
            altitude_of_periapsis,
            altitude_of_apoapsis,
            inclination,
            argument_of_periapsis,
            right_ascension_of_ascending_node,
            true_anomaly,
            spacecraft_masses,
        )

        propagator_num, generator, logger = setup_propagation_objects(
            initialOrbit=initialOrbit, initialState=initialState
        )

        state_at_separation = propagate_to_separation_apsis(
            propagator_num, initialDate, initialOrbit, impulse_at_which_apsis
        )

    except Exception as global_exc:
        print(f"Global unhandled exception caught: {global_exc}")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
