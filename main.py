import argparse
import logging
from math import degrees, radians
from datetime import datetime as dt
import os
from typing import Dict
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

import orekit
from orekit.pyhelpers import (
    setup_orekit_curdir,
    absolutedate_to_datetime,
    datetime_to_absolutedate,
)

from org.hipparchus.geometry.euclidean.threed import Vector3D
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

from src.ephemeris import save_final_states, save_shooter_ephemeris
from src.plots import plot_separation_as_function_of_dv
from src.utils import get_true_anomaly_from_equinoctial, wrap_to_2_pi


minStep = 0.001
maxstep = 1000.0
initStep = 60.0
positionTolerance = 1.0
step = 60

APSIS_CATEGORIZATION_TRUE_ANOMALY_TOLERANCE_DEG = float(
    os.getenv("APSIS_CATEGORIZATION_TRUE_ANOMALY_TOLERANCE_DEG", 10.0)
)


logging.basicConfig(level=logging.DEBUG)
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

    return propagator_num, logger


def setup_propagation_context(initialDate: AbsoluteDate, initialOrbit: KeplerianOrbit):
    propagation_time = initialOrbit.getKeplerianPeriod()
    samples = int(propagation_time / step)
    linspace = np.linspace(
        absolutedate_to_datetime(initialDate),
        absolutedate_to_datetime(initialDate.shiftedBy(float(propagation_time))),
        samples + 1,
    )
    return linspace


def propagate_to_chosen_apsis(
    propagator: NumericalPropagator,
    initialDate: AbsoluteDate,
    propagate_by_duration: float,
    target_apsis: str,
):
    target_true_anomaly = 0.0 if target_apsis == "periapsis" else 180.0
    _date = initialDate
    while True:
        state = propagator.propagate(_date.shiftedBy(propagate_by_duration))

        true_anomaly_from_equinoctial = degrees(
            wrap_to_2_pi(get_true_anomaly_from_equinoctial(state))
        )
        module_logger.debug(f"Propagated to {true_anomaly_from_equinoctial=}")
        if (
            abs(true_anomaly_from_equinoctial - target_true_anomaly)
            < APSIS_CATEGORIZATION_TRUE_ANOMALY_TOLERANCE_DEG
        ):
            module_logger.info(
                f"Reached {target_apsis=} at {state.getDate()}: state = {state}"
            )
            break
        else:
            module_logger.debug(
                f"Propagation reached an ApsideEvent at {state.getDate()}: state = {state}. Updating the date and continuing..."
            )
            _date = state.getDate()

    return state


def create_prograde_impulse(state, dv):

    impulse = SmallManeuverAnalyticalModel(
        state, Vector3D(float(dv), 0.0, 0.0), 340.0
    )  # A completely made up Isp, but reasonable; besides, I am not calculating mass decrements

    return impulse


def shooter(
    initial_state: SpacecraftState,
    propagator: NumericalPropagator,
    dv_initial_guess: float,
    shoot_to_apsis: str,
    perturbation=0.1,
    number_of_parallel_trajectories=10,
):

    module_logger.warning("Starting shooting procedure...")
    propagated_states_map = {}
    ephem_generators_map = {}
    dv_distribution = norm(loc=dv_initial_guess, scale=perturbation)
    dvs_array = dv_distribution.rvs(size=number_of_parallel_trajectories)
    module_logger.info(
        f"Starting with an initial guess for the DV of {dv_initial_guess}, "
        f"and performing {number_of_parallel_trajectories} perturbatios with STD {perturbation}"
    )
    for cnt, dv in enumerate(dvs_array):

        module_logger.info(f"Performing perturbation {cnt+1} of {len(dvs_array)}")
        impulse = create_prograde_impulse(initial_state, dv)
        perturbed_state = impulse.apply(initial_state.shiftedBy(1.0))
        generator = propagator.getEphemerisGenerator()
        module_logger.debug(
            f"Resetting propagator at perturbed_state date of {perturbed_state.getDate().toString()}"
        )
        propagator.resetInitialState(perturbed_state)
        propagated_states_map.update(
            {
                dv: propagate_to_chosen_apsis(
                    propagator,
                    perturbed_state.getDate(),
                    perturbed_state.getKeplerianPeriod(),
                    shoot_to_apsis,
                )
            }
        )
        ephem_generators_map.update({dv: generator.getGeneratedEphemeris()})

    module_logger.info("Finished shooting")

    return propagated_states_map, ephem_generators_map


def calc_separation(
    primary_state: SpacecraftState, secondary_state: SpacecraftState
) -> float:

    # First make sure the two states are indeed at the same time, if not shift the secondary.
    # Orekit can shift the secondary by small margins directly via Keplerian + quadratic effects, which is okay over small steps
    _primary_date = absolutedate_to_datetime(primary_state.getDate())
    _secondary_date = absolutedate_to_datetime(secondary_state.getDate())
    # if _primary_date != _secondary_date:
    #     module_logger.debug(
    #         f"Processing RIC. Primary state epoch: {_primary_date}; Secondary state epoch: {_secondary_date}"
    #     )
    #     secondary_state = secondary_state.shiftedBy(
    #         (_primary_date - _secondary_date).total_seconds()
    #     )

    x_distance = (
        primary_state.getPosition().getX() - secondary_state.getPosition().getX()
    )
    y_distance = (
        primary_state.getPosition().getY() - secondary_state.getPosition().getY()
    )
    z_distance = (
        primary_state.getPosition().getZ() - secondary_state.getPosition().getZ()
    )

    return Vector3D(x_distance, y_distance, z_distance).getNorm()


def process_rics_at_final_epoch(
    reference_state: SpacecraftState,
    parallel_trajectories_map: Dict[float, SpacecraftState],
):

    dvs_ric_map = dict.fromkeys(parallel_trajectories_map.keys())

    for dv in dvs_ric_map.keys():

        dvs_ric_map[dv] = calc_separation(
            reference_state, parallel_trajectories_map[dv]
        )

    return dvs_ric_map


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

        propagator_num, logger = setup_propagation_objects(
            initialOrbit=initialOrbit, initialState=initialState
        )

        state_at_separation = propagate_to_chosen_apsis(
            propagator_num,
            initialDate,
            initialOrbit.getKeplerianPeriod(),
            impulse_at_which_apsis,
        )
        module_logger.info(
            f"Simulating separation at epoch {absolutedate_to_datetime(state_at_separation.getDate()).isoformat()}"
        )
        state_after_one_rev = propagate_to_chosen_apsis(
            propagator_num,
            state_at_separation.getDate(),
            state_at_separation.getKeplerianPeriod(),
            impulse_at_which_apsis,
        )

        dv_states_map, ephem_generators_map = shooter(
            state_at_separation,
            propagator_num,
            1.0,
            impulse_at_which_apsis,
            perturbation=0.1,
        )

        dv_separations_map = process_rics_at_final_epoch(
            state_after_one_rev, dv_states_map
        )

        save_shooter_ephemeris(
            ephem_generators_map.values(),
            "shooting",
        )
        plot_separation_as_function_of_dv(impulse_at_which_apsis, dv_separations_map)

    except Exception as global_exc:
        print(f"Global unhandled exception caught: {global_exc}")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
