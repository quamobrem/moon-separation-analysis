import argparse
import logging
from math import degrees, radians
from datetime import datetime as dt
import os
from typing import Callable, Dict, List, Tuple
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit, newton
import matplotlib.pyplot as plt

import orekit
from orekit.pyhelpers import (
    setup_orekit_curdir,
    absolutedate_to_datetime,
    datetime_to_absolutedate,
)

from org.hipparchus.geometry.euclidean.threed import Vector3D, Rotation
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants
from org.orekit.frames import FramesFactory, Transform, Frame

from org.orekit.propagation import BoundedPropagator
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import SpacecraftState
from org.orekit.orbits import OrbitType
from org.orekit.forces.gravity import NewtonianAttraction, ThirdBodyAttraction
from org.orekit.propagation.events import ApsideDetector, EventsLogger
from org.orekit.propagation.events.handlers import StopOnEvent
from org.orekit.forces.maneuvers import SmallManeuverAnalyticalModel


from orekit import JArray_double

from src.ephemeris import (
    save_final_states,
    save_shooter_ephemeris,
    unravel_bounded_propagator,
)
from src.reporting import (
    plot_fitted_curve,
    plot_root,
    report_results,
    scatter_separation_as_function_of_dv,
)
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
    parser.add_argument("--dV_initial_guess", default=1.0)
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
    starting_state: SpacecraftState,
    propagate_by_duration: float,
    target_apsis: str,
) -> Tuple[SpacecraftState, BoundedPropagator]:
    target_true_anomaly = 0.0 if target_apsis == "periapsis" else 180.0
    _date = starting_state.getDate()
    generator = propagator.getEphemerisGenerator()
    ephemeris = []
    module_logger.debug(
        f"Resetting propagator at perturbed_state date of {starting_state.getDate().toString()}"
    )
    propagator.resetInitialState(starting_state)
    condition = True
    while condition:
        end_state = propagator.propagate(_date.shiftedBy(propagate_by_duration))

        true_anomaly_from_equinoctial = degrees(
            wrap_to_2_pi(get_true_anomaly_from_equinoctial(end_state))
        )
        module_logger.debug(f"Propagated to {true_anomaly_from_equinoctial=}")
        if (
            abs(true_anomaly_from_equinoctial - target_true_anomaly)
            < APSIS_CATEGORIZATION_TRUE_ANOMALY_TOLERANCE_DEG
        ):
            module_logger.info(
                f"Reached {target_apsis=} at {end_state.getDate()}: state = {end_state}"
            )
            condition = False
        else:
            module_logger.debug(
                f"Propagation reached an ApsideEvent at {end_state.getDate()}: state = {end_state}. Updating the date and continuing..."
            )
            _date = end_state.getDate()

        ephemeris.extend(unravel_bounded_propagator(generator.getGeneratedEphemeris()))

    return end_state, ephemeris


def create_prograde_impulse(state, dv):

    # Creating the impulse via SmallManeuverAnalyticalModel. Without specifying the Frame,
    # the dV vector is assumed in the spacecraft body frame, instead of in RTN. So have to set the X to -1
    # A completely made up Isp, but reasonable; besides, I am not calculating mass decrements
    impulse = SmallManeuverAnalyticalModel(
        state, Vector3D(-1 * float(dv), 0.0, 0.0), 340.0
    )
    return impulse


def final_propagation(
    propagator: NumericalPropagator,
    initial_state: SpacecraftState,
    solved_dv: float,
    impulse_at_which_apsis: str,
):
    impulse = create_prograde_impulse(initial_state, solved_dv)
    perturbed_state = impulse.apply(initial_state.shiftedBy(1.0))
    return propagate_to_chosen_apsis(
        propagator,
        perturbed_state,
        perturbed_state.getKeplerianPeriod(),
        impulse_at_which_apsis,
    )


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
    ephem_map = {}
    dv_distribution = norm(loc=dv_initial_guess, scale=perturbation)
    dvs_array = dv_distribution.rvs(size=number_of_parallel_trajectories)
    module_logger.info(
        f"Starting with an initial guess for the DV of {dv_initial_guess}, "
        f"and performing {number_of_parallel_trajectories} perturbatios with STD {perturbation}"
    )
    dvs_array = np.insert(dvs_array, 0, [0.0])
    for cnt, dv in enumerate(dvs_array):

        if cnt == 0:
            module_logger.info("Propagating the unperturbed trajectory")
        else:
            module_logger.info(
                f"Performing perturbation {cnt} of {len(dvs_array)}, imposing a dV of {dv} m/s"
            )
        impulse = create_prograde_impulse(initial_state, dv)
        perturbed_state = impulse.apply(
            initial_state.shiftedBy(1.0)
        )  # Must shift the time by at least 1 second for the impulse to be applied

        state_at_chosen_apsis, ephemeris = propagate_to_chosen_apsis(
            propagator,
            perturbed_state,
            perturbed_state.getKeplerianPeriod(),
            shoot_to_apsis,
        )
        propagated_states_map.update({dv: state_at_chosen_apsis})
        ephem_map.update({cnt: ephemeris})

    module_logger.info("Finished shooting")

    unperturbed_final_state = propagated_states_map.pop(0.0)

    return unperturbed_final_state, propagated_states_map, ephem_map


def calc_separation(
    primary_state: SpacecraftState, secondary_state: SpacecraftState
) -> float:

    # Interpolate the secondary state at the same date of the primary
    # This can easily be achieved by Orekit via Keplerian + quadratic propagation,
    # which is accurate enough for small shifts. Both forward and backward propagation are allowed
    primary_epoch = absolutedate_to_datetime(primary_state.getDate())
    secondary_epoch = absolutedate_to_datetime(secondary_state.getDate())
    shift = (primary_epoch - secondary_epoch).total_seconds()
    secondary_state_shifted = secondary_state.shiftedBy(shift)

    # Forming the RIC frame transform
    primary_pos_inertial = np.array(
        primary_state.getPVCoordinates().getPosition().toArray()
    )
    secondary_pos_inertial = np.array(
        secondary_state_shifted.getPVCoordinates().getPosition().toArray()
    )
    primary_vel_inertial = np.array(
        primary_state.getPVCoordinates().getVelocity().toArray()
    )
    secondary_vel_inertial = np.array(
        secondary_state_shifted.getPVCoordinates().getVelocity().toArray()
    )

    R = primary_pos_inertial / np.linalg.vector_norm(primary_pos_inertial)
    r_cross_v = np.cross(primary_pos_inertial, primary_vel_inertial)
    r_cross_v_norm = np.linalg.vector_norm(r_cross_v)
    I = r_cross_v / r_cross_v_norm
    C = np.cross(I, R)

    DCM_inertial_to_RIC = [R, I, C]

    x_distance = (
        primary_state.getPVCoordinates().getPosition().getX()
        - secondary_state_shifted.getPVCoordinates().getPosition().getX()
    )
    y_distance = (
        primary_state.getPVCoordinates().getPosition().getY()
        - secondary_state_shifted.getPVCoordinates().getPosition().getY()
    )
    z_distance = (
        primary_state.getPVCoordinates().getPosition().getZ()
        - secondary_state_shifted.getPVCoordinates().getPosition().getZ()
    )

    distance_vector_inertial = np.array([x_distance, y_distance, z_distance])

    distance_vector_RIC = np.matmul(DCM_inertial_to_RIC, distance_vector_inertial)

    module_logger.debug(
        f"Calculated RIC between reference and perturbed trajectory. "
        f"Epoch: {primary_state.getDate().toString()}. "
        f"Primary position [km]: {primary_pos_inertial / 1000}. "
        f"Primary velocity [km/sec]: {primary_vel_inertial / 1000}. "
        f"Secondary position [km]: {secondary_pos_inertial / 1000}. "
        f"Secondary position [km/sec]: {secondary_vel_inertial / 1000}. "
        f"RIC: {distance_vector_RIC}"
    )

    return np.sqrt(np.sum([i**2 for i in distance_vector_RIC]))


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


def get_curve_fit(x, y):

    # A 3rd degree polynomial is probably overkill because the relationship is very linear
    def polynomial(x, a, b, c, d):
        return a * x**3 + b * x**2 + c * x + d

    popt, _ = curve_fit(polynomial, x, y)

    module_logger.debug(
        f"The curve fit of the dV data produced the following coefficients a * x**3 + b * x**2 + c * x + d = {popt}"
    )

    return polynomial, popt


def root_finder(
    func: Callable,
    coefficients: List[float],
    initial_guess: float,
    desired_value: float,
) -> float:

    updated_coefficients = coefficients.copy()
    updated_coefficients[-1] = coefficients[-1] - desired_value

    root, results = newton(
        func, initial_guess, args=updated_coefficients, full_output=True
    )

    module_logger.debug(f"Newton root finding: convergence={results.converged}")

    return root


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
        dV_initial_guess = args.dV_initial_guess
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

        state_at_separation, _ = propagate_to_chosen_apsis(
            propagator_num,
            initialState,
            initialOrbit.getKeplerianPeriod(),
            impulse_at_which_apsis,
        )
        module_logger.info(
            f"Simulating separation at epoch {absolutedate_to_datetime(state_at_separation.getDate()).isoformat()}"
        )

        reference_final_state, dv_states_map, ephem_generators_map = shooter(
            state_at_separation,
            propagator_num,
            dV_initial_guess,
            impulse_at_which_apsis,
            perturbation=dV_initial_guess / 10,
        )

        dv_separations_map = process_rics_at_final_epoch(
            reference_final_state, dv_states_map
        )

        func, coefficients = get_curve_fit(
            list(dv_separations_map.keys()), list(dv_separations_map.values())
        )

        root = root_finder(func, coefficients, dV_initial_guess, separation_to_achieve)

        final_state, _ = final_propagation(
            propagator_num, state_at_separation, root, impulse_at_which_apsis
        )

        report_results(
            root,
            separation_to_achieve,
            impulse_at_which_apsis,
            state_at_separation,
            final_state,
        )

        fig, ax = plt.subplots()
        ax = scatter_separation_as_function_of_dv(
            ax, impulse_at_which_apsis, dv_separations_map
        )

        ax = plot_fitted_curve(
            ax, np.array(list(dv_separations_map.keys())), func, coefficients
        )

        ax = plot_root(ax, root, separation_to_achieve)

        ax.legend()
        fig.savefig("final_plot.png")

        save_shooter_ephemeris(
            ephem_generators_map.values(),
            "shooting",
        )
    except Exception as global_exc:
        print(f"Global unhandled exception caught: {global_exc}")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
