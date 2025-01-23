# Moon Separation Analysis Tool

This few scripts are a tool to estimate the minimum dV necessary to be applied at either periapsis or apoapsis to achieve a certain separation after one revolution in Lunar Orbit. 

--- 

## How to run this 

### Requirements

This app requires:

* `python3.13`
* `orekit`, which in turns requires an `openjdk` instance
* `numpy`, `scipy` and `matplotlib`

### Setup

The app can be run as a python module directly or via a Docker container for ease of packaging. 

#### Python Setup 

To setup for running as a Python module:

1. install `conda` from https://github.com/conda-forge/miniforge?tab=readme-ov-file#install
2. create a new `conda` environment with `conda create -n <choose a name for the environment> -f environment.yml`
    1. this will create the virtual environment and install the necessary dependencies
3. initialize the environment by running `conda init`
4. activate the environment by running `conda activate`
5. verify the activation of the environment by running `python -m main -h`. This will load up the main dependencies and 
   print the command line options available. 

Running it directly in Python has the advantage of direct access to the resulting ephemeris, 
which can be found in the `shooting` folder that results after a successful run. 

#### Docker Setup

For ease of setup, a containerized version of the scripts is provided. 

A installed version of `Docker` is necessary to run the following operations.

Two options are available:

1. Build the Docker image locally by running in the root of the directory `docker build -t moon-separation-analysis:latest .` (mind the dot!)
2. Pull the Docker image from my docker hub: `docker pull luciobruma/moon-separation-analysis:0.0.1`

The scripts can then be run via `docker run` command; a command line argument must be provided to override the default,
which simply prints the help menu again. Thus, for a basic run: `docker run luciobruma/moon-separation-analysis:0.0.1 --impulse_at_which_apsis apoapsis`

To access the resulting ephemeris, either mount a Docker volume to the `shooting` directory 
(`docker run -v <absolute_path_to_volume>:/shooting luciobruma/moon-separation-analysis:0.0.1 --impulse_at_which_apsis apoapsis`) 
or perform a `docker cp` command after the run. 

## Design Choices

I chose numerical propagation for this problem because I wanted to be able to model 
the CR3BP forces accurately. 

While performing the exercise I realised that the given conditions of periapsis and apoapsis 
(100 km and 10000 km, respectively) fall well within the Sphere of Influence of the Moon (~66000 km of radius). 
This would have allowed a simpler approximation via 2-Body Keplerian motion, which would have made the formulation 
of the equations of motion much simpler (and easily differentiable for minimization).
The choice of numerical propagation made it hard to establish an analytical relationship between 
the imparted dV impulse and the separation after 1 revolution. 

To solve this, I chose the approach of perturbing the state at separation by a normal distribution of dV values 
around an initial guess and propagating to the final epoch and record the RIC against the reference trajectory, then use curve fitting 
to create a simple univariate polynomial which could then be easily solved for via a root finding algorithm. The relationship between `dV` 
and separation proved to be quite linear, which seemed sensible. This approach is akin in spirit to Monte Carlo simulation. 

So the choice of numerical propagation is dubious, but still, during validation (see below), I noticed that the 3rd body effects of the Earth and the Sun do 
have a large effect in accuracy. 

I chose `Orekit`, in its python wrapper, as the main astrodynamics library, `numpy` for some linear algebra and 
`scipy` for the curve fitting and root finding algorithms.


## Validation

The trajectory was validated against the same scenario modelled in STK via Astrogator. The separation at final epoch agrees to around 10 km. 
The difference can probably be ascribed to better force modelling in STK, namely: 
- the better gravity modelling of the Moon potential
- the application of Solar Radiation Pressure

## Potential Improvements

There are several improvements that could be done to this application. 

1. Python was assigned as the programming language to be used, but without this constraint Orekit in its native Java implementation would be superior and more 
   performant. Otherwise, a different native python library would be also a good choice, but I was not familiar with other ones.
1. Orekit provides `CR3BP` classes directly. I couldn't make the force model work but I think this could improve the accuracy. 
1. SRP and Lunar gravity potential could be used to improve the propagation accuracy
1. The impulse was modelled via the `SmallManeuverAnalyticalModel` class of Orekit, using a constructor which abstracts a bit the direction vector with which the `dV` is imparted. Because of time constraints I didn't perform the necessary frame transformations needed to explicitly calculate the thrust vector in inertial frame. 
1. Although potentially unnecessary, the Monte Carlo approach could be useful for similar scenarios. Orekit provides `Field` classes specifically made for this purpose and parallel propagation capabilities that I think would be worth investigating.
1. `conda` is unfortunately needed for `orekit` python wrapper. The `miniforge` Docker images are very bloated. Skipping the wrapper could solve this directly. Otherwise, a multistage 
   Docker build could help get away with just the environment and not `conda`. 

