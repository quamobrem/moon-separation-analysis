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

1. install `conda` from 