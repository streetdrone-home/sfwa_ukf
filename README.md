# Sfwa Ukf
[![Build Status](https://github.com/PickNikRobotics/sfwa_ukf/actions/workflows/build_and_test.yaml/badge.svg)](https://github.com/PickNikRobotics/sfwa_ukf/actions/workflows/build_and_test.yaml)

This is a ROS2 port of [sfwa/ukf](https://github.com/sfwa/ukf)

Unscented Kalman filter library. Several different UKF implementations are
provided:

  * Standard Unscented Kalman Filter for state estimation, as originally
    described in [[1]](#ref1), with extensions for quaternions as described
    in [[2]](#ref2).
  * Square-root Unscented Kalman Filter for state estimation, implemented as
    described in [[3]](#ref3).
  * Optimised form of square-root Unscented Kalman filter for parameter
    estimation, implemented as described in [[3]](#ref3).

This library makes use of the [Eigen](https://eigen.tuxfamily.org) library
for linear algebra routines and matrix and vector operations. Heavy use is
made of C++11 and C++14 features in an attempt to avoid any dynamic memory
allocations and maximise opportunities for compile-time optimisations.

A primary goal of the library is to provide efficient UKF implementations for
use on embedded systems, so there is a strong focus on having minimal
dependencies and avoiding non-deterministic operations.

The filter can be compiled using either single or double precision by
choosing one of the following preprocessor definitions:

  * UKF_SINGLE_PRECISION
  * UKF_DOUBLE_PRECISION

## Usage

The library contains a number of class templates which are to be specialised
for the particular application. There are three main classes which need to be
specialised to make up an implementation:

  * State vector
  * Measurement vector
  * Core

The following examples are all derived from the [unit tests](test/), so have
a look at them for more detail.

### State vector

The state vector and measurement vector are made up of a number of fields,
each of which contains a key and a type. The reason for this is to allow
the filter to handle quaternions in the state vector transparently, so
something like the following is allowed:

```C++
enum MyFields {
    LatLon,
    Altitude,
    Velocity,
    Attitude
};

using MyStateVector = sfwa_ukf::StateVector<
    sfwa_ukf::Field<LatLon, sfwa_ukf::Vector<2>>,
    sfwa_ukf::Field<Altitude, real_t>,
    sfwa_ukf::Field<Velocity, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Attitude, sfwa_ukf::Quaternion>
>;
```

Internally, these fields are all stored together as one contiguous Eigen
column vector, and all key lookups are done at compile time.

UKF scaling parameters can be adjusted in the following way:

```C++
template <> constexpr real_t sfwa_ukf::Parameters::AlphaSquared<MyStateVector> = 1e-6;
```

For a description of what the scaling parameters do, see [[2]](#ref2) or read
the comments in the [code](include/UKF/StateVector.h).

### Measurement vector

The measurement vector can be specialised in a similar way, but with the
choice of a fixed or dynamic measurement vector:

```C++
enum MyFields {
    StaticPressure,
    DynamicPressure,
    Accelerometer,
    Gyroscope
};

using MyMeasurementVector = sfwa_ukf::FixedMeasurementVector<
    sfwa_ukf::Field<Accelerometer, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Gyroscope, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<StaticPressure, real_t>,
    sfwa_ukf::Field<DynamicPressure, real_t>
>;
```

or:

```C++
enum MyFields {
    StaticPressure,
    DynamicPressure,
    Accelerometer,
    Gyroscope
};

using MyMeasurementVector = sfwa_ukf::DynamicMeasurementVector<
    sfwa_ukf::Field<Accelerometer, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Gyroscope, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<StaticPressure, real_t>,
    sfwa_ukf::Field<DynamicPressure, real_t>
>;
```

For the fixed measurement vector, all measurements have to be provided every
filter iteration. The dynamic measurement vector allows for a filter where
some measurements are not available at every iteration, and so should only be
fused when they are available.

There is a small performance overhead for the dynamic measurement vector,
but it does not do any dynamic memory allocation.

### Core

The Core class contains the filter state and step function. It can be
specialised as follows:

```C++
using MyCore = sfwa_ukf::Core<
    MyStateVector,
    MyMeasurementVector,
    sfwa_ukf::IntegratorRK4
>;
```

Here, the user-specialised state vector and measurement vector classes are
provided as template parameters, along with the integrator to use for the
process model. For a list of available integrators, see
[Integrator.h](include/UKF/Integrator.h).

The square-root state estimation filter can be used instead like this:

```C++
using MyCore = sfwa_ukf::SquareRootCore<
    MyStateVector,
    MyMeasurementVector,
    sfwa_ukf::IntegratorRK4
>;
```

Specialisations of the process model (for state estimation filters) and
measurement model must also be provided.

### Process model

The process model is implemented as an ODE, with a user-provided function to
calculate the derivative. The ODE is then solved using the integrator method
specified in the Core class specialisation.

Here is an example of a process model for a simple state vector:

```C++
using ProcessModelTestStateVector = sfwa_ukf::StateVector<
    sfwa_ukf::Field<Position, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Velocity, sfwa_ukf::Vector<3>>
>;

template <> template <>
ProcessModelTestStateVector ProcessModelTestStateVector::derivative<>() const {
    ProcessModelTestStateVector temp;
    /* Position derivative. */
    temp.set_field<Position>(get_field<Velocity>());

    /* Velocity derivative. */
    temp.set_field<Velocity>(sfwa_ukf::Vector<3>(0, 0, 0));

    return temp;
}
```

Also, the process model can take an arbitrary number of user-specified
inputs, like this:

```C++
template <> template <>
ProcessModelTestStateVector ProcessModelTestStateVector::derivative<sfwa_ukf::Vector<3>>(
        const sfwa_ukf::Vector<3>& acceleration) const {
    ProcessModelTestStateVector temp;
    /* Position derivative. */
    temp.set_field<Position>(get_field<Velocity>());

    /* Velocity derivative. */
    temp.set_field<Velocity>(acceleration);

    return temp;
}
```

### Measurement model

The measurement model is specified per field, in order to allow the expected
measurement vector to be constructed for the dynamic measurement vector where
not all measurements may be available each iteration. Each measurement model
function takes a state vector as an input.

The state vector type is provided to the measurement model specialisation as
a template parameter; this allows a measurement vector class to be shared
across multiple state vectors, with difference process model defined for
each.

Here is an example of a measurement model for a simple measurement vector and
state vector:

```C++
using MyMeasurementVector = sfwa_ukf::DynamicMeasurementVector<
    sfwa_ukf::Field<StaticPressure, real_t>,
    sfwa_ukf::Field<DynamicPressure, real_t>,
    sfwa_ukf::Field<Accelerometer, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Gyroscope, sfwa_ukf::Vector<3>>
>;

using MyStateVector = sfwa_ukf::StateVector<
    sfwa_ukf::Field<Velocity, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<AngularVelocity, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Attitude, sfwa_ukf::Quaternion>,
    sfwa_ukf::Field<Altitude, real_t>
>;

template <> template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * sfwa_ukf::Vector<3>(0, 0, -9.8);
}

template <> template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope>(const MyStateVector& state) {
    return state.get_field<AngularVelocity>();
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, StaticPressure>(const MyStateVector& state) {
    return 101.3 - 1.2*(state.get_field<Altitude>() / 100.0);
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, DynamicPressure>(const MyStateVector& state) {
    return 0.5 * 1.225 * state.get_field<Velocity>().squaredNorm();
}
```

As with the process model, the measurement model can take an arbitrary number
of user-specified inputs:

```C++
template <> template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer, sfwa_ukf::Vector<3>>(const MyStateVector& state, const sfwa_ukf::Vector<3>& input) {
    return state.get_field<Attitude>() * sfwa_ukf::Vector<3>(0, 0, -9.8) + input;
}

template <> template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope, sfwa_ukf::Vector<3>>(const MyStateVector& state, const sfwa_ukf::Vector<3>& input) {
    return state.get_field<AngularVelocity>();
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, StaticPressure, sfwa_ukf::Vector<3>>(const MyStateVector& state, const sfwa_ukf::Vector<3>& input) {
    return 101.3 - 1.2*(state.get_field<Altitude>() / 100.0);
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, DynamicPressure, sfwa_ukf::Vector<3>>(const MyStateVector& state, const sfwa_ukf::Vector<3>& input) {
    return 0.5 * 1.225 * state.get_field<Velocity>().squaredNorm();
}
```

### Initialisation

The filter state, covariance, process noise covariance and measurement noise
covariance should be initialised to appropriate values, e.g.:

```C++
MyCore test_filter;
test_filter.state.set_field<Position>(sfwa_ukf::Vector<3>(0, 0, 0));
test_filter.state.set_field<Velocity>(sfwa_ukf::Vector<3>(0, 0, 0));
test_filter.state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
test_filter.state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(0, 0, 0));
test_filter.covariance = MyStateVector::CovarianceMatrix::Zero();
test_filter.covariance.diagonal() << 10000, 10000, 10000, 100, 100, 100, 1, 1, 5, 10, 10, 10;
test_filter.process_noise_covariance = MyStateVector::CovarianceMatrix::Identity()*1e-5;
test_filter.measurement_covariance << 10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;

```

Or, for the SR-UKF:
```C++
MyCore test_filter;
test_filter.state.set_field<Position>(sfwa_ukf::Vector<3>(0, 0, 0));
test_filter.state.set_field<Velocity>(sfwa_ukf::Vector<3>(0, 0, 0));
test_filter.state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
test_filter.state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(0, 0, 0));
test_filter.root_covariance = MyStateVector::CovarianceMatrix::Zero();
test_filter.root_covariance.diagonal() << 100, 100, 100, 10, 10, 10, 1, 1, 2.2, 3.2, 3.2, 3.2;
test_filter.process_noise_root_covariance = MyStateVector::CovarianceMatrix::Identity()*3e-2;
test_filter.measurement_root_covariance << 10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;
test_filter.measurement_root_covariance = test_filter.measurement_root_covariance.array().sqrt();
```

Currently, only a diagonal measurement noise covariance matrix is supported.

### Iteration

The general steps for carrying out a filter iteration are something like:

```C++
MyMeasurementVector meas;
meas.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0.0, 0.0, 9.8));
meas.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.0, 0.0, 0.0));
meas.set_field<StaticPressure>(101300.0);
meas.set_field<DynamicPressure>(101300.0);

test_filter.step(0.01, meas);
```

Or, if it's necessary to do things with the internal filter state (e.g.
filter health monitoring), then iteration can be split up into three steps:

```C++
MyMeasurementVector meas;
meas.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0.0, 0.0, 9.8));
meas.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.0, 0.0, 0.0));
meas.set_field<StaticPressure>(101300.0);
meas.set_field<DynamicPressure>(101300.0);

test_filter.a_priori_step(0.01);
/*
At this point, the state and covariance (or root_covariance) variables
reflect the a priori state and (root_)covariance.
*/
test_filter.innovation_step(meas);
/* Innovation and innovation_(root_)covariance variables are now set. */
test_filter.a_posteriori_step();
/* State and (root_)covariance variables are set to the a priori values.
```


## References

<a name="ref1">[1]</a> "A New Extension of the Kalman Filter to Nonlinear Systems:,
S. J. Julier and J. K. Uhlmann,
https://www.cs.unc.edu/~welch/kalman/media/pdf/Julier1997_SPIE_KF.pdf

<a name="ref2">[2]</a> "Unscented Filtering for Spacecraft Attitude Estimation", John L.
Crassidis and F. Landis Markley, http://www.acsu.buffalo.edu/~johnc/uf_att.pdf

<a name="ref3">[3]</a> "The Square-Root Unscented Kalman Filter for State and Parameter-Estimation",
Rudolph van der Merwe and Eric A. Want,


## Installation

### Build from Source

These instructions assume you are running on Ubuntu 20.04:

1. [Install ROS2 foxy](https://docs.ros.org/en/foxy/Installation/Ubuntu-Install-Debians.html). You can stop following along with the tutorial after you complete the section titled: [Environment setup](https://docs.ros.org/en/foxy/Installation/Ubuntu-Install-Debians.html#environment-setup). Make sure you setup your environment with:

        source /opt/ros/foxy/setup.bash

   > *NOTE:* You may want to add that line to your `~/.bashrc`

2. [Install ROS2 Build Tools](https://docs.ros.org/en/foxy/Installation/Ubuntu-Development-Setup.html#install-development-tools-and-ros-tools)

   > *NOTE:* If installing on a fresh OS, run `sudo rosdep init` and `rosdep update` before the install script

3. Create a colcon workspace (*Note:* Feel free to change `~/ws_ros2` to whatever absolute path you want):

        export COLCON_WS=~/ws_ros2/
        mkdir -p $COLCON_WS/src

4. Get the repo:

        cd $COLCON_WS/src
        git clone git@github.com:PickNikRobotics/sfwa_ukf.git

5. Download the required repositories and install any dependencies:

        cd $COLCON_WS/src
        vcs import < sfwa_ukf/sfwa_ukf.repos
        rosdep install --ignore-src --from-paths . -y

        # Pick a ROS_DOMAIN_ID that doesn't clash with others
        echo 'export ROS_DOMAIN_ID='<YOUR-NUMBER> >> ~/.bashrc

7. Configure and build the workspace:

        cd $COLCON_WS
        colcon build --symlink-install --event-handlers log-

8. Source the workspace.

        source $COLCON_WS/install/setup.bash

> *Note*: Whenever you open a new terminal be sure to run `source /opt/ros/foxy/setup.bash`, `export COLCON_WS=~/ws_ros2`, and `source $COLCON_WS/install/setup.bash`. You may want to add those commands to your `~/.bashrc`

## For Developers

### Quickly update code repositories

To make sure you have the latest repos:

      cd $COLCON_WS/src/sfwa_ukf
      git checkout main
      git pull origin main
      cd $COLCON_WS/src
      vcs import < sfwa_ukf/sfwa_ukf.repos
      rosdep install --from-paths . --ignore-src -y

### Setup pre-commit

pre-commit is a tool to automatically run formatting checks on each commit, which saves you from manually running clang-format (or, crucially, from forgetting to run them!).

Install pre-commit like this:

```
pip3 install pre-commit
```

Run this in the top directory of the repo to set up the git hooks:

```
pre-commit install
```

### Testing and Linting

To test the packages in sfwa_ukf, use the following command with [colcon](https://colcon.readthedocs.io/en/released/).

    export TEST_PACKAGES="PROJECT_PACKAGE_NAMES"
    colcon build --packages-up-to ${TEST_PACKAGES}
    colcon test --packages-select ${TEST_PACKAGES}
    colcon test-result

### Using ccache

ccache is a useful tool to speed up compilation times with GCC or any other sufficiently similar compiler.

To install ccache on Linux:

    sudo apt-get install ccache

For other OS, search the package manager or software store for ccache, or refer to the [ccache website](https://ccache.dev/)

#### Setup

To use ccache after installing it there are two methods. you can add it to your PATH or you can configure it for more specific uses.

ccache must be in front of your regular compiler or it won't be called. It is recommended that you add this line to your `.bashrc`:

    export PATH=/usr/lib/ccache:$PATH

To configure ccache for more particular uses, set the CC and CXX environment variables before invoking make, cmake, catkin_make or catkin build.

For more information visit the [ccache website](https://ccache.dev/).
