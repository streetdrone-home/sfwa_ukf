#include <Eigen/Core>
#include <Eigen/Geometry>

#include <benchmark/benchmark.h>

#include "sfwa_ukf/Core.h"
#include "sfwa_ukf/Integrator.h"
#include "sfwa_ukf/MeasurementVector.h"
#include "sfwa_ukf/StateVector.h"
#include "sfwa_ukf/Types.h"

enum MyStateFields
{
  Position,
  Velocity,
  Attitude,
  AngularVelocity
};

using MyStateVector =
    sfwa_ukf::StateVector<sfwa_ukf::Field<Position, sfwa_ukf::Vector<3>>, sfwa_ukf::Field<Velocity, sfwa_ukf::Vector<3>>,
                          sfwa_ukf::Field<Attitude, sfwa_ukf::Quaternion>,
                          sfwa_ukf::Field<AngularVelocity, sfwa_ukf::Vector<3>>>;

namespace sfwa_ukf
{
namespace Parameters
{
template <>
constexpr real_t AlphaSquared<MyStateVector> = 1e-6;
}

/*
State vector process model. One version takes body frame kinematic
acceleration and angular acceleration as inputs, the other doesn't (assumes
zero accelerations).
*/
template <>
template <>
MyStateVector MyStateVector::derivative<sfwa_ukf::Vector<3>, sfwa_ukf::Vector<3>>(
    const sfwa_ukf::Vector<3>& acceleration, const sfwa_ukf::Vector<3>& angular_acceleration) const
{
  MyStateVector temp;

  /* Position derivative. */
  temp.set_field<Position>(get_field<Velocity>());

  /* Velocity derivative. */
  temp.set_field<Velocity>(get_field<Attitude>().conjugate() * acceleration);

  /* Attitude derivative. */
  sfwa_ukf::Quaternion temp_q;
  temp_q.vec() = get_field<AngularVelocity>();
  temp_q.w() = 0;
  temp.set_field<Attitude>(temp_q);

  /* Angular velocity derivative. */
  temp.set_field<AngularVelocity>(angular_acceleration);

  return temp;
}

template <>
template <>
MyStateVector MyStateVector::derivative<>() const
{
  return derivative(sfwa_ukf::Vector<3>(0, 0, 0), sfwa_ukf::Vector<3>(0, 0, 0));
}

}  // namespace sfwa_ukf

/* Set up measurement vector class. */
enum MyMeasurementFields
{
  GPS_Position,
  GPS_Velocity,
  Accelerometer,
  Magnetometer,
  Gyroscope
};

using MyMeasurementVector = sfwa_ukf::DynamicMeasurementVector<
    sfwa_ukf::Field<GPS_Position, sfwa_ukf::Vector<3>>, sfwa_ukf::Field<GPS_Velocity, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Accelerometer, sfwa_ukf::Vector<3>>, sfwa_ukf::Field<Magnetometer, sfwa_ukf::FieldVector>,
    sfwa_ukf::Field<Gyroscope, sfwa_ukf::Vector<3>>>;

using MyCore = sfwa_ukf::Core<MyStateVector, MyMeasurementVector, sfwa_ukf::IntegratorRK4>;

namespace sfwa_ukf
{
/*
Define measurement model to be used in tests. NOTE: These are just for
testing, don't expect them to make any physical sense whatsoever.
*/
template <>
template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement<MyStateVector, GPS_Position>(const MyStateVector& state)
{
  return state.get_field<Position>();
}

template <>
template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement<MyStateVector, GPS_Velocity>(const MyStateVector& state)
{
  return state.get_field<Velocity>();
}

template <>
template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement<MyStateVector, Accelerometer>(const MyStateVector& state)
{
  return state.get_field<Attitude>() * sfwa_ukf::Vector<3>(0, 0, -9.8);
}

template <>
template <>
sfwa_ukf::FieldVector MyMeasurementVector::expected_measurement<MyStateVector, Magnetometer>(const MyStateVector& state)
{
  return state.get_field<Attitude>() * sfwa_ukf::FieldVector(1, 0, 0);
}

template <>
template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement<MyStateVector, Gyroscope>(const MyStateVector& state)
{
  return state.get_field<AngularVelocity>();
}

}  // namespace sfwa_ukf

MyCore create_initialised_test_filter()
{
  MyCore test_filter;
  test_filter.state.set_field<Position>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_filter.state.set_field<Velocity>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_filter.state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
  test_filter.state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_filter.covariance = MyStateVector::CovarianceMatrix::Zero();
  test_filter.covariance.diagonal() << 10000, 10000, 10000, 100, 100, 100, 1, 1, 5, 10, 10, 10;
  test_filter.measurement_covariance << 10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;

  real_t a, b;
  real_t dt = 0.01;
  a = std::sqrt(0.1 * dt * dt);
  b = std::sqrt(0.1 * dt);
  test_filter.process_noise_covariance << a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, b, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, b;

  return test_filter;
}

void Core_APrioriStep(benchmark::State& state)
{
  MyCore test_filter = create_initialised_test_filter();

  while (state.KeepRunning())
  {
    test_filter.a_priori_step(0.01);
  }
}

BENCHMARK(Core_APrioriStep);

void Core_InnovationStep(benchmark::State& state)
{
  MyCore test_filter = create_initialised_test_filter();
  MyMeasurementVector m;

  m.set_field<GPS_Position>(sfwa_ukf::Vector<3>(100, 10, -50));
  m.set_field<GPS_Velocity>(sfwa_ukf::Vector<3>(20, 0, 0));
  m.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, -9.8));
  m.set_field<Magnetometer>(sfwa_ukf::FieldVector(0, -1, 0));
  m.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.5, 0, 0));

  test_filter.a_priori_step(0.01);

  while (state.KeepRunning())
  {
    test_filter.innovation_step(m);
  }
}

BENCHMARK(Core_InnovationStep);

void Core_APosterioriStep(benchmark::State& state)
{
  MyCore test_filter = create_initialised_test_filter();
  MyMeasurementVector m;

  m.set_field<GPS_Position>(sfwa_ukf::Vector<3>(100, 10, -50));
  m.set_field<GPS_Velocity>(sfwa_ukf::Vector<3>(20, 0, 0));
  m.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, -9.8));
  m.set_field<Magnetometer>(sfwa_ukf::FieldVector(0, -1, 0));
  m.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.5, 0, 0));

  test_filter.a_priori_step(0.01);
  test_filter.innovation_step(m);

  MyStateVector::CovarianceMatrix initial_cov = test_filter.covariance;
  MyStateVector initial_state = test_filter.state;

  while (state.KeepRunning())
  {
    test_filter.covariance = initial_cov;
    test_filter.state = initial_state;
    test_filter.a_posteriori_step();
  }
}

BENCHMARK(Core_APosterioriStep);
