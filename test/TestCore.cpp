#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <gtest/gtest.h>

#include "comparisons.h"
#include "sfwa_ukf/Core.h"
#include "sfwa_ukf/Integrator.h"
#include "sfwa_ukf/MeasurementVector.h"
#include "sfwa_ukf/StateVector.h"
#include "sfwa_ukf/Types.h"

/* Set up state vector class. */
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

/*
These versions of the predicted measurement functions take kinematic
acceleration and angular acceleration as inputs. Note that in reality, the
inputs would probably be a control vector and the accelerations would be
calculated using the state vector and a dynamics model.
*/
template <>
template <>
sfwa_ukf::Vector<3>
MyMeasurementVector::expected_measurement<MyStateVector, GPS_Position>(const MyStateVector& state,
                                                                       const sfwa_ukf::Vector<3>& acceleration,
                                                                       const sfwa_ukf::Vector<3>& angular_acceleration)
{
  return state.get_field<Position>();
}

template <>
template <>
sfwa_ukf::Vector<3>
MyMeasurementVector::expected_measurement<MyStateVector, GPS_Velocity>(const MyStateVector& state,
                                                                       const sfwa_ukf::Vector<3>& acceleration,
                                                                       const sfwa_ukf::Vector<3>& angular_acceleration)
{
  return state.get_field<Velocity>();
}

template <>
template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement<MyStateVector, Accelerometer, sfwa_ukf::Vector<3>>(
    const MyStateVector& state, const sfwa_ukf::Vector<3>& acceleration, const sfwa_ukf::Vector<3>& angular_acceleration)
{
  return state.get_field<Attitude>() * sfwa_ukf::Vector<3>(0, 0, -9.8) + acceleration;
}

template <>
template <>
sfwa_ukf::FieldVector MyMeasurementVector::expected_measurement<MyStateVector, Magnetometer, sfwa_ukf::Vector<3>>(
    const MyStateVector& state, const sfwa_ukf::Vector<3>& acceleration, const sfwa_ukf::Vector<3>& angular_acceleration)
{
  return state.get_field<Attitude>() * sfwa_ukf::FieldVector(1, 0, 0);
}

template <>
template <>
sfwa_ukf::Vector<3> MyMeasurementVector::expected_measurement<MyStateVector, Gyroscope, sfwa_ukf::Vector<3>>(
    const MyStateVector& state, const sfwa_ukf::Vector<3>& acceleration, const sfwa_ukf::Vector<3>& angular_acceleration)
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

TEST(CoreTest, Initialisation)
{
  MyCore test_filter = create_initialised_test_filter();
}

/*
All these tests check that the estimated state matches the 'actual' state to
within 2-sigma.
*/
TEST(CoreTest, APrioriStep)
{
  MyCore test_filter = create_initialised_test_filter();

  test_filter.a_priori_step(0.01);

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}

TEST(CoreTest, APrioriStepWithInputs)
{
  MyCore test_filter = create_initialised_test_filter();

  test_filter.a_priori_step(0.01, sfwa_ukf::Vector<3>(0, 0, -5), sfwa_ukf::Vector<3>(1, 0, 0));

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}

TEST(CoreTest, InnovationStep)
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

  /*
  With the field vector, we expect the determinant to be approximately zero,
  so allow for it to be slightly negative due to numerical precision.
  */
  EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

  test_filter.a_posteriori_step();

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}

TEST(CoreTest, InnovationStepPartialMeasurement)
{
  MyCore test_filter = create_initialised_test_filter();
  MyMeasurementVector m;

  m.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, -9.8));
  m.set_field<Magnetometer>(sfwa_ukf::FieldVector(0, -1, 0));
  m.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.5, 0, 0));

  test_filter.a_priori_step(0.01);
  test_filter.innovation_step(m);

  EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

  test_filter.a_posteriori_step();

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}

TEST(CoreTest, InnovationStepWithInputs)
{
  MyCore test_filter = create_initialised_test_filter();
  MyMeasurementVector m;

  m.set_field<GPS_Position>(sfwa_ukf::Vector<3>(100, 10, -50));
  m.set_field<GPS_Velocity>(sfwa_ukf::Vector<3>(20, 0, 0));
  m.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, -15));
  m.set_field<Magnetometer>(sfwa_ukf::FieldVector(0, -1, 0));
  m.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.5, 0, 0));

  test_filter.a_priori_step(0.01, sfwa_ukf::Vector<3>(0, 0, -5), sfwa_ukf::Vector<3>(1, 0, 0));
  test_filter.innovation_step(m, sfwa_ukf::Vector<3>(0, 0, -5), sfwa_ukf::Vector<3>(1, 0, 0));

  EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

  test_filter.a_posteriori_step();

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}

TEST(CoreTest, InnovationStepPartialMeasurementWithInputs)
{
  MyCore test_filter = create_initialised_test_filter();
  MyMeasurementVector m;

  m.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, -15));
  m.set_field<Magnetometer>(sfwa_ukf::FieldVector(0, -1, 0));
  m.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.5, 0, 0));

  test_filter.a_priori_step(0.01, sfwa_ukf::Vector<3>(0, 0, -5), sfwa_ukf::Vector<3>(1, 0, 0));
  test_filter.innovation_step(m, sfwa_ukf::Vector<3>(0, 0, -5), sfwa_ukf::Vector<3>(1, 0, 0));

  EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

  test_filter.a_posteriori_step();

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}

TEST(CoreTest, APosterioriStep)
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
  test_filter.a_posteriori_step();

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}

TEST(CoreTest, FullStep)
{
  MyCore test_filter = create_initialised_test_filter();
  MyMeasurementVector m;

  m.set_field<GPS_Position>(sfwa_ukf::Vector<3>(100, 10, -50));
  m.set_field<GPS_Velocity>(sfwa_ukf::Vector<3>(20, 0, 0));
  m.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, -14.8));
  m.set_field<Magnetometer>(sfwa_ukf::FieldVector(0, -1, 0));
  m.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0.5, 0, 0));

  test_filter.step(0.01, m, sfwa_ukf::Vector<3>(0, 0, -5), sfwa_ukf::Vector<3>(1, 0, 0));

  EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
  EXPECT_LT((sfwa_ukf::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm()) * 2);
  EXPECT_LT(
      2 * std::acos(std::abs(sfwa_ukf::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
      std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm()) * 2);
  EXPECT_LT((sfwa_ukf::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
            std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm()) * 2);
}
