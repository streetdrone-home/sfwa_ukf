#include <Eigen/Core>
#include <Eigen/Geometry>

#include <benchmark/benchmark.h>

#include "sfwa_ukf/MeasurementVector.h"
#include "sfwa_ukf/StateVector.h"
#include "sfwa_ukf/Types.h"

enum MyStateVectorFields
{
  AngularVelocity,
  Altitude,
  Velocity,
  Attitude
};

using MyStateVector =
    sfwa_ukf::StateVector<sfwa_ukf::Field<Velocity, sfwa_ukf::Vector<3>>,
                          sfwa_ukf::Field<AngularVelocity, sfwa_ukf::Vector<3>>,
                          sfwa_ukf::Field<Attitude, sfwa_ukf::Quaternion>, sfwa_ukf::Field<Altitude, real_t>>;

enum MyFields
{
  StaticPressure,
  DynamicPressure,
  Accelerometer,
  Gyroscope
};

using MV_Fixed =
    sfwa_ukf::FixedMeasurementVector<sfwa_ukf::Field<Accelerometer, sfwa_ukf::Vector<3>>,
                                     sfwa_ukf::Field<Gyroscope, sfwa_ukf::Vector<3>>,
                                     sfwa_ukf::Field<StaticPressure, real_t>, sfwa_ukf::Field<DynamicPressure, real_t>>;

namespace sfwa_ukf
{
template <>
template <>
sfwa_ukf::Vector<3> MV_Fixed::expected_measurement<MyStateVector, Accelerometer>(const MyStateVector& state)
{
  return state.get_field<Attitude>() * sfwa_ukf::Vector<3>(0, 0, -9.8);
}

template <>
template <>
sfwa_ukf::Vector<3> MV_Fixed::expected_measurement<MyStateVector, Gyroscope>(const MyStateVector& state)
{
  return state.get_field<AngularVelocity>();
}

template <>
template <>
real_t MV_Fixed::expected_measurement<MyStateVector, StaticPressure>(const MyStateVector& state)
{
  return 101.3 - 1.2 * (state.get_field<Altitude>() / 100.0);
}

template <>
template <>
real_t MV_Fixed::expected_measurement<MyStateVector, DynamicPressure>(const MyStateVector& state)
{
  return 0.5 * 1.225 * state.get_field<Velocity>().squaredNorm();
}

}  // namespace sfwa_ukf

using MV_Dynamic = sfwa_ukf::DynamicMeasurementVector<
    sfwa_ukf::Field<Accelerometer, sfwa_ukf::Vector<3>>, sfwa_ukf::Field<Gyroscope, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<StaticPressure, real_t>, sfwa_ukf::Field<DynamicPressure, real_t>>;

namespace sfwa_ukf
{
template <>
template <>
sfwa_ukf::Vector<3> MV_Dynamic::expected_measurement<MyStateVector, Accelerometer>(const MyStateVector& state)
{
  return state.get_field<Attitude>() * sfwa_ukf::Vector<3>(0, 0, -9.8);
}

template <>
template <>
sfwa_ukf::Vector<3> MV_Dynamic::expected_measurement<MyStateVector, Gyroscope>(const MyStateVector& state)
{
  return state.get_field<AngularVelocity>();
}

template <>
template <>
real_t MV_Dynamic::expected_measurement<MyStateVector, StaticPressure>(const MyStateVector& state)
{
  return 101.3 - 1.2 * (state.get_field<Altitude>() / 100.0);
}

template <>
template <>
real_t MV_Dynamic::expected_measurement<MyStateVector, DynamicPressure>(const MyStateVector& state)
{
  return 0.5 * 1.225 * state.get_field<Velocity>().squaredNorm();
}

}  // namespace sfwa_ukf

/*
Tests to compare set/get performance between fixed and dynamic measurement vectors.
*/
void MeasurementVectorFixed_SetGetField(benchmark::State& state)
{
  MV_Fixed test_measurement;
  while (state.KeepRunning())
  {
    test_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(1, 2, 3));
    benchmark::DoNotOptimize(test_measurement.get_field<Accelerometer>());
  }
}

BENCHMARK(MeasurementVectorFixed_SetGetField);

void MeasurementVectorDynamic_SetGetField(benchmark::State& state)
{
  MV_Dynamic test_measurement;
  while (state.KeepRunning())
  {
    test_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(1, 2, 3));
    benchmark::DoNotOptimize(test_measurement.get_field<Accelerometer>());
  }
}

BENCHMARK(MeasurementVectorDynamic_SetGetField);

void MeasurementVectorFixed_SigmaPointGeneration(benchmark::State& state)
{
  MyStateVector test_state;
  MV_Fixed test_measurement;

  test_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<StaticPressure>(0);
  test_measurement.set_field<DynamicPressure>(0);

  test_state.set_field<Velocity>(sfwa_ukf::Vector<3>(1, 2, 3));
  test_state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(1, 0, 0));
  test_state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(1000);

  MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

  while (state.KeepRunning())
  {
    benchmark::DoNotOptimize(test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points));
  }
}

BENCHMARK(MeasurementVectorFixed_SigmaPointGeneration);

void MeasurementVectorDynamic_SigmaPointGeneration(benchmark::State& state)
{
  MyStateVector test_state;
  MV_Dynamic test_measurement;

  test_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<StaticPressure>(0);
  test_measurement.set_field<DynamicPressure>(0);

  test_state.set_field<Velocity>(sfwa_ukf::Vector<3>(1, 2, 3));
  test_state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(1, 0, 0));
  test_state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(1000);

  MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

  while (state.KeepRunning())
  {
    benchmark::DoNotOptimize(test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points));
  }
}

BENCHMARK(MeasurementVectorDynamic_SigmaPointGeneration);

void MeasurementVectorFixed_FullMeasurementCalculation(benchmark::State& state)
{
  MyStateVector test_state;
  MV_Fixed test_measurement;

  test_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<StaticPressure>(0);
  test_measurement.set_field<DynamicPressure>(0);

  test_state.set_field<Velocity>(sfwa_ukf::Vector<3>(1, 2, 3));
  test_state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(1, 0, 0));
  test_state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(1000);

  MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

  MV_Fixed::SigmaPointDistribution<MyStateVector> measurement_sigma_points;
  MV_Fixed mean_measurement;
  MV_Fixed::SigmaPointDeltas<MyStateVector> sigma_point_deltas;
  MV_Fixed::CovarianceMatrix calculated_covariance;

  while (state.KeepRunning())
  {
    measurement_sigma_points = test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);
    mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    sigma_point_deltas = mean_measurement.calculate_sigma_point_deltas<MyStateVector>(measurement_sigma_points);
    calculated_covariance = mean_measurement.calculate_sigma_point_covariance<MyStateVector>(sigma_point_deltas);
  }
}

BENCHMARK(MeasurementVectorFixed_FullMeasurementCalculation);

void MeasurementVectorDynamic_FullMeasurementCalculation(benchmark::State& state)
{
  MyStateVector test_state;
  MV_Dynamic test_measurement;

  test_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<Gyroscope>(sfwa_ukf::Vector<3>(0, 0, 0));
  test_measurement.set_field<StaticPressure>(0);
  test_measurement.set_field<DynamicPressure>(0);

  test_state.set_field<Velocity>(sfwa_ukf::Vector<3>(1, 2, 3));
  test_state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(1, 0, 0));
  test_state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(1000);

  MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

  MV_Dynamic::SigmaPointDistribution<MyStateVector> measurement_sigma_points(test_measurement.size(),
                                                                             MyStateVector::num_sigma());
  MV_Dynamic mean_measurement(test_measurement.size());
  MV_Dynamic::SigmaPointDeltas<MyStateVector> sigma_point_deltas(test_measurement.size(), MyStateVector::num_sigma());
  MV_Dynamic::CovarianceMatrix calculated_covariance(test_measurement.size(), test_measurement.size());

  while (state.KeepRunning())
  {
    measurement_sigma_points = test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);
    mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    sigma_point_deltas = mean_measurement.calculate_sigma_point_deltas<MyStateVector>(measurement_sigma_points);
    calculated_covariance = mean_measurement.calculate_sigma_point_covariance<MyStateVector>(sigma_point_deltas);
  }
}

BENCHMARK(MeasurementVectorDynamic_FullMeasurementCalculation);

void MeasurementVectorFixed_MeasurementCovariance(benchmark::State& state)
{
  MV_Fixed test_measurement, expected_measurement;
  MV_Fixed::CovarianceVector measurement_covariance;

  measurement_covariance.set_field<Accelerometer>(sfwa_ukf::Vector<3>(1, 2, 3));
  measurement_covariance.set_field<Gyroscope>(sfwa_ukf::Vector<3>(4, 5, 6));
  measurement_covariance.set_field<StaticPressure>(7);
  measurement_covariance.set_field<DynamicPressure>(8);

  expected_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(1, 2, 3));
  expected_measurement.set_field<Gyroscope>(sfwa_ukf::Vector<3>(4, 5, 6));
  expected_measurement.set_field<StaticPressure>(7);
  expected_measurement.set_field<DynamicPressure>(8);

  while (state.KeepRunning())
  {
    benchmark::DoNotOptimize(
        test_measurement.calculate_measurement_covariance(measurement_covariance, expected_measurement));
  }
}

BENCHMARK(MeasurementVectorFixed_MeasurementCovariance);

void MeasurementVectorDynamic_MeasurementCovariance(benchmark::State& state)
{
  MV_Dynamic test_measurement, expected_measurement;
  MV_Dynamic::CovarianceVector measurement_covariance;

  measurement_covariance.set_field<Accelerometer>(sfwa_ukf::Vector<3>(1, 2, 3));
  measurement_covariance.set_field<Gyroscope>(sfwa_ukf::Vector<3>(4, 5, 6));
  measurement_covariance.set_field<StaticPressure>(7);
  measurement_covariance.set_field<DynamicPressure>(8);

  expected_measurement.set_field<Accelerometer>(sfwa_ukf::Vector<3>(1, 2, 3));
  expected_measurement.set_field<Gyroscope>(sfwa_ukf::Vector<3>(4, 5, 6));
  expected_measurement.set_field<StaticPressure>(7);
  expected_measurement.set_field<DynamicPressure>(8);

  while (state.KeepRunning())
  {
    benchmark::DoNotOptimize(
        test_measurement.calculate_measurement_covariance(measurement_covariance, expected_measurement));
  }
}

BENCHMARK(MeasurementVectorDynamic_MeasurementCovariance);
