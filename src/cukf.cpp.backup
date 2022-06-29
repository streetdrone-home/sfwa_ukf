#include "sfwa_ukf/cukf.h"

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "sfwa_ukf/Core.h"
#include "sfwa_ukf/Integrator.h"
#include "sfwa_ukf/MeasurementVector.h"
#include "sfwa_ukf/StateVector.h"
#include "sfwa_ukf/Types.h"

/*
This file implements the same UKF as was used for the 2014 UAV Challenge.
There are a number of efficiency and performance improvements that could be
made with the new library structure, but the main purpose of this is for
direct comparison of performance of the new library against the old
implementation.
*/

enum SFWA_States
{
  LatLon,
  Altitude,
  Velocity,
  Acceleration,
  Attitude,
  AngularVelocity,
  AngularAcceleration,
  WindVelocity,
  GyroBias
};

using SFWA_StateVector = sfwa_ukf::StateVector<
    sfwa_ukf::Field<LatLon, sfwa_ukf::Vector<2>>,              /* Latitude and longitude (rad) */
    sfwa_ukf::Field<Altitude, real_t>,                         /* Altitude above the WGS84 ellipsoid (m) */
    sfwa_ukf::Field<Velocity, sfwa_ukf::Vector<3>>,            /* Velocity (NEED frame, m/s) */
    sfwa_ukf::Field<Acceleration, sfwa_ukf::Vector<3>>,        /* Acceleration (body frame, m/s^2) */
    sfwa_ukf::Field<Attitude, sfwa_ukf::Quaternion>,           /* Attitude as a quaternion (NEED frame to body frame) */
    sfwa_ukf::Field<AngularVelocity, sfwa_ukf::Vector<3>>,     /* Angular velocity (body frame, rad/s) */
    sfwa_ukf::Field<AngularAcceleration, sfwa_ukf::Vector<3>>, /* Angular acceleration (body frame, rad/s^2) */
    sfwa_ukf::Field<WindVelocity, sfwa_ukf::Vector<3>>,        /* Wind velocity (NEED frame, m/s) */
    sfwa_ukf::Field<GyroBias, sfwa_ukf::Vector<3>>             /* Gyro bias (body frame, rad/s) */
    >;

/* WGS84 reference ellipsoid constants. */
#define WGS84_A (6378137.0)
#define WGS84_B (6356752.314245)
#define WGS84_A2 (WGS84_A * WGS84_A)
#define WGS84_B2 (WGS84_B * WGS84_B)
#define WGS84_AB2 (WGS84_A2 * WGS84_B2)

#define G_ACCEL (9.80665)
#define RHO (1.225)

/* SFWA vehicle dynamics model. Not used for this comparison. */
sfwa_ukf::Vector<6> x8_dynamics_model(const SFWA_StateVector& state, const sfwa_ukf::Vector<3>& control)
{
  sfwa_ukf::Vector<6> output;

  /* Cache state data for convenience */
  sfwa_ukf::Quaternion attitude = state.get_field<Attitude>();
  real_t yaw_rate = state.get_field<AngularVelocity>()[2], pitch_rate = state.get_field<AngularVelocity>()[1],
         roll_rate = state.get_field<AngularVelocity>()[0];

  /* External axes */
  sfwa_ukf::Vector<3> airflow = attitude * (state.get_field<WindVelocity>() - state.get_field<Velocity>());

  /*
  Calculate axial airflow
  */
  real_t airflow_x2, airflow_y2, airflow_z2, airflow_v2;
  airflow_x2 = airflow[0] * airflow[0];
  airflow_y2 = airflow[1] * airflow[1];
  airflow_z2 = airflow[2] * airflow[2];
  airflow_v2 = airflow_x2 + airflow_y2 + airflow_z2;

  /*
  Determine motor thrust and torque.
  */
  real_t rpm = control[0] * 12000.0, thrust, ve2 = (0.0025 * 0.0025) * rpm * rpm;
  /* 1 / 3.8kg times area * density of air */
  thrust = (ve2 - airflow_v2) * (0.26315789473684 * 0.5 * RHO * 0.02);

  /*
  Calculate airflow in the horizontal and vertical planes, as well as
  pressure
  */
  real_t v_inv, vertical_v, vertical_v_inv, qbar;

  qbar = (RHO * 0.5) * airflow_v2;
  v_inv = 1.0 / std::sqrt(std::max(1.0, airflow_v2));

  vertical_v = std::sqrt(airflow_x2 + airflow_z2);
  vertical_v_inv = 1.0 / std::max(1.0, vertical_v);

  /* Work out sin/cos of alpha and beta */
  real_t sin_alpha, cos_alpha, sin_beta, cos_beta, sin_cos_alpha;

  sin_beta = airflow[1] * v_inv;
  cos_beta = vertical_v * v_inv;

  sin_alpha = -airflow[2] * vertical_v_inv;
  cos_alpha = -airflow[0] * vertical_v_inv;

  sin_cos_alpha = sin_alpha * cos_alpha;

  /* Work out aerodynamic forces in wind frame */
  real_t lift, drag, side_force;

  /* 0.26315789473684 is the reciprocal of mass (3.8kg) */
  lift = (qbar * 0.26315789473684) * (0.8 * sin_cos_alpha + 0.18);
  drag = (qbar * 0.26315789473684) * (0.05 + 0.7 * sin_alpha * sin_alpha);
  side_force = (qbar * 0.26315789473684) * 0.2 * sin_beta * cos_beta;

  /* Convert aerodynamic forces from wind frame to body frame */
  real_t x_aero_f = lift * sin_alpha - drag * cos_alpha - side_force * sin_beta,
         z_aero_f = lift * cos_alpha + drag * sin_alpha, y_aero_f = side_force * cos_beta;

  output.segment<3>(0)
      << sfwa_ukf::Vector<3>(x_aero_f + thrust, y_aero_f, -z_aero_f) + (attitude * sfwa_ukf::Vector<3>(0, 0, G_ACCEL));

  /* Determine moments */
  real_t pitch_moment, yaw_moment, roll_moment, left_aileron = control[1] - 0.5, right_aileron = control[2] - 0.5;
  pitch_moment = 0.0 - 0.0 * sin_alpha - 0.0 * pitch_rate - 0.1 * (left_aileron + right_aileron) * vertical_v * 0.1;
  roll_moment = 0.05 * sin_beta - 0.1 * roll_rate + 0.15 * (left_aileron - right_aileron) * vertical_v * 0.1;
  yaw_moment =
      -0.02 * sin_beta - 0.05 * yaw_rate - 0.02 * (std::abs(left_aileron) + std::abs(right_aileron)) * vertical_v * 0.1;
  pitch_moment *= qbar;
  roll_moment *= qbar;
  yaw_moment *= qbar;

  /*
  Calculate angular acceleration (tau / inertia tensor).
  Inertia tensor is:
      0.3 0 -0.0334
      0 0.17 0
      -0.0334 0 0.405
  So inverse is:
      3.36422 0 0.277444
      0 5.88235 0
      0.277444 0 2.49202
  */
  output.segment<3>(3) << sfwa_ukf::Vector<3>((3.364222 * roll_moment + 0.27744448 * yaw_moment),
                                              10.8823528 * pitch_moment,
                                              (0.27744448 * roll_moment + 2.4920163 * yaw_moment));

  return output;
}

namespace sfwa_ukf
{
/* SFWA state vector process model. */
template <>
template <>
SFWA_StateVector SFWA_StateVector::derivative<>() const
{
  SFWA_StateVector output;

  /* Calculate the normal and meridional radii of curvature. */
  real_t lat = get_field<LatLon>()[0];
  real_t tempA = WGS84_A * std::cos(lat), tempB = WGS84_B * std::sin(lat), temp = tempA * tempA + tempB * tempB,
         temp_sqrt = std::sqrt(temp);
  real_t M = WGS84_AB2 / (temp_sqrt * temp);
  real_t N = WGS84_A2 / temp_sqrt;

  /*
  Calculate change in position. Using the small angle approximation, this
  becomes very simple – no trig required for latitude derivative, and one
  cosine function for longitude derivative.
  */
  sfwa_ukf::Vector<3> vel = get_field<Velocity>();
  output.set_field<LatLon>(
      Vector<2>(vel[0] / (M + get_field<Altitude>()), (vel[1] / (N + get_field<Altitude>())) * std::cos(lat)));
  output.set_field<Altitude>(-vel[2]);

  /* Calculate change in velocity. */
  output.set_field<Velocity>(get_field<Attitude>().conjugate() * get_field<Acceleration>());

  /* Change in linear acceleration is zero. */
  output.set_field<Acceleration>(sfwa_ukf::Vector<3>(0, 0, 0));

  /* Calculate change in attitude. */
  sfwa_ukf::Quaternion omega_q;
  omega_q.vec() = get_field<AngularVelocity>() * 0.5;
  omega_q.w() = 0;
  output.set_field<Attitude>(omega_q.conjugate() * get_field<Attitude>());

  /* Calculate change in angular velocity (just angular acceleration). */
  output.set_field<AngularVelocity>(get_field<AngularAcceleration>());

  /* Change in angular acceleration is zero. */
  output.set_field<AngularAcceleration>(sfwa_ukf::Vector<3>(0, 0, 0));

  /* Change in wind velocity is zero. */
  output.set_field<WindVelocity>(sfwa_ukf::Vector<3>(0, 0, 0));

  /* Change in gyro bias is zero. */
  output.set_field<GyroBias>(sfwa_ukf::Vector<3>(0, 0, 0));

  return output;
}

}  // namespace sfwa_ukf

enum SFWA_Measurements
{
  Accelerometer,
  Gyroscope,
  Magnetometer,
  GPS_Position,
  GPS_Velocity,
  Airspeed,
  PressureAltitude
};

using SFWA_MeasurementVector = sfwa_ukf::DynamicMeasurementVector<
    sfwa_ukf::Field<Accelerometer, sfwa_ukf::Vector<3>>, sfwa_ukf::Field<Gyroscope, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<Magnetometer, sfwa_ukf::Vector<3>>, sfwa_ukf::Field<GPS_Position, sfwa_ukf::Vector<3>>,
    sfwa_ukf::Field<GPS_Velocity, sfwa_ukf::Vector<3>>, sfwa_ukf::Field<Airspeed, real_t>,
    sfwa_ukf::Field<PressureAltitude, real_t>>;

/*
Hard code the magnetic field vector to the WMM value for
-37.954690, 145.237575 for the purposes of this comparison. Value is in
microtesla.
*/
static sfwa_ukf::Vector<3> local_mag_field = sfwa_ukf::Vector<3>(21.2584, 4.4306, -55.9677);

using SFWA_UKF = sfwa_ukf::Core<SFWA_StateVector, SFWA_MeasurementVector, sfwa_ukf::IntegratorRK4>;

namespace sfwa_ukf
{
/* SFWA measurement model. */
template <>
template <>
sfwa_ukf::Vector<3>
SFWA_MeasurementVector::expected_measurement<SFWA_StateVector, Accelerometer>(const SFWA_StateVector& state)
{
  return state.get_field<Acceleration>() + state.get_field<Attitude>() * sfwa_ukf::Vector<3>(0, 0, -G_ACCEL);
}

template <>
template <>
sfwa_ukf::Vector<3>
SFWA_MeasurementVector::expected_measurement<SFWA_StateVector, Gyroscope>(const SFWA_StateVector& state)
{
  return state.get_field<AngularVelocity>() + state.get_field<GyroBias>();
}

template <>
template <>
sfwa_ukf::Vector<3>
SFWA_MeasurementVector::expected_measurement<SFWA_StateVector, Magnetometer>(const SFWA_StateVector& state)
{
  return state.get_field<Attitude>() * local_mag_field;
}

template <>
template <>
sfwa_ukf::Vector<3>
SFWA_MeasurementVector::expected_measurement<SFWA_StateVector, GPS_Position>(const SFWA_StateVector& state)
{
  return sfwa_ukf::Vector<3>(state.get_field<LatLon>()[0], state.get_field<LatLon>()[1], state.get_field<Altitude>());
}

template <>
template <>
sfwa_ukf::Vector<3>
SFWA_MeasurementVector::expected_measurement<SFWA_StateVector, GPS_Velocity>(const SFWA_StateVector& state)
{
  return state.get_field<Velocity>();
}

template <>
template <>
real_t SFWA_MeasurementVector::expected_measurement<SFWA_StateVector, Airspeed>(const SFWA_StateVector& state)
{
  return (state.get_field<Attitude>() * (state.get_field<Velocity>() - state.get_field<WindVelocity>()))[0];
}

template <>
template <>
real_t SFWA_MeasurementVector::expected_measurement<SFWA_StateVector, PressureAltitude>(const SFWA_StateVector& state)
{
  return state.get_field<Altitude>();
}

}  // namespace sfwa_ukf
static SFWA_UKF ukf;
static SFWA_MeasurementVector meas;

/*
The following functions provide a ctypes-compatible interface for ease of
testing.
*/

void ukf_init()
{
  ukf.state.set_field<LatLon>(sfwa_ukf::Vector<2>(-0.662434, 2.534874));
  ukf.state.set_field<Altitude>(0);
  ukf.state.set_field<Velocity>(sfwa_ukf::Vector<3>(0, 0, 0));
  ukf.state.set_field<Acceleration>(sfwa_ukf::Vector<3>(0, 0, 0));
  ukf.state.set_field<Attitude>(sfwa_ukf::Quaternion(1, 0, 0, 0));
  ukf.state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(0, 0, 0));
  ukf.state.set_field<AngularAcceleration>(sfwa_ukf::Vector<3>(0, 0, 0));
  ukf.state.set_field<WindVelocity>(sfwa_ukf::Vector<3>(0, 0, 0));
  ukf.state.set_field<GyroBias>(sfwa_ukf::Vector<3>(0, 0, 0));
  ukf.covariance = SFWA_StateVector::CovarianceMatrix::Zero();
  ukf.covariance.diagonal() << 1, 1, 10000, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 10, 10, 10, 1e-6, 1e-6, 1e-6, 1e-6,
      1e-6, 1e-6, 100, 100, 100, 0.01, 0.01, 0.01;
  ukf.measurement_covariance << 1, 1, 1, 0.01, 0.01, 0.01, 10, 10, 10, 1e-20, 1e-20, 1, 1, 1, 1, 1, 1;

  local_mag_field = sfwa_ukf::Vector<3>(21.2584, 4.4306, -55.9677);

  ukf.process_noise_covariance = SFWA_StateVector::CovarianceMatrix::Zero();
  ukf.process_noise_covariance.diagonal() << 1e-17, 1e-17, 1e-4, /* lat, lon, alt */
      2e-3, 2e-3, 2e-3,                                          /* velocity N, E, D */
      2e-2, 2e-2, 2e-2,                                          /* acceleration x, y, z */
      7e-8, 7e-8, 7e-8,                                          /* attitude roll, pitch, yaw */
      1e-3, 1e-3, 1e-3,                                          /* angular velocity roll, pitch, yaw */
      1e-3, 1e-3, 1e-3,                                          /* angular acceleration roll, pitch, yaw */
      1e-7, 1e-7, 1e-7,    /* wind velocity N, E, D -- NOTE: from FCS armed mode */
      1e-12, 1e-12, 1e-12; /* gyro bias x, y, z -- NOTE: from FCS armed mode */
}

void ukf_set_position(real_t lat, real_t lon, real_t alt)
{
  ukf.state.set_field<LatLon>(sfwa_ukf::Vector<2>(lat, lon));
  ukf.state.set_field<Altitude>(alt);
}

void ukf_set_velocity(real_t x, real_t y, real_t z)
{
  ukf.state.set_field<Velocity>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_set_acceleration(real_t x, real_t y, real_t z)
{
  ukf.state.set_field<Acceleration>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_set_attitude(real_t w, real_t x, real_t y, real_t z)
{
  ukf.state.set_field<Attitude>(sfwa_ukf::Quaternion(w, x, y, z));
}

void ukf_set_angular_velocity(real_t x, real_t y, real_t z)
{
  ukf.state.set_field<AngularVelocity>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_set_angular_acceleration(real_t x, real_t y, real_t z)
{
  ukf.state.set_field<AngularAcceleration>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_set_wind_velocity(real_t x, real_t y, real_t z)
{
  ukf.state.set_field<WindVelocity>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_set_gyro_bias(real_t x, real_t y, real_t z)
{
  ukf.state.set_field<GyroBias>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_get_state(struct ukf_state_t* in)
{
  in->position[0] = ukf.state.get_field<LatLon>()[0];
  in->position[1] = ukf.state.get_field<LatLon>()[1];
  in->position[2] = ukf.state.get_field<Altitude>();
  in->velocity[0] = ukf.state.get_field<Velocity>()[0];
  in->velocity[1] = ukf.state.get_field<Velocity>()[1];
  in->velocity[2] = ukf.state.get_field<Velocity>()[2];
  in->acceleration[0] = ukf.state.get_field<Acceleration>()[0];
  in->acceleration[1] = ukf.state.get_field<Acceleration>()[1];
  in->acceleration[2] = ukf.state.get_field<Acceleration>()[2];
  in->attitude[0] = ukf.state.get_field<Attitude>().x();
  in->attitude[1] = ukf.state.get_field<Attitude>().y();
  in->attitude[2] = ukf.state.get_field<Attitude>().z();
  in->attitude[3] = ukf.state.get_field<Attitude>().w();
  in->angular_velocity[0] = ukf.state.get_field<AngularVelocity>()[0];
  in->angular_velocity[1] = ukf.state.get_field<AngularVelocity>()[1];
  in->angular_velocity[2] = ukf.state.get_field<AngularVelocity>()[2];
  in->angular_acceleration[0] = ukf.state.get_field<AngularAcceleration>()[0];
  in->angular_acceleration[1] = ukf.state.get_field<AngularAcceleration>()[1];
  in->angular_acceleration[2] = ukf.state.get_field<AngularAcceleration>()[2];
  in->wind_velocity[0] = ukf.state.get_field<WindVelocity>()[0];
  in->wind_velocity[1] = ukf.state.get_field<WindVelocity>()[1];
  in->wind_velocity[2] = ukf.state.get_field<WindVelocity>()[2];
  in->gyro_bias[0] = ukf.state.get_field<GyroBias>()[0];
  in->gyro_bias[1] = ukf.state.get_field<GyroBias>()[1];
  in->gyro_bias[2] = ukf.state.get_field<GyroBias>()[2];
}

void ukf_set_state(struct ukf_state_t* in)
{
  ukf.state.set_field<LatLon>(sfwa_ukf::Vector<2>(in->position[0], in->position[1]));
  ukf.state.set_field<Altitude>(in->position[2]);
  ukf.state.set_field<Velocity>(sfwa_ukf::Vector<3>(in->velocity[0], in->velocity[1], in->velocity[2]));
  ukf.state.set_field<Acceleration>(sfwa_ukf::Vector<3>(in->acceleration[0], in->acceleration[1], in->acceleration[2]));
  ukf.state.set_field<Attitude>(
      sfwa_ukf::Quaternion(in->attitude[3], in->attitude[0], in->attitude[1], in->attitude[2]));
  ukf.state.set_field<AngularVelocity>(
      sfwa_ukf::Vector<3>(in->angular_velocity[0], in->angular_velocity[1], in->angular_velocity[2]));
  ukf.state.set_field<AngularAcceleration>(
      sfwa_ukf::Vector<3>(in->angular_acceleration[0], in->angular_acceleration[1], in->angular_acceleration[2]));
  ukf.state.set_field<WindVelocity>(
      sfwa_ukf::Vector<3>(in->wind_velocity[0], in->wind_velocity[1], in->wind_velocity[2]));
  ukf.state.set_field<GyroBias>(sfwa_ukf::Vector<3>(in->gyro_bias[0], in->gyro_bias[1], in->gyro_bias[2]));
}

void ukf_get_state_covariance(
    real_t state_covariance[SFWA_StateVector::covariance_size() * SFWA_StateVector::covariance_size()])
{
  Eigen::Map<typename SFWA_StateVector::CovarianceMatrix> covariance_map(state_covariance);
  covariance_map = ukf.covariance;
}

void ukf_get_state_covariance_diagonal(real_t state_covariance_diagonal[SFWA_StateVector::covariance_size()])
{
  Eigen::Map<sfwa_ukf::Vector<SFWA_StateVector::covariance_size()>> covariance_map(state_covariance_diagonal);
  covariance_map = ukf.covariance.diagonal();
}

void ukf_get_state_error(real_t state_error[SFWA_StateVector::covariance_size()])
{
  Eigen::Map<typename SFWA_StateVector::StateVectorDelta> error_map(state_error);
  error_map = ukf.covariance.cwiseAbs().rowwise().sum().cwiseSqrt();
}

void ukf_sensor_clear()
{
  meas = SFWA_MeasurementVector();
}

void ukf_sensor_set_accelerometer(real_t x, real_t y, real_t z)
{
  meas.set_field<Accelerometer>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_sensor_set_gyroscope(real_t x, real_t y, real_t z)
{
  meas.set_field<Gyroscope>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_sensor_set_magnetometer(real_t x, real_t y, real_t z)
{
  meas.set_field<Magnetometer>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_sensor_set_gps_position(real_t lat, real_t lon, real_t alt)
{
  meas.set_field<GPS_Position>(sfwa_ukf::Vector<3>(lat, lon, alt));
}

void ukf_sensor_set_gps_velocity(real_t x, real_t y, real_t z)
{
  meas.set_field<GPS_Velocity>(sfwa_ukf::Vector<3>(x, y, z));
}

void ukf_sensor_set_pitot_tas(real_t tas)
{
  meas.set_field<Airspeed>(tas);
}

void ukf_sensor_set_barometer_amsl(real_t amsl)
{
  meas.set_field<PressureAltitude>(amsl);
}

void ukf_set_params(struct ukf_ioboard_params_t* in)
{
  local_mag_field << in->mag_field[0], in->mag_field[1], in->mag_field[2];
  ukf.measurement_covariance << in->accel_covariance[0], in->accel_covariance[1], in->accel_covariance[2],
      in->gyro_covariance[0], in->gyro_covariance[1], in->gyro_covariance[2], in->mag_covariance[0],
      in->mag_covariance[1], in->mag_covariance[2], in->gps_position_covariance[0], in->gps_position_covariance[1],
      in->gps_position_covariance[2], in->gps_velocity_covariance[0], in->gps_velocity_covariance[1],
      in->gps_velocity_covariance[2], in->pitot_covariance, in->barometer_amsl_covariance;
}

void ukf_choose_dynamics(enum ukf_model_t t)
{
}

void ukf_set_custom_dynamics_model(ukf_model_function_t func)
{
}

void ukf_iterate(float dt, real_t control_vector[UKF_CONTROL_DIM])
{
  ukf.step(dt, meas);
}

void ukf_set_process_noise(real_t process_noise_covariance[SFWA_StateVector::covariance_size()])
{
  Eigen::Map<typename SFWA_StateVector::StateVectorDelta> covariance_map(process_noise_covariance);
  ukf.process_noise_covariance = SFWA_StateVector::CovarianceMatrix::Zero();
  ukf.process_noise_covariance.diagonal() << covariance_map;
}

uint32_t ukf_config_get_state_dim()
{
  return SFWA_StateVector::covariance_size();
}

uint32_t ukf_config_get_measurement_dim()
{
  return SFWA_MeasurementVector::max_size();
}

uint32_t ukf_config_get_control_dim()
{
  return UKF_CONTROL_DIM;
}

enum ukf_precision_t ukf_config_get_precision()
{
  return UKF_PRECISION_DOUBLE;
}
