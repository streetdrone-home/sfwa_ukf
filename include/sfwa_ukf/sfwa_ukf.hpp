/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2022, PickNik Inc
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of PickNik Inc nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

#pragma once

// C++
#include <memory>
#include <string>

// ROS
#include <eigen3/Eigen/Eigen>

#include <rclcpp/logging.hpp>
#include <rclcpp/rclcpp.hpp>

#include <sensor_msgs/msg/imu.hpp>

#include "ahrs.h"

namespace sfwa_ukf
{
struct ImuState
{
  builtin_interfaces::msg::Time time;
  Eigen::Vector3d acc;
  Eigen::Vector3d angvel;
  Eigen::Quaterniond orientation;
  ImuState& operator=(const ImuState& rhs)
  {
    acc.noalias() = rhs.acc;
    angvel.noalias() = rhs.angvel;
    orientation = rhs.orientation;
    time = rhs.time;
    return *this;
  }
};

class sfwa_ukf : public rclcpp::Node
{
public:
  /** \brief Constructor */
  sfwa_ukf();

  ImuState imu_;
  ImuState prev_imu_;

private:
  void InitFilter();
  void UpdateMeasurements();
  void ImuCallback(const sensor_msgs::msg::Imu::SharedPtr msg)
  {
    prev_imu_ = imu_;
    imu_.time = msg->header.stamp;
    imu_.acc = Eigen::Vector3d(msg->linear_acceleration.x, msg->linear_acceleration.y, msg->linear_acceleration.z);
    imu_.angvel = Eigen::Vector3d(msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z);
    imu_.orientation =
        Eigen::Quaterniond(msg->orientation.w, msg->orientation.x, msg->orientation.y, msg->orientation.z);

    // update ukf measurements
    UpdateMeasurements();

    // ukf update
    float timediff = imu_.time.sec - prev_imu_.time.sec;
    std::cout << timediff << std::endl;
    ukf_iterate(timediff);

    // print estimation result
    ukf_state_t curr_est_result;
    ukf_get_state(&curr_est_result);
    Eigen::Quaterniond curr_attitude(curr_est_result.attitude[3], curr_est_result.attitude[0],
                                     curr_est_result.attitude[1], curr_est_result.attitude[2]);
    std::cout << curr_attitude.toRotationMatrix().eulerAngles(0, 1, 2).transpose() << std::endl;
  }

  rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_subscription_;
};  // end class sfwa_ukf

// Create std pointers for this class
typedef std::shared_ptr<sfwa_ukf> sfwa_ukfPtr;
typedef std::shared_ptr<const sfwa_ukf> sfwa_ukfConstPtr;

}  // namespace sfwa_ukf
