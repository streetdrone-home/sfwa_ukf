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

/* Author: Sina Aghli
   Desc: TODO(GITHUB_NAME):
*/
#include <sfwa_ukf/sfwa_ukf.hpp>

namespace sfwa_ukf
{
sfwa_ukf::sfwa_ukf() : Node("sfwa_ukf")
{
  RCLCPP_INFO(get_logger(), "Initializing sfwa_ukf...");
  imu_subscription_ = this->create_subscription<sensor_msgs::msg::Imu>(
      "/imu_xsens_mti_ros", 1, std::bind(&sfwa_ukf::ImuCallback, this, std::placeholders::_1));

  InitFilter();
}

void sfwa_ukf::InitFilter()
{
  std::cout << "Initializing UKF!" << std::endl;
  RCLCPP_INFO(this->get_logger(), "Initializing UKF!");

  ukf_init();

  // set initial state
  ukf_state_t init_state;
  init_state.attitude[0] = imu_.orientation.x();
  init_state.attitude[1] = imu_.orientation.y();
  init_state.attitude[2] = imu_.orientation.z();
  init_state.attitude[3] = imu_.orientation.w();
  init_state.acceleration[0] = 0;
  init_state.acceleration[1] = 0;
  init_state.acceleration[2] = 0;
  init_state.angular_velocity[0] = 0;
  init_state.angular_velocity[1] = 0;
  init_state.angular_velocity[2] = 0;

  ukf_set_state(&init_state);
}

void sfwa_ukf::UpdateMeasurements()
{
  ukf_sensor_set_gyroscope(imu_.angvel[0], imu_.angvel[1], imu_.angvel[2]);
  ukf_sensor_set_accelerometer(imu_.acc[0], imu_.acc[1], imu_.acc[2]);
}

}  // end namespace sfwa_ukf
