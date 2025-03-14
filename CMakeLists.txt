cmake_minimum_required(VERSION 3.16)
project(sfwa_ukf CXX)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wfatal-errors)

  # -Wsign-conversion -Wall -Werror -Wextra -Wpedantic -Wshadow -Winit-self -Wredundant-decls
endif()

set(THIS_PACKAGE_INCLUDE_DEPENDS
  rclcpp
  eigen3_cmake_module
  Eigen3
  sensor_msgs
)

# Define floating point precision
add_definitions(-DUKF_DOUBLE_PRECISION)

find_package(ament_cmake REQUIRED)

foreach(Dependency IN ITEMS ${THIS_PACKAGE_INCLUDE_DEPENDS})
  find_package(${Dependency} REQUIRED)
endforeach()

# ##########
# # Build ##
# ##########

# Declare a C++ library
add_library(
  ${PROJECT_NAME}
  src/sfwa_ukf.cpp
  src/ahrs.cpp
)
ament_target_dependencies(${PROJECT_NAME} ${THIS_PACKAGE_INCLUDE_DEPENDS})
target_include_directories(${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

# Declare a C++ executable
add_executable(
  ${PROJECT_NAME}_node
  src/main.cpp
)
target_link_libraries(${PROJECT_NAME}_node ${PROJECT_NAME})
ament_target_dependencies(${PROJECT_NAME}_node ${THIS_PACKAGE_INCLUDE_DEPENDS})

# ############
# # Install ##
# ############
install(
  DIRECTORY include/
  DESTINATION include
)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib/${PROJECT_NAME}
)

install(TARGETS ${PROJECT_NAME}_node
  DESTINATION lib/${PROJECT_NAME})

# ############
# # Testing ##
# ############
if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)
  find_package(ament_lint_auto REQUIRED)
  find_package(ament_cmake_pytest REQUIRED)

  # ROS2 linters, but disable copyright test. PickNik's copyright's may not conform
  # to this test
  # set(ament_cmake_copyright_FOUND TRUE)
  # ament_lint_auto_find_test_dependencies()

  # GTests
  ament_add_gtest(sfwa_ukf_test
    test/TestStateVector.cpp
    test/TestFixedMeasurementVector.cpp
    test/TestDynamicMeasurementVector.cpp
    test/TestCore.cpp
    test/TestSquareRootCore.cpp)
  target_link_libraries(sfwa_ukf_test
    ${PROJECT_NAME})
  ament_target_dependencies(sfwa_ukf_test ${THIS_PACKAGE_INCLUDE_DEPENDS})

  # Benchmarks
  find_package(ament_cmake_google_benchmark REQUIRED)
  ament_add_google_benchmark(sfwa_ukf_benchmark
    benchmark/BenchmarkMain.cpp
    benchmark/CoreBenchmark.cpp
    benchmark/MeasurementVectorBenchmark.cpp
    benchmark/SquareRootCoreBenchmark.cpp
    benchmark/StateVectorBenchmark.cpp)
  target_link_libraries(sfwa_ukf_benchmark
    ${PROJECT_NAME})
  ament_target_dependencies(sfwa_ukf_benchmark ${THIS_PACKAGE_INCLUDE_DEPENDS})
endif()

ament_export_targets(${PROJECT_NAME}Targets HAS_LIBRARY_TARGET)
ament_export_dependencies(${THIS_PACKAGE_INCLUDE_DEPENDS})
ament_package()
