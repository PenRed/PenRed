add_executable(penred_measurements_tests
  test_main.cpp
  test_grids.cpp
  test_results.cpp
  test_measurement.cpp
  test_profiles.cpp
  test_io.cpp
)
target_include_directories(penred_measurements_tests PRIVATE ${PROJECT_SOURCE_DIR}/kernel/includes)
target_compile_features(penred_measurements_tests PRIVATE cxx_std_14)

install(
  TARGETS penred_measurements_tests
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/measurement
  )
