########################
# PROJECT CONFIGURATION
########################
message("Set up project name and cmake version")
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
project(SchwingerModel)

########################
# MAIN
########################
message("Add main executable")
add_executable(main main.cpp)
# add_compile_options(-Wall -Wextra)

########################
# TESTS
########################

# UNIT TESTS
add_executable(unit_test    tests/unit_test.cpp)
add_test(NAME Test          COMMAND unit_test)


enable_testing()
message("Done!")