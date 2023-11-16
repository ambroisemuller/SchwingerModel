########################
# PROJECT CONFIGURATION
########################
message("Set up project name and cmake version")
set(CMAKE_CXX_STANDARD 17)
project(SchwingerModel)

########################
# MAIN
########################
message("Add main executable")
add_executable(main main.cpp)

########################
# TESTS
########################

# UNIT TESTS
add_executable(test         tests/test.cpp)
add_test(NAME Test          COMMAND test)


enable_testing()
message("Done!")