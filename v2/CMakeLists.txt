########################
# PROJECT CONFIGURATION
########################
message("Set up project name and cmake version")
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
project(SchwingerModel)

########################
# FLAGS
########################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

########################
# MAIN
########################
message("Add main executable")
find_package(Threads REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main Threads::Threads)
# add_compile_options(-Wall -Wextra)

########################
# TESTS
########################

# UNIT TESTS
add_executable(unit_test    tests/unit_test.cpp)
target_link_libraries(unit_test Threads::Threads)
add_test(NAME Test          COMMAND unit_test)


enable_testing()
message("Done!")
