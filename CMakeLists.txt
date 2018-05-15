cmake_minimum_required(VERSION 2.8)
project(quantum-fdtd)

if (APPLE)
    cmake_policy(SET CMP0042 NEW)
endif ()

# --- Metadata ---
set(qfdtd_VERSION_MAJOR 0)
set(qfdtd_VERSION_MINOR 1)
set(PACKAGE "quantum-fdtd")
set(PACKAGE_BUGREPORT "true.merrill@gtri.gatech.edu")
set(PACKAGE_NAME ${PACKAGE})
set(PACKAGE_VERSION ${iontrap_VERSION_MAJOR}.${iontrap_VERSION_MINOR})
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Options ---
option(BUILD_SHARED_LIBS "Build a shared library" OFF)
option(VERBOSE "Enable verbose compilation" OFF)
option(ENABLE_TESTING "Build and run the unit tests" OFF)

# --- Targets ---

add_library(qfdtd
    src/core.h
    src/core.cc
    src/wavefunction.h
    src/wavefunction.cc
)

if (${ENABLE_TESTING})
    add_subdirectory(test)
endif ()
