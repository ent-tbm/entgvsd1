cmake_minimum_required(VERSION 3.0)
project (entgvsd
    LANGUAGES Fortran)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include (setup_rpath)

# We need this separately because we use the command ncgen
find_package(NetCDF4_C REQUIRED)
include_directories(${NETCDF4_C_INCLUDE_DIR})

find_package(NetCDF4_Fortran REQUIRED)
include_directories(${NETCDF4_FORTRAN_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBS ${NETCDF4_FORTRAN_LIBRARY})

add_subdirectory(slib)
