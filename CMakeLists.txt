cmake_minimum_required(VERSION 3.0)
project (entgvsd
    LANGUAGES Fortran)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include (setup_rpath)

find_package(NetCDF4_Fortran REQUIRED)
include_directories(${NETCDF_FORTRAN_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBS ${NETCDF4_FORTRAN_LIBRARY})

find_package(Ncl REQUIRED)
include_directories(${PATH_INC_NCL})
list(APPEND EXTERNAL_LIBS ${NCL_LIB_FILE})

add_subdirectory(slib)
