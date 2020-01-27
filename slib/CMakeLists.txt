include_directories(
    ${PROJECT_SOURCE_DIR}/slib
    ${CMAKE_CURRENT_BINARY_DIR})

set(entgvsd_SOURCES
    arrayutil.f
    ent_params.f90
    chunker.f90
    hntr.f90
    ent_labels.F90
    assign_laimax.F90
    cropmerge_laisparse_splitbare.f90
    conversions.f
    carrer.F90
    gridfill.f90
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd_config.f90
)

set(AT "@")

# ----- Determine lib directories from found libraries
get_filename_component(NETCDF4_FORTRAN_LIB_DIR ${NETCDF4_FORTRAN_LIBRARY} DIRECTORY)
get_filename_component(NETCDF4_C_LIB_DIR ${NETCDF4_C_LIBRARY} DIRECTORY)


# ----- Write CMake-derived paths into place accessible to C++/Fortran
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/entgvsd_config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd_config.hpp)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/entgvsd_config.f90.in
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd_config.f90)


execute_process( COMMAND
    cp --no-clobber 
    ${CMAKE_CURRENT_SOURCE_DIR}/ent_params.f90.in
    ${CMAKE_CURRENT_SOURCE_DIR}/ent_params.f90)

# If @ONLY is specified, only variables of the form @VAR@ will be
# replaced and ${ VAR} will be ignored. This is useful for configuring
# scripts that use ${VAR}.
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/entgvsd.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd.sh @ONLY)



add_library(entgvsd SHARED ${entgvsd_SOURCES})
target_link_libraries(entgvsd ${EXTERNAL_LIBS})
install(TARGETS entgvsd DESTINATION lib)

install(DIRECTORY ${PROJECT_BINARY_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.mod")
install(DIRECTORY ${PROJECT_BINARY_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/entgvsd.sh DESTINATION bin RENAME entgvsd)
install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/xent.py DESTINATION bin RENAME xent)

