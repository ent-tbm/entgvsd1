include_directories(
    ${PROJECT_SOURCE_DIR}/slib
    ${CMAKE_CURRENT_BINARY_DIR})

set(entgvsd_SOURCES
    arrayutil.f
#    convertnc_util.f90
#    entgvsd_util.f90
    paths.f90
    chunkparams.f90
    chunker.f90
    hntr.f90
    ent_labels.F90
    assign_laimax.F90
    cropmerge_laisparse_splitbare.f90
    conversions.f
    carrer.F90
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
install(PROGRAMS entgvsd_create_nc.sh DESTINATION bin RENAME entgvsd_create_nc)
install(PROGRAMS entgvsd_link_input.py DESTINATION bin RENAME entgvsd_link_input)
install(PROGRAMS cdlparser.py DESTINATION bin)

