include_directories(
    ${PROJECT_SOURCE_DIR}/slib
    ${CMAKE_CURRENT_BINARY_DIR})

set(entgvsd_SOURCES
    arrayutil.f
    convertnc_util.f90
    entgvsd_util.f90
    paths_mod.f90
    chunkparams.f90
    #chunker_mod.f90
    #hntr_mod.f90
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd_config_mod.f90
)

# ----- Write CMake-derived paths into place accessible to C++/Fortran
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/entgvsd_config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd_config.hpp)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/entgvsd_config_mod.f90.in
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd_config_mod.f90)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/entgvsd.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/entgvsd.sh)



add_library(entgvsd SHARED ${entgvsd_SOURCES})
target_link_libraries(entgvsd ${EXTERNAL_LIBS})
install(TARGETS entgvsd DESTINATION lib)

install(DIRECTORY ${PROJECT_BINARY_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.mod")
install(DIRECTORY ${PROJECT_BINARY_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/entgvsd.sh DESTINATION bin RENAME entgvsd)
install(PROGRAMS entgvsd_create_nc.sh DESTINATION bin RENAME entgvsd_create_nc)
install(PROGRAMS entgvsd_link_input.py DESTINATION bin RENAME entgvsd_link_input)
install(PROGRAMS cdlparser.py DESTINATION bin)

