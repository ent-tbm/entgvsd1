include_directories(${PROJECT_SOURCE_DIR}/slib)

set(entgvsd_SOURCES
    arrayutil.f
)

add_library(entgvsd SHARED ${entgvsd_SOURCES})
target_link_libraries(entgvsd ${EXTERNAL_LIBS})
install(TARGETS entgvsd DESTINATION lib)

install(DIRECTORY ${PROJECT_BINARY_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.mod")
