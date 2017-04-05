# - Publish headers

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/src/libs/libfecom)
include_directories(${PROJECT_SOURCE_DIR}/src/libs/libfecom)

set(FeCom_PROGRAMS
  main_decoder_serializer.cxx
  main_histogram_builder.cxx
  )

find_package(Bayeux 3.0.0
  REQUIRED
  )

find_package(Boost 1.60.0
  REQUIRED
  date_time
  )

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

option(FECOM_ENABLE_PROGRAMS "Build unit testing system for libfecom" ON)

if(FECOM_ENABLE_PROGRAMS)
  foreach(_programsource ${FeCom_PROGRAMS})
    get_filename_component(_programname "${_programsource}" NAME_WE)
    set(_programname "fecom-${_programname}")
    add_executable(${_programname}
      ${_programsource}
      )
    target_link_libraries(${_programname} fecom
      # ${FeCom_OTHER_LIBS}
      )

    # - Ensure link to internal and external deps
    set_target_properties(${_programname} PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

    if(UNIX AND NOT APPLE)
      set_target_properties(${_programname}
        PROPERTIES INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
        )
    elseif(APPLE)
      # Temporary setting - needs testing
      set_target_properties(${_programname}
        PROPERTIES
        INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}"
        )
    endif()

    install(TARGETS ${_programname}
      # EXPORT FalaiseTargets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
  endforeach()
endif()