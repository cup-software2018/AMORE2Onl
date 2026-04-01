# ------------------------------------------------------------
# Tests for AMOREHDF5 module
# ------------------------------------------------------------
set ( TEST_SOURCES
    test_amore_event.cc
    test_amore_hit.cc

    # test_amore_subrun.cc
    # amore_edge_test.cc
)

foreach ( src ${TEST_SOURCES} )
    get_filename_component ( test_name ${src} NAME_WE )

    add_executable ( ${test_name} ${src} )

    # Linking HDF5Utils automatically provides its PUBLIC include directories
    # and transative dependencies (like ROOT and HDF5).
    # Therefore, explicit target_include_directories and repetitive linking are removed.
    target_link_libraries ( ${test_name}
        PRIVATE
            HDF5Utils
            AMOREHDF5
    )

    install ( TARGETS ${test_name}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endforeach ()