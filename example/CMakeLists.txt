# test 1
add_executable(test_stvk_tet test_stvk_tet.cc)
target_link_libraries(test_stvk_tet
    elastic
)

add_executable(test_simulator test_simulator.cc)
target_link_libraries(test_simulator
    elastic
)

add_executable(test_sofa test_sofa.cc)
target_link_libraries(test_sofa
    elastic
)

add_executable(test_multiobjects test_multiobjects.cc)
target_link_libraries(test_multiobjects
    elastic
)

add_executable(test_arpack test_arpack.cc)
target_link_libraries(test_arpack
    ${ARPACK_LIBRARIES}
)

add_executable(test_reduced test_reduced.cc)
target_link_libraries(test_reduced
    elastic
    ${Boost_LIBRARIES}
)

add_executable(test_eigen_solver test_eigen_solver.cc)
target_link_libraries(test_eigen_solver
    elastic
)

add_executable(test_unit test_unit.cc)
target_link_libraries(test_unit
    elastic
    ${Boost_LIBRARIES}
)

add_executable(test_subspace test_subspace.cc)
target_link_libraries(test_subspace
    elastic
    ${Boost_LIBRARIES}
)

add_executable(test_reduced_stvk test_reduced_stvk.cc)
target_link_libraries(test_reduced_stvk
    elastic
    ${Boost_LIBRARIES}
)

if(OPENCL_FOUND)
    add_executable(test_cl test_cl.cc)
    target_link_libraries(test_cl
    ${OPENCL_LIBRARIES}
    )
endif(OPENCL_FOUND)

add_executable(test_stvk_full test_stvk_full.cc)
target_link_libraries(test_stvk_full
    elastic
    ${Boost_LIBRARIES}
)
