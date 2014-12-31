# test 1
add_executable(test_stvk_tet test_stvk_tet.cc)
target_link_libraries(test_stvk_tet
elastic
hj-util
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
)

add_executable(test_simulator test_simulator.cc)
target_link_libraries(test_simulator
elastic
jtf-mesh
hj-util
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
)

add_executable(test_sofa test_sofa.cc)
target_link_libraries(test_sofa
elastic
jtf-mesh
hj-util
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
)

add_executable(test_multiobjects test_multiobjects.cc)
target_link_libraries(test_multiobjects
elastic
jtf-mesh
hj-util
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
)

add_executable(test_arpack test_arpack.cc)
target_link_libraries(test_arpack
${ARPACK_LIBRARIES}
)

add_executable(test_reduced test_reduced.cc)
target_link_libraries(test_reduced
elastic
jtf-mesh
hj-util
${UMFPACK_LIBRARIES}
${Boost_LIBRARIES}
)

add_executable(test_eigen_solver test_eigen_solver.cc)
target_link_libraries(test_eigen_solver
elastic
)

add_executable(test_unit test_unit.cc)
target_link_libraries(test_unit
elastic
hj-util
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
${Boost_LIBRARIES}
)

add_executable(test_subspace test_subspace.cc)
target_link_libraries(test_subspace
elastic
jtf-mesh
hj-util
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
${Boost_LIBRARIES}
)

