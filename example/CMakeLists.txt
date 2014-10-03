# test 1
add_executable(test_stvk_tet test_stvk_tet.cc)
target_link_libraries(test_stvk_tet
elastic
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
)

add_executable(test_simulator test_simulator.cc)
target_link_libraries(test_simulator
elastic
jtf-mesh
${UMFPACK_LIBRARIES}
${LAPACK_LIBRARIES}
)
