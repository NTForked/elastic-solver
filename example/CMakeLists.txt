# test 1
add_executable(test_stvk_tet test_stvk_tet.cc)
target_link_libraries(test_stvk_tet
stvk
${LAPACK_LIBRARIES}
)
