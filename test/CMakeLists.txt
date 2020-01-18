file(GLOB math_test_srcs "*.cpp")

enable_testing()

add_executable(test.math ${math_test_srcs})
target_link_libraries(test.math common_3rd_lib)
target_link_libraries(test.math pthread rt)
