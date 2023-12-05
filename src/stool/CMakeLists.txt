find_package(Arrow REQUIRED)

add_executable( daily_y daily_y.cpp )
target_link_libraries( daily_y common_3rd_lib ${OpenMP_CXX_FLAGS}
        ${FST_LIB} ${LZ4_LIB} ${ZSTD_LIB} arrow arrow_dataset pthread rt dl)
target_compile_options(daily_y PRIVATE "${OpenMP_CXX_FLAGS}")
