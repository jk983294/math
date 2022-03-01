file(GLOB CommonSrc "*.cpp")
add_library(ornate SHARED ${CommonSrc})
target_link_libraries( ornate blas pthread )
set_property(TARGET ornate PROPERTY POSITION_INDEPENDENT_CODE ON)
