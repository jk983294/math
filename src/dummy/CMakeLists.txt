file(GLOB mathDummySrc "*.cpp")
add_library(mathdummystatic STATIC ${mathDummySrc})
install(TARGETS mathdummystatic ARCHIVE DESTINATION lib)

add_library(mathdummydyn SHARED ${mathDummySrc})
install(TARGETS mathdummydyn LIBRARY DESTINATION lib)
