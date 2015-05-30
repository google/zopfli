cmake_minimum_required(VERSION 2.8)

project(Zopfli)

option(BUILD_SHARED_LIBS "Build Zopfli with shared libraries" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif()
 
set(zopflilib_src
  src/zopfli/blocksplitter.c
  src/zopfli/cache.c
  src/zopfli/deflate.c
  src/zopfli/gzip_container.c
  src/zopfli/hash.c
  src/zopfli/katajainen.c
  src/zopfli/lz77.c
  src/zopfli/squeeze.c
  src/zopfli/tree.c
  src/zopfli/util.c
  src/zopfli/zlib_container.c
  src/zopfli/zopfli_lib.c
)

set(zopflipnglib_src
  src/zopflipng/zopflipng_lib.cc
)

set (lodepng_src
  src/zopflipng/lodepng/lodepng.cpp
  src/zopflipng/lodepng/lodepng_util.cpp
)

#
# libzopfli object files shared by both libraries
#
add_library(zopflilib_obj OBJECT
  ${zopflilib_src}
)
if(BUILD_SHARED_LIBS)
  set_property(TARGET zopflilib_obj PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

#
# libzopfli
#
add_library(libzopfli
  $<TARGET_OBJECTS:zopflilib_obj>
)
set_target_properties(libzopfli PROPERTIES
  OUTPUT_NAME zopfli
  VERSION 1.0.1
  SOVERSION 1
)
if(UNIX)
  target_link_libraries(libzopfli m)
endif()

#
# libzopflipng
#
add_library(libzopflipng
  ${zopflipnglib_src}
  ${lodepng_src}
  $<TARGET_OBJECTS:zopflilib_obj>
)
set_target_properties(libzopflipng PROPERTIES
  OUTPUT_NAME zopflipng
  VERSION 1.0.0
  SOVERSION 1
)

#
# zopfli
#
add_executable(zopfli src/zopfli/zopfli_bin.c)
target_link_libraries(zopfli libzopfli)

#
# zopflipng
#
add_executable(zopflipng src/zopflipng/zopflipng_bin.cc)
target_link_libraries(zopflipng libzopflipng)
