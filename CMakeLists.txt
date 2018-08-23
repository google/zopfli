cmake_minimum_required(VERSION 2.8.11)

project(Zopfli)

#
# Library version
#
set(ZOPFLI_VERSION_MAJOR 1)
set(ZOPFLI_VERSION_MINOR 0)
set(ZOPFLI_VERSION_PATCH 2)
set(ZOPFLI_VERSION "${ZOPFLI_VERSION_MAJOR}.${ZOPFLI_VERSION_MINOR}.${ZOPFLI_VERSION_PATCH}")

option(BUILD_SHARED_LIBS "Build Zopfli with shared libraries" OFF)
option(ZOPFLI_BUILD_INSTALL "Add Zopfli install target" ON)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Zopfli build type is default (CMAKE_BUILD_TYPE empty)")
  endif()
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
target_include_directories(libzopfli INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src/zopfli")
set_target_properties(libzopfli PROPERTIES
  OUTPUT_NAME zopfli
  VERSION ${ZOPFLI_VERSION}
  SOVERSION ${ZOPFLI_VERSION_MAJOR}
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
target_include_directories(libzopflipng INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src/zopflipng")
set_target_properties(libzopflipng PROPERTIES
  OUTPUT_NAME zopflipng
  VERSION ${ZOPFLI_VERSION}
  SOVERSION ${ZOPFLI_VERSION_MAJOR}
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

#
# Install
#
if(ZOPFLI_BUILD_INSTALL)
  include(GNUInstallDirs)
  install(TARGETS libzopfli libzopflipng zopfli zopflipng
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
  install(FILES include/zopfli.h include/zopflipng_lib.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
endif()
