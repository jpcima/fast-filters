cmake_minimum_required(VERSION 3.0)

project(fast-filters
  VERSION 0.1
  LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(WarningsEnabled)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

add_subdirectory("${PROJECT_SOURCE_DIR}/../coreutil" coreutil.dir)

set(fast-filters_SOURCES
  sources/coredsp/filter.h
  sources/coredsp/noise.h)
add_library(fast-filters STATIC ${fast-filters_SOURCES})
target_include_directories(fast-filters PUBLIC sources)
target_link_libraries(fast-filters coreutil)
