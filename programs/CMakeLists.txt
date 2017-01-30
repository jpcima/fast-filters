cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

add_subdirectory(.. fast-filters.dir)

add_executable(fir-benchmark fir-benchmark.cc)
target_link_libraries(fir-benchmark fast-filters)

add_executable(iir-benchmark iir-benchmark.cc)
target_link_libraries(iir-benchmark fast-filters)
