cmake_minimum_required(VERSION 3.14)
project(btools)

set(CMAKE_CXX_STANDARD 11)

add_executable(btools btools_main.cpp btools_intersect.cpp btools.cpp tools.cpp jobs.cpp btools_annotate2.cpp)
target_link_libraries(btools pthread)
