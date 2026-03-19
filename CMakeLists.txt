cmake_minimum_required(VERSION 4.0)

project(pic-engine
    VERSION 0.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

add_executable(pic-sim
    src/main.cpp
)

target_link_libraries(pic-sim
    PRIVATE
        core
        solvers
        domain_objects
)

add_subdirectory(src)

