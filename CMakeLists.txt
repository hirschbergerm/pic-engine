cmake_minimum_required(VERSION 4.0)

project(pic-engine)

add_executable(pic src/main.cpp)

target_include_directories(pic PRIVATE ${PROJECT_SOURCE_DIR}/external/eigen-5.0.0)

