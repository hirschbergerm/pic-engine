add_subdirectory(core)
add_subdirectory(domain_objects)
add_subdirectory(solvers)
add_subdirectory(output)

add_executable(pic-sim
    main.cpp
)

target_include_directories(pic-sim
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(eigen INTERFACE)

target_include_directories(eigen
    INTERFACE ${PROJECT_SOURCE_DIR}/external/eigen
)

target_link_libraries(pic-sim
    PRIVATE
        core
        eigen
        domain_objects
        solvers
        output
)
