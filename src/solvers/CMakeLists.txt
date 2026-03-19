add_library(solvers
    potential_solver.cpp
)

target_include_directories(solvers
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(solvers
    PUBLIC
        core
)  