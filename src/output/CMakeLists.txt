add_library(output
    output.cpp
)

target_include_directories(output
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(output
    PUBLIC 
        core
        domain_objects
)