add_library(domain_objects
    world.cpp
    field.cpp
    species.cpp
)

target_include_directories(domain_objects
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(domain_objects
    PUBLIC
        core
        eigen
)