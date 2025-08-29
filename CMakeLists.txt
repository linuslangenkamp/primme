cmake_minimum_required(VERSION 3.15)
project(PRIMME LANGUAGES C)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared library")

if(NOT CMAKE_C_FLAGS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC")
endif()

file(GLOB_RECURSE PRIMME_SOURCES
    "src/*.c"
    "src/*/*.c"
)

file(GLOB_RECURSE PRIMME_HEADERS
    "include/*.h"
    "include/*.inc"
)

add_library(primme ${PRIMME_SOURCES} ${PRIMME_HEADERS})

target_include_directories(primme PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

find_package(LAPACK)
if(LAPACK_FOUND)
    target_link_libraries(primme PUBLIC ${LAPACK_LIBRARIES})
endif()

set_target_properties(primme PROPERTIES
    OUTPUT_NAME "primme"
    VERSION 1.1.0
    SOVERSION 1
)

install(TARGETS primme
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
)

install(DIRECTORY include/ DESTINATION include)
