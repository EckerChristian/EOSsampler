cmake_minimum_required(VERSION 3.16)

project(eos_inference LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_WARNINGS "Enable common compiler warnings" ON)

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)

set(COMMON_SOURCES
    src/config.cpp
    src/hdf5_io.cpp
)

add_executable(EOS
    src/main_eos.cpp
    src/eos_sampler.cpp
    ${COMMON_SOURCES}
)

add_executable(TOV
    src/main_tov.cpp
    src/tov_solver.cpp
    src/likelihoods.cpp
    ${COMMON_SOURCES}
)

foreach(target EOS TOV)
    target_include_directories(${target}
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${HDF5_INCLUDE_DIRS}
    )

    target_link_libraries(${target}
        PRIVATE
            MPI::MPI_CXX
            ${HDF5_LIBRARIES}
    )

    target_compile_options(${target}
        PRIVATE
            $<$<CONFIG:Release>:-O3>
    )

    if(ENABLE_WARNINGS)
        target_compile_options(${target}
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
        )
    endif()
endforeach()
