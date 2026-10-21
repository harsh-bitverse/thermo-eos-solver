cmake_minimum_required(VERSION 3.20)
project(ThermoSoftware VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warning flags
if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Core thermodynamic library
add_library(thermo_core
    src/Component.cpp
    src/ComponentDatabase.cpp
    src/CubicSolver.cpp
    src/PengRobinsonEOS.cpp
)

target_include_directories(thermo_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# CLI application executable
add_executable(thermo_cli
    src/cli_main.cpp
    src/cli/CLIController.cpp
)
target_include_directories(thermo_cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(thermo_cli PRIVATE thermo_core)

# Automated testing
enable_testing()
add_subdirectory(tests)
