cmake_minimum_required(VERSION 2.9)

project(goja)

# CFLAGS:

add_definitions(-std=c++14 -Wall -Wunused-parameter)

# Searching for dependencies:

list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT 6.24 REQUIRED COMPONENTS Tree)
if(ROOT_USE_FILE)
  include(${ROOT_USE_FILE})
endif()
message(STATUS "ROOT version: ${ROOT_VERSION}")
add_definitions(${ROOT_DEFINITIONS})

find_package(Boost 1.71.0 REQUIRED COMPONENTS program_options)
add_definitions(${Boost_DEFINITIONS})

include_directories(${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cpp ${PROJECT_SOURCE_DIR}/*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.h)

# Building:

add_executable(goja goja.cpp ${sources} ${headers})

# Linking:

target_link_libraries(goja ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Installing:

set_target_properties(goja PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS goja DESTINATION bin)
install (FILES ${PROJECT_SOURCE_DIR}/goja_manager.py
               ${PROJECT_SOURCE_DIR}/verify_gate_results.py
               ${PROJECT_SOURCE_DIR}/verify_goja_results.py
               ${PROJECT_SOURCE_DIR}/../examples/python_scripts/nema_common.py
               ${PROJECT_SOURCE_DIR}/../examples/python_scripts/analyze_coincidences.py
         DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                                     GROUP_READ GROUP_EXECUTE
                                     WORLD_READ WORLD_EXECUTE)
