cmake_minimum_required(VERSION 2.6)

project(smear-events)

# CFLAGS:

add_definitions(-std=c++11 -Wall -Wunused-parameter)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/*.h)

# Building:

add_executable(SmearEvents ${sources} ${headers})

# Installing:

set_target_properties(SmearEvents PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS SmearEvents DESTINATION bin)
