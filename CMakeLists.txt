cmake_minimum_required(VERSION 3.0.2)

project(IK_Simulation)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.55.0 REQUIRED program_options)

include_directories(include ${Boost_INCLUDE_DIRS})

file(GLOB HEADERS "include/*.hpp")
file(GLOB SOURCES "src/*.cpp")

add_executable(IK_Simulation ${SOURCES} ${HEADERS})
target_link_libraries(IK_Simulation ${Boost_LIBRARIES})