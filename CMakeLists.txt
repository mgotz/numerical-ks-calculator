cmake_minimum_required(VERSION 3.0.2)

project(ks_calculator)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.55.0 REQUIRED program_options)

#offer some cmake options, that will then be passed to the program compilation
#via header inclusion
set(LOGLEVEL 0 CACHE STRING "determine the amount of output 0 is minimal, 1 is status every 1000 steps and 2 logs information every step")

option(CONSTFIELD "help string describing option" OFF)
if(CONSTFIELD)
	set(CONSTFIELD 1)
else()
	set(CONSTFIELD 0)
endif()


option(CHARGETIMELINE "output files containing E-Field, electron and ion distributions, every 1000 steps" OFF)
if(CHARGETIMELINE)
        set(CHARGETIMELINE 1)
else()
        set(CHARGETIMELINE 0)
endif()


option(ELECTRONCOMBO "OFF saves computations compared to just setting electron recombination coefficient=0" ON)
if(ELECTRONCOMBO)
        set(ELECTRONCOMBO 1)
else()
        set(ELECTRONCOMBO 0)
endif()


option(CONSTIONMOB "set to on if ion mobility does not depend on E-field: saves computation time" OFF)
if(CONSTIONMOB)
        set(CONSTIONMOB 1)
else()
        set(CONSTIONMOB 0)
endif()


option(CONSTRECOMB "set to on if ion recombination does not depend on E-field: saves computation time" OFF)
if(CONSTRECOMB)
        set(CONSTRECOMB 1)
else()
        set(CONSTRECOMB 0)
endif()


option(CONSTELECTRONCOMBO "set to on if electron-ion recombination does not depend on E-field: saves computation time" OFF)
if(CONSTELECTRONCOMBO)
        set(CONSTELECTRONCOMBO 1)
else()
        set(CONSTELECTRONCOMBO 0)
endif()


# configure a header file to pass those CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/config.hpp.in"
  "${PROJECT_BINARY_DIR}/config.hpp"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")


include_directories(include ${Boost_INCLUDE_DIRS})

file(GLOB HEADERS "include/*.hpp")
file(GLOB SOURCES "src/*.cpp")

add_executable(ks_calculator ${SOURCES} ${HEADERS})
target_link_libraries(ks_calculator ${Boost_LIBRARIES})
