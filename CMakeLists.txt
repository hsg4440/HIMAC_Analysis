cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(CSAnalysis)

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/source)
include_directories(${SOURCE_DIR})

file(GLOB SOURCES "${SOURCE_DIR}/*.cc")
file(GLOB HEADERS "${SOURCE_DIR}/*.hh")

add_executable(CSAnalysis CS_analysis.C ${SOURCES} ${HEADERS})
target_link_libraries(CSAnalysis ${ROOT_LIBRARIES})

# add_executable(stoppedCSAnalysis stoppedCS_analysis.C ${SOURCES} ${HEADERS})
# target_link_libraries(stoppedCSAnalysis ${ROOT_LIBRARIES})

