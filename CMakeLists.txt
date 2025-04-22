cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(Simulation)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

#find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)

file(GLOB MACRO_FILES "*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

file(GLOB DATA_FILES "*.dat")
file(COPY ${DATA_FILES} DESTINATION ${PROJECT_BINARY_DIR})

add_executable(sim sim.cc ${sources} ${headers})


# Link to Geant4 and ROOT
target_link_libraries(sim ${Geant4_LIBRARIES})
target_link_libraries(sim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

add_custom_target(Simulation DEPENDS sim)
