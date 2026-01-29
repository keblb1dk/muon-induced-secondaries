cmake_minimum_required(VERSION 3.5)
project(MuSecondaries)

find_package(Geant4 REQUIRED)
#find_package(Geant4 REQUIRED COMPONENTS analysis)

include(${Geant4_USE_FILE})

add_executable(MuSecondaries
	main.cc 
	src/DetectorConstruction.cc
	src/PrimaryGeneratorAction.cc
	src/RunAction.cc
	src/SteppingAction.cc
	src/EventAction.cc
	src/PrimaryGeneratorMessenger.cc
	src/DetectorMessenger.cc)

target_include_directories(MuSecondaries PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(MuSecondaries ${Geant4_LIBRARIES})

install(TARGETS MuSecondaries DESTINATION bin)
