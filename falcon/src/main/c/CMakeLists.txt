include(ExternalProject)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

cmake_minimum_required(VERSION 2.8.12)

project(falcon-genomics)

find_package(Boost 1.53.0 COMPONENTS
	  system thread iostreams filesystem REQUIRED)
find_package(JNI REQUIRED)
find_package(Threads)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Blaze REQUIRED)
find_package(Google REQUIRED)
find_package(KSight REQUIRED)
find_package(IntelGKL REQUIRED)
find_package(PairhmmClient REQUIRED)
find_package(FalconLM REQUIRED)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# get version str from git
execute_process(
    COMMAND git describe --tags --always
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -m64 -DBOOST_NO_CXX11_SCOPED_ENUMS -fPIC -DUSELICENSE")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-z,defs")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -DNDEBUG")

if(DEPLOYMENT) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEPLOY_${DEPLOYMENT}")
    set(BIN_VERSION "${GIT_VERSION}-${DEPLOYMENT}")
    message(STATUS "destination deployment platform is ${DEPLOYMENT}")
else()
    set(BIN_VERSION "${GIT_VERSION}")
endif()

if(NO_PROFILE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_PROFILE")
else()
  message(STATUS "Profiling is enabled")
endif()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(BIN_VERSION "${BIN_VERSION}-dev")
endif()

message(STATUS "version: ${BIN_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=\\\"${BIN_VERSION}\\\"")

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${JNI_INCLUDE_DIRS}
    ${Blaze_INCLUDE_DIRS}
    ${Google_INCLUDE_DIRS}
    ${KSight_INCLUDE_DIRS}
    ${PMMClient_INCLUDE_DIRS}
    ${IntelGKL_INCLUDE_DIRS}
    ${FalconLM_INCLUDE_DIRS})

link_directories(
    ${Boost_LIBRARY_DIRS}
    ${Blaze_LIBRARY_DIRS}
    ${Google_LIBRARY_DIRS}
    ${KSight_LIBRARY_DIRS}
    ${PMMClient_LIBRARY_DIRS}
    ${IntelGKL_LIBRARY_DIRS}
    ${FalconLM_LIBRARY_DIRS})

# find source files
file(GLOB_RECURSE SRC_LIST src/*.cpp)

add_library(falcon_genomics SHARED ${SRC_LIST})

add_dependencies(falcon_genomics 
    Blaze 
    FalconLM
    Google 
    PMMClient
    IntelGKL 
    JsonCPP
    KSight)

target_link_libraries(falcon_genomics
    ${Boost_LIBRARIES} 
    ${PMMClient_LIBRARIES}
    ${Blaze_LIBRARIES}
    ${KSight_LIBRARIES}
    ${IntelGKL_LIBRARIES}
    ${FalconLM_LIBRARIES}
    ${Google_LIBRARIES}
    ${CMAKE_DL_LIBS} 
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})

#enable_testing()
#add_subdirectory(test)
