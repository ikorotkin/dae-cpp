cmake_minimum_required(VERSION 3.0)

project(dae-cpp)

enable_language(CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(DAE_FORTRAN_STYLE "Matrices will be defined using FORTRAN style (one-based indexing of columns and rows)" OFF)
option(DAE_LONG_INT "Use long integer representation for huge systems (more than ~10^7 equations)" OFF)
option(DAE_SINGLE "Use single precision in the solver)" OFF)
option(DAE_TEST "Build automatic solver test" ON)

if(DAE_TEST)
  unset(DAE_BUILD_EXAMPLES CACHE)
endif(DAE_TEST)

option(DAE_BUILD_EXAMPLES "Build all examples" ON)


#### Locate Intel MKL ####

if(DEFINED ENV{MKLROOT})
  message(STATUS "Found MKLROOT environment variable...")
  set(MKL_ROOT_DIR $ENV{MKLROOT})
else(DEFINED ENV{MKLROOT})
  message(STATUS "MKLROOT environment variable not found. Will use default MKL paths...")
  if(WIN32)
    set(MKL_ROOT_DIR "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl")
  else(WIN32)
    set(MKL_ROOT_DIR "/opt/intel/mkl")
  endif(WIN32)
endif(DEFINED ENV{MKLROOT})

set(DAE_MKL_DIR "${MKL_ROOT_DIR}" CACHE STRING "Path to Intel MKL root directory")

find_path(DAE_MKL_INCLUDE NAMES mkl.h HINTS ${MKL_ROOT_DIR}/include)

if(WIN32)
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    # 64 bits
    set(DAE_MKL_LIB "${DAE_MKL_DIR}/lib/intel64_win")
    set(DAE_INTEL_LIB "${DAE_MKL_DIR}/../compiler/lib/intel64_win")
  elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # 32 bits
    set(DAE_MKL_LIB "${DAE_MKL_DIR}/lib/ia32_win")
    set(DAE_INTEL_LIB "${DAE_MKL_DIR}/../compiler/lib/ia32_win")
  endif()
else(WIN32)
  # MKL library
  if(EXISTS ${DAE_MKL_DIR}/lib/intel64)
    set(DAE_MKL_LIB "${DAE_MKL_DIR}/lib/intel64")
  else(EXISTS ${DAE_MKL_DIR}/lib/intel64)
    set(DAE_MKL_LIB "${DAE_MKL_DIR}/lib")
  endif(EXISTS ${DAE_MKL_DIR}/lib/intel64)
  # Intel OMP library
  find_path(DAE_INTEL_LIB NAMES libiomp5.so libiomp5.a PATHS
    ${DAE_MKL_DIR}/../lib/intel64
    ${DAE_MKL_DIR}/../compiler/lib/intel64
    ${DAE_MKL_DIR}/../lib
    DOC "Path to libiomp5 (Intel OMP library)")
endif(WIN32)


#### Locate Python if needed ####

option(PLOTTING "Use Python interface (matplotlib) for plotting" OFF)

if(PLOTTING)
  find_package(PythonLibs)
  if(PYTHONLIBS_FOUND)
    set(PYTHON_INCLUDE ${PYTHON_INCLUDE_DIRS} CACHE STRING "Path to Python include file (Python.h) for plotting")
    set(PYTHON_LIB ${PYTHON_LIBRARIES} CACHE STRING "Python library for plotting")
  else(PYTHONLIBS_FOUND)
    set(PYTHON_INCLUDE "" CACHE STRING "Path to Python include file (Python.h) for plotting")
    set(PYTHON_LIB "" CACHE STRING "Python library for plotting")
  endif(PYTHONLIBS_FOUND)
  # A very rude attempt to find numpy library (numpy/arrayobject.h) for plotting
  find_path(PYTHON_NUMPY_INCLUDE numpy/arrayobject.h PATHS
    /usr/local/lib/python3.5/dist-packages/numpy/core/include
    /usr/local/lib/python3.6/dist-packages/numpy/core/include
    /usr/local/lib/python3.7/dist-packages/numpy/core/include
    /usr/local/lib/python3.8/dist-packages/numpy/core/include
    /usr/local/lib/python3.5/site-packages/numpy/core/include
    /usr/local/lib/python3.6/site-packages/numpy/core/include
    /usr/local/lib/python3.7/site-packages/numpy/core/include
    /usr/local/lib/python3.8/site-packages/numpy/core/include
    ~/.local/lib/python3.5/site-packages/numpy/core/include
    ~/.local/lib/python3.6/site-packages/numpy/core/include
    ~/.local/lib/python3.7/site-packages/numpy/core/include
    ~/.local/lib/python3.8/site-packages/numpy/core/include
    ${PYTHON_INCLUDE}/../Lib/site-packages/numpy/core/include
    ${PYTHON_LIB}/../Lib/site-packages/numpy/core/include
    DOC "Path to Python numpy include file (numpy/arrayobject.h) for plotting")
  include_directories(${PYTHON_INCLUDE} ${PYTHON_NUMPY_INCLUDE})
endif(PLOTTING)


#### Set compiler options ####

unset(CMAKE_CXX_FLAGS CACHE)

if(UNIX AND NOT APPLE)
  message(STATUS "Unix (not Apple) OS detected...")
  if(DAE_LONG_INT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DMKL_ILP64 -m64 -fopenmp -Wall -Wl,--no-as-needed")
  else(DAE_LONG_INT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -m64 -fopenmp -Wall -Wl,--no-as-needed")
  endif(DAE_LONG_INT)

elseif(APPLE)
  message(STATUS "Apple OS detected...")
  if(DAE_LONG_INT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DMKL_ILP64 -m64 -fopenmp -Wall -Wl,-rpath,${DAE_MKL_LIB}")
  else(DAE_LONG_INT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -m64 -fopenmp -Wall -Wl,-rpath,${DAE_MKL_LIB}")
  endif(DAE_LONG_INT)

elseif(WIN32)
  message(STATUS "Windows OS detected...")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp /EHsc")

else(UNIX AND NOT APPLE)
  message(STATUS "Unknown OS detected... Cannot set up compiler options.")

endif(UNIX AND NOT APPLE)


#### Add include and library paths ####

configure_file(${PROJECT_SOURCE_DIR}/src/cmake_config.h.in ${PROJECT_SOURCE_DIR}/src/cmake_config.h)

set(EXT_INCLUDE "${PROJECT_SOURCE_DIR}/src/external")

include_directories(${DAE_MKL_INCLUDE} ${EXT_INCLUDE})
include_directories("${PROJECT_SOURCE_DIR}/src")

link_directories(${DAE_MKL_LIB})
link_directories(${DAE_INTEL_LIB})

add_subdirectory(src)


#### Build examples if needed ####

if(DAE_BUILD_EXAMPLES)

  set(EXAMPLE_LIST "perovskite" "diffusion_2d" "robertson" "simple_dae" "two_bodies")

  if(WIN32)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      # 64 bits
      find_library(LIB_MKL_INTEL_LP64 mkl_intel_lp64.lib HINTS ${DAE_MKL_LIB})
    elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
      # 32 bits
      find_library(LIB_MKL_INTEL_LP64 mkl_intel_c.lib HINTS ${DAE_MKL_LIB})
    endif()
    find_library(LIB_MKL_THREAD mkl_intel_thread.lib HINTS ${DAE_MKL_LIB})
    find_library(LIB_MKL_CORE mkl_core.lib HINTS ${DAE_MKL_LIB})
    find_library(LIB_INTEL_OPENMP libiomp5md.lib HINTS ${DAE_INTEL_LIB})
    set(LIB_PTHREAD "")
    set(LIB_MATH "")
    set(LIB_DL "")

  else(WIN32)
    if(DAE_LONG_INT)
      unset(LIB_MKL_INTEL_LP64 CACHE)
      find_library(LIB_MKL_INTEL_LP64 mkl_intel_ilp64 HINTS ${DAE_MKL_LIB})
    else(DAE_LONG_INT)
      find_library(LIB_MKL_INTEL_LP64 mkl_intel_lp64 HINTS ${DAE_MKL_LIB})
    endif(DAE_LONG_INT)
    find_library(LIB_MKL_THREAD mkl_intel_thread HINTS ${DAE_MKL_LIB})
    find_library(LIB_MKL_CORE mkl_core HINTS ${DAE_MKL_LIB})
    find_library(LIB_INTEL_OPENMP iomp5 HINTS ${DAE_INTEL_LIB})
    find_library(LIB_PTHREAD pthread)
    find_library(LIB_MATH m)
    find_library(LIB_DL dl)

  endif(WIN32)

  foreach(EXAMPLE_NAME ${EXAMPLE_LIST})

    FILE(GLOB SOURCES ${PROJECT_SOURCE_DIR}/examples/${EXAMPLE_NAME}/*.cpp)

    add_executable(${EXAMPLE_NAME} ${SOURCES})

    target_link_libraries(${EXAMPLE_NAME} daecpp_static)
    target_link_libraries(${EXAMPLE_NAME} ${LIB_MKL_INTEL_LP64} ${LIB_MKL_THREAD} ${LIB_MKL_CORE} ${LIB_INTEL_OPENMP} ${LIB_PTHREAD} ${LIB_MATH} ${LIB_DL})

    if(PLOTTING)
      target_link_libraries(${EXAMPLE_NAME} ${PYTHON_LIB})
    endif(PLOTTING)

    install(TARGETS ${EXAMPLE_NAME} DESTINATION bin)

  endforeach()

  if(DAE_TEST)

    include(CTest)

    foreach(EXAMPLE_NAME ${EXAMPLE_LIST})
      add_test(${EXAMPLE_NAME}_TEST ${EXAMPLE_NAME})
    endforeach()

  endif(DAE_TEST)

endif(DAE_BUILD_EXAMPLES)
