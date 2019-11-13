## Project initialization ######################################################

# the most advanced feature we use is the ProcessorCount module, which is only
# available since CMake 2.8.5
cmake_minimum_required(VERSION 2.8.5)

# Set the default CMAKE_BUILD_TYPE
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE INTERNAL
        "Set the default CMake build type to Release" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

project(CosTuuM)

# Enable unit testing using CTest
enable_testing()

## Macro definitions ###########################################################

# Add minimal cpp file for flag testing
execute_process(COMMAND ${CMAKE_COMMAND} -E echo
                "int main(int argc, char **argv){ return 0; }"
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/mintest.cpp)

# Macro used to safely add compiler flags (to the C/CXX compiler):
# we try to add the compiler flag and then compile a minimal cpp file
# if this succeeds, we keep the flag
# if it fails, we reset the flags, and throw an error/warning depending on the
# value of the FATAL parameter
macro(add_compiler_flag FLAG FATAL)
  set(BACKUP_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
  try_compile(FLAG_WORKS ${PROJECT_BINARY_DIR}
                         ${PROJECT_BINARY_DIR}/mintest.cpp)
  if(NOT FLAG_WORKS)
    if(${FATAL} STREQUAL "REQUIRED")
      message(FATAL_ERROR "The compiler does not support the ${FLAG} flag!")
    else(${FATAL} STREQUAL "REQUIRED")
      message(STATUS "Not using unsupported compiler flag ${FLAG}.")
      set(CMAKE_CXX_FLAGS ${BACKUP_CXX_FLAGS})
    endif(${FATAL} STREQUAL "REQUIRED")
  endif(NOT FLAG_WORKS)
endmacro(add_compiler_flag)

# Variable used to count the number of configuration options (CMake does not
# actually define counters, so instead we just add +1 to a variable and let the
# CXX preprocessor do the actual addition)
set(CONFIGURATION_OPTIONS_NUMBER "0")

# Macro used to add a configuration option
# We set the given variable name to the given value (as a CMake variable), and
# also add the name and value to two separate lists that will be configured in
# a cpp file. We count the number of configuration variables by increasing the
# counter by 1.
macro(add_configuration_option NAME VALUE)
  set(${NAME} ${VALUE})
  set(CONFIGURATION_OPTIONS_NUMBER "${CONFIGURATION_OPTIONS_NUMBER}+1")
  set(CONFIGURATION_OPTIONS_KEYS "${CONFIGURATION_OPTIONS_KEYS}\"${NAME}\",")
  set(CONFIGURATION_OPTIONS_VALUES
      "${CONFIGURATION_OPTIONS_VALUES}\"${VALUE}\",")
endmacro(add_configuration_option)

## System configuration checks #################################################

# Check if we are compiling on Windows
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  add_configuration_option(HAVE_WINDOWS True)
else(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  add_configuration_option(HAVE_WINDOWS False)
endif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

# If we are on Windows, make sure to add the _USE_MATH_DEFINES define to the
# compilation command
# And the NOMINMAX define to disable use of the min and max macros defined
# by windows.h
# And the _CRT_SECURE_NO_WARNINGS define to disable warnings about unsafe
# sprintf usages
if(HAVE_WINDOWS)
  add_compiler_flag("-D_USE_MATH_DEFINES" REQUIRED)
  add_compiler_flag("-DNOMINMAX" REQUIRED)
  add_compiler_flag("-D_CRT_SECURE_NO_WARNINGS" REQUIRED)
endif(HAVE_WINDOWS)

# Check if we are on a POSIX system. If not, some functionality will not work
execute_process(COMMAND ${CMAKE_COMMAND} -E echo
                "#include <unistd.h>\nint main(int, char**){ return 0; }"
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/posixtest.cpp)
try_compile(DETECT_POSIX ${PROJECT_BINARY_DIR}
                       ${PROJECT_BINARY_DIR}/posixtest.cpp)
if(NOT DETECT_POSIX)
  add_configuration_option(HAVE_POSIX False)
  message(WARNING "This is not a POSIX system! Some features will not work!")
else(NOT DETECT_POSIX)
  add_configuration_option(HAVE_POSIX True)
endif(NOT DETECT_POSIX)

# Find Git
find_package(Git)

# The version number.
execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always
                OUTPUT_VARIABLE GIT_BUILD_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
string(REGEX REPLACE "v([0-9]*).*" "\\1" GIT_VERSION_MAJOR ${GIT_BUILD_STRING})
string(REGEX REPLACE "v[0-9]*.([0-9]*).*" "\\1" GIT_VERSION_MINOR
       ${GIT_BUILD_STRING})
message(STATUS
        "This is CosTuuM version ${GIT_VERSION_MAJOR}.${GIT_VERSION_MINOR}")
set(CosTuuM_VERSION_MAJOR ${GIT_VERSION_MAJOR})
set(CosTuuM_VERSION_MINOR ${GIT_VERSION_MINOR})

# Enable C++11, since we need this for delegating constructors
add_compiler_flag("-std=c++11" REQUIRED)

# Enable all standard compiler warnings and enforce them
add_compiler_flag("-Wall -Werror" OPTIONAL)

# Enable the address sanitizer in debug builds
# (to symbolize the code, run
#   export ASAN_SYMBOLIZER_PATH=<path to llvm-symbolizer>
#   export ASAN_OPTIONS=symbolize=1
#  before running the code)
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  add_compiler_flag("-fsanitize=address -fno-omit-frame-pointer" OPTIONAL)
endif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")

# Add a Profile CMAKE_BUILD_TYPE
set(CMAKE_CXX_FLAGS_PROFILE "-pg -g")

# If not explicitly set, try to find the number of cores on the system
if(MAX_NUMBER_OF_THREADS)
  add_configuration_option(MAX_NUM_THREADS ${MAX_NUMBER_OF_THREADS})
else(MAX_NUMBER_OF_THREADS)
  include(ProcessorCount)
  ProcessorCount(MAX_NUMBER_OF_THREADS)
  if(MAX_NUMBER_OF_THREADS EQUAL 0)
    # If we cannot determine the number of cores, we set it to a safe large
    # value
    message(WARNING
            "Cannot determine core count on this system, using value 128")
    set(MAX_NUMBER_OF_THREADS 128)
  endif(MAX_NUMBER_OF_THREADS EQUAL 0)
  add_configuration_option(MAX_NUM_THREADS ${MAX_NUMBER_OF_THREADS})
endif(MAX_NUMBER_OF_THREADS)

# If we have multiple cores available: find OpenMP.
if(MAX_NUM_THREADS GREATER 1)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_configuration_option(HAVE_OPENMP True)
    message(STATUS "OpenMP found. Shared memory parallelization will work.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  else(OPENMP_FOUND)
    add_configuration_option(HAVE_OPENMP False)
    message(WARNING
            "OpenMP not found. Shared memory parallelization will not work.")
  endif(OPENMP_FOUND)
else(MAX_NUM_THREADS GREATER 1)
  add_configuration_option(HAVE_OPENMP False)
  message(WARNING "Only 1 core available, so not enabling OpenMP support.")
endif(MAX_NUM_THREADS GREATER 1)

find_package(Boost 1.53.0)
if(Boost_FOUND)
  message(STATUS "Found Boost multiprecision.")
  include_directories(${Boost_INCLUDE_DIRS})
  add_configuration_option(HAVE_MULTIPRECISION True)
else(Boost_FOUND)
  message(WARNING
          "Boost multiprecision not found. Exact arithmetics will not work!")
  add_configuration_option(HAVE_MULTIPRECISION False)
endif(Boost_FOUND)

# Check if we need to activate fast math
if(ACTIVATE_FAST_MATH)
  message(STATUS "Activating Fast Math...")
  add_compiler_flag("-ffast-math" OPTIONAL)
  add_configuration_option(HAVE_FAST_MATH True)
else(ACTIVATE_FAST_MATH)
  add_configuration_option(HAVE_FAST_MATH False)
endif(ACTIVATE_FAST_MATH)

# Check if we need to activate vectorization
if(ACTIVATE_VECTORIZATION)
  message(STATUS "Activating vectorization...")
  add_compiler_flag("-ftree-vectorize" OPTIONAL)
  add_configuration_option(HAVE_VECTORIZATION True)
else(ACTIVATE_VECTORIZATION)
  add_configuration_option(HAVE_VECTORIZATION False)
endif(ACTIVATE_VECTORIZATION)

# Check if we need to activate assertions
if(ACTIVATE_ASSERTIONS)
  add_configuration_option(HAVE_ASSERTIONS True)
  message(STATUS "Assertions are activated.")
else(ACTIVATE_ASSERTIONS)
  add_configuration_option(HAVE_ASSERTIONS False)
  message(STATUS "Assertions not activated.")
endif(ACTIVATE_ASSERTIONS)

# Check if we need to activate quad precision
if(ACTIVATE_QUAD_PRECISION)
  add_configuration_option(HAVE_QUAD_PRECISION True)
  message(STATUS "Using quad precision.")
else(ACTIVATE_QUAD_PRECISION)
  add_configuration_option(HAVE_QUAD_PRECISION False)
  message(STATUS "Using double precision.")
endif(ACTIVATE_QUAD_PRECISION)

if(ACTIVATE_TASKS)
  add_configuration_option(HAVE_TASKS True)
  set(LIBQUICKSCHED quicksched)
  include_directories(${PROJECT_SOURCE_DIR}/quicksched)
  message(STATUS "Using task-based parallel algorithm.")
else(ACTIVATE_TASKS)
  add_configuration_option(HAVE_TASKS False)
  message(STATUS "Not using task-based parallel algorithm.")
  set(LIBQUICKSCHED "")
endif(ACTIVATE_TASKS)

## Code configuration ##########################################################

# Tell CMake that headers are in one of the src folders.
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src)

# Set the path where the generated executables are stored.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Set the path where the generated libraries are stored
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# Add the command and target that creates the CompilerInfo.cpp file
# Note that we need to add the 'dummy_file_that_is_not_created' to the outputs
# to force rerunnning this command every time the code is compiled
add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                          dummy_file_that_is_not_created
                   COMMENT "Generating CompilerInfo.cpp"
                   COMMAND ${CMAKE_COMMAND}
                           -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
                           -DCOMPILER_NAME=${CMAKE_CXX_COMPILER_ID}
                           -DCOMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
                           -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                           -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                           -DHAVE_WINDOWS=${HAVE_WINDOWS}
                           -P ${PROJECT_SOURCE_DIR}/write_compiler_info.cmake
                   DEPENDS ${PROJECT_SOURCE_DIR}/.git/HEAD
                           ${PROJECT_SOURCE_DIR}/src/CompilerInfo.cpp.in
                   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_custom_target(CompilerInfo
                  DEPENDS ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp)

# Prepare the Python module setup script
find_package(PythonLibs 3.0)
find_package(PythonInterp 3.0)
if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                              "import numpy; print(numpy.get_include())"
                  RESULT_VARIABLE NUMPY_RESULT
                  OUTPUT_VARIABLE NUMPY_OUTPUT
                  ERROR_VARIABLE NUMPY_ERROR)
  if(NUMPY_RESULT STREQUAL "0")
    message(STATUS "Python and NumPy found. Building Python module.")
    configure_file(${PROJECT_SOURCE_DIR}/setup.py.in
                   ${PROJECT_BINARY_DIR}/setup.py @ONLY)
    configure_file(${PROJECT_SOURCE_DIR}/install_python_module.cmake.in
                   ${PROJECT_BINARY_DIR}/install_python_module.cmake @ONLY)
    add_configuration_option(HAVE_PYTHON True)
  else(NUMPY_RESULT STREQUAL "0")
    message(WARNING "Could not find NumPy. Not building Python module.")
    add_configuration_option(HAVE_PYTHON False)
  endif(NUMPY_RESULT STREQUAL "0")
else(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
  message(WARNING "Python not found. Not building Python module.")
  add_configuration_option(HAVE_PYTHON False)
endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)

# Enter the folder containing the source files and configure the executables.
add_subdirectory(src)

# Enter the folder containing the unit tests
add_subdirectory(test)

if(HAVE_TASKS)
  add_subdirectory(quicksched)
endif(HAVE_TASKS)

## Documentation configuration #################################################

# Generate documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  message(STATUS
          "Doxygen found! You can generate documentation using 'make doc'")
  if(DOXYGEN_DOT_FOUND)
    set(HAVE_DOT YES)
  else(DOXYGEN_DOT_FOUND)
    set(HAVE_DOT NO)
  endif(DOXYGEN_DOT_FOUND)
  # Configure the Doxyfile (sets the correct output path and Dot path)
  configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in
                 ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
  # Add a command that will generate the doxygen documentation
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/doc/html/index.html
                     COMMAND ${DOXYGEN_EXECUTABLE}
                             ${PROJECT_BINARY_DIR}/Doxyfile
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Generating API documentation with Doxygen"
                     VERBATIM)
  # Add a custom target that calls the command defined above
  add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html)
else(DOXYGEN_FOUND)
# Could not find doxygen. Inform the user that documentation will not be
# available.
message(WARNING "Doxygen not found. 'make doc' will not work!")
endif(DOXYGEN_FOUND)
