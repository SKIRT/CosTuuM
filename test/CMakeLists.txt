### Convenient macros to automate unit test generation #########################
### Do not touch the code below unless you know what you're doing! #############

# Add a new unit test
# A new target with the test sources is constructed, and a CTest test with the
# same name is created. The new test is also added to the global list of test
# contained in the check target
macro(add_unit_test)
    set(options PARALLEL)
    set(oneValueArgs NAME)
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating ${TEST_NAME}")
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/bin/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})

    if(TEST_PARALLEL AND HAVE_MPI)
      set(TESTCOMMAND ${MPIEXEC})
      set(TESTARGS ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
                   "./${TEST_NAME}" ${MPIEXEC_POSTFLAGS})
      set(TESTCOMMAND ${TESTCOMMAND} ${TESTARGS})
    else(TEST_PARALLEL AND HAVE_MPI)
      set(TESTCOMMAND ${TEST_NAME})
    endif(TEST_PARALLEL AND HAVE_MPI)
    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
             COMMAND ${TESTCOMMAND})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})

    # create a custom target for this test alone
    add_custom_target(check_${TEST_NAME}
                      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                              -R ${TEST_NAME}
                      DEPENDS ${TEST_NAME})
endmacro(add_unit_test)

### Actual unit test generation ################################################
### Add new unit tests below ###################################################

## Unit test for SpecialFunctions
set(TESTSPECIALFUNCTIONS_SOURCES
    testSpecialFunctions.cpp
)
add_unit_test(NAME testSpecialFunctions
              SOURCES ${TESTSPECIALFUNCTIONS_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/plot_test_special_functions.py
               ${PROJECT_BINARY_DIR}/bin/test/plot_test_special_functions.py
               COPYONLY)

## Unit test for Matrix
set(TESTMATRIX_SOURCES
    testMatrix.cpp
)
add_unit_test(NAME testMatrix
              SOURCES ${TESTMATRIX_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/test_matrix.py
               ${PROJECT_BINARY_DIR}/bin/test/test_matrix.py
               COPYONLY)

## Unit test for TMatrixCalculator
set(TESTTMATRIXCALCULATOR_SOURCES
    testTMatrixCalculator.cpp
)
add_unit_test(NAME testTMatrixCalculator
              SOURCES ${TESTTMATRIXCALCULATOR_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/test_tmatrixcalculator.txt
               ${PROJECT_BINARY_DIR}/bin/test/test_tmatrixcalculator.txt
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/plot_tmatrixcalculator_result.py
               ${PROJECT_BINARY_DIR}/bin/test/plot_tmatrixcalculator_result.py
               COPYONLY)

## Unit test for OrientationDistribution
set(TESTORIENTATIONDISTRIBUTION_SOURCES
    testOrientationDistribution.cpp
)
add_unit_test(NAME testOrientationDistribution
              SOURCES ${TESTORIENTATIONDISTRIBUTION_SOURCES})
configure_file(
  ${PROJECT_SOURCE_DIR}/test/plot_test_orientationdistribution.py
  ${PROJECT_BINARY_DIR}/bin/test/plot_test_orientationdistribution.py
  COPYONLY)

### Done adding unit tests. Create the 'make check' target #####################
### Do not touch these lines unless you know what you're doing! ################
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                        DEPENDS ${TESTNAMES})
