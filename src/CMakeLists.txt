# Rules for files that need to be configured.

# Configuration files
configure_file(${PROJECT_SOURCE_DIR}/src/Configuration.hpp.in
               ${PROJECT_BINARY_DIR}/src/Configuration.hpp @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/ConfigurationInfo.cpp.in
               ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.cpp @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/ConfigurationInfo.hpp.in
               ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.hpp @ONLY)

set_source_files_properties(${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                            PROPERTIES GENERATED TRUE)

# Configured data files
configure_file(${PROJECT_SOURCE_DIR}/data/callindex.out_CpaD03_0.01
               ${PROJECT_BINARY_DIR}/data/callindex.out_CpaD03_0.01
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/callindex.out_CpaD03_0.10
               ${PROJECT_BINARY_DIR}/data/callindex.out_CpaD03_0.10
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/callindex.out_CpeD03_0.01
               ${PROJECT_BINARY_DIR}/data/callindex.out_CpeD03_0.01
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/callindex.out_CpeD03_0.10
               ${PROJECT_BINARY_DIR}/data/callindex.out_CpeD03_0.10
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/callindex.out_silD03
               ${PROJECT_BINARY_DIR}/data/callindex.out_silD03
               COPYONLY)
set(DRAINEDUSTPROPERTIESDATALOCATION ${PROJECT_BINARY_DIR}/data/)
message(STATUS
  "DraineDustPropertiesDataLocation.hpp: ${DRAINEDUSTPROPERTIESDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/DraineDustPropertiesDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/DraineDustPropertiesDataLocation.hpp
               @ONLY)

set(CosTuuM_SOURCES
    CommandLineOption.cpp
    CommandLineParser.cpp
    CosTuuM.cpp
    ParameterFile.cpp

    ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
    ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.cpp
)

add_executable(CosTuuM ${CosTuuM_SOURCES})
target_link_libraries(CosTuuM ${LIBQUICKSCHED})
add_dependencies(CosTuuM CompilerInfo)

if(HAVE_PYTHON)
  # Add a target for the Python module (this simply calls the Python build
  # script; this script figures out dependencies itself)
  # we overwrite the default build folder to make sure the library file ends up
  # in the bin/ folder
  add_custom_target(CTMmodule ALL
                    COMMENT "Building Python module"
                    COMMAND ${PYTHON_EXECUTABLE}
                              ${PROJECT_BINARY_DIR}/setup.py build
                              --build-lib ${PROJECT_BINARY_DIR}/bin)

  # add an installation script, so that 'make install' installs the Python
  # module (it will then be available for the current user from anywhere in
  # the system)
  install(SCRIPT ${PROJECT_BINARY_DIR}/install_python_module.cmake)
endif(HAVE_PYTHON)
