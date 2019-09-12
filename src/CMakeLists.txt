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

set(CTM_SOURCES
    CompilerInfo.cpp
    ConfigurationInfo.cpp
    CTM.cpp
)

add_executable(CTM ${CTM_SOURCES})
add_dependencies(CTM CompilerInfo)
