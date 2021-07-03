#
# Build configuration for the SmoothSystematics package
#
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(CommonSmoothingTool)
  set(HAS_PARENT 0)
else()
  set(HAS_PARENT 1)
endif()
cmake_minimum_required(VERSION 3.0)

# Set the name of the package:
find_package(AnalysisBase QUIET)

IF(${HAS_PARENT})
  set(HAS_CommonSmoothingTool 1 PARENT_SCOPE)
ENDIF()

# External(s) used by the package:
find_package( ROOT COMPONENTS Core Hist RIO Tree Gpad )

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )

if( ${AnalysisBase_FOUND} )
  IF(NOT ${HAS_PARENT})
    atlas_ctest_setup() # Set up the project:
    atlas_project( CommonSmoothingTool 1.0.0
      USE AnalysisBase ${AnalysisBase_VERSION} )

    # Generate an environment setup script:
    lcg_generate_env( SH_FILE ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh )
    install( FILES ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh DESTINATION . )

    # Set up the usage of CPack:
    set(CMAKE_INSTALL_PREFIX /InstallArea/${ATLAS_PLATFORM})
    atlas_cpack_setup()
  ENDIF()
  
  file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/lib\n")  
  
  atlas_subdir( CommonSystSmoothingTool )
  
  # Set up which packages this package depends on:
  atlas_depends_on_subdirs(
    PUBLIC
    PRIVATE
    )
  
  # # Generate a dictionary for the library:
  atlas_add_root_dictionary( CommonSystSmoothingToolLib CommonSystSmoothingToolLibCintDict
    ROOT_HEADERS  SmoothSystematics/*.h Root/LinkDef.h 
    EXTERNAL_PACKAGES ROOT 
    )
  
  # Build the shared library of the package:
  atlas_add_library( CommonSystSmoothingToolLib SmoothSystematics/*.h Root/*.cxx Root/LinkDef.h ${CommonSystSmoothingToolLibCintDict}
    PUBLIC_HEADERS SmoothSystematics 
    INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} 
    LINK_LIBRARIES ${ROOT_LIBRARIES} 
    )
  
  
  # Build the executables of the package:
  atlas_add_executable(testSmooth util/testSmooth.cxx
    LINK_LIBRARIES CommonSystSmoothingToolLib )
  
  atlas_add_executable(testInputFile util/testInputFile.cxx 
    LINK_LIBRARIES CommonSystSmoothingToolLib  )
  # Install the script(s) of the package:
  atlas_install_scripts( python/*.py scripts/* )
else()
  file(GLOB SmoothingToolSources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "Root/[A-Za-z]*.cxx")
  file(GLOB SmoothingToolHeaders RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "CommonSystSmoothing/[A-Za-z]*.h")  

  include(${ROOT_USE_FILE})
  include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")  
  ROOT_GENERATE_DICTIONARY(G__CommonSmoothingTool ${SmoothingToolHeaders} LINKDEF Root/LinkDef.h)
  add_library( CommonSmoothingTool SHARED ${SmoothingToolSources} G__CommonSmoothingTool.cxx )

  add_executable(testSmooth util/testSmooth.cxx)
  target_link_libraries ( testSmooth CommonSmoothingTool ${ROOT_LIBRARIES} )  
  
  add_executable(testInputFile util/testInputFile.cxx)
  target_link_libraries ( testInputFile CommonSmoothingTool ${ROOT_LIBRARIES} )      

  file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_CURRENT_BINARY_DIR}\n")

  # Export the package for use from the build-tree
  # (this registers the build-tree with a global CMake-registry)
  export(PACKAGE CommonSmoothingTool)

  set(CommonSmoothingTool_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set(CommonSmoothingTool_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set(CommonSmoothingTool_LIBRARIES CommonSmooothingTool PARENT_SCOPE)
  set(CommonSmoothingTool_LIBRARY_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)
  
  set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
  set(CONF_LIBRARY_DIRS "${CMAKE_CURRENT_BINARY_DIR}")
  set(CONF_LIBRARIES    CommonSmoothingTool)
  configure_file(CommonSmoothingToolConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/CommonSmoothingToolConfig.cmake" @ONLY)

  # Install the CommonSmoothingToolsConfig.cmake
  install(FILES
    "${PROJECT_BINARY_DIR}/CommonSmoothingToolsConfig.cmake"
    DESTINATION "${PROJECT_SOURCE_DIR}" COMPONENT dev
    )
  
   
endif()

add_test(NAME testInput WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/share COMMAND  testInputFile)
add_test(NAME testSmooth WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/share COMMAND  testSmooth)

enable_testing()





