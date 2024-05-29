#------------------------------------------------------------------------------#
# Copyright (c) 2023, Triad National Security, LLC.
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.19)

#------------------------------------------------------------------------------#
# CMake policies.
#------------------------------------------------------------------------------#

cmake_policy(SET CMP0144 NEW) # ``PACKAGE_NAME''_ROOT

#------------------------------------------------------------------------------#
# Project.
#------------------------------------------------------------------------------#

project(FleCSI-SP-FVM LANGUAGES C CXX)
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

set(default_build_type "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS
    "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------#
# Custom unit test function. This needs to be defined unconditionally, i.e.,
# outside of the top-level project condition.
#------------------------------------------------------------------------------#

function(fvm_add_unit name)
  flecsi_add_test(${name} ${ARGN})
  if(TARGET ${name})
    target_link_libraries(${name} PUBLIC FleCSI-SP-FVM::FleCSI-SP-FVM)
  endif()
endfunction()

if(PROJECT_IS_TOP_LEVEL)

  #----------------------------------------------------------------------------#
  # Add local module path.
  #----------------------------------------------------------------------------#

  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

  #----------------------------------------------------------------------------#
  # Find the FleCSI package. This must be done early so that the FleCSI cmake
  # packages will be available to the project, e.g., 'include(FleCSI/format)'.
  #
  # CMake documentation for FleCSI client applications is included in the
  # FleCSI userguide.
  #----------------------------------------------------------------------------#

  if(NOT FleCSI_FOUND)
    find_package(FleCSI 2 REQUIRED)
  endif()

  #----------------------------------------------------------------------------#
  # Formatting
  #----------------------------------------------------------------------------#

  option(ENABLE_FORMAT "Enable format target" OFF)
  mark_as_advanced(ENABLE_FORMAT)

  if(ENABLE_FORMAT)
    include(FleCSI/format)

    set(CLANG_FORMAT_VERSION "13...<14" CACHE STRING
        "Set the required version (major[.minor[.patch]]) of clang-format")
    mark_as_advanced(CLANG_FORMAT_VERSION)

    flecsi_add_format_target(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}
      "${CLANG_FORMAT_VERSION}")
  endif()

  #----------------------------------------------------------------------------#
  # Documentation
  #----------------------------------------------------------------------------#

  option(ENABLE_DOCUMENTATION "Enable documentation" OFF)

  if(ENABLE_DOCUMENTATION)
    include(FleCSI/documentation)
    include(FleCSI/doxygen)
    include(FleCSI/sphinx)

    # Avoid race condition by proactively creating output directory
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/.doxygen/xml)
    flecsi_add_doxygen_target(fvm
      CONFIGS
        ${CMAKE_SOURCE_DIR}/doc/doxygen/fvm.conf.in
        ${CMAKE_SOURCE_DIR}/doc/doxygen/conf.in
    )

    file(COPY ${CMAKE_SOURCE_DIR}/doc/sphinx
      DESTINATION ${CMAKE_BINARY_DIR}
    )
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/.sphinx)
    file(RENAME ${CMAKE_BINARY_DIR}/sphinx
      ${CMAKE_BINARY_DIR}/.sphinx
    )
    # Create output directory
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/.sphinx/fvm/src)
    execute_process(COMMAND
      ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/spec
      ${CMAKE_BINARY_DIR}/.sphinx/fvm/src/spec
    )
    execute_process(COMMAND
      ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/examples
      ${CMAKE_BINARY_DIR}/.sphinx/fvm/src/examples
    )

    flecsi_add_sphinx_target(main
      CONFIG ${CMAKE_BINARY_DIR}/.sphinx
      OUTPUT ${CMAKE_BINARY_DIR}/doc
    )
  endif()

  #----------------------------------------------------------------------------#
  # Unit Tests
  #----------------------------------------------------------------------------#

  option(ENABLE_UNIT_TESTS "Enable unit testing" OFF)
  mark_as_advanced(ENABLE_UNIT_TESTS)
  if(ENABLE_UNIT_TESTS)
    include(FleCSI/unit)
    flecsi_enable_testing()
    add_subdirectory(spec)
  endif()
endif()

#------------------------------------------------------------------------------#
# Add library.
#------------------------------------------------------------------------------#

add_library(FleCSI-SP-FVM INTERFACE)
add_library(FleCSI-SP-FVM::FleCSI-SP-FVM ALIAS FleCSI-SP-FVM)
target_include_directories(FleCSI-SP-FVM INTERFACE
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

#------------------------------------------------------------------------------#
# Add examples.
#------------------------------------------------------------------------------#

add_subdirectory(examples/muscl)
