#------------------------------------------------------------------------------#
# Copyright (c) 2023, Triad National Security, LLC.
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.19)

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

  set(PAGES_REPO "git@re-git.lanl.gov:flecsi-sp/fvm.git"
    CACHE STRING "Pages repository for pages branch")
  option(ENABLE_DOCUMENTATION "Enable documentation" OFF)

  if(ENABLE_DOCUMENTATION)
    include(FleCSI/documentation)
    include(deploy)

    option(ENABLE_DOXYGEN "Enable Doxygen documentation" ON)
    option(ENABLE_SPHINX  "Enable Sphinx documentation"  ON)
    mark_as_advanced(ENABLE_DOXYGEN)
    mark_as_advanced(ENABLE_SPHINX)

    if(ENABLE_DOXYGEN)
      include(FleCSI/doxygen)

      # Create Doxygen output directory eagerly to avoid race condition.
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/api)

      flecsi_add_doxygen_target(user
        CONFIGS
          ${CMAKE_SOURCE_DIR}/doc/doxygen/user.conf.in
          ${CMAKE_SOURCE_DIR}/doc/doxygen/conf.in
      )
    endif()

    if(ENABLE_SPHINX)
      include(FleCSI/sphinx)

      flecsi_add_sphinx_target(main
        CONFIG ${CMAKE_SOURCE_DIR}/doc/sphinx
        OUTPUT ${CMAKE_BINARY_DIR}/doc
      )

      # This is for local install of the documentation. To deploy the
      # documentation to the flecsi website, use the "deploy-documentation"
      # target.
      install(
        DIRECTORY
          ${CMAKE_BINARY_DIR}/doc/
        DESTINATION
          ${CMAKE_INSTALL_PREFIX}/share/doc/fvm
        OPTIONAL
      )
    endif()

    if((NOT PAGES_REPO STREQUAL "") AND ENABLE_SPHINX AND ENABLE_DOXYGEN)
      add_doc_deployment(deploy-docs PAGES_REPO ${PAGES_REPO})
    endif()
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
add_subdirectory(examples/test)
