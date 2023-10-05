#------------------------------------------------------------------------------#
# Add application headers.
#------------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

#------------------------------------------------------------------------------#
# Yaml
#------------------------------------------------------------------------------#

find_package(yaml-cpp REQUIRED)

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#

add_executable(muscl
  src/advance.cc
  src/initialize.cc
  src/muscl.cc
  src/tasks/advance.cc
  src/tasks/initialize.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#

target_link_libraries(muscl FleCSI::FleCSI
  FleCSI-SP-FVM::FleCSI-SP-FVM
  yaml-cpp
)
