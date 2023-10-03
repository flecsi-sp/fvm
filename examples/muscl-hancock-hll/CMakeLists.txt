#------------------------------------------------------------------------------#
# Add application headers.
#------------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#

add_executable(muscl
  src/initialize.cc
  src/muscl.cc
  src/tasks/initialize.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#

target_link_libraries(muscl FleCSI::FleCSI
  FleCSI-SP-FVM::FleCSI-SP-FVM
)
