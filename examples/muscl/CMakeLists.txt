#------------------------------------------------------------------------------#
# Yaml
#------------------------------------------------------------------------------#

find_package(yaml-cpp REQUIRED)

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#

add_executable(muscl
  src/advance.cc
  src/analyze.cc
  src/init.cc
  src/muscl.cc
  src/tasks/hydro.cc
  src/tasks/init.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#

target_link_libraries(muscl FleCSI::FleCSI
  FleCSI-SP-FVM::FleCSI-SP-FVM
  yaml-cpp
)
