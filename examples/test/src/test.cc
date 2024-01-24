#include "test.hh"

#include <flecsi/runtime.hh>

using namespace flecsi;

int
main(int argc, char ** argv) {
  const run::dependencies_guard dg;
  return flecsi::runtime().control<flecsi::run::call>(test::execute);
} // main
