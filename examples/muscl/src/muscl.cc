#include "advance.hh"
#include "analyze.hh"
#include "control.hh"
#include "initialize.hh"

#include <flecsi/runtime.hh>

using namespace flecsi;

int
main(int argc, char ** argv) {
  run::arguments args(argc, argv);
  const run::dependencies_guard dg(args.dep);
  const runtime run(args.cfg);
  flecsi::flog::add_output_stream("flog", std::clog, true);
  return run.main<muscl::control>(args.act);
} // main
