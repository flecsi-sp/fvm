#include "analyze.hh"
#include "io.hh"
#include "tasks/io.hh"

#include <flecsi/flog.hh>

using namespace flecsi;

void
muscl::action::analyze(control_policy & cp) {
  flog(info) << "analyze action" << std::endl;
  execute<tasks::io::raw, mpi>(
    io::name{"output-"} << cp.step(), m, r(m), ru(m), rE(m));
} // analyze
