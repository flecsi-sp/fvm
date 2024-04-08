#include "analyze.hh"
#include "io.hh"
#include "tasks/io.hh"

#include <flecsi/flog.hh>

using namespace flecsi;

void
muscl::action::analyze(control_policy & cp) {
  auto lm = data::launch::make(m);
  execute<tasks::io::raw, mpi>(
    io::name{"output-"} << cp.step(), lm, r(lm), ru(lm), rE(lm));
} // analyze
