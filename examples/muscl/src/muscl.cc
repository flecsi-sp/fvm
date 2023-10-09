#include "advance.hh"
#include "analyze.hh"
#include "control.hh"
#include "initialize.hh"

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

int
main(int argc, char ** argv) {
  auto status = flecsi::initialize(argc, argv);
  status = muscl::control::check_status(status);

  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }

  flecsi::flog::add_output_stream("clog", std::clog, true);

  status = flecsi::start(muscl::control::execute);

  flecsi::finalize();

  return status;
} // main
