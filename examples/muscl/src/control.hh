#ifndef POISSON_CONTROL_HH
#define POISSON_CONTROL_HH

#include "state.hh"
#include "tasks/advance.hh"

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>
#include <flecsi/run/control.hh>

namespace muscl {

enum class cp { initialize, advance, analyze, finalize };

inline const char *
operator*(cp control_point) {
  switch(control_point) {
    case cp::initialize:
      return "initialize";
    case cp::advance:
      return "advance";
    case cp::analyze:
      return "analyze";
    case cp::finalize:
      return "finalize";
  }
  flog_fatal("invalid control point");
}

struct control_policy : flecsi::run::control_base {

  using control_points_enum = cp;
  struct node_policy {};
  using control = flecsi::run::control<control_policy>;

  void init(double t0, double tf, double max_steps, double cfl) {
    t0_ = t0;
    tf_ = tf;
    max_steps_ = max_steps;
    cfl_ = cfl;
  } // init

  auto & dt_future() {
    return dt_future_;
  }

  static bool cycle_control(control_policy & cp) {
    cp.dt_future_.get();
    flog(info) << "cycle control: " << cp.step_++ << std::endl;
    return cp.step_ < cp.max_steps_;
  } // cycle_control

  using evolve = cycle<cycle_control, point<cp::advance>, point<cp::analyze>>;

  using control_points =
    list<point<cp::initialize>, evolve, point<cp::finalize>>;

private:
  std::size_t step_{0};
  double t0_;
  double tf_;
  double max_steps_{5};
  double cfl_;
  flecsi::future<double, flecsi::exec::launch_type_t::index> dt_future_;
}; // struct control_policy

using control = flecsi::run::control<control_policy>;

} // namespace muscl

#endif // POISSON_CONTROL_HH
