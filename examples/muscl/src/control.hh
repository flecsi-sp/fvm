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

  void
  init(double t0, double tf, std::size_t max_steps, double cfl, double max_dt) {
    t_ = t0;
    tf_ = tf;
    max_steps_ = max_steps;
    cfl_ = cfl;
    max_dt_ = max_dt;
  } // init

  auto & dtmin() {
    return dtmin_;
  }

  double dt() {
    return dt_;
  }

  static bool cycle_control(control_policy & cp) {
    cp.dt_ = cp.cfl_ * cp.dtmin_.get();
    cp.t_ += cp.dt_;
    ++cp.step_;
    flog(info) << "step: " << cp.step_ << " dt: " << cp.dt_ << std::endl;
    return cp.step_ < cp.max_steps_ && cp.dt_ < cp.max_dt_;
  } // cycle_control

  using evolve = cycle<cycle_control, point<cp::advance>, point<cp::analyze>>;

  using control_points =
    list<point<cp::initialize>, evolve, point<cp::finalize>>;

private:
  std::size_t step_{0};
  double t_;
  double tf_;
  std::size_t max_steps_;
  double cfl_;
  double dt_;
  double max_dt_;
  flecsi::future<double> dtmin_;
}; // struct control_policy

using control = flecsi::run::control<control_policy>;

} // namespace muscl

#endif // POISSON_CONTROL_HH
