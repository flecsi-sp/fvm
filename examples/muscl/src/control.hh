#ifndef POISSON_CONTROL_HH
#define POISSON_CONTROL_HH

#include "state.hh"
#include "tasks/hydro.hh"

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

  void init(double t0,
    double tf,
    std::size_t max_steps,
    double cfl,
    double max_dt,
    std::size_t log_frequency) {
    t_ = t0;
    tf_ = tf;
    max_steps_ = max_steps;
    cfl_ = cfl;
    max_dt_ = max_dt;
    log_frequency_ = log_frequency;
  } // init

  auto & dtmin() {
    return dtmin_;
  }

  double dt() {
    return dt_;
  }

  static bool cycle_control(control_policy & cp) {
    bool exec_cycle = cp.step_ < cp.max_steps_ && cp.t_ < cp.tf_;

#if 1 // FIXME: Debug
    flog(warn) << "dtmin: " << cp.dtmin_.get() << std::endl;
#endif

    cp.dt_ = cp.cfl_ * cp.dtmin_.get();
    cp.dt_ = cp.t_ + cp.dt_ > cp.tf_ ? cp.tf_ - cp.t_ : cp.dt_;
    cp.dt_ = std::min(cp.dt_, cp.max_dt_);

    if((cp.step_ % cp.log_frequency_) == 0 || cp.step_ == cp.max_steps_ ||
       cp.t_ == cp.tf_) {
      flog(info) << "step: " << cp.step_ << " time: " << cp.t_
                 << " dt: " << cp.dt_ << std::endl;
      flecsi::flog::flush();
    } // if

    cp.t_ += cp.dt_;
    ++cp.step_;

    return exec_cycle;
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
  std::size_t log_frequency_;
  flecsi::future<double> dtmin_;
}; // struct control_policy

using control = flecsi::run::control<control_policy>;

} // namespace muscl

#endif // POISSON_CONTROL_HH
