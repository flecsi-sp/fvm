#include "initialize.hh"
#include "options.hh"
#include "state.hh"
#include "tasks/boundary.hh"
#include "tasks/initialize.hh"

#include <flecsi/flog.hh>
#include <yaml-cpp/yaml.h>

using namespace flecsi;
using namespace muscl;

void
action::initialize(control_policy & cp) {
  YAML::Node config = YAML::LoadFile(opt::config.value());

  /*--------------------------------------------------------------------------*
    Control model initialization.
   *--------------------------------------------------------------------------*/

  cp.init(config["t0"].as<double>(),
    config["tf"].as<double>(),
    config["max_steps"].as<std::size_t>(),
    config["cfl"].as<double>(),
    config["max_dt"].as<double>());

  /*--------------------------------------------------------------------------*
    Set mesh resolution.
   *--------------------------------------------------------------------------*/

  mesh::gcoord axis_extents{
    opt::x_extents.value(), opt::y_extents.value(), opt::z_extents.value()};

  /*--------------------------------------------------------------------------*
    Gamma.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::gamma>(gamma(gt), config["gamma"].as<double>());

  /*--------------------------------------------------------------------------*
    Set boundaries.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::boundaries>(bmap(gt),
    utils::mesh_boundary(config["boundary"]["xlow"].as<std::string>()),
    utils::mesh_boundary(config["boundary"]["xhigh"].as<std::string>()),
    utils::mesh_boundary(config["boundary"]["ylow"].as<std::string>()),
    utils::mesh_boundary(config["boundary"]["yhigh"].as<std::string>()),
    utils::mesh_boundary(config["boundary"]["zlow"].as<std::string>()),
    utils::mesh_boundary(config["boundary"]["zhigh"].as<std::string>()));

  /*--------------------------------------------------------------------------*
    Topology allocations.
   *--------------------------------------------------------------------------*/

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();
  ct.allocate(num_colors);

  {
    mesh::cslot coloring;
    coloring.allocate(num_colors, axis_extents);

    mesh::grect geom;
    geom[0][0] = config["coords"][0][0].as<double>();
    geom[0][1] = config["coords"][1][0].as<double>();
    geom[1][0] = config["coords"][0][1].as<double>();
    geom[1][1] = config["coords"][1][1].as<double>();
    geom[2][0] = config["coords"][0][2].as<double>();
    geom[2][1] = config["coords"][1][2].as<double>();

    m.allocate(coloring.get(), geom);
  } // scope

  /*--------------------------------------------------------------------------*
    Initialize problem state.
   *--------------------------------------------------------------------------*/

  if(config["problem"].as<std::string>() == "sod") {
    execute<tasks::init::sod>(m, r(m), ru(m), rE(m), gamma(gt));
  }
  else {
    flog_fatal(
      "unsupported problem(" << config["problem"].as<std::string>() << ")");
  } // if

  execute<tasks::apply_boundaries>(m, bmap(gt), r(m), ru(m), rE(m));
  execute<tasks::init::primitives>(
    m, r(m), ru(m), rE(m), u(m), p(m), lmax(ct), gamma(gt));
  cp.dtmin() = reduce<tasks::dtmin, exec::fold::min>(m, lmax(ct));
} // action::initialize
