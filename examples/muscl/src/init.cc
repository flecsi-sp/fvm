#include "init.hh"
#include "options.hh"
#include "state.hh"
#include "tasks/boundary.hh"
#include "tasks/hydro.hh"
#include "tasks/init.hh"
#include "tasks/utils.hh"

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
    config["max_dt"].as<double>(),
    config["log_frequency"].as<std::size_t>());

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

  auto bf = execute<tasks::init::boundaries>(bmap(gt),
    utils::mesh_boundary(config["boundaries"]["xlow"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["xhigh"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["ylow"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["yhigh"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["zlow"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["zhigh"].as<std::string>()));

  /*--------------------------------------------------------------------------*
    Topology allocations.
   *--------------------------------------------------------------------------*/

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();
  ct.allocate(num_colors);

  {
    mesh::grect geom;
    geom[0][0] = config["coords"][0][0].as<double>();
    geom[0][1] = config["coords"][1][0].as<double>();
    geom[1][0] = config["coords"][0][1].as<double>();
    geom[1][1] = config["coords"][1][1].as<double>();
    geom[2][0] = config["coords"][0][2].as<double>();
    geom[2][1] = config["coords"][1][2].as<double>();

    m.allocate(mesh::mpi_coloring{num_colors, axis_extents, bf.get()}, geom);
  } // scope

  /*--------------------------------------------------------------------------*
    Fake initialization to avoid legion warnings.
   *--------------------------------------------------------------------------*/

  // clang-format off
  execute<tasks::init::touch>(m,
    r(m), ru(m), rE(m),
    u(m), p(m),
    q(m), qu(m), qE(m),
    dr_ds(m), du_ds(m), dp_ds(m),
    rTail(m), ruTail(m), rETail(m), uTail(m), pTail(m),
    rHead(m), ruHead(m), rEHead(m), uHead(m), pHead(m),
    rF(m), ruF(m), rEF(m));
  execute<tasks::init::touch1>(m,
    r(m), ru(m), rE(m),
    u(m), p(m),
    q(m), qu(m), qE(m),
    dr_ds(m), du_ds(m), dp_ds(m),
    rTail(m), ruTail(m), rETail(m), uTail(m), pTail(m),
    rHead(m), ruHead(m), rEHead(m), uHead(m), pHead(m),
    rF(m), ruF(m), rEF(m));
  // clang-format on

  /*--------------------------------------------------------------------------*
    Initialize problem state.
   *--------------------------------------------------------------------------*/

  if(config["problem"].as<std::string>() == "sod") {
    execute<tasks::init::shock<sod>>(m, r(m), ru(m), rE(m), gamma(gt));
  }
  else if(config["problem"].as<std::string>() == "rankine-hugoniot") {
    execute<tasks::init::shock<rankine_hugoniot>>(
      m, r(m), ru(m), rE(m), gamma(gt));
  }
  else {
    flog_fatal(
      "unsupported problem(" << config["problem"].as<std::string>() << ")");
  } // if

  /*--------------------------------------------------------------------------*
    Initialize time advance.
   *--------------------------------------------------------------------------*/

  execute<tasks::apply_boundaries>(m, bmap(gt), r(m), ru(m), rE(m));
  execute<tasks::hydro::update_primitives>(
    m, r(m), ru(m), rE(m), u(m), p(m), gamma(gt));
  execute<tasks::hydro::update_eigenvalues>(
    m, r(m), u(m), p(m), lmax(ct), gamma(gt));
  cp.dtmin() = reduce<tasks::hydro::update_dtmin, exec::fold::min>(m, lmax(ct));
} // action::initialize
