#include "initialize.hh"
#include "options.hh"
#include "state.hh"
#include "tasks/initialize.hh"

#include <flecsi/flog.hh>
//#include <yaml-cpp/yaml.h>

using namespace flecsi;
using namespace muscl;

void
action::initialize(control_policy & cp) {
  //flog(info) << "config: " << opt::config.value() << std::endl;
  //YAML::Node config = YAML::LoadFile(opt::config.value());

  //flog(info) << "config: " << config["problem"].as<std::string>() << std::endl;
  mesh::gcoord axis_extents{
    opt::x_extents.value(), opt::y_extents.value(), opt::z_extents.value()};

  // FIXME: Boundaries need to be set by problem type.
  // Add command-line option to specify problem type and do proper
  // initialization.
  mesh::bmap boundaries;
  boundaries[mesh::x_axis][mesh::low] = mesh::inflow;
  boundaries[mesh::x_axis][mesh::high] = mesh::outflow;
  boundaries[mesh::y_axis][mesh::low] = mesh::inflow;
  boundaries[mesh::y_axis][mesh::high] = mesh::outflow;
  boundaries[mesh::z_axis][mesh::low] = mesh::inflow;
  boundaries[mesh::z_axis][mesh::high] = mesh::outflow;

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();
  //ct.allocate(num_colors);

  {
    mesh::cslot coloring;
    coloring.allocate(num_colors, axis_extents, boundaries);

    // FIXME: problem size must be consistent with initialization, e.g., sod.
    mesh::grect geom;
    geom[0][0] = 0.0;
    geom[0][1] = 1.0;
    geom[1] = geom[0];
    geom[2] = geom[0];

    m.allocate(coloring.get(), geom);
  } // scope

  //execute<tasks::check>(m);
 // execute<tasks::sod>(m, r(m), ru(m), rE(m), 1.4);
  //execute<tasks::init>(m, r(m), ru(m), rE(m), u(m), p(m), lmax(ct), 1.4);
} // action::initialize
