#ifndef TEST_TEST_HH
#define TEST_TEST_HH

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

#include "options.hh"
#include "state.hh"
#include "tasks/initialize.hh"
#include "tasks/utils.hh"

#include <yaml-cpp/yaml.h>

namespace test {

namespace utils {

inline fvm::mesh::boundary_type
mesh_boundary(std::string const & b) {
  if(b == "flow") {
    return fvm::mesh::boundary_type::flow;
  }
  if(b == "reflecting") {
    return fvm::mesh::boundary_type::reflecting;
  }
  if(b == "periodic") {
    return fvm::mesh::boundary_type::periodic;
  }
  else {
    flog_fatal("invalid boundary type(" << b << ")");
  } // if
} // mesh_boundary

} // namespace utils

inline int
execute() {
  YAML::Node config = YAML::LoadFile(opt::config.value());

  mesh::gcoord axis_extents{
    opt::x_extents.value(), opt::y_extents.value(), opt::z_extents.value()};

  auto bf = flecsi::execute<tasks::init::boundaries>(bmap(gt),
    utils::mesh_boundary(config["boundaries"]["xlow"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["xhigh"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["ylow"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["yhigh"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["zlow"].as<std::string>()),
    utils::mesh_boundary(config["boundaries"]["zhigh"].as<std::string>()));

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();

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

  flecsi::execute<tasks::init::monotonic>(m, f(m));
  flecsi::execute<tasks::init::boundary<mesh::x_axis, mesh::low>>(
    m, f(m), bmap(gt));
  flecsi::execute<tasks::init::boundary<mesh::x_axis, mesh::high>>(
    m, f(m), bmap(gt));
  flecsi::execute<tasks::init::boundary<mesh::y_axis, mesh::low>>(
    m, f(m), bmap(gt));
  flecsi::execute<tasks::init::boundary<mesh::y_axis, mesh::high>>(
    m, f(m), bmap(gt));
  flecsi::execute<tasks::init::boundary<mesh::z_axis, mesh::low>>(
    m, f(m), bmap(gt));
  flecsi::execute<tasks::init::boundary<mesh::z_axis, mesh::high>>(
    m, f(m), bmap(gt));

  flecsi::execute<tasks::util::mesh_info<mesh::x_axis>>(m);
  flecsi::execute<tasks::util::cell_info<mesh::domain::quantities>>(m);
  flecsi::execute<tasks::util::field_info<mesh::domain::all>>(m, f(m), 2);

  return 0;
} // execute

} // namespace test

#endif // TEST_TEST_HH
