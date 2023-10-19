#ifndef TEST_TEST_HH
#define TEST_TEST_HH

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

#include "options.hh"
#include "state.hh"
#include "tasks/utils.hh"
#include "tasks/initialize.hh"

namespace test {
inline int execute() {
  mesh::gcoord axis_extents{
    opt::x_extents.value(), opt::y_extents.value(), opt::z_extents.value()};

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();

  {
    mesh::cslot coloring;
    coloring.allocate(num_colors, axis_extents);

    mesh::grect geom;
    geom[0][0] = 0.0;
    geom[0][1] = 1.0;
    geom[1][0] = 0.0;
    geom[1][1] = 1.0;
    geom[2][0] = 0.0;
    geom[2][1] = 1.0;

    m.allocate(coloring.get(), geom);
  } // scope

  flecsi::execute<tasks::init::monotonic>(m, f(m));
  flecsi::execute<tasks::init::boundary<mesh::x_axis, mesh::low>>(m, f(m));
  flecsi::execute<tasks::init::boundary<mesh::x_axis, mesh::high>>(m, f(m));
  flecsi::execute<tasks::init::boundary<mesh::y_axis, mesh::low>>(m, f(m));
  flecsi::execute<tasks::init::boundary<mesh::y_axis, mesh::high>>(m, f(m));
  //flecsi::execute<tasks::init::boundary<mesh::z_axis, mesh::low>>(m, f(m));
  //flecsi::execute<tasks::init::boundary<mesh::z_axis, mesh::high>>(m, f(m));

  flecsi::execute<tasks::util::mesh_info<mesh::x_axis>>(m);
  //flecsi::execute<tasks::util::cell_info<mesh::domain::quantities>>(m);
  flecsi::execute<tasks::util::field_info<mesh::domain::all>>(m, f(m), 2);

  return 0;
} // execute

} // namespace test

#endif // TEST_TEST_HH
