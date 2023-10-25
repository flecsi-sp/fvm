#ifndef MUSCL_TASKS_IO_HH
#define MUSCL_TASKS_IO_HH

#include "../io.hh"
#include "../state.hh"

#include <fstream>

namespace muscl::tasks::io {

void
raw(muscl::io::name const & base,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto ru = m.mdspan<mesh::cells>(ru_a);
  auto rE = m.mdspan<mesh::cells>(rE_a);

  std::ofstream file(
    base.str() + "-" + std::to_string(flecsi::process()) + ".raw");

  {
    auto ccoords = m.color_indeces();
    file << ccoords[mesh::x_axis] << " " << ccoords[mesh::y_axis] << " "
         << ccoords[mesh::z_axis] << std::endl;
    auto ccolors = m.axis_colors();
    file << ccolors[mesh::x_axis] << " " << ccolors[mesh::y_axis] << " "
         << ccolors[mesh::z_axis] << std::endl;
  } // scope

  // Density
  for(auto k : m.cells<mesh::z_axis, mesh::domain::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::domain::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::domain::quantities>()) {
        file << r[k][j][i] << std::endl;
      } // for
    } // for
  } // for

  // Momentum
  for(auto k : m.cells<mesh::z_axis, mesh::domain::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::domain::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::domain::quantities>()) {
        file << ru[k][j][i].x << " " << ru[k][j][i].y << " " << ru[k][j][i].z
             << std::endl;
      } // for
    } // for
  } // for

  // Total Energy
  for(auto k : m.cells<mesh::z_axis, mesh::domain::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::domain::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::domain::quantities>()) {
        file << rE[k][j][i] << std::endl;
      } // for
    } // for
  } // for
} // raw

} // namespace muscl::tasks::io

#endif // MUSCL_TASKS_IO_HH
