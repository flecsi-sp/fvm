#ifndef MUSCL_TASKS_IO_HH
#define MUSCL_TASKS_IO_HH

#include "../io.hh"
#include "../state.hh"

#include <fstream>

namespace muscl::tasks::io {

void inline raw(muscl::io::name const & base,
  multi<mesh::accessor<ro>> mm,
  multi<field<double>::accessor<ro, ro>> r_ma,
  multi<field<vec3>::accessor<ro, ro>> ru_ma,
  multi<field<double>::accessor<ro, ro>> rE_ma) {

  std::size_t i{0};
  for(auto const [c, m] : mm.components()) {
    auto r_a = r_ma.accessors()[i];
    auto ru_a = ru_ma.accessors()[i];
    auto rE_a = rE_ma.accessors()[i];
    auto r = m.mdcolex<mesh::cells>(r_a);
    auto ru = m.mdcolex<mesh::cells>(ru_a);
    auto rE = m.mdcolex<mesh::cells>(rE_a);

    std::ofstream file(base.str() + "-" + std::to_string(c) + ".raw");

    file << m.size<mesh::x_axis, mesh::domain::quantities>() << " "
         << m.size<mesh::y_axis, mesh::domain::quantities>() << " "
         << m.size<mesh::z_axis, mesh::domain::quantities>() << std::endl;
    file << m.size<mesh::x_axis, mesh::domain::global>() << " "
         << m.size<mesh::y_axis, mesh::domain::global>() << " "
         << m.size<mesh::z_axis, mesh::domain::global>() << std::endl;

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
          file << r(i, j, k) << std::endl;
        } // for
      } // for
    } // for

    // Momentum
    for(auto k : m.cells<mesh::z_axis, mesh::domain::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::domain::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::domain::quantities>()) {
          file << ru(i, j, k).x << " " << ru(i, j, k).y << " " << ru(i, j, k).z
               << std::endl;
        } // for
      } // for
    } // for

    // Total Energy
    for(auto k : m.cells<mesh::z_axis, mesh::domain::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::domain::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::domain::quantities>()) {
          file << rE(i, j, k) << std::endl;
        } // for
      } // for
    } // for
    ++i;
  } // for
} // raw

} // namespace muscl::tasks::io

#endif // MUSCL_TASKS_IO_HH
