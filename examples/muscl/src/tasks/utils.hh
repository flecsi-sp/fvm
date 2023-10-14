#ifndef MUSCL_TASKS_UTIL_HH
#define MUSCL_TASKS_UTIL_HH

#include "../state.hh"

#include <sstream>

namespace muscl::tasks::util {

template<mesh::axis A>
inline void
mesh_info(mesh::accessor<ro> m) {
  std::stringstream ss;

  ss << "Mesh Info:" << std::endl;
  ss << "  sizes:" << std::endl;
  ss << "    quantities: " << m.size<A, mesh::domain::quantities>()
     << std::endl;
  ss << "    predictor: " << m.size<A, mesh::domain::predictor>() << std::endl;
  ss << "    corrector: " << m.size<A, mesh::domain::corrector>() << std::endl;
  ss << "    all: " << m.size<A, mesh::domain::all>() << std::endl;
  flog(info) << ss.str() << std::endl;
} // mesh_info

template<mesh::domain DM>
inline void
cell_info(mesh::accessor<ro> m) {
  std::stringstream ss;
  for(auto k : m.cells<mesh::z_axis, DM>()) {
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << "cell: (" << i << ", " << j << ", " << k << ")" << std::endl;
        ss << "\tcenter: (" << m.center<mesh::x_axis>(i) << ", "
           << m.center<mesh::x_axis>(j) << ", " << m.center<mesh::x_axis>(k)
           << ")" << std::endl;
        ss << "\ttail: (" << m.tail<mesh::x_axis>(i) << ", "
           << m.tail<mesh::x_axis>(j) << ", " << m.tail<mesh::x_axis>(k) << ")"
           << std::endl;
        ss << "\thead: (" << m.head<mesh::x_axis>(i) << ", "
           << m.head<mesh::x_axis>(j) << ", " << m.head<mesh::x_axis>(k) << ")"
           << std::endl;
      } // for
    } // for
  } // for

  flog(info) << ss.str() << std::endl;
} // cell_info

template<mesh::domain DM>
inline void
print_conserved(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  flecsi::util::id zslice) {
  {
    auto r = m.mdspan<mesh::cells>(r_a);
    std::stringstream ss;
    ss << "DENSITY:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM, true>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << r[zslice][j][i] << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
#if 0
  {
    auto ru = m.mdspan<mesh::cells>(ru_a);
    std::stringstream ss;
    ss << "MOMENTUM:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM, true>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << ru[zslice][j][i] << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto rE = m.mdspan<mesh::cells>(rE_a);
    std::stringstream ss;
    ss << "TOTAL ENERGY:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM, true>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << rE[zslice][j][i] << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
#endif
} // print_conserved

template<mesh::domain DM>
inline void
print_primitives(mesh::accessor<ro> m,
  field<vec3>::accessor<ro, ro> u_a,
  field<double>::accessor<ro, ro> p_a,
  flecsi::util::id zslice) {
  {
    auto u = m.mdspan<mesh::cells>(u_a);
    std::stringstream ss;
    ss << "VELOCITY:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << u[zslice][j][i] << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto p = m.mdspan<mesh::cells>(p_a);
    std::stringstream ss;
    ss << "PRESSURE:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << p[zslice][j][i] << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
} // print_conserved
} // namespace muscl::tasks::util

#endif // MUSCL_TASKS_UTIL_HH
