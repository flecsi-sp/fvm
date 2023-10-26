#ifndef MUSCL_TASKS_UTIL_HH
#define MUSCL_TASKS_UTIL_HH

#include "../state.hh"

#include <sstream>

namespace muscl::tasks::util {

template<mesh::domain DM>
inline void
print_conserved(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  flecsi::util::id zslice) {
  {
    auto r = m.mdcolex<mesh::cells>(r_a);
    std::stringstream ss;
    ss << "DENSITY:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM, true>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << r(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
#if 0
  {
    auto ru = m.mdcolex<mesh::cells>(ru_a);
    std::stringstream ss;
    ss << "MOMENTUM:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM, true>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << ru(i,j,zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto rE = m.mdcolex<mesh::cells>(rE_a);
    std::stringstream ss;
    ss << "TOTAL ENERGY:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM, true>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << rE(i,j,zslice) << " ";
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
    auto u = m.mdcolex<mesh::cells>(u_a);
    std::stringstream ss;
    ss << "VELOCITY:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << u(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto p = m.mdcolex<mesh::cells>(p_a);
    std::stringstream ss;
    ss << "PRESSURE:" << std::endl;
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << p(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
} // print_conserved
} // namespace muscl::tasks::util

#endif // MUSCL_TASKS_UTIL_HH
