#ifndef MUSCL_TASKS_UTIL_HH
#define MUSCL_TASKS_UTIL_HH

#include "../state.hh"

#include <sstream>

namespace muscl::tasks::util {

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
print(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  flecsi::util::id zslice) {
  {
    auto r = m.mdspan<mesh::cells>(r_a);
    std::stringstream ss;
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << r[zslice][j][i] << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto ru = m.mdspan<mesh::cells>(ru_a);
    std::stringstream ss;
    for(auto j : m.cells<mesh::y_axis, DM>()) {
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
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        ss << rE[zslice][j][i] << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
} // print

} // namespace muscl::tasks::util

#endif // MUSCL_TASKS_UTIL_HH
