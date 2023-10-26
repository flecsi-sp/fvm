#ifndef TEST_TASKS_UTILS_HH
#define TEST_TASKS_UTILS_HH

#include "../state.hh"

namespace test::tasks::util {

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
  for(auto k : m.cells<mesh::z_axis, DM>()) {
    for(auto j : m.cells<mesh::y_axis, DM>()) {
      for(auto i : m.cells<mesh::x_axis, DM>()) {
        std::stringstream ss;
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
        flog(info) << ss.str() << std::endl;
      } // for
    } // for
  } // for

} // cell_info

template<mesh::domain DM>
inline void
field_info(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> f_a,
  flecsi::util::id zslice) {
  auto f = m.mdcolex<mesh::cells>(f_a);
  std::stringstream ss;
  for(auto j : m.cells<mesh::y_axis, DM, true>()) {
    for(auto i : m.cells<mesh::x_axis, DM>()) {
      ss << f(i, j, zslice) << " ";
    } // for
    ss << std::endl;
  } // for

  flog(info) << ss.str() << std::endl;
} // field_info

} // namespace test::tasks::util

#endif // TEST_TASKS_UTILS_HH
