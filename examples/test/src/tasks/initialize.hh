#ifndef TEST_TASKS_INITIALIZE_HH
#define TEST_TASKS_INITIALIZE_HH

#include "../state.hh"

namespace test::tasks::init {
void monotonic(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a);

void color(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a);

template<mesh::axis A, mesh::boundary B>
void
boundary(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a) {
  auto f = m.mdspan<mesh::cells>(f_a);

  if constexpr(A == mesh::axis::x_axis) {
    if (B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          f[k][j][0] = -1.0;
          f[k][j][1] = -1.0;
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t i = m.size<mesh::x_axis, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          f[k][j][i - 1] = -1.0;
          f[k][j][i - 2] = -1.0;
        } // for
      } // for
    } // if
  }
  else if constexpr(A == mesh::axis::y_axis) {
    if (B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f[k][0][i] = -2.0;
          f[k][1][i] = -2.0;
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t j = m.size<A, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f[k][j - 1][i] = -2.0;
          f[k][j - 2][i] = -2.0;
        } // for
      } // for
    } // if
  }
  else {
    if (B == mesh::boundary::low && m.is_low<A>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f[0][j][i] = -3.0;
          f[1][j][i] = -3.0;
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t k = m.size<A, mesh::all>();
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f[k - 1][j][i] = -3.0;
          f[k - 2][j][i] = -3.0;
        } // for
      } // for
    } // if
  } // if
} // boundary

} // namespace test::tasks::init

#endif // TEST_TASKS_INITIALIZE_HH
