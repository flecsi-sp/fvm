#ifndef TEST_TASKS_INIT_HH
#define TEST_TASKS_INIT_HH

#include "../state.hh"

namespace test::tasks::init {
void monotonic(mesh::accessor<ro> m, field<double>::accessor<rw, ro> f_a);

void color(mesh::accessor<ro> m, field<double>::accessor<rw, ro> f_a);

mesh::periodic_axes boundaries(single<mesh::bmap>::accessor<wo> bmap_a,
  mesh::boundary_type xlow,
  mesh::boundary_type xhigh,
  mesh::boundary_type ylow,
  mesh::boundary_type yhigh,
  mesh::boundary_type zlow,
  mesh::boundary_type zhigh);

template<mesh::axis A, mesh::boundary B>
inline void
boundary(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a,
  single<mesh::bmap>::accessor<ro> bmap_a) {
  auto f = m.mdcolex<mesh::cells>(f_a);
  auto & bmap = *bmap_a;

  if(A == mesh::axis::x_axis &&
     !(bmap[mesh::x_axis][mesh::low] == mesh::boundary_type::periodic)) {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          f(0, j, k) = -1.0;
          f(1, j, k) = -1.0;
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t i = m.size<mesh::x_axis, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          f(i - 1, j, k) = -2.0;
          f(i - 2, j, k) = -2.0;
        } // for
      } // for
    } // if
  }
  else if(A == mesh::axis::y_axis &&
          !(bmap[mesh::y_axis][mesh::low] == mesh::boundary_type::periodic)) {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f(i, 0, k) = -3.0;
          f(i, 1, k) = -3.0;
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t j = m.size<A, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f(i, j - 1, k) = -4.0;
          f(i, j - 2, k) = -4.0;
        } // for
      } // for
    } // if
  }
  else if(!(bmap[mesh::z_axis][mesh::low] == mesh::boundary_type::periodic)) {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f(i, j, 0) = -5.0;
          f(i, j, 1) = -5.0;
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t k = m.size<A, mesh::all>();
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          f(i, j, k - 1) = -6.0;
          f(i, j, k - 2) = -6.0;
        } // for
      } // for
    } // if
  } // if
} // boundary

} // namespace test::tasks::init

#endif // TEST_TASKS_INIT_HH
