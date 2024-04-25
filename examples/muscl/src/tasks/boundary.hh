#ifndef MUSCL_TASKS_BOUNDARY_HH
#define MUSCL_TASKS_BOUNDARY_HH

#include "../state.hh"

namespace muscl::tasks {

template<mesh::axis A, mesh::boundary B>
void
flow(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a) {
  auto r = m.mdcolex<mesh::cells>(r_a);
  auto ru = m.mdcolex<mesh::cells>(ru_a);
  auto rE = m.mdcolex<mesh::cells>(rE_a);

  if constexpr(A == mesh::axis::x_axis) {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          r(0, j, k) = r(2, j, k);
          r(1, j, k) = r(2, j, k);
          ru(0, j, k).x = ru(2, j, k).x;
          ru(1, j, k).x = ru(2, j, k).x;
          ru(0, j, k).y = ru(2, j, k).y;
          ru(1, j, k).y = ru(2, j, k).y;
          ru(0, j, k).z = ru(2, j, k).z;
          ru(1, j, k).z = ru(2, j, k).z;
          rE(0, j, k) = rE(2, j, k);
          rE(1, j, k) = rE(2, j, k);
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t i = m.size<A, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          r(i - 1, j, k) = r(i - 3, j, k);
          r(i - 2, j, k) = r(i - 3, j, k);
          ru(i - 1, j, k).x = ru(i - 3, j, k).x;
          ru(i - 2, j, k).x = ru(i - 3, j, k).x;
          ru(i - 1, j, k).y = ru(i - 3, j, k).y;
          ru(i - 2, j, k).y = ru(i - 3, j, k).y;
          ru(i - 1, j, k).z = ru(i - 3, j, k).z;
          ru(i - 2, j, k).z = ru(i - 3, j, k).z;
          rE(i - 1, j, k) = rE(i - 3, j, k);
          rE(i - 2, j, k) = rE(i - 3, j, k);
        } // for
      } // for
    } // if
  }
  else if constexpr(A == mesh::axis::y_axis) {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, 0, k) = r(i, 2, k);
          r(i, 1, k) = r(i, 2, k);
          ru(i, 0, k).x = ru(i, 2, k).x;
          ru(i, 1, k).x = ru(i, 2, k).x;
          ru(i, 0, k).y = ru(i, 2, k).y;
          ru(i, 1, k).y = ru(i, 2, k).y;
          ru(i, 0, k).z = ru(i, 2, k).z;
          ru(i, 1, k).z = ru(i, 2, k).z;
          rE(i, 0, k) = rE(i, 2, k);
          rE(i, 1, k) = rE(i, 2, k);
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t j = m.size<A, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, j - 1, k) = r(i, j - 3, k);
          r(i, j - 2, k) = r(i, j - 3, k);
          ru(i, j - 1, k).x = ru(i, j - 3, k).x;
          ru(i, j - 2, k).x = ru(i, j - 3, k).x;
          ru(i, j - 1, k).y = ru(i, j - 3, k).y;
          ru(i, j - 2, k).y = ru(i, j - 3, k).y;
          ru(i, j - 1, k).z = ru(i, j - 3, k).z;
          ru(i, j - 2, k).z = ru(i, j - 3, k).z;
          rE(i, j - 1, k) = rE(i, j - 3, k);
          rE(i, j - 2, k) = rE(i, j - 3, k);
        } // for
      } // for
    } // if
  }
  else {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, j, 0) = r(i, j, 2);
          r(i, j, 1) = r(i, j, 2);
          ru(i, j, 0).x = ru(i, j, 2).x;
          ru(i, j, 1).x = ru(i, j, 2).x;
          ru(i, j, 0).y = ru(i, j, 2).y;
          ru(i, j, 1).y = ru(i, j, 2).y;
          ru(i, j, 0).z = ru(i, j, 2).z;
          ru(i, j, 1).z = ru(i, j, 2).z;
          rE(i, j, 0) = rE(i, j, 2);
          rE(i, j, 1) = rE(i, j, 2);
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t k = m.size<A, mesh::all>();
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, j, k - 1) = r(i, j, k - 3);
          r(i, j, k - 2) = r(i, j, k - 3);
          ru(i, j, k - 1).x = ru(i, j, k - 3).x;
          ru(i, j, k - 2).x = ru(i, j, k - 3).x;
          ru(i, j, k - 1).y = ru(i, j, k - 3).y;
          ru(i, j, k - 2).y = ru(i, j, k - 3).y;
          ru(i, j, k - 1).z = ru(i, j, k - 3).z;
          ru(i, j, k - 2).z = ru(i, j, k - 3).z;
          rE(i, j, k - 1) = rE(i, j, k - 3);
          rE(i, j, k - 2) = rE(i, j, k - 3);
        } // for
      } // for
    } // if
  } // if
} // flow

template<mesh::axis A, mesh::boundary B>
void
reflecting(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a) {
  auto r = m.mdcolex<mesh::cells>(r_a);
  auto ru = m.mdcolex<mesh::cells>(ru_a);
  auto rE = m.mdcolex<mesh::cells>(rE_a);

  if constexpr(A == mesh::axis::x_axis) {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          r(0, j, k) = r(3, j, k);
          r(1, j, k) = r(2, j, k);
          ru(0, j, k).x = ru(3, j, k).x;
          ru(1, j, k).x = ru(2, j, k).x;
          ru(0, j, k).y = ru(3, j, k).y;
          ru(1, j, k).y = ru(2, j, k).y;
          ru(0, j, k).z = ru(3, j, k).z;
          ru(1, j, k).z = ru(2, j, k).z;
          rE(0, j, k) = rE(3, j, k);
          rE(1, j, k) = rE(2, j, k);
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t i = m.size<A, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
          r(i - 1, j, k) = r(i - 4, j, k);
          r(i - 2, j, k) = r(i - 3, j, k);
          ru(i - 1, j, k).x = ru(i - 4, j, k).x;
          ru(i - 2, j, k).x = ru(i - 3, j, k).x;
          ru(i - 1, j, k).y = ru(i - 4, j, k).y;
          ru(i - 2, j, k).y = ru(i - 3, j, k).y;
          ru(i - 1, j, k).z = ru(i - 4, j, k).z;
          ru(i - 2, j, k).z = ru(i - 3, j, k).z;
          rE(i - 1, j, k) = rE(i - 4, j, k);
          rE(i - 2, j, k) = rE(i - 3, j, k);
        } // for
      } // for
    } // if
  }
  else if constexpr(A == mesh::axis::y_axis) {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, 0, k) = r(i, 3, k);
          r(i, 1, k) = r(i, 2, k);
          ru(i, 0, k).x = ru(i, 3, k).x;
          ru(i, 1, k).x = ru(i, 2, k).x;
          ru(i, 0, k).y = ru(i, 3, k).y;
          ru(i, 1, k).y = ru(i, 2, k).y;
          ru(i, 0, k).z = ru(i, 3, k).z;
          ru(i, 1, k).z = ru(i, 2, k).z;
          rE(i, 0, k) = rE(i, 3, k);
          rE(i, 1, k) = rE(i, 2, k);
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t j = m.size<A, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, j - 1, k) = r(i, j - 4, k);
          r(i, j - 2, k) = r(i, j - 3, k);
          ru(i, j - 1, k).x = ru(i, j - 4, k).x;
          ru(i, j - 2, k).x = ru(i, j - 3, k).x;
          ru(i, j - 1, k).y = ru(i, j - 4, k).y;
          ru(i, j - 2, k).y = ru(i, j - 3, k).y;
          ru(i, j - 1, k).z = ru(i, j - 4, k).z;
          ru(i, j - 2, k).z = ru(i, j - 3, k).z;
          rE(i, j - 1, k) = rE(i, j - 4, k);
          rE(i, j - 2, k) = rE(i, j - 3, k);
        } // for
      } // for
    } // if
  }
  else {
    if(B == mesh::boundary::low && m.is_low<A>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, j, 0) = r(i, j, 3);
          r(i, j, 1) = r(i, j, 2);
          ru(i, j, 0).x = ru(i, j, 3).x;
          ru(i, j, 1).x = ru(i, j, 2).x;
          ru(i, j, 0).y = ru(i, j, 3).y;
          ru(i, j, 1).y = ru(i, j, 2).y;
          ru(i, j, 0).z = ru(i, j, 3).z;
          ru(i, j, 1).z = ru(i, j, 2).z;
          rE(i, j, 0) = rE(i, j, 3);
          rE(i, j, 1) = rE(i, j, 2);
        } // for
      } // for
    }
    else if(B == mesh::boundary::high && m.is_high<A>()) {
      const std::size_t k = m.size<A, mesh::all>();
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r(i, j, k - 1) = r(i, j, k - 4);
          r(i, j, k - 2) = r(i, j, k - 3);
          ru(i, j, k - 1).x = ru(i, j, k - 4).x;
          ru(i, j, k - 2).x = ru(i, j, k - 3).x;
          ru(i, j, k - 1).y = ru(i, j, k - 4).y;
          ru(i, j, k - 2).y = ru(i, j, k - 3).y;
          ru(i, j, k - 1).z = ru(i, j, k - 4).z;
          ru(i, j, k - 2).z = ru(i, j, k - 3).z;
          rE(i, j, k - 1) = rE(i, j, k - 4);
          rE(i, j, k - 2) = rE(i, j, k - 3);
        } // for
      } // for
    } // if
  } // if
} // reflecting

inline void
apply_boundaries(mesh::accessor<ro> m,
  single<mesh::bmap>::accessor<ro> bmap_a,
  field<double>::accessor<rw, ro> ra,
  field<vec3>::accessor<rw, ro> rua,
  field<double>::accessor<rw, ro> rEa) {
  const mesh::bmap & bm = *bmap_a;

  {
    const auto xlow = bm[0][0];
    if(xlow == mesh::boundary_type::flow) {
      flow<mesh::x_axis, mesh::low>(m, ra, rua, rEa);
    }
    else if(xlow == mesh::boundary_type::reflecting) {
      reflecting<mesh::x_axis, mesh::low>(m, ra, rua, rEa);
    }
    else { /* periodic */
    } // if
  } // scope

  {
    const auto xhigh = bm[0][1];
    if(xhigh == mesh::boundary_type::flow) {
      flow<mesh::x_axis, mesh::high>(m, ra, rua, rEa);
    }
    else if(xhigh == mesh::boundary_type::reflecting) {
      reflecting<mesh::x_axis, mesh::high>(m, ra, rua, rEa);
    }
    else { /* periodic */
    } // if
  } // scope

  {
    const auto ylow = bm[1][0];
    if(ylow == mesh::boundary_type::flow) {
      flow<mesh::y_axis, mesh::low>(m, ra, rua, rEa);
    }
    else if(ylow == mesh::boundary_type::reflecting) {
      reflecting<mesh::y_axis, mesh::low>(m, ra, rua, rEa);
    }
    else { /* periodic */
    } // if
  } // scope

  {
    const auto yhigh = bm[1][1];
    if(yhigh == mesh::boundary_type::flow) {
      flow<mesh::y_axis, mesh::high>(m, ra, rua, rEa);
    }
    else if(yhigh == mesh::boundary_type::reflecting) {
      reflecting<mesh::y_axis, mesh::high>(m, ra, rua, rEa);
    }
    else { /* periodic */
    } // if
  } // scope

  {
    const auto zlow = bm[2][0];
    if(zlow == mesh::boundary_type::flow) {
      flow<mesh::z_axis, mesh::low>(m, ra, rua, rEa);
    }
    else if(zlow == mesh::boundary_type::reflecting) {
      reflecting<mesh::z_axis, mesh::low>(m, ra, rua, rEa);
    }
    else { /* periodic */
    } // if
  } // scope

  {
    const auto zhigh = bm[2][1];
    if(zhigh == mesh::boundary_type::flow) {
      flow<mesh::z_axis, mesh::high>(m, ra, rua, rEa);
    }
    else if(zhigh == mesh::boundary_type::reflecting) {
      reflecting<mesh::z_axis, mesh::high>(m, ra, rua, rEa);
    }
    else { /* periodic */
    } // if
  } // scope
} // apply_boundaries

} // namespace muscl::tasks

#endif // MUSCL_TASKS_BOUNDARY_HH
