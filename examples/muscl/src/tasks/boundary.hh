#ifndef MUSCL_TASKS_BOUNDARY_HH
#define MUSCL_TASKS_BOUNDARY_HH

#include "../boundary.hh"

namespace muscl::tasks {

inline void
apply_boundaries(mesh::accessor<ro> m,
  single<mesh::bmap>::accessor<ro> bmap_a,
  field<double>::accessor<rw, ro> ra,
  field<vec3>::accessor<rw, ro> rua,
  field<double>::accessor<rw, ro> rEa) {
  const mesh::bmap & bm = *bmap_a;

  {
    const auto xlow = bm[0][0];
    if(xlow == mesh::boundary_type::inflow ||
       xlow == mesh::boundary_type::outflow) {
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
    if(xhigh == mesh::boundary_type::inflow ||
       xhigh == mesh::boundary_type::outflow) {
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
    if(ylow == mesh::boundary_type::inflow ||
       ylow == mesh::boundary_type::outflow) {
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
    if(yhigh == mesh::boundary_type::inflow ||
       yhigh == mesh::boundary_type::outflow) {
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
    if(zlow == mesh::boundary_type::inflow ||
       zlow == mesh::boundary_type::outflow) {
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
    if(zhigh == mesh::boundary_type::inflow ||
       zhigh == mesh::boundary_type::outflow) {
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
