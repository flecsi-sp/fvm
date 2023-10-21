#include "initialize.hh"

/*----------------------------------------------------------------------------*
  Gamma.
 *----------------------------------------------------------------------------*/

void
muscl::tasks::init::gamma(single<double>::accessor<wo> gamma_a, double g) {
  (*gamma_a) = g;
} // gamma

/*----------------------------------------------------------------------------*
  Boundary initialization.
 *----------------------------------------------------------------------------*/

void
muscl::tasks::init::boundaries(single<mesh::bmap>::accessor<wo> bmap_a,
  mesh::boundary_type xlow,
  mesh::boundary_type xhigh,
  mesh::boundary_type ylow,
  mesh::boundary_type yhigh,
  mesh::boundary_type zlow,
  mesh::boundary_type zhigh) {
  auto & bmap = *bmap_a;
  bmap[mesh::x_axis][mesh::low] = xlow;
  bmap[mesh::x_axis][mesh::high] = xhigh;
  bmap[mesh::y_axis][mesh::low] = ylow;
  bmap[mesh::y_axis][mesh::high] = yhigh;
  bmap[mesh::z_axis][mesh::low] = zlow;
  bmap[mesh::z_axis][mesh::high] = zhigh;
} // boundaries
