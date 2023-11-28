#include "init.hh"

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

muscl::mesh::periodic_axes
muscl::tasks::init::boundaries(single<mesh::bmap>::accessor<wo> bmap_a,
  mesh::boundary_type xlow,
  mesh::boundary_type xhigh,
  mesh::boundary_type ylow,
  mesh::boundary_type yhigh,
  mesh::boundary_type zlow,
  mesh::boundary_type zhigh) {
  auto & bmap = *bmap_a;
  mesh::periodic_axes p{false, false, false};

  if(xlow == mesh::boundary_type::periodic ||
     xhigh == mesh::boundary_type::periodic) {
    p[mesh::x_axis] = true;
    bmap[mesh::x_axis][mesh::low] = mesh::boundary_type::periodic;
    bmap[mesh::x_axis][mesh::high] = mesh::boundary_type::periodic;
  }
  else {
    bmap[mesh::x_axis][mesh::low] = xlow;
    bmap[mesh::x_axis][mesh::high] = xhigh;
  } // if

  if(ylow == mesh::boundary_type::periodic ||
     yhigh == mesh::boundary_type::periodic) {
    p[mesh::y_axis] = true;
    bmap[mesh::y_axis][mesh::low] = mesh::boundary_type::periodic;
    bmap[mesh::y_axis][mesh::high] = mesh::boundary_type::periodic;
  }
  else {
    bmap[mesh::y_axis][mesh::low] = ylow;
    bmap[mesh::y_axis][mesh::high] = yhigh;
  } // if

  if(zlow == mesh::boundary_type::periodic ||
     zhigh == mesh::boundary_type::periodic) {
    p[mesh::z_axis] = true;
    bmap[mesh::z_axis][mesh::low] = mesh::boundary_type::periodic;
    bmap[mesh::z_axis][mesh::high] = mesh::boundary_type::periodic;
  }
  else {
    bmap[mesh::z_axis][mesh::low] = zlow;
    bmap[mesh::z_axis][mesh::high] = zhigh;
  } // if

  return p;
} // boundaries
