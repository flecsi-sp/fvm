#include "initialize.hh"

using namespace test;

mesh::periodic_axes
tasks::init::boundaries(single<mesh::bmap>::accessor<wo> bmap_a,
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

void
test::tasks::init::monotonic(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a) {
  auto f = m.mdspan<mesh::cells>(f_a);

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    auto kg = m.global_id<mesh::z_axis>(k);
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      auto jg = m.global_id<mesh::y_axis>(j);
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        auto ig = m.global_id<mesh::y_axis>(j);
        f[k][j][i] = m.global_id<mesh::x_axis>(i) +
                     jg * m.size<mesh::x_axis, mesh::domain::global>() +
                     kg * m.size<mesh::x_axis, mesh::domain::global>() *
                       m.size<mesh::y_axis, mesh::domain::global>();
      } // for
    } // for
  } // for
} // monotone

void
test::tasks::init::color(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a) {
  auto f = m.mdspan<mesh::cells>(f_a);

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        f[k][j][i] = flecsi::color();
      } // for
    } // for
  } // for
} // color
