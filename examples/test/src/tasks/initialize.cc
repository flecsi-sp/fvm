#include "initialize.hh"

using namespace test;

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
