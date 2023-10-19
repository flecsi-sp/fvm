#include "hydro.hh"

#include <limits>

void
muscl::tasks::hydro::update_primitives(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  field<vec3>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<double>::accessor<ro> gamma_a) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto ru = m.mdspan<mesh::cells>(ru_a);
  auto rE = m.mdspan<mesh::cells>(rE_a);
  auto u = m.mdspan<mesh::cells>(u_a);
  auto p = m.mdspan<mesh::cells>(p_a);
  auto const gamma = *gamma_a;

  // Initialize primitive quantities.
  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        u[k][j][i].x = ru[k][j][i].x / r[k][j][i];
        u[k][j][i].y = ru[k][j][i].y / r[k][j][i];
        u[k][j][i].z = ru[k][j][i].z / r[k][j][i];
        p[k][j][i] =
          (gamma - 1.0) * (rE[k][j][i] - 0.5 * r[k][j][i] *
                                           (utils::sqr(u[k][j][i].x) +
                                             utils::sqr(u[k][j][i].y) +
                                             utils::sqr(u[k][j][i].z)));
      } // for
    } // for
  } // for
} // update_primitives

void
muscl::tasks::hydro::update_eigenvalues(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<vec3>::accessor<wo> lmax_a,
  single<double>::accessor<ro> gamma_a) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto u = m.mdspan<mesh::cells>(u_a);
  auto p = m.mdspan<mesh::cells>(p_a);
  auto & lmax = *lmax_a;
  auto const gamma = *gamma_a;

  // Iinitialize max eigenvalues for dt.
  lmax.x = std::numeric_limits<double>::min();
  lmax.y = std::numeric_limits<double>::min();
  lmax.z = std::numeric_limits<double>::min();
  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        const double c = std::sqrt(gamma * p[k][j][i] / r[k][j][i]);
        lmax.x = std::max(std::abs(u[k][j][i].x) + c, lmax.x);
        lmax.y = std::max(std::abs(u[k][j][i].y) + c, lmax.y);
        lmax.z = std::max(std::abs(u[k][j][i].z) + c, lmax.z);
      } // for
    } // for
  } // for
} // update_eigenvalues

double
muscl::tasks::hydro::update_dtmin(mesh::accessor<ro> m,
  single<vec3>::accessor<ro> lmax) {
  return std::min(m.delta<mesh::x_axis>() / lmax->x,
    std::min(
      m.delta<mesh::y_axis>() / lmax->y, m.delta<mesh::z_axis>() / lmax->z));
} // update_dtmin
