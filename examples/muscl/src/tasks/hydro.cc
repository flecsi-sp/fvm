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
  auto r = m.mdcolex<mesh::cells>(r_a);
  auto ru = m.mdcolex<mesh::cells>(ru_a);
  auto rE = m.mdcolex<mesh::cells>(rE_a);
  auto u = m.mdcolex<mesh::cells>(u_a);
  auto p = m.mdcolex<mesh::cells>(p_a);
  auto const gamma = *gamma_a;

  // Initialize primitive quantities.
  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        u(i, j, k).x = ru(i, j, k).x / r(i, j, k);
        u(i, j, k).y = ru(i, j, k).y / r(i, j, k);
        u(i, j, k).z = ru(i, j, k).z / r(i, j, k);
        p(i, j, k) =
          (gamma - 1.0) * (rE(i, j, k) - 0.5 * r(i, j, k) *
                                           (utils::sqr(u(i, j, k).x) +
                                             utils::sqr(u(i, j, k).y) +
                                             utils::sqr(u(i, j, k).z)));
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
  auto r = m.mdcolex<mesh::cells>(r_a);
  auto u = m.mdcolex<mesh::cells>(u_a);
  auto p = m.mdcolex<mesh::cells>(p_a);
  auto & lmax = *lmax_a;
  auto const gamma = *gamma_a;

  // Iinitialize max eigenvalues for dt.
  lmax.x = std::numeric_limits<double>::min();
  lmax.y = std::numeric_limits<double>::min();
  lmax.z = std::numeric_limits<double>::min();
  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        const double c = std::sqrt(gamma * p(i, j, k) / r(i, j, k));
        lmax.x = std::max(std::abs(u(i, j, k).x) + c, lmax.x);
        lmax.y = std::max(std::abs(u(i, j, k).y) + c, lmax.y);
        lmax.z = std::max(std::abs(u(i, j, k).z) + c, lmax.z);
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
