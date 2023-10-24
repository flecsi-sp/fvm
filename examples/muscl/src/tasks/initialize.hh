#ifndef MUSCL_TASKS_INITIALIZE_HH
#define MUSCL_TASKS_INITIALIZE_HH

#include "../types.hh"
#include "../utils.hh"

namespace muscl {

/*----------------------------------------------------------------------------*
  Shock tube initialization types.
 *----------------------------------------------------------------------------*/

struct sod {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double vL = 0.0;
  static constexpr double wL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.125;
  static constexpr double uR = 0.0;
  static constexpr double vR = 0.0;
  static constexpr double wR = 0.0;
  static constexpr double pR = 0.1;

  static constexpr double x0 = 0.5;
}; // struct sod

struct rankine_hugoniot {
  static constexpr double rL = 2.299156e-01;
  static constexpr double uL = 1.270775e+00;
  static constexpr double vL = 0.0;
  static constexpr double wL = 0.0;
  static constexpr double pL = 6.700238e-01;

  static constexpr double rR = 1.250000e-01;
  static constexpr double uR = 0.000000e+00;
  static constexpr double vR = 0.0;
  static constexpr double wR = 0.0;
  static constexpr double pR = 2.276638e-01;

  static constexpr double x0 = 0.5;
}; // struct rankine_hugoniot

namespace tasks::init {

/*----------------------------------------------------------------------------*
  Ideal gas parameter.
 *----------------------------------------------------------------------------*/

void gamma(single<double>::accessor<wo> gamma_a, double g);

/*----------------------------------------------------------------------------*
  Fake initialization tasks to avoid legion errors.
 *----------------------------------------------------------------------------*/

inline void
touch(mesh::accessor<ro> m,
  field<double>::accessor<wo, wo> r_a,
  field<vec3>::accessor<wo, wo> ru_a,
  field<double>::accessor<wo, wo> rE_a,
  field<vec3>::accessor<wo, wo> u_a,
  field<double>::accessor<wo, wo> p_a,
  field<double>::accessor<wo, wo> q_a,
  field<vec3>::accessor<wo, wo> qu_a,
  field<double>::accessor<wo, wo> qE_a,
  field<double>::accessor<wo, wo> dr_ds_a,
  field<vec3>::accessor<wo, wo> du_ds_a,
  field<double>::accessor<wo, wo> dp_ds_a,
  field<double>::accessor<wo, wo> rTail_a,
  field<vec3>::accessor<wo, wo> ruTail_a,
  field<double>::accessor<wo, wo> rETail_a,
  field<vec3>::accessor<wo, wo> uTail_a,
  field<double>::accessor<wo, wo> pTail_a,
  field<double>::accessor<wo, wo> rHead_a,
  field<vec3>::accessor<wo, wo> ruHead_a,
  field<double>::accessor<wo, wo> rEHead_a,
  field<vec3>::accessor<wo, wo> uHead_a,
  field<double>::accessor<wo, wo> pHead_a,
  field<double>::accessor<wo, wo> rF_a,
  field<vec3>::accessor<wo, wo> ruF_a,
  field<double>::accessor<wo, wo> rEF_a) {} // touch

inline void
touch1(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  field<vec3>::accessor<rw, ro> u_a,
  field<double>::accessor<rw, ro> p_a,
  field<double>::accessor<rw, ro> q_a,
  field<vec3>::accessor<rw, ro> qu_a,
  field<double>::accessor<rw, ro> qE_a,
  field<double>::accessor<rw, ro> dr_ds_a,
  field<vec3>::accessor<rw, ro> du_ds_a,
  field<double>::accessor<rw, ro> dp_ds_a,
  field<double>::accessor<rw, ro> rTail_a,
  field<vec3>::accessor<rw, ro> ruTail_a,
  field<double>::accessor<rw, ro> rETail_a,
  field<vec3>::accessor<rw, ro> uTail_a,
  field<double>::accessor<rw, ro> pTail_a,
  field<double>::accessor<rw, ro> rHead_a,
  field<vec3>::accessor<rw, ro> ruHead_a,
  field<double>::accessor<rw, ro> rEHead_a,
  field<vec3>::accessor<rw, ro> uHead_a,
  field<double>::accessor<rw, ro> pHead_a,
  field<double>::accessor<rw, ro> rF_a,
  field<vec3>::accessor<rw, ro> ruF_a,
  field<double>::accessor<rw, ro> rEF_a) {} // touch

/*----------------------------------------------------------------------------*
  Boundary input translation.
 *----------------------------------------------------------------------------*/

mesh::periodic_axes boundaries(single<mesh::bmap>::accessor<wo> bmap_a,
  mesh::boundary_type xlow,
  mesh::boundary_type xhigh,
  mesh::boundary_type ylow,
  mesh::boundary_type yhigh,
  mesh::boundary_type zlow,
  mesh::boundary_type zhigh);

/*----------------------------------------------------------------------------*
  Shock Tube.
 *----------------------------------------------------------------------------*/

template<typename T>
void
shock(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  single<double>::accessor<ro> gamma_a) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto ru = m.mdspan<mesh::cells>(ru_a);
  auto rE = m.mdspan<mesh::cells>(rE_a);
  auto const gamma = *gamma_a;
  const double mult = 1.0 / (gamma - 1.0);

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        const auto x = m.head<mesh::x_axis>(i);

        if(x < T::x0) {
          r[k][j][i] = T::rL;
          ru[k][j][i].x = T::rL * T::uL;
          ru[k][j][i].y = T::rL * T::vL;
          ru[k][j][i].z = T::rL * T::wL;
          rE[k][j][i] =
            mult * T::pL +
            0.5 * T::rL *
              (utils::sqr(T::uL) + utils::sqr(T::vL) + utils::sqr(T::wL));
        }
        else {
          r[k][j][i] = T::rR;
          ru[k][j][i].x = T::rR * T::uR;
          ru[k][j][i].y = T::rR * T::vR;
          ru[k][j][i].z = T::rR * T::wR;
          rE[k][j][i] =
            mult * T::pR +
            0.5 * T::rR *
              (utils::sqr(T::uR) + utils::sqr(T::vR) + utils::sqr(T::wR));
        } // if
      } // for
    } // for
  } // for
} // shock

} // namespace tasks::init
} // namespace muscl

#endif // MUSCL_TASKS_INITIALIZE_HH
