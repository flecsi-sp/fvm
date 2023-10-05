#ifndef MUSCL_TASKS_ADVANCE_HH
#define MUSCL_TASKS_ADVANCE_HH

#include "../state.hh"
#include "../utils.hh"

namespace muscl::tasks {
double dtmin(mesh::accessor<ro> m, single<velocity>::accessor<ro> lmax);

template<typename L>
void
advance(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<velocity>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  field<velocity>::accessor<rw, ro> u_a,
  field<double>::accessor<rw, ro> p_a,
  field<double>::accessor<rw, ro> q_a,
  field<velocity>::accessor<rw, ro> qu_a,
  field<double>::accessor<rw, ro> qE_a,
  field<double>::accessor<rw, ro> dr_ds_a,
  field<velocity>::accessor<rw, ro> du_ds_a,
  field<double>::accessor<rw, ro> dp_ds_a,
  field<double>::accessor<rw, ro> rTail_a,
  field<velocity>::accessor<rw, ro> ruTail_a,
  field<double>::accessor<rw, ro> rETail_a,
  field<velocity>::accessor<rw, ro> uTail_a,
  field<double>::accessor<rw, ro> pTail_a,
  field<double>::accessor<rw, ro> rHead_a,
  field<velocity>::accessor<rw, ro> ruHead_a,
  field<double>::accessor<rw, ro> rEHead_a,
  field<velocity>::accessor<rw, ro> uHead_a,
  field<double>::accessor<rw, ro> pHead_a,
  field<double>::accessor<rw, ro> rF_a,
  field<velocity>::accessor<rw, ro> ruF_a,
  field<double>::accessor<rw, ro> rEF_a,
  single<double>::accessor<ro> gamma_a,
  double dt) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto ru = m.mdspan<mesh::cells>(ru_a);
  auto rE = m.mdspan<mesh::cells>(rE_a);
  auto u = m.mdspan<mesh::cells>(u_a);
  auto p = m.mdspan<mesh::cells>(p_a);
  auto q = m.mdspan<mesh::cells>(q_a);
  auto qu = m.mdspan<mesh::cells>(qu_a);
  auto qE = m.mdspan<mesh::cells>(qE_a);
  auto dr_ds = m.mdspan<mesh::cells>(dr_ds_a);
  auto du_ds = m.mdspan<mesh::cells>(du_ds_a);
  auto dp_ds = m.mdspan<mesh::cells>(dp_ds_a);
  auto rTail = m.mdspan<mesh::cells>(rTail_a);
  auto ruTail = m.mdspan<mesh::cells>(ruTail_a);
  auto rETail = m.mdspan<mesh::cells>(rETail_a);
  auto uTail = m.mdspan<mesh::cells>(uTail_a);
  auto pTail = m.mdspan<mesh::cells>(pTail_a);
  auto rHead = m.mdspan<mesh::cells>(rHead_a);
  auto ruHead = m.mdspan<mesh::cells>(ruHead_a);
  auto rEHead = m.mdspan<mesh::cells>(rEHead_a);
  auto uHead = m.mdspan<mesh::cells>(uHead_a);
  auto pHead = m.mdspan<mesh::cells>(pHead_a);
  auto rF = m.mdspan<mesh::cells>(rF_a);
  auto ruF = m.mdspan<mesh::cells>(ruF_a);
  auto rEF = m.mdspan<mesh::cells>(rEF_a);
  auto const gamma = *gamma_a;

  /*--------------------------------------------------------------------------*
    Compute X slopes.
   *--------------------------------------------------------------------------*/

  const double xslope_factor{1.0 / m.delta<mesh::x_axis>()};

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {
        dr_ds[k][j][i] =
          xslope_factor * L::limit(r[k][j][i - 1], r[k][j][i], r[k][j][i + 1]);
        du_ds[k][j][i].x =
          xslope_factor *
          L::limit(u[k][j][i - 1].x, u[k][j][i].x, u[k][j][i + 1].x);
        du_ds[k][j][i].y =
          xslope_factor *
          L::limit(u[k][j][i - 1].y, u[k][j][i].y, u[k][j][i + 1].y);
        du_ds[k][j][i].z =
          xslope_factor *
          L::limit(u[k][j][i - 1].z, u[k][j][i].z, u[k][j][i + 1].z);
        dp_ds[k][j][i] =
          xslope_factor * L::limit(p[k][j][i - 1], p[k][j][i], p[k][j][i + 1]);
      } // for
    } // for
  } // for

  /*--------------------------------------------------------------------------*
    Reconstruct faces.
   *--------------------------------------------------------------------------*/

  double xextrap_factor{m.delta<mesh::x_axis>() / 2.0};
  double p2rEFactor{1.0 / (gamma - 1.0)};

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {
        // Primitive quantities
        rTail[k][j][i] = r[k][j][i] + xextrap_factor * dr_ds[k][j][i];
        rHead[k][j][i] = r[k][j][i] - xextrap_factor * dr_ds[k][j][i];
        uTail[k][j][i].x = u[k][j][i].x + xextrap_factor * du_ds[k][j][i].x;
        uHead[k][j][i].x = u[k][j][i].x - xextrap_factor * du_ds[k][j][i].x;
        uTail[k][j][i].y = u[k][j][i].y + xextrap_factor * du_ds[k][j][i].y;
        uHead[k][j][i].y = u[k][j][i].y - xextrap_factor * du_ds[k][j][i].y;
        uTail[k][j][i].z = u[k][j][i].z + xextrap_factor * du_ds[k][j][i].z;
        uHead[k][j][i].z = u[k][j][i].z - xextrap_factor * du_ds[k][j][i].z;
        pTail[k][j][i] = p[k][j][i] + xextrap_factor * dp_ds[k][j][i];
        pHead[k][j][i] = p[k][j][i] - xextrap_factor * dp_ds[k][j][i];

        // Conserved quantities
        ruTail[k][j][i].x = rTail[k][j][i] * uTail[k][j][i].x;
        ruHead[k][j][i].x = rHead[k][j][i] * uHead[k][j][i].x;
        ruTail[k][j][i].y = rTail[k][j][i] * uTail[k][j][i].y;
        ruHead[k][j][i].y = rHead[k][j][i] * uHead[k][j][i].y;
        ruTail[k][j][i].z = rTail[k][j][i] * uTail[k][j][i].z;
        ruHead[k][j][i].z = rHead[k][j][i] * uHead[k][j][i].z;
        rETail[k][j][i] =
          p2rEFactor * pTail[k][j][i] +
          0.5 * rTail[k][j][i] *
            (utils::sqr(uTail[k][j][i].x) + utils::sqr(uTail[k][j][i].y) +
              utils::sqr(uTail[k][j][i].z));
        rEHead[k][j][i] =
          p2rEFactor * pHead[k][j][i] +
          0.5 * rHead[k][j][i] *
            (utils::sqr(uHead[k][j][i].x) + utils::sqr(uHead[k][j][i].y) +
              utils::sqr(uHead[k][j][i].z));
      } // for
    } // for
  } // for

  /*--------------------------------------------------------------------------*
    Predictor step: approximate flux by averaging the governing flux functions
    at either side of the cell.
   *--------------------------------------------------------------------------*/

  double xcourant{dt / (2.0 * m.delta<mesh::x_axis>())};

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {
        // Density
        qu[k][j][i].x =
          r[k][j][i] - xcourant * (ruTail[k][j][i].x - ruHead[k][j][i].x);

        // X Momentum
        const double f_qu_Tail =
          ruTail[k][j][i].x * uTail[k][j][i].x + pTail[k][j][i];
        const double f_qu_Head =
          ruHead[k][j][i].x * uHead[k][j][i].x + pHead[k][j][i];
        qu[k][j][i].x = ru[k][j][i].x - xcourant * (f_qu_Tail - f_qu_Head);

        // Y Momentum
        const double f_qv_Tail =
          rTail[k][j][i] * uTail[k][j][i].x * uTail[k][j][i].y;
        const double f_qv_Head =
          rHead[k][j][i] * uHead[k][j][i].x * uHead[k][j][i].y;
        qu[k][j][i].y = ru[k][j][i].y - xcourant * (f_qv_Tail - f_qv_Head);

        // Z Momentum
        const double f_qw_Tail =
          rTail[k][j][i] * uTail[k][j][i].z * uTail[k][j][i].x;
        const double f_qw_Head =
          rHead[k][j][i] * uHead[k][j][i].z * uHead[k][j][i].x;
        qu[k][j][i].z = ru[k][j][i].z - xcourant * (f_qw_Tail - f_qw_Head);

        // Total Energy
        const double f_qE_Tail =
          (rETail[k][j][i] + pTail[k][j][i]) * uTail[k][j][i].x;
        const double f_qE_Head =
          (rEHead[k][j][i] + pHead[k][j][i]) * uHead[k][j][i].x;
        qE[k][j][i] = rE[k][j][i] - xcourant * (f_qE_Tail - f_qE_Head);

        // Primitives
        u[k][j][i].x = qu[k][j][i].x / q[k][j][i];
        u[k][j][i].y = qu[k][j][i].y / q[k][j][i];
        u[k][j][i].z = qu[k][j][i].z / q[k][j][i];
        p[k][j][i] =
          (gamma - 1.0) * (qE[k][j][i] - 0.5 * q[k][j][i] *
                                           (utils::sqr(u[k][j][i].x) +
                                             utils::sqr(u[k][j][i].y) +
                                             utils::sqr(u[k][j][i].z)));
      } // for
    } // for
  } // for

  /*--------------------------------------------------------------------------*
    Reconstruct faces from intermediates.
   *--------------------------------------------------------------------------*/

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::corrector>()) {
      } // for
    } // for
  } // for

} // advance

} // namespace muscl::tasks

#endif // MUSCL_TASKS_ADVANCE_HH
