#ifndef MUSCL_TASKS_HYDRO_HH
#define MUSCL_TASKS_HYDRO_HH

#include "../state.hh"
#include "../utils.hh"

namespace muscl::tasks::hydro {

void update_primitives(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<velocity>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  field<velocity>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<double>::accessor<ro> gamma_a);

void update_eigenvalues(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<velocity>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<velocity>::accessor<wo> lmax_a,
  single<double>::accessor<ro> gamma_a);

double update_dtmin(mesh::accessor<ro> m, single<velocity>::accessor<ro> lmax);

template<mesh::axis A, typename L>
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

  if constexpr(A == mesh::axis::x_axis) {
    /*------------------------------------------------------------------------*
      Compute slopes.
     *------------------------------------------------------------------------*/

    const double slope_factor{1.0 / m.delta<A>()};

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {
          dr_ds[k][j][i] =
            slope_factor * L::limit(r[k][j][i - 1], r[k][j][i], r[k][j][i + 1]);
          du_ds[k][j][i].x =
            slope_factor *
            L::limit(u[k][j][i - 1].x, u[k][j][i].x, u[k][j][i + 1].x);
          du_ds[k][j][i].y =
            slope_factor *
            L::limit(u[k][j][i - 1].y, u[k][j][i].y, u[k][j][i + 1].y);
          du_ds[k][j][i].z =
            slope_factor *
            L::limit(u[k][j][i - 1].z, u[k][j][i].z, u[k][j][i + 1].z);
          dp_ds[k][j][i] =
            slope_factor * L::limit(p[k][j][i - 1], p[k][j][i], p[k][j][i + 1]);
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Reconstruct faces.
     *------------------------------------------------------------------------*/

    double extrap_factor{m.delta<A>() / 2.0};
    double p2rEFactor{1.0 / (gamma - 1.0)};

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {
          // Primitive quantities
          rTail[k][j][i] = r[k][j][i] + extrap_factor * dr_ds[k][j][i];
          rHead[k][j][i] = r[k][j][i] - extrap_factor * dr_ds[k][j][i];
          uTail[k][j][i].x = u[k][j][i].x + extrap_factor * du_ds[k][j][i].x;
          uHead[k][j][i].x = u[k][j][i].x - extrap_factor * du_ds[k][j][i].x;
          uTail[k][j][i].y = u[k][j][i].y + extrap_factor * du_ds[k][j][i].y;
          uHead[k][j][i].y = u[k][j][i].y - extrap_factor * du_ds[k][j][i].y;
          uTail[k][j][i].z = u[k][j][i].z + extrap_factor * du_ds[k][j][i].z;
          uHead[k][j][i].z = u[k][j][i].z - extrap_factor * du_ds[k][j][i].z;
          pTail[k][j][i] = p[k][j][i] + extrap_factor * dp_ds[k][j][i];
          pHead[k][j][i] = p[k][j][i] - extrap_factor * dp_ds[k][j][i];

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

    /*------------------------------------------------------------------------*
      Predictor step: approximate flux by averaging the governing flux
      functions at either side of the cell.
     *------------------------------------------------------------------------*/

    double courant{dt / (2.0 * m.delta<A>())};

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {
          // Density
          qu[k][j][i].x =
            r[k][j][i] - courant * (ruTail[k][j][i].x - ruHead[k][j][i].x);

          // Momentum
          // ru^2 + p
          // ruv
          // ruw
          const velocity f_qu_Tail{
            ruTail[k][j][i].x * uTail[k][j][i].x + pTail[k][j][i],
            rTail[k][j][i] * uTail[k][j][i].x * uTail[k][j][i].y,
            rTail[k][j][i] * uTail[k][j][i].x * uTail[k][j][i].z};
          const velocity f_qu_Head{
            ruHead[k][j][i].x * uHead[k][j][i].x + pHead[k][j][i],
            rHead[k][j][i] * uHead[k][j][i].x * uHead[k][j][i].y,
            rHead[k][j][i] * uHead[k][j][i].x * uHead[k][j][i].z};
          qu[k][j][i].x = ru[k][j][i].x - courant * (f_qu_Tail.x - f_qu_Head.x);
          qu[k][j][i].y = ru[k][j][i].y - courant * (f_qu_Tail.y - f_qu_Head.y);
          qu[k][j][i].z = ru[k][j][i].z - courant * (f_qu_Tail.z - f_qu_Head.z);

          // Total Energy
          const double f_qE_Tail =
            (rETail[k][j][i] + pTail[k][j][i]) * uTail[k][j][i].x;
          const double f_qE_Head =
            (rEHead[k][j][i] + pHead[k][j][i]) * uHead[k][j][i].x;
          qE[k][j][i] = rE[k][j][i] - courant * (f_qE_Tail - f_qE_Head);

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

    /*------------------------------------------------------------------------*
      Reconstruct faces from intermediates.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::corrector>()) {
          // Primitive quantities
          rTail[k][j][i] = q[k][j][i - 1] + extrap_factor * dr_ds[k][j][i - 1];
          rHead[k][j][i] = q[k][j][i] - extrap_factor * dr_ds[k][j][i];
          uTail[k][j][i].x =
            u[k][j][i - 1].x + extrap_factor * du_ds[k][j][i - 1].x;
          uHead[k][j][i].x = u[k][j][i].x - extrap_factor * du_ds[k][j][i].x;
          uTail[k][j][i].y =
            u[k][j][i - 1].y + extrap_factor * du_ds[k][j][i - 1].y;
          uHead[k][j][i].y = u[k][j][i].y - extrap_factor * du_ds[k][j][i].y;
          uTail[k][j][i].z =
            u[k][j][i - 1].z + extrap_factor * du_ds[k][j][i - 1].z;
          uHead[k][j][i].z = u[k][j][i].z - extrap_factor * du_ds[k][j][i].z;
          pTail[k][j][i] = p[k][j][i - 1] + extrap_factor * dp_ds[k][j][i - 1];
          pHead[k][j][i] = p[k][j][i] - extrap_factor * dp_ds[k][j][i];

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

    /*------------------------------------------------------------------------*
      Compute corrector setp: HLL Riemann solver.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::corrector>()) {
          const double cT = std::sqrt(gamma * pTail[k][j][i] / rTail[k][j][i]);
          const double cH = std::sqrt(gamma * pHead[k][j][i] / rHead[k][j][i]);

          // Update min/max eigenvalues on tail face.
          const double LminT = uTail[k][j][i].x - cT;
          const double LmaxT = uTail[k][j][i].x + cT;

          // Update min/max eigenvalues on head face.
          const double LminH = uHead[k][j][i].x - cH;
          const double LmaxH = uHead[k][j][i].x + cH;

          // Values for HLL.
          const double Lminus{std::min(LminH, std::min(LminT, double{0.0}))};
          const double Lplus{std::max(LmaxH, std::max(LmaxT, double{0.0}))};
          const double Ldiv{1.0 / (Lplus - Lminus)};
          const double Lmult{Lplus * Lminus};
          const double delta_r{rHead[k][j][i] - rTail[k][j][i]};
          const velocity delta_ru{ruHead[k][j][i].x - ruTail[k][j][i].x,
            ruHead[k][j][i].y - ruTail[k][j][i].y,
            ruHead[k][j][i].z - ruTail[k][j][i].z};
          const double delta_rE{rEHead[k][j][i] - rETail[k][j][i]};

          const double f_r_T{ruTail[k][j][i].x};
          const double f_r_H{ruHead[k][j][i].x};
          const velocity f_ru_T{
            ruTail[k][j][i].x * uTail[k][j][i].x + pTail[k][j][i],
            ruTail[k][j][i].x * uTail[k][j][i].y,
            ruTail[k][j][i].x * uTail[k][j][i].z};
          const velocity f_ru_H{
            ruHead[k][j][i].x * uHead[k][j][i].x + pHead[k][j][i],
            ruHead[k][j][i].x * uHead[k][j][i].y,
            ruHead[k][j][i].x * uHead[k][j][i].z};
          const double f_rE_T{
            (rETail[k][j][i] + pTail[k][j][i]) * uTail[k][j][i].x};
          const double f_rE_H{
            (rEHead[k][j][i] + pHead[k][j][i]) * uHead[k][j][i].x};

          rF[k][j][i] =
            Ldiv * (Lplus * f_r_T - Lminus * f_r_H + Lmult * delta_r);
          ruF[k][j][i].x =
            Ldiv * (Lplus * f_ru_T.x - Lminus * f_ru_H.x + Lmult * delta_ru.x);
          ruF[k][j][i].y =
            Ldiv * (Lplus * f_ru_T.y - Lminus * f_ru_H.y + Lmult * delta_ru.y);
          ruF[k][j][i].z =
            Ldiv * (Lplus * f_ru_T.z - Lminus * f_ru_H.z + Lmult * delta_ru.z);
          rEF[k][j][i] =
            Ldiv * (Lplus * f_rE_T - Lminus * f_rE_H + Lmult * delta_rE);
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Update averages (in-place).
     *------------------------------------------------------------------------*/

    const double update_factor = dt / m.delta<A>();
    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r[k][j][i] =
            r[k][j][i] - update_factor * (rF[k][j][i + 1] - rF[k][j][i]);
          ru[k][j][i].x = ru[k][j][i].x -
                          update_factor * (ruF[k][j][i + 1].x - ruF[k][j][i].x);
          ru[k][j][i].y = ru[k][j][i].y -
                          update_factor * (ruF[k][j][i + 1].y - ruF[k][j][i].y);
          ru[k][j][i].z = ru[k][j][i].z -
                          update_factor * (ruF[k][j][i + 1].z - ruF[k][j][i].z);
          rE[k][j][i] =
            rE[k][j][i] - update_factor * (rEF[k][j][i + 1] - rEF[k][j][i]);
        } // for
      } // for
    } // for
  }
  else if constexpr(A == mesh::axis::y_axis) {
    /*------------------------------------------------------------------------*
      Compute slopes.
     *------------------------------------------------------------------------*/

    const double slope_factor{1.0 / m.delta<A>()};

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::predictor>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          dr_ds[k][j][i] =
            slope_factor * L::limit(r[k][j - 1][i], r[k][j][i], r[k][j + 1][i]);
          du_ds[k][j][i].x =
            slope_factor *
            L::limit(u[k][j - 1][i].x, u[k][j][i].x, u[k][j + 1][i].x);
          du_ds[k][j][i].y =
            slope_factor *
            L::limit(u[k][j - 1][i].y, u[k][j][i].y, u[k][j + 1][i].y);
          du_ds[k][j][i].z =
            slope_factor *
            L::limit(u[k][j - 1][i].z, u[k][j][i].z, u[k][j + 1][i].z);
          dp_ds[k][j][i] =
            slope_factor * L::limit(p[k][j - 1][i], p[k][j][i], p[k][j + 1][i]);
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Reconstruct faces.
     *------------------------------------------------------------------------*/

    double extrap_factor{m.delta<A>() / 2.0};
    double p2rEFactor{1.0 / (gamma - 1.0)};

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::predictor>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          // Primitive quantities
          rTail[k][j][i] = r[k][j][i] + extrap_factor * dr_ds[k][j][i];
          rHead[k][j][i] = r[k][j][i] - extrap_factor * dr_ds[k][j][i];
          uTail[k][j][i].x = u[k][j][i].x + extrap_factor * du_ds[k][j][i].x;
          uHead[k][j][i].x = u[k][j][i].x - extrap_factor * du_ds[k][j][i].x;
          uTail[k][j][i].y = u[k][j][i].y + extrap_factor * du_ds[k][j][i].y;
          uHead[k][j][i].y = u[k][j][i].y - extrap_factor * du_ds[k][j][i].y;
          uTail[k][j][i].z = u[k][j][i].z + extrap_factor * du_ds[k][j][i].z;
          uHead[k][j][i].z = u[k][j][i].z - extrap_factor * du_ds[k][j][i].z;
          pTail[k][j][i] = p[k][j][i] + extrap_factor * dp_ds[k][j][i];
          pHead[k][j][i] = p[k][j][i] - extrap_factor * dp_ds[k][j][i];

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

    /*------------------------------------------------------------------------*
      Predictor step: approximate flux by averaging the governing flux
      functions at either side of the cell.
     *------------------------------------------------------------------------*/

    double courant{dt / (2.0 * m.delta<A>())};

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::predictor>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          // Density
          qu[k][j][i].y =
            r[k][j][i] - courant * (ruTail[k][j][i].y - ruHead[k][j][i].y);

          // Momentum
          // ruv
          // rv^2+p
          // rvw
          const velocity f_qu_Tail{
            rTail[k][j][i] * uTail[k][j][i].x * uTail[k][j][i].y,
            ruTail[k][j][i].y * uTail[k][j][i].y + pTail[k][j][i],
            rTail[k][j][i] * uTail[k][j][i].y * uTail[k][j][i].z};
          const velocity f_qu_Head{
            rHead[k][j][i] * uHead[k][j][i].x * uHead[k][j][i].y,
            ruHead[k][j][i].y * uHead[k][j][i].y + pHead[k][j][i],
            rHead[k][j][i] * uHead[k][j][i].y * uHead[k][j][i].z};
          qu[k][j][i].x = ru[k][j][i].x - courant * (f_qu_Tail.x - f_qu_Head.x);
          qu[k][j][i].y = ru[k][j][i].y - courant * (f_qu_Tail.y - f_qu_Head.y);
          qu[k][j][i].z = ru[k][j][i].z - courant * (f_qu_Tail.z - f_qu_Head.z);

          // Total Energy
          const double f_qE_Tail =
            (rETail[k][j][i] + pTail[k][j][i]) * uTail[k][j][i].y;
          const double f_qE_Head =
            (rEHead[k][j][i] + pHead[k][j][i]) * uHead[k][j][i].y;
          qE[k][j][i] = rE[k][j][i] - courant * (f_qE_Tail - f_qE_Head);

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

    /*------------------------------------------------------------------------*
      Reconstruct faces from intermediates.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::corrector>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          // Primitive quantities
          rTail[k][j][i] = q[k][j - 1][i] + extrap_factor * dr_ds[k][j - 1][i];
          rHead[k][j][i] = q[k][j][i] - extrap_factor * dr_ds[k][j][i];
          uTail[k][j][i].x =
            u[k][j - 1][i].x + extrap_factor * du_ds[k][j - 1][i].x;
          uHead[k][j][i].x = u[k][j][i].x - extrap_factor * du_ds[k][j][i].x;
          uTail[k][j][i].y =
            u[k][j - 1][i].y + extrap_factor * du_ds[k][j - 1][i].y;
          uHead[k][j][i].y = u[k][j][i].y - extrap_factor * du_ds[k][j][i].y;
          uTail[k][j][i].z =
            u[k][j - 1][i].z + extrap_factor * du_ds[k][j - 1][i].z;
          uHead[k][j][i].z = u[k][j][i].z - extrap_factor * du_ds[k][j][i].z;
          pTail[k][j][i] = p[k][j - 1][i] + extrap_factor * dp_ds[k][j - 1][i];
          pHead[k][j][i] = p[k][j][i] - extrap_factor * dp_ds[k][j][i];

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

    /*------------------------------------------------------------------------*
      Compute corrector setp: HLL Riemann solver.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::corrector>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          const double cT = std::sqrt(gamma * pTail[k][j][i] / rTail[k][j][i]);
          const double cH = std::sqrt(gamma * pHead[k][j][i] / rHead[k][j][i]);

          // Update min/max eigenvalues on tail face.
          const double LminT = uTail[k][j][i].y - cT;
          const double LmaxT = uTail[k][j][i].y + cT;

          // Update min/max eigenvalues on head face.
          const double LminH = uHead[k][j][i].y - cH;
          const double LmaxH = uHead[k][j][i].y + cH;

          // Values for HLL.
          const double Lminus{std::min(LminH, std::min(LminT, double{0.0}))};
          const double Lplus{std::max(LmaxH, std::max(LmaxT, double{0.0}))};
          const double Ldiv{1.0 / (Lplus - Lminus)};
          const double Lmult{Lplus * Lminus};
          const double delta_r{rHead[k][j][i] - rTail[k][j][i]};
          const velocity delta_ru{ruHead[k][j][i].x - ruTail[k][j][i].x,
            ruHead[k][j][i].y - ruTail[k][j][i].y,
            ruHead[k][j][i].z - ruTail[k][j][i].z};
          const double delta_rE{rEHead[k][j][i] - rETail[k][j][i]};

          const double f_r_T{ruTail[k][j][i].y};
          const double f_r_H{ruHead[k][j][i].y};
          const velocity f_ru_T{
            rTail[k][j][i] * uTail[k][j][i].x * uTail[k][j][i].y,
            ruTail[k][j][i].y * uTail[k][j][i].y + pTail[k][j][i],
            rTail[k][j][i] * uTail[k][j][i].y * uTail[k][j][i].z};
          const velocity f_ru_H{
            rHead[k][j][i] * uHead[k][j][i].x * uHead[k][j][i].y,
            ruHead[k][j][i].y * uHead[k][j][i].y + pHead[k][j][i],
            rHead[k][j][i] * uHead[k][j][i].y * uHead[k][j][i].z};
          const double f_rE_T{
            (rETail[k][j][i] + pTail[k][j][i]) * uTail[k][j][i].y};
          const double f_rE_H{
            (rEHead[k][j][i] + pHead[k][j][i]) * uHead[k][j][i].y};

          rF[k][j][i] =
            Ldiv * (Lplus * f_r_T - Lminus * f_r_H + Lmult * delta_r);
          ruF[k][j][i].x =
            Ldiv * (Lplus * f_ru_T.x - Lminus * f_ru_H.x + Lmult * delta_ru.x);
          ruF[k][j][i].y =
            Ldiv * (Lplus * f_ru_T.y - Lminus * f_ru_H.y + Lmult * delta_ru.y);
          ruF[k][j][i].z =
            Ldiv * (Lplus * f_ru_T.z - Lminus * f_ru_H.z + Lmult * delta_ru.z);
          rEF[k][j][i] =
            Ldiv * (Lplus * f_rE_T - Lminus * f_rE_H + Lmult * delta_rE);
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Update averages (in-place).
     *------------------------------------------------------------------------*/

    const double update_factor = dt / m.delta<A>();
    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r[k][j][i] =
            r[k][j][i] - update_factor * (rF[k][j + 1][i] - rF[k][j][i]);
          ru[k][j][i].x = ru[k][j][i].x -
                          update_factor * (ruF[k][j + 1][i].x - ruF[k][j][i].x);
          ru[k][j][i].y = ru[k][j][i].y -
                          update_factor * (ruF[k][j + 1][i].y - ruF[k][j][i].y);
          ru[k][j][i].z = ru[k][j][i].z -
                          update_factor * (ruF[k][j + 1][i].z - ruF[k][j][i].z);
          rE[k][j][i] =
            rE[k][j][i] - update_factor * (rEF[k][j + 1][i] - rEF[k][j][i]);
        } // for
      } // for
    } // for
  }
  else if constexpr(A == mesh::axis::z_axis) {
    /*------------------------------------------------------------------------*
      Compute slopes.
     *------------------------------------------------------------------------*/

    const double slope_factor{1.0 / m.delta<A>()};

    for(auto k : m.cells<mesh::z_axis, mesh::predictor>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          dr_ds[k][j][i] =
            slope_factor * L::limit(r[k-1][j][i], r[k][j][i], r[k+1][j][i]);
          du_ds[k][j][i].x =
            slope_factor *
            L::limit(u[k-1][j][i].x, u[k][j][i].x, u[k+1][j][i].x);
          du_ds[k][j][i].y =
            slope_factor *
            L::limit(u[k-1][j][i].y, u[k][j][i].y, u[k+1][j][i].y);
          du_ds[k][j][i].z =
            slope_factor *
            L::limit(u[k-1][j][i].z, u[k][j][i].z, u[k+1][j][i].z);
          dp_ds[k][j][i] =
            slope_factor * L::limit(p[k-1][j][i], p[k][j][i], p[k+1][j][i]);
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Reconstruct faces.
     *------------------------------------------------------------------------*/

    double extrap_factor{m.delta<A>() / 2.0};
    double p2rEFactor{1.0 / (gamma - 1.0)};

    for(auto k : m.cells<mesh::z_axis, mesh::predictor>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          // Primitive quantities
          rTail[k][j][i] = r[k][j][i] + extrap_factor * dr_ds[k][j][i];
          rHead[k][j][i] = r[k][j][i] - extrap_factor * dr_ds[k][j][i];
          uTail[k][j][i].x = u[k][j][i].x + extrap_factor * du_ds[k][j][i].x;
          uHead[k][j][i].x = u[k][j][i].x - extrap_factor * du_ds[k][j][i].x;
          uTail[k][j][i].y = u[k][j][i].y + extrap_factor * du_ds[k][j][i].y;
          uHead[k][j][i].y = u[k][j][i].y - extrap_factor * du_ds[k][j][i].y;
          uTail[k][j][i].z = u[k][j][i].z + extrap_factor * du_ds[k][j][i].z;
          uHead[k][j][i].z = u[k][j][i].z - extrap_factor * du_ds[k][j][i].z;
          pTail[k][j][i] = p[k][j][i] + extrap_factor * dp_ds[k][j][i];
          pHead[k][j][i] = p[k][j][i] - extrap_factor * dp_ds[k][j][i];

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

    /*------------------------------------------------------------------------*
      Predictor step: approximate flux by averaging the governing flux
      functions at either side of the cell.
     *------------------------------------------------------------------------*/

    double courant{dt / (2.0 * m.delta<A>())};

    for(auto k : m.cells<mesh::z_axis, mesh::predictor>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          // Density
          qu[k][j][i].z =
            r[k][j][i] - courant * (ruTail[k][j][i].z - ruHead[k][j][i].z);

          // Momentum
          // ruw
          // rvw
          // rw^2+p
          const velocity f_qu_Tail{
            rTail[k][j][i] * uTail[k][j][i].x * uTail[k][j][i].z,
            rTail[k][j][i] * uTail[k][j][i].y * uTail[k][j][i].z,
            ruTail[k][j][i].z * uTail[k][j][i].z + pTail[k][j][i]};
          const velocity f_qu_Head{
            rHead[k][j][i] * uHead[k][j][i].x * uHead[k][j][i].z,
            rHead[k][j][i] * uHead[k][j][i].y * uHead[k][j][i].z,
            ruHead[k][j][i].z * uHead[k][j][i].z + pHead[k][j][i]};
          qu[k][j][i].x = ru[k][j][i].x - courant * (f_qu_Tail.x - f_qu_Head.x);
          qu[k][j][i].y = ru[k][j][i].y - courant * (f_qu_Tail.y - f_qu_Head.y);
          qu[k][j][i].z = ru[k][j][i].z - courant * (f_qu_Tail.z - f_qu_Head.z);

          // Total Energy
          const double f_qE_Tail =
            (rETail[k][j][i] + pTail[k][j][i]) * uTail[k][j][i].z;
          const double f_qE_Head =
            (rEHead[k][j][i] + pHead[k][j][i]) * uHead[k][j][i].z;
          qE[k][j][i] = rE[k][j][i] - courant * (f_qE_Tail - f_qE_Head);

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

    /*------------------------------------------------------------------------*
      Reconstruct faces from intermediates.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::corrector>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          // Primitive quantities
          rTail[k][j][i] = q[k-1][j][i] + extrap_factor * dr_ds[k-1][j][i];
          rHead[k][j][i] = q[k][j][i] - extrap_factor * dr_ds[k][j][i];
          uTail[k][j][i].x =
            u[k-1][j][i].x + extrap_factor * du_ds[k-1][j][i].x;
          uHead[k][j][i].x = u[k][j][i].x - extrap_factor * du_ds[k][j][i].x;
          uTail[k][j][i].y =
            u[k-1][j][i].y + extrap_factor * du_ds[k-1][j][i].y;
          uHead[k][j][i].y = u[k][j][i].y - extrap_factor * du_ds[k][j][i].y;
          uTail[k][j][i].z =
            u[k-1][j][i].z + extrap_factor * du_ds[k-1][j][i].z;
          uHead[k][j][i].z = u[k][j][i].z - extrap_factor * du_ds[k][j][i].z;
          pTail[k][j][i] = p[k-1][j][i] + extrap_factor * dp_ds[k-1][j][i];
          pHead[k][j][i] = p[k][j][i] - extrap_factor * dp_ds[k][j][i];

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

    /*------------------------------------------------------------------------*
      Compute corrector setp: HLL Riemann solver.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::corrector>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          const double cT = std::sqrt(gamma * pTail[k][j][i] / rTail[k][j][i]);
          const double cH = std::sqrt(gamma * pHead[k][j][i] / rHead[k][j][i]);

          // Update min/max eigenvalues on tail face.
          const double LminT = uTail[k][j][i].z - cT;
          const double LmaxT = uTail[k][j][i].z + cT;

          // Update min/max eigenvalues on head face.
          const double LminH = uHead[k][j][i].z - cH;
          const double LmaxH = uHead[k][j][i].z + cH;

          // Values for HLL.
          const double Lminus{std::min(LminH, std::min(LminT, double{0.0}))};
          const double Lplus{std::max(LmaxH, std::max(LmaxT, double{0.0}))};
          const double Ldiv{1.0 / (Lplus - Lminus)};
          const double Lmult{Lplus * Lminus};
          const double delta_r{rHead[k][j][i] - rTail[k][j][i]};
          const velocity delta_ru{ruHead[k][j][i].x - ruTail[k][j][i].x,
            ruHead[k][j][i].y - ruTail[k][j][i].y,
            ruHead[k][j][i].z - ruTail[k][j][i].z};
          const double delta_rE{rEHead[k][j][i] - rETail[k][j][i]};

          const double f_r_T{ruTail[k][j][i].z};
          const double f_r_H{ruHead[k][j][i].z};
          const velocity f_ru_T{
            rTail[k][j][i] * uTail[k][j][i].x * uTail[k][j][i].z,
            rTail[k][j][i] * uTail[k][j][i].y * uTail[k][j][i].z,
            ruTail[k][j][i].z * uTail[k][j][i].z + pTail[k][j][i]};
          const velocity f_ru_H{
            rHead[k][j][i] * uHead[k][j][i].x * uHead[k][j][i].z,
            rHead[k][j][i] * uHead[k][j][i].y * uHead[k][j][i].z,
            ruHead[k][j][i].z * uHead[k][j][i].z + pHead[k][j][i]};
          const double f_rE_T{
            (rETail[k][j][i] + pTail[k][j][i]) * uTail[k][j][i].z};
          const double f_rE_H{
            (rEHead[k][j][i] + pHead[k][j][i]) * uHead[k][j][i].z};

          rF[k][j][i] =
            Ldiv * (Lplus * f_r_T - Lminus * f_r_H + Lmult * delta_r);
          ruF[k][j][i].x =
            Ldiv * (Lplus * f_ru_T.x - Lminus * f_ru_H.x + Lmult * delta_ru.x);
          ruF[k][j][i].y =
            Ldiv * (Lplus * f_ru_T.y - Lminus * f_ru_H.y + Lmult * delta_ru.y);
          ruF[k][j][i].z =
            Ldiv * (Lplus * f_ru_T.z - Lminus * f_ru_H.z + Lmult * delta_ru.z);
          rEF[k][j][i] =
            Ldiv * (Lplus * f_rE_T - Lminus * f_rE_H + Lmult * delta_rE);
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Update averages (in-place).
     *------------------------------------------------------------------------*/

    const double update_factor = dt / m.delta<A>();
    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          r[k][j][i] =
            r[k][j][i] - update_factor * (rF[k+1][j][i] - rF[k][j][i]);
          ru[k][j][i].x = ru[k][j][i].x -
                          update_factor * (ruF[k+1][j][i].x - ruF[k][j][i].x);
          ru[k][j][i].y = ru[k][j][i].y -
                          update_factor * (ruF[k+1][j][i].y - ruF[k][j][i].y);
          ru[k][j][i].z = ru[k][j][i].z -
                          update_factor * (ruF[k+1][j][i].z - ruF[k][j][i].z);
          rE[k][j][i] =
            rE[k][j][i] - update_factor * (rEF[k+1][j][i] - rEF[k][j][i]);
        } // for
      } // for
    } // for
  } // if
} // advance

} // namespace muscl::tasks::hydro

#endif // MUSCL_TASKS_HYDRO_HH
