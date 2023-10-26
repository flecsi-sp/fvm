#ifndef MUSCL_TASKS_HYDRO_HH
#define MUSCL_TASKS_HYDRO_HH

#include "../state.hh"
#include "../utils.hh"
#include "utils.hh"

namespace muscl::tasks::hydro {

void update_primitives(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  field<vec3>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<double>::accessor<ro> gamma_a);

void update_eigenvalues(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<vec3>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<vec3>::accessor<wo> lmax_a,
  single<double>::accessor<ro> gamma_a);

double update_dtmin(mesh::accessor<ro> m, single<vec3>::accessor<ro> lmax);

#define DEBUG_PRINT(Ms, Mm, Mr, Mru, MrE)                                      \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Ms << std::endl;                                                     \
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {                  \
      for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {                 \
        ss << Mr(i, j, 2) << " ";                                              \
      }                                                                        \
      ss << std::endl;                                                         \
    }                                                                          \
    flog(info) << ss.str() << std::endl;                                       \
  }                                                                            \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Ms << std::endl;                                                     \
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {                  \
      for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {                 \
        ss << Mru(i, j, 2) << " ";                                             \
      }                                                                        \
      ss << std::endl;                                                         \
    }                                                                          \
    flog(info) << ss.str() << std::endl;                                       \
  }                                                                            \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Ms << std::endl;                                                     \
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {                  \
      for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {                 \
        ss << MrE(i, j, 2) << " ";                                             \
      }                                                                        \
      ss << std::endl;                                                         \
    }                                                                          \
    flog(info) << ss.str() << std::endl;                                       \
  }

template<mesh::axis A, typename L>
void
advance(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  field<vec3>::accessor<rw, ro> u_a,
  field<double>::accessor<rw, ro> p_a,
  field<double>::accessor<wo, ro> q_a,
  field<vec3>::accessor<wo, ro> qu_a,
  field<double>::accessor<wo, ro> qE_a,
  field<double>::accessor<wo, ro> dr_ds_a,
  field<vec3>::accessor<wo, ro> du_ds_a,
  field<double>::accessor<wo, ro> dp_ds_a,
  field<double>::accessor<wo, ro> rTail_a,
  field<vec3>::accessor<wo, ro> ruTail_a,
  field<double>::accessor<wo, ro> rETail_a,
  field<vec3>::accessor<wo, ro> uTail_a,
  field<double>::accessor<wo, ro> pTail_a,
  field<double>::accessor<wo, ro> rHead_a,
  field<vec3>::accessor<wo, ro> ruHead_a,
  field<double>::accessor<wo, ro> rEHead_a,
  field<vec3>::accessor<wo, ro> uHead_a,
  field<double>::accessor<wo, ro> pHead_a,
  field<double>::accessor<wo, ro> rF_a,
  field<vec3>::accessor<wo, ro> ruF_a,
  field<double>::accessor<wo, ro> rEF_a,
  single<double>::accessor<ro> gamma_a,
  double dt) {
  auto r = m.mdcolex<mesh::cells>(r_a);
  auto ru = m.mdcolex<mesh::cells>(ru_a);
  auto rE = m.mdcolex<mesh::cells>(rE_a);
  auto u = m.mdcolex<mesh::cells>(u_a);
  auto p = m.mdcolex<mesh::cells>(p_a);
  auto q = m.mdcolex<mesh::cells>(q_a);
  auto qu = m.mdcolex<mesh::cells>(qu_a);
  auto qE = m.mdcolex<mesh::cells>(qE_a);
  auto dr_ds = m.mdcolex<mesh::cells>(dr_ds_a);
  auto du_ds = m.mdcolex<mesh::cells>(du_ds_a);
  auto dp_ds = m.mdcolex<mesh::cells>(dp_ds_a);
  auto rTail = m.mdcolex<mesh::cells>(rTail_a);
  auto ruTail = m.mdcolex<mesh::cells>(ruTail_a);
  auto rETail = m.mdcolex<mesh::cells>(rETail_a);
  auto uTail = m.mdcolex<mesh::cells>(uTail_a);
  auto pTail = m.mdcolex<mesh::cells>(pTail_a);
  auto rHead = m.mdcolex<mesh::cells>(rHead_a);
  auto ruHead = m.mdcolex<mesh::cells>(ruHead_a);
  auto rEHead = m.mdcolex<mesh::cells>(rEHead_a);
  auto uHead = m.mdcolex<mesh::cells>(uHead_a);
  auto pHead = m.mdcolex<mesh::cells>(pHead_a);
  auto rF = m.mdcolex<mesh::cells>(rF_a);
  auto ruF = m.mdcolex<mesh::cells>(ruF_a);
  auto rEF = m.mdcolex<mesh::cells>(rEF_a);
  auto const gamma = *gamma_a;

  if constexpr(A == mesh::axis::x_axis) {
    /*------------------------------------------------------------------------*
      Compute slopes.
     *------------------------------------------------------------------------*/

    const double slope_factor{1.0 / m.delta<A>()};

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::predictor>()) {
          dr_ds(i, j, k) =
            slope_factor * L::limit(r(i - 1, j, k), r(i, j, k), r(i + 1, j, k));
          du_ds(i, j, k).x =
            slope_factor *
            L::limit(u(i - 1, j, k).x, u(i, j, k).x, u(i + 1, j, k).x);
          du_ds(i, j, k).y =
            slope_factor *
            L::limit(u(i - 1, j, k).y, u(i, j, k).y, u(i + 1, j, k).y);
          du_ds(i, j, k).z =
            slope_factor *
            L::limit(u(i - 1, j, k).z, u(i, j, k).z, u(i + 1, j, k).z);
          dp_ds(i, j, k) =
            slope_factor * L::limit(p(i - 1, j, k), p(i, j, k), p(i + 1, j, k));
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
          rTail(i, j, k) = r(i, j, k) + extrap_factor * dr_ds(i, j, k);
          rHead(i, j, k) = r(i, j, k) - extrap_factor * dr_ds(i, j, k);
          uTail(i, j, k).x = u(i, j, k).x + extrap_factor * du_ds(i, j, k).x;
          uHead(i, j, k).x = u(i, j, k).x - extrap_factor * du_ds(i, j, k).x;
          uTail(i, j, k).y = u(i, j, k).y + extrap_factor * du_ds(i, j, k).y;
          uHead(i, j, k).y = u(i, j, k).y - extrap_factor * du_ds(i, j, k).y;
          uTail(i, j, k).z = u(i, j, k).z + extrap_factor * du_ds(i, j, k).z;
          uHead(i, j, k).z = u(i, j, k).z - extrap_factor * du_ds(i, j, k).z;
          pTail(i, j, k) = p(i, j, k) + extrap_factor * dp_ds(i, j, k);
          pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

          // Conserved quantities
          ruTail(i, j, k).x = rTail(i, j, k) * uTail(i, j, k).x;
          ruHead(i, j, k).x = rHead(i, j, k) * uHead(i, j, k).x;
          ruTail(i, j, k).y = rTail(i, j, k) * uTail(i, j, k).y;
          ruHead(i, j, k).y = rHead(i, j, k) * uHead(i, j, k).y;
          ruTail(i, j, k).z = rTail(i, j, k) * uTail(i, j, k).z;
          ruHead(i, j, k).z = rHead(i, j, k) * uHead(i, j, k).z;
          rETail(i, j, k) =
            p2rEFactor * pTail(i, j, k) +
            0.5 * rTail(i, j, k) *
              (utils::sqr(uTail(i, j, k).x) + utils::sqr(uTail(i, j, k).y) +
                utils::sqr(uTail(i, j, k).z));
          rEHead(i, j, k) =
            p2rEFactor * pHead(i, j, k) +
            0.5 * rHead(i, j, k) *
              (utils::sqr(uHead(i, j, k).x) + utils::sqr(uHead(i, j, k).y) +
                utils::sqr(uHead(i, j, k).z));
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
          q(i, j, k) =
            r(i, j, k) - courant * (ruTail(i, j, k).x - ruHead(i, j, k).x);

          // Momentum
          // ru^2 + p
          // ruv
          // ruw
          const vec3 f_qu_Tail{
            ruTail(i, j, k).x * uTail(i, j, k).x + pTail(i, j, k),
            rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).y,
            rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).z};
          const vec3 f_qu_Head{
            ruHead(i, j, k).x * uHead(i, j, k).x + pHead(i, j, k),
            rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).y,
            rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).z};
          qu(i, j, k).x = ru(i, j, k).x - courant * (f_qu_Tail.x - f_qu_Head.x);
          qu(i, j, k).y = ru(i, j, k).y - courant * (f_qu_Tail.y - f_qu_Head.y);
          qu(i, j, k).z = ru(i, j, k).z - courant * (f_qu_Tail.z - f_qu_Head.z);

          // Total Energy
          const double f_qE_Tail =
            (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).x;
          const double f_qE_Head =
            (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).x;
          qE(i, j, k) = rE(i, j, k) - courant * (f_qE_Tail - f_qE_Head);

          // Primitives
          u(i, j, k).x = qu(i, j, k).x / q(i, j, k);
          u(i, j, k).y = qu(i, j, k).y / q(i, j, k);
          u(i, j, k).z = qu(i, j, k).z / q(i, j, k);
          p(i, j, k) =
            (gamma - 1.0) * (qE(i, j, k) - 0.5 * q(i, j, k) *
                                             (utils::sqr(u(i, j, k).x) +
                                               utils::sqr(u(i, j, k).y) +
                                               utils::sqr(u(i, j, k).z)));
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
          rTail(i, j, k) = q(i - 1, j, k) + extrap_factor * dr_ds(i - 1, j, k);
          rHead(i, j, k) = q(i, j, k) - extrap_factor * dr_ds(i, j, k);
          uTail(i, j, k).x =
            u(i - 1, j, k).x + extrap_factor * du_ds(i - 1, j, k).x;
          uHead(i, j, k).x = u(i, j, k).x - extrap_factor * du_ds(i, j, k).x;
          uTail(i, j, k).y =
            u(i - 1, j, k).y + extrap_factor * du_ds(i - 1, j, k).y;
          uHead(i, j, k).y = u(i, j, k).y - extrap_factor * du_ds(i, j, k).y;
          uTail(i, j, k).z =
            u(i - 1, j, k).z + extrap_factor * du_ds(i - 1, j, k).z;
          uHead(i, j, k).z = u(i, j, k).z - extrap_factor * du_ds(i, j, k).z;
          pTail(i, j, k) = p(i - 1, j, k) + extrap_factor * dp_ds(i - 1, j, k);
          pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

          // Conserved quantities
          ruTail(i, j, k).x = rTail(i, j, k) * uTail(i, j, k).x;
          ruHead(i, j, k).x = rHead(i, j, k) * uHead(i, j, k).x;
          ruTail(i, j, k).y = rTail(i, j, k) * uTail(i, j, k).y;
          ruHead(i, j, k).y = rHead(i, j, k) * uHead(i, j, k).y;
          ruTail(i, j, k).z = rTail(i, j, k) * uTail(i, j, k).z;
          ruHead(i, j, k).z = rHead(i, j, k) * uHead(i, j, k).z;
          rETail(i, j, k) =
            p2rEFactor * pTail(i, j, k) +
            0.5 * rTail(i, j, k) *
              (utils::sqr(uTail(i, j, k).x) + utils::sqr(uTail(i, j, k).y) +
                utils::sqr(uTail(i, j, k).z));
          rEHead(i, j, k) =
            p2rEFactor * pHead(i, j, k) +
            0.5 * rHead(i, j, k) *
              (utils::sqr(uHead(i, j, k).x) + utils::sqr(uHead(i, j, k).y) +
                utils::sqr(uHead(i, j, k).z));
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Compute corrector setp: HLL Riemann solver.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::corrector>()) {
          const double cT = std::sqrt(gamma * pTail(i, j, k) / rTail(i, j, k));
          const double cH = std::sqrt(gamma * pHead(i, j, k) / rHead(i, j, k));

          // Update min/max eigenvalues on tail face.
          const double LminT = uTail(i, j, k).x - cT;
          const double LmaxT = uTail(i, j, k).x + cT;

          // Update min/max eigenvalues on head face.
          const double LminH = uHead(i, j, k).x - cH;
          const double LmaxH = uHead(i, j, k).x + cH;

          // Values for HLL.
          const double Lminus{std::min(LminH, std::min(LminT, double{0.0}))};
          const double Lplus{std::max(LmaxH, std::max(LmaxT, double{0.0}))};
          const double Ldiv{1.0 / (Lplus - Lminus)};
          const double Lmult{Lplus * Lminus};
          const double delta_r{rHead(i, j, k) - rTail(i, j, k)};
          const vec3 delta_ru{ruHead(i, j, k).x - ruTail(i, j, k).x,
            ruHead(i, j, k).y - ruTail(i, j, k).y,
            ruHead(i, j, k).z - ruTail(i, j, k).z};
          const double delta_rE{rEHead(i, j, k) - rETail(i, j, k)};

          const double f_r_T{ruTail(i, j, k).x};
          const double f_r_H{ruHead(i, j, k).x};
          const vec3 f_ru_T{
            ruTail(i, j, k).x * uTail(i, j, k).x + pTail(i, j, k),
            ruTail(i, j, k).x * uTail(i, j, k).y,
            ruTail(i, j, k).x * uTail(i, j, k).z};
          const vec3 f_ru_H{
            ruHead(i, j, k).x * uHead(i, j, k).x + pHead(i, j, k),
            ruHead(i, j, k).x * uHead(i, j, k).y,
            ruHead(i, j, k).x * uHead(i, j, k).z};
          const double f_rE_T{
            (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).x};
          const double f_rE_H{
            (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).x};

          rF(i, j, k) =
            Ldiv * (Lplus * f_r_T - Lminus * f_r_H + Lmult * delta_r);
          ruF(i, j, k).x =
            Ldiv * (Lplus * f_ru_T.x - Lminus * f_ru_H.x + Lmult * delta_ru.x);
          ruF(i, j, k).y =
            Ldiv * (Lplus * f_ru_T.y - Lminus * f_ru_H.y + Lmult * delta_ru.y);
          ruF(i, j, k).z =
            Ldiv * (Lplus * f_ru_T.z - Lminus * f_ru_H.z + Lmult * delta_ru.z);
          rEF(i, j, k) =
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
          r(i, j, k) =
            r(i, j, k) - update_factor * (rF(i + 1, j, k) - rF(i, j, k));
          ru(i, j, k).x = ru(i, j, k).x -
                          update_factor * (ruF(i + 1, j, k).x - ruF(i, j, k).x);
          ru(i, j, k).y = ru(i, j, k).y -
                          update_factor * (ruF(i + 1, j, k).y - ruF(i, j, k).y);
          ru(i, j, k).z = ru(i, j, k).z -
                          update_factor * (ruF(i + 1, j, k).z - ruF(i, j, k).z);
          rE(i, j, k) =
            rE(i, j, k) - update_factor * (rEF(i + 1, j, k) - rEF(i, j, k));
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
          dr_ds(i, j, k) =
            slope_factor * L::limit(r(i, j - 1, k), r(i, j, k), r(i, j + 1, k));
          du_ds(i, j, k).x =
            slope_factor *
            L::limit(u(i, j - 1, k).x, u(i, j, k).x, u(i, j + 1, k).x);
          du_ds(i, j, k).y =
            slope_factor *
            L::limit(u(i, j - 1, k).y, u(i, j, k).y, u(i, j + 1, k).y);
          du_ds(i, j, k).z =
            slope_factor *
            L::limit(u(i, j - 1, k).z, u(i, j, k).z, u(i, j + 1, k).z);
          dp_ds(i, j, k) =
            slope_factor * L::limit(p(i, j - 1, k), p(i, j, k), p(i, j + 1, k));
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
          rTail(i, j, k) = r(i, j, k) + extrap_factor * dr_ds(i, j, k);
          rHead(i, j, k) = r(i, j, k) - extrap_factor * dr_ds(i, j, k);
          uTail(i, j, k).x = u(i, j, k).x + extrap_factor * du_ds(i, j, k).x;
          uHead(i, j, k).x = u(i, j, k).x - extrap_factor * du_ds(i, j, k).x;
          uTail(i, j, k).y = u(i, j, k).y + extrap_factor * du_ds(i, j, k).y;
          uHead(i, j, k).y = u(i, j, k).y - extrap_factor * du_ds(i, j, k).y;
          uTail(i, j, k).z = u(i, j, k).z + extrap_factor * du_ds(i, j, k).z;
          uHead(i, j, k).z = u(i, j, k).z - extrap_factor * du_ds(i, j, k).z;
          pTail(i, j, k) = p(i, j, k) + extrap_factor * dp_ds(i, j, k);
          pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

          // Conserved quantities
          ruTail(i, j, k).x = rTail(i, j, k) * uTail(i, j, k).x;
          ruHead(i, j, k).x = rHead(i, j, k) * uHead(i, j, k).x;
          ruTail(i, j, k).y = rTail(i, j, k) * uTail(i, j, k).y;
          ruHead(i, j, k).y = rHead(i, j, k) * uHead(i, j, k).y;
          ruTail(i, j, k).z = rTail(i, j, k) * uTail(i, j, k).z;
          ruHead(i, j, k).z = rHead(i, j, k) * uHead(i, j, k).z;
          rETail(i, j, k) =
            p2rEFactor * pTail(i, j, k) +
            0.5 * rTail(i, j, k) *
              (utils::sqr(uTail(i, j, k).x) + utils::sqr(uTail(i, j, k).y) +
                utils::sqr(uTail(i, j, k).z));
          rEHead(i, j, k) =
            p2rEFactor * pHead(i, j, k) +
            0.5 * rHead(i, j, k) *
              (utils::sqr(uHead(i, j, k).x) + utils::sqr(uHead(i, j, k).y) +
                utils::sqr(uHead(i, j, k).z));
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
          q(i, j, k) =
            r(i, j, k) - courant * (ruTail(i, j, k).y - ruHead(i, j, k).y);

          // Momentum
          // ruv
          // rv^2+p
          // rvw
          const vec3 f_qu_Tail{
            rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).y,
            ruTail(i, j, k).y * uTail(i, j, k).y + pTail(i, j, k),
            rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z};
          const vec3 f_qu_Head{
            rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).y,
            ruHead(i, j, k).y * uHead(i, j, k).y + pHead(i, j, k),
            rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z};
          qu(i, j, k).x = ru(i, j, k).x - courant * (f_qu_Tail.x - f_qu_Head.x);
          qu(i, j, k).y = ru(i, j, k).y - courant * (f_qu_Tail.y - f_qu_Head.y);
          qu(i, j, k).z = ru(i, j, k).z - courant * (f_qu_Tail.z - f_qu_Head.z);

          // Total Energy
          const double f_qE_Tail =
            (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).y;
          const double f_qE_Head =
            (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).y;
          qE(i, j, k) = rE(i, j, k) - courant * (f_qE_Tail - f_qE_Head);

          // Primitives
          u(i, j, k).x = qu(i, j, k).x / q(i, j, k);
          u(i, j, k).y = qu(i, j, k).y / q(i, j, k);
          u(i, j, k).z = qu(i, j, k).z / q(i, j, k);
          p(i, j, k) =
            (gamma - 1.0) * (qE(i, j, k) - 0.5 * q(i, j, k) *
                                             (utils::sqr(u(i, j, k).x) +
                                               utils::sqr(u(i, j, k).y) +
                                               utils::sqr(u(i, j, k).z)));
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
          rTail(i, j, k) = q(i, j - 1, k) + extrap_factor * dr_ds(i, j - 1, k);
          rHead(i, j, k) = q(i, j, k) - extrap_factor * dr_ds(i, j, k);
          uTail(i, j, k).x =
            u(i, j - 1, k).x + extrap_factor * du_ds(i, j - 1, k).x;
          uHead(i, j, k).x = u(i, j, k).x - extrap_factor * du_ds(i, j, k).x;
          uTail(i, j, k).y =
            u(i, j - 1, k).y + extrap_factor * du_ds(i, j - 1, k).y;
          uHead(i, j, k).y = u(i, j, k).y - extrap_factor * du_ds(i, j, k).y;
          uTail(i, j, k).z =
            u(i, j - 1, k).z + extrap_factor * du_ds(i, j - 1, k).z;
          uHead(i, j, k).z = u(i, j, k).z - extrap_factor * du_ds(i, j, k).z;
          pTail(i, j, k) = p(i, j - 1, k) + extrap_factor * dp_ds(i, j - 1, k);
          pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

          // Conserved quantities
          ruTail(i, j, k).x = rTail(i, j, k) * uTail(i, j, k).x;
          ruHead(i, j, k).x = rHead(i, j, k) * uHead(i, j, k).x;
          ruTail(i, j, k).y = rTail(i, j, k) * uTail(i, j, k).y;
          ruHead(i, j, k).y = rHead(i, j, k) * uHead(i, j, k).y;
          ruTail(i, j, k).z = rTail(i, j, k) * uTail(i, j, k).z;
          ruHead(i, j, k).z = rHead(i, j, k) * uHead(i, j, k).z;
          rETail(i, j, k) =
            p2rEFactor * pTail(i, j, k) +
            0.5 * rTail(i, j, k) *
              (utils::sqr(uTail(i, j, k).x) + utils::sqr(uTail(i, j, k).y) +
                utils::sqr(uTail(i, j, k).z));
          rEHead(i, j, k) =
            p2rEFactor * pHead(i, j, k) +
            0.5 * rHead(i, j, k) *
              (utils::sqr(uHead(i, j, k).x) + utils::sqr(uHead(i, j, k).y) +
                utils::sqr(uHead(i, j, k).z));
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Compute corrector setp: HLL Riemann solver.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::corrector>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          const double cT = std::sqrt(gamma * pTail(i, j, k) / rTail(i, j, k));
          const double cH = std::sqrt(gamma * pHead(i, j, k) / rHead(i, j, k));

          // Update min/max eigenvalues on tail face.
          const double LminT = uTail(i, j, k).y - cT;
          const double LmaxT = uTail(i, j, k).y + cT;

          // Update min/max eigenvalues on head face.
          const double LminH = uHead(i, j, k).y - cH;
          const double LmaxH = uHead(i, j, k).y + cH;

          // Values for HLL.
          const double Lminus{std::min(LminH, std::min(LminT, double{0.0}))};
          const double Lplus{std::max(LmaxH, std::max(LmaxT, double{0.0}))};
          const double Ldiv{1.0 / (Lplus - Lminus)};
          const double Lmult{Lplus * Lminus};
          const double delta_r{rHead(i, j, k) - rTail(i, j, k)};
          const vec3 delta_ru{ruHead(i, j, k).x - ruTail(i, j, k).x,
            ruHead(i, j, k).y - ruTail(i, j, k).y,
            ruHead(i, j, k).z - ruTail(i, j, k).z};
          const double delta_rE{rEHead(i, j, k) - rETail(i, j, k)};

          const double f_r_T{ruTail(i, j, k).y};
          const double f_r_H{ruHead(i, j, k).y};
          const vec3 f_ru_T{
            rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).y,
            ruTail(i, j, k).y * uTail(i, j, k).y + pTail(i, j, k),
            rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z};
          const vec3 f_ru_H{
            rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).y,
            ruHead(i, j, k).y * uHead(i, j, k).y + pHead(i, j, k),
            rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z};
          const double f_rE_T{
            (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).y};
          const double f_rE_H{
            (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).y};

          rF(i, j, k) =
            Ldiv * (Lplus * f_r_T - Lminus * f_r_H + Lmult * delta_r);
          ruF(i, j, k).x =
            Ldiv * (Lplus * f_ru_T.x - Lminus * f_ru_H.x + Lmult * delta_ru.x);
          ruF(i, j, k).y =
            Ldiv * (Lplus * f_ru_T.y - Lminus * f_ru_H.y + Lmult * delta_ru.y);
          ruF(i, j, k).z =
            Ldiv * (Lplus * f_ru_T.z - Lminus * f_ru_H.z + Lmult * delta_ru.z);
          rEF(i, j, k) =
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
          r(i, j, k) =
            r(i, j, k) - update_factor * (rF(i, j + 1, k) - rF(i, j, k));
          ru(i, j, k).x = ru(i, j, k).x -
                          update_factor * (ruF(i, j + 1, k).x - ruF(i, j, k).x);
          ru(i, j, k).y = ru(i, j, k).y -
                          update_factor * (ruF(i, j + 1, k).y - ruF(i, j, k).y);
          ru(i, j, k).z = ru(i, j, k).z -
                          update_factor * (ruF(i, j + 1, k).z - ruF(i, j, k).z);
          rE(i, j, k) =
            rE(i, j, k) - update_factor * (rEF(i, j + 1, k) - rEF(i, j, k));
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
          dr_ds(i, j, k) =
            slope_factor * L::limit(r(i, j, k - 1), r(i, j, k), r(i, j, k + 1));
          du_ds(i, j, k).x =
            slope_factor *
            L::limit(u(i, j, k - 1).x, u(i, j, k).x, u(i, j, k + 1).x);
          du_ds(i, j, k).y =
            slope_factor *
            L::limit(u(i, j, k - 1).y, u(i, j, k).y, u(i, j, k + 1).y);
          du_ds(i, j, k).z =
            slope_factor *
            L::limit(u(i, j, k - 1).z, u(i, j, k).z, u(i, j, k + 1).z);
          dp_ds(i, j, k) =
            slope_factor * L::limit(p(i, j, k - 1), p(i, j, k), p(i, j, k + 1));
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
          rTail(i, j, k) = r(i, j, k) + extrap_factor * dr_ds(i, j, k);
          rHead(i, j, k) = r(i, j, k) - extrap_factor * dr_ds(i, j, k);
          uTail(i, j, k).x = u(i, j, k).x + extrap_factor * du_ds(i, j, k).x;
          uHead(i, j, k).x = u(i, j, k).x - extrap_factor * du_ds(i, j, k).x;
          uTail(i, j, k).y = u(i, j, k).y + extrap_factor * du_ds(i, j, k).y;
          uHead(i, j, k).y = u(i, j, k).y - extrap_factor * du_ds(i, j, k).y;
          uTail(i, j, k).z = u(i, j, k).z + extrap_factor * du_ds(i, j, k).z;
          uHead(i, j, k).z = u(i, j, k).z - extrap_factor * du_ds(i, j, k).z;
          pTail(i, j, k) = p(i, j, k) + extrap_factor * dp_ds(i, j, k);
          pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

          // Conserved quantities
          ruTail(i, j, k).x = rTail(i, j, k) * uTail(i, j, k).x;
          ruHead(i, j, k).x = rHead(i, j, k) * uHead(i, j, k).x;
          ruTail(i, j, k).y = rTail(i, j, k) * uTail(i, j, k).y;
          ruHead(i, j, k).y = rHead(i, j, k) * uHead(i, j, k).y;
          ruTail(i, j, k).z = rTail(i, j, k) * uTail(i, j, k).z;
          ruHead(i, j, k).z = rHead(i, j, k) * uHead(i, j, k).z;
          rETail(i, j, k) =
            p2rEFactor * pTail(i, j, k) +
            0.5 * rTail(i, j, k) *
              (utils::sqr(uTail(i, j, k).x) + utils::sqr(uTail(i, j, k).y) +
                utils::sqr(uTail(i, j, k).z));
          rEHead(i, j, k) =
            p2rEFactor * pHead(i, j, k) +
            0.5 * rHead(i, j, k) *
              (utils::sqr(uHead(i, j, k).x) + utils::sqr(uHead(i, j, k).y) +
                utils::sqr(uHead(i, j, k).z));
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
          q(i, j, k) =
            r(i, j, k) - courant * (ruTail(i, j, k).z - ruHead(i, j, k).z);

          // Momentum
          // ruw
          // rvw
          // rw^2+p
          const vec3 f_qu_Tail{
            rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).z,
            rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z,
            ruTail(i, j, k).z * uTail(i, j, k).z + pTail(i, j, k)};
          const vec3 f_qu_Head{
            rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).z,
            rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z,
            ruHead(i, j, k).z * uHead(i, j, k).z + pHead(i, j, k)};
          qu(i, j, k).x = ru(i, j, k).x - courant * (f_qu_Tail.x - f_qu_Head.x);
          qu(i, j, k).y = ru(i, j, k).y - courant * (f_qu_Tail.y - f_qu_Head.y);
          qu(i, j, k).z = ru(i, j, k).z - courant * (f_qu_Tail.z - f_qu_Head.z);

          // Total Energy
          const double f_qE_Tail =
            (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).z;
          const double f_qE_Head =
            (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).z;
          qE(i, j, k) = rE(i, j, k) - courant * (f_qE_Tail - f_qE_Head);

          // Primitives
          u(i, j, k).x = qu(i, j, k).x / q(i, j, k);
          u(i, j, k).y = qu(i, j, k).y / q(i, j, k);
          u(i, j, k).z = qu(i, j, k).z / q(i, j, k);
          p(i, j, k) =
            (gamma - 1.0) * (qE(i, j, k) - 0.5 * q(i, j, k) *
                                             (utils::sqr(u(i, j, k).x) +
                                               utils::sqr(u(i, j, k).y) +
                                               utils::sqr(u(i, j, k).z)));
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
          rTail(i, j, k) = q(i, j, k - 1) + extrap_factor * dr_ds(i, j, k - 1);
          rHead(i, j, k) = q(i, j, k) - extrap_factor * dr_ds(i, j, k);
          uTail(i, j, k).x =
            u(i, j, k - 1).x + extrap_factor * du_ds(i, j, k - 1).x;
          uHead(i, j, k).x = u(i, j, k).x - extrap_factor * du_ds(i, j, k).x;
          uTail(i, j, k).y =
            u(i, j, k - 1).y + extrap_factor * du_ds(i, j, k - 1).y;
          uHead(i, j, k).y = u(i, j, k).y - extrap_factor * du_ds(i, j, k).y;
          uTail(i, j, k).z =
            u(i, j, k - 1).z + extrap_factor * du_ds(i, j, k - 1).z;
          uHead(i, j, k).z = u(i, j, k).z - extrap_factor * du_ds(i, j, k).z;
          pTail(i, j, k) = p(i, j, k - 1) + extrap_factor * dp_ds(i, j, k - 1);
          pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

          // Conserved quantities
          ruTail(i, j, k).x = rTail(i, j, k) * uTail(i, j, k).x;
          ruHead(i, j, k).x = rHead(i, j, k) * uHead(i, j, k).x;
          ruTail(i, j, k).y = rTail(i, j, k) * uTail(i, j, k).y;
          ruHead(i, j, k).y = rHead(i, j, k) * uHead(i, j, k).y;
          ruTail(i, j, k).z = rTail(i, j, k) * uTail(i, j, k).z;
          ruHead(i, j, k).z = rHead(i, j, k) * uHead(i, j, k).z;
          rETail(i, j, k) =
            p2rEFactor * pTail(i, j, k) +
            0.5 * rTail(i, j, k) *
              (utils::sqr(uTail(i, j, k).x) + utils::sqr(uTail(i, j, k).y) +
                utils::sqr(uTail(i, j, k).z));
          rEHead(i, j, k) =
            p2rEFactor * pHead(i, j, k) +
            0.5 * rHead(i, j, k) *
              (utils::sqr(uHead(i, j, k).x) + utils::sqr(uHead(i, j, k).y) +
                utils::sqr(uHead(i, j, k).z));
        } // for
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Compute corrector setp: HLL Riemann solver.
     *------------------------------------------------------------------------*/

    for(auto k : m.cells<mesh::z_axis, mesh::corrector>()) {
      for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
          const double cT = std::sqrt(gamma * pTail(i, j, k) / rTail(i, j, k));
          const double cH = std::sqrt(gamma * pHead(i, j, k) / rHead(i, j, k));

          // Update min/max eigenvalues on tail face.
          const double LminT = uTail(i, j, k).z - cT;
          const double LmaxT = uTail(i, j, k).z + cT;

          // Update min/max eigenvalues on head face.
          const double LminH = uHead(i, j, k).z - cH;
          const double LmaxH = uHead(i, j, k).z + cH;

          // Values for HLL.
          const double Lminus{std::min(LminH, std::min(LminT, double{0.0}))};
          const double Lplus{std::max(LmaxH, std::max(LmaxT, double{0.0}))};
          const double Ldiv{1.0 / (Lplus - Lminus)};
          const double Lmult{Lplus * Lminus};
          const double delta_r{rHead(i, j, k) - rTail(i, j, k)};
          const vec3 delta_ru{ruHead(i, j, k).x - ruTail(i, j, k).x,
            ruHead(i, j, k).y - ruTail(i, j, k).y,
            ruHead(i, j, k).z - ruTail(i, j, k).z};
          const double delta_rE{rEHead(i, j, k) - rETail(i, j, k)};

          const double f_r_T{ruTail(i, j, k).z};
          const double f_r_H{ruHead(i, j, k).z};
          const vec3 f_ru_T{
            rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).z,
            rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z,
            ruTail(i, j, k).z * uTail(i, j, k).z + pTail(i, j, k)};
          const vec3 f_ru_H{
            rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).z,
            rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z,
            ruHead(i, j, k).z * uHead(i, j, k).z + pHead(i, j, k)};
          const double f_rE_T{
            (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).z};
          const double f_rE_H{
            (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).z};

          rF(i, j, k) =
            Ldiv * (Lplus * f_r_T - Lminus * f_r_H + Lmult * delta_r);
          ruF(i, j, k).x =
            Ldiv * (Lplus * f_ru_T.x - Lminus * f_ru_H.x + Lmult * delta_ru.x);
          ruF(i, j, k).y =
            Ldiv * (Lplus * f_ru_T.y - Lminus * f_ru_H.y + Lmult * delta_ru.y);
          ruF(i, j, k).z =
            Ldiv * (Lplus * f_ru_T.z - Lminus * f_ru_H.z + Lmult * delta_ru.z);
          rEF(i, j, k) =
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
          r(i, j, k) =
            r(i, j, k) - update_factor * (rF(i, j, k + 1) - rF(i, j, k));
          ru(i, j, k).x = ru(i, j, k).x -
                          update_factor * (ruF(i, j, k + 1).x - ruF(i, j, k).x);
          ru(i, j, k).y = ru(i, j, k).y -
                          update_factor * (ruF(i, j, k + 1).y - ruF(i, j, k).y);
          ru(i, j, k).z = ru(i, j, k).z -
                          update_factor * (ruF(i, j, k + 1).z - ruF(i, j, k).z);
          rE(i, j, k) =
            rE(i, j, k) - update_factor * (rEF(i, j, k + 1) - rEF(i, j, k));
        } // for
      } // for
    } // for
  } // if
} // advance

} // namespace muscl::tasks::hydro

#endif // MUSCL_TASKS_HYDRO_HH
