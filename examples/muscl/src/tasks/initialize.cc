#include "initialize.hh"
#include "../utils.hh"

#include <flecsi/flog.hh>

/*----------------------------------------------------------------------------*
  Sanity check.
 *----------------------------------------------------------------------------*/

void
muscl::tasks::init::check(mesh::accessor<ro> m) {

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    const auto z = m.center<mesh::z_axis>(k);
    const auto ztail = m.tail<mesh::z_axis>(k);
    const auto zhead = m.head<mesh::z_axis>(k);
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      const auto y = m.center<mesh::y_axis>(j);
      const auto ytail = m.tail<mesh::y_axis>(j);
      const auto yhead = m.head<mesh::y_axis>(j);
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        const auto x = m.center<mesh::x_axis>(i);
        const auto xtail = m.tail<mesh::x_axis>(i);
        const auto xhead = m.head<mesh::x_axis>(i);

        flog(info) << "(" << i << "," << j << "," << k << ")" << std::endl;
        flog(info) << "x: " << xhead << " -> " << x << " -> " << xtail
                   << std::endl;
        flog(info) << "y: " << yhead << " -> " << y << " -> " << ytail
                   << std::endl;
        flog(info) << "z: " << zhead << " -> " << z << " -> " << ztail
                   << std::endl;
        flog(info) << std::endl;
      } // for
    } // for
  } // for
} // check

/*----------------------------------------------------------------------------*
  Gamma.
 *----------------------------------------------------------------------------*/

void
muscl::tasks::init::gamma(single<double>::accessor<wo> gamma_a, double g) {
  (*gamma_a) = g;
}

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

/*----------------------------------------------------------------------------*
  Sod Shock Tube.
 *----------------------------------------------------------------------------*/

static const double sodrL = 1.0;
static const double soduL = 0.0;
static const double sodvL = 0.0;
static const double sodwL = 0.0;
static const double sodpL = 1.0;

static const double sodrR = 0.125;
static const double soduR = 0.0;
static const double sodvR = 0.0;
static const double sodwR = 0.0;
static const double sodpR = 0.1;

static const double sodx0 = 0.5;

void
muscl::tasks::init::sod(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<velocity>::accessor<rw, ro> ru_a,
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
        const auto x = m.tail<mesh::x_axis>(i);

        if(x < sodx0) {
          r[k][j][i] = sodrL;
          ru[k][j][i].x = soduL;
          ru[k][j][i].y = sodvL;
          ru[k][j][i].z = sodwL;
          rE[k][j][i] =
            mult * sodpL +
            0.5 * sodrL * (soduL * soduL + sodvL * sodvL + sodwL * sodwL);
        }
        else {
          r[k][j][i] = sodrR;
          ru[k][j][i].x = soduR;
          ru[k][j][i].y = sodvR;
          ru[k][j][i].z = sodwR;
          rE[k][j][i] =
            mult * sodpR +
            0.5 * sodrR * (soduR * soduR + sodvR * sodvR + sodwR * sodwR);
        } // if
      } // for
    } // for
  } // for
} // sod

/*----------------------------------------------------------------------------*
  Hydro initialization.
 *----------------------------------------------------------------------------*/

void
muscl::tasks::init::primitives(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<velocity>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  field<velocity>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<velocity>::accessor<wo> lmax_a,
  single<double>::accessor<ro> gamma_a) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto ru = m.mdspan<mesh::cells>(ru_a);
  auto rE = m.mdspan<mesh::cells>(rE_a);
  auto u = m.mdspan<mesh::cells>(u_a);
  auto p = m.mdspan<mesh::cells>(p_a);
  auto & lmax = *lmax_a;
  auto const gamma = *gamma_a;

  // Initialize primitive quantities.
  for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
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

  // Iinitialize max eigenvalues for dt.
  lmax.x = std::numeric_limits<double>::min();
  lmax.y = std::numeric_limits<double>::min();
  lmax.z = std::numeric_limits<double>::min();
  for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
        const double c = std::sqrt(gamma * p[k][j][i] / r[k][j][i]);
        lmax.x = std::max(std::abs(u[k][j][i].x) + c, lmax.x);
        lmax.y = std::max(std::abs(u[k][j][i].y) + c, lmax.y);
        lmax.z = std::max(std::abs(u[k][j][i].z) + c, lmax.z);
      } // for
    } // for
  } // for
} // init
