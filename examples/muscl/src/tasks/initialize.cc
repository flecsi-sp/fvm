#include "initialize.hh"
#include "../utils.hh"

#include <flecsi/flog.hh>

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

        if(x < sodx0) {
          r[k][j][i] = sodrL;
          ru[k][j][i].x = sodrL * soduL;
          ru[k][j][i].y = sodrL * sodvL;
          ru[k][j][i].z = sodrL * sodwL;
          rE[k][j][i] =
            mult * sodpL +
            0.5 * sodrL *
              (utils::sqr(soduL) + utils::sqr(sodvL) + utils::sqr(sodwL));
        }
        else {
          r[k][j][i] = sodrR;
          ru[k][j][i].x = sodrR * soduR;
          ru[k][j][i].y = sodrR * sodvR;
          ru[k][j][i].z = sodrR * sodwR;
          rE[k][j][i] =
            mult * sodpR +
            0.5 * sodrR *
              (utils::sqr(soduR) + utils::sqr(sodvR) + utils::sqr(sodwR));
        } // if
      } // for
    } // for
  } // for
} // sod
