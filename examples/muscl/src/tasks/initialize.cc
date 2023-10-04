#include "initialize.hh"

#include <flecsi/flog.hh>

/*----------------------------------------------------------------------------*
  Sanity check.
 *----------------------------------------------------------------------------*/

void
muscl::tasks::check(mesh::accessor<ro> m) {

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
}

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
muscl::tasks::sod(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ra,
  field<velocity>::accessor<rw, ro> rua,
  field<double>::accessor<rw, ro> rEa,
  double gamma) {
  auto r = m.mdspan<mesh::cells>(ra);
  auto ru = m.mdspan<mesh::cells>(rua);
  auto rE = m.mdspan<mesh::cells>(rEa);

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
