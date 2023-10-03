#include "initialize.hh"

#include <flecsi/flog.hh>

void muscl::tasks::check(mesh::accessor<ro> m) {

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    const auto z = m.value<mesh::z_axis>(k);
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      const auto y = m.value<mesh::y_axis>(j);
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        const auto x = m.value<mesh::x_axis>(i);

        flog(info) << "(" << i << "," << j << "," << k << ")" << std::endl;
      } // for
    } // for
  } // for
}

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

void muscl::tasks::sod(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ra,
  field<double>::accessor<rw, ro> rua,
  field<double>::accessor<rw, ro> rva,
  field<double>::accessor<rw, ro> rwa,
  field<double>::accessor<rw, ro> rEa,
  double gamma) {
  auto r = m.mdspan<mesh::cells>(ra);
  auto ru = m.mdspan<mesh::cells>(rua);
  auto rv = m.mdspan<mesh::cells>(rva);
  auto rw = m.mdspan<mesh::cells>(rwa);
  auto rE = m.mdspan<mesh::cells>(rEa);

  const double mult = 1.0/(gamma - 1.0);

  for(auto k : m.cells<mesh::z_axis, mesh::quantities>()) {
    for(auto j : m.cells<mesh::y_axis, mesh::quantities>()) {
      for(auto i : m.cells<mesh::x_axis, mesh::quantities>()) {
        const auto x = m.value<mesh::x_axis>(i);

        if(x<sodx0) {
          r[k][j][i] = sodrL;
          ru[k][j][i] = soduL;
          rv[k][j][i] = sodvL;
          rw[k][j][i] = sodwL;
          rE[k][j][i] =
            mult*sodpL + 0.5*sodrL*(soduL*soduL+sodvL*sodvL+sodwL*sodwL);
        }
        else {
          r[k][j][i] = sodrR;
          ru[k][j][i] = soduR;
          rv[k][j][i] = sodvR;
          rw[k][j][i] = sodwR;
          rE[k][j][i] =
            mult*sodpR + 0.5*sodrR*(soduR*soduR+sodvR*sodvR+sodwR*sodwR);
        } // if
      } // for
    } // for
  } // for
} // sod
