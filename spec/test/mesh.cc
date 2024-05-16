#include "../mesh.hh"

#include <flecsi/util/unit.hh>

using namespace flecsi;
using namespace fvm;

const flecsi::field<double>::definition<mesh, mesh::cells> rho;

void
init_mesh(mesh::accessor<ro> m, field<double>::accessor<wo, wo> f_a) {
} // init_mesh

int
verify_mesh(mesh::accessor<ro> m, field<double>::accessor<wo, wo> f_a) {
  UNIT("TASK") {
  };
} // verify_mesh

int
fvm_mesh() {
  UNIT("DRIVER") {
  }; // UNIT
} // mesh

flecsi::util::unit::driver<fvm_mesh> driver;
