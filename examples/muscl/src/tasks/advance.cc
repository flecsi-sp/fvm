#include "advance.hh"

#include <limits>

using namespace muscl;
using namespace muscl::tasks;

double
dtmin(mesh::accessor<ro> m, single<velocity>::accessor<ro> lmax) {
  return std::min(m.delta<mesh::x_axis>() / lmax->x,
    std::min(
      m.delta<mesh::y_axis>() / lmax->y, m.delta<mesh::z_axis>() / lmax->z));
} // dtmin
