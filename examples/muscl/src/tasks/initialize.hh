#ifndef MUSCL_TASKS_INITIALIZE_HH
#define MUSCL_TASKS_INITIALIZE_HH

#include "../types.hh"

namespace muscl::tasks::init {

void gamma(single<double>::accessor<wo> gamma_a, double g);

void boundaries(single<mesh::bmap>::accessor<wo> bmap_a,
  mesh::boundary_type xlow,
  mesh::boundary_type xhigh,
  mesh::boundary_type ylow,
  mesh::boundary_type yhigh,
  mesh::boundary_type zlow,
  mesh::boundary_type zhigh);

void sod(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<velocity>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  single<double>::accessor<ro> gamma_a);

} // namespace muscl::tasks::init

#endif // MUSCL_TASKS_INITIALIZE_HH
