#ifndef MUSCL_TASKS_INITIALIZE_HH
#define MUSCL_TASKS_INITIALIZE_HH

#include "../types.hh"

namespace muscl::tasks::init {

void check(mesh::accessor<ro> m);

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

void primitives(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<velocity>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  field<velocity>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<velocity>::accessor<wo> lmax_a,
  single<double>::accessor<ro> gamma_a);

} // namespace muscl::tasks::init

#endif // MUSCL_TASKS_INITIALIZE_HH
