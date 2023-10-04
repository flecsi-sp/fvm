#ifndef MUSCL_TASKS_INITIALIZE_HH
#define MUSCL_TASKS_INITIALIZE_HH

#include "../types.hh"

namespace muscl::tasks {

void check(mesh::accessor<ro> m);

void sod(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ra,
  field<velocity>::accessor<rw, ro> rua,
  field<double>::accessor<rw, ro> rEa,
  double gamma);

void init(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  field<velocity>::accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  field<velocity>::accessor<wo, ro> u_a,
  field<double>::accessor<wo, ro> p_a,
  single<velocity>::accessor<wo> lmax,
  double gamma);

} // namespace muscl::tasks

#endif // MUSCL_TASKS_INITIALIZE_HH
