#ifndef MUSCL_TASKS_INITIALIZE_HH
#define MUSCL_TASKS_INITIALIZE_HH

#include "../types.hh"

namespace muscl::tasks::init {

void gamma(single<double>::accessor<wo> gamma_a, double g);

inline void
touch(mesh::accessor<ro> m,
  field<double>::accessor<wo, wo> r_a,
  field<vec3>::accessor<wo, wo> ru_a,
  field<double>::accessor<wo, wo> rE_a,
  field<vec3>::accessor<wo, wo> u_a,
  field<double>::accessor<wo, wo> p_a,
  field<double>::accessor<wo, wo> q_a,
  field<vec3>::accessor<wo, wo> qu_a,
  field<double>::accessor<wo, wo> qE_a,
  field<double>::accessor<wo, wo> dr_ds_a,
  field<vec3>::accessor<wo, wo> du_ds_a,
  field<double>::accessor<wo, wo> dp_ds_a,
  field<double>::accessor<wo, wo> rTail_a,
  field<vec3>::accessor<wo, wo> ruTail_a,
  field<double>::accessor<wo, wo> rETail_a,
  field<vec3>::accessor<wo, wo> uTail_a,
  field<double>::accessor<wo, wo> pTail_a,
  field<double>::accessor<wo, wo> rHead_a,
  field<vec3>::accessor<wo, wo> ruHead_a,
  field<double>::accessor<wo, wo> rEHead_a,
  field<vec3>::accessor<wo, wo> uHead_a,
  field<double>::accessor<wo, wo> pHead_a,
  field<double>::accessor<wo, wo> rF_a,
  field<vec3>::accessor<wo, wo> ruF_a,
  field<double>::accessor<wo, wo> rEF_a) {} // touch

inline void
touch1(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  field<vec3>::accessor<rw, ro> u_a,
  field<double>::accessor<rw, ro> p_a,
  field<double>::accessor<rw, ro> q_a,
  field<vec3>::accessor<rw, ro> qu_a,
  field<double>::accessor<rw, ro> qE_a,
  field<double>::accessor<rw, ro> dr_ds_a,
  field<vec3>::accessor<rw, ro> du_ds_a,
  field<double>::accessor<rw, ro> dp_ds_a,
  field<double>::accessor<rw, ro> rTail_a,
  field<vec3>::accessor<rw, ro> ruTail_a,
  field<double>::accessor<rw, ro> rETail_a,
  field<vec3>::accessor<rw, ro> uTail_a,
  field<double>::accessor<rw, ro> pTail_a,
  field<double>::accessor<rw, ro> rHead_a,
  field<vec3>::accessor<rw, ro> ruHead_a,
  field<double>::accessor<rw, ro> rEHead_a,
  field<vec3>::accessor<rw, ro> uHead_a,
  field<double>::accessor<rw, ro> pHead_a,
  field<double>::accessor<rw, ro> rF_a,
  field<vec3>::accessor<rw, ro> ruF_a,
  field<double>::accessor<rw, ro> rEF_a) {} // touch

void boundaries(single<mesh::bmap>::accessor<wo> bmap_a,
  mesh::boundary_type xlow,
  mesh::boundary_type xhigh,
  mesh::boundary_type ylow,
  mesh::boundary_type yhigh,
  mesh::boundary_type zlow,
  mesh::boundary_type zhigh);

void sod(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  single<double>::accessor<ro> gamma_a);

void monotonic(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  single<double>::accessor<ro> gamma_a);

void color(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  single<double>::accessor<ro> gamma_a);

} // namespace muscl::tasks::init

#endif // MUSCL_TASKS_INITIALIZE_HH
