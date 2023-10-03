#ifndef MUSCL_TASKS_INITIALIZE_HH
#define MUSCL_TASKS_INITIALIZE_HH

#include "../types.hh"

namespace muscl::tasks {

void check(mesh::accessor<ro> m);

void sod(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ra,
  field<double>::accessor<rw, ro> rua,
  field<double>::accessor<rw, ro> rva,
  field<double>::accessor<rw, ro> rwa,
  field<double>::accessor<rw, ro> rEa,
  double gamma);

} // namespace muscl::tasks

#endif // MUSCL_TASKS_INITIALIZE_HH
