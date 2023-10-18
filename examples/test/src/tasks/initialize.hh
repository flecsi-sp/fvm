#ifndef TEST_TASKS_INITIALIZE_HH
#define TEST_TASKS_INITIALIZE_HH

#include "../state.hh"

namespace test::tasks::init {
void monotonic(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a);

void color(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> f_a);

} // namespace test::tasks::init

#endif // TEST_TASKS_INITIALIZE_HH
