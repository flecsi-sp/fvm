#include "finalize.hh"
#include "state.hh"

using namespace muscl;

void action::finalize(muscl::control_policy & cp) {
  m.deallocate();
} // action::finalize
