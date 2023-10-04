#ifndef MUSCL_ADVANCE_HH
#define MUSCL_ADVANCE_HH

#include "control.hh"

namespace muscl::action {

void advance(control_policy & cp);
inline control::action<advance, cp::advance> advance_action;

} // namespace muscl::action

#endif // MUSCL_ADVANCE_HH
