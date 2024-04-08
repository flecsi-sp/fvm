#ifndef MUSCL_FINALIZE_HH
#define MUSCL_FINALIZE_HH

#include "control.hh"

namespace muscl::action {

void finalize(control_policy & cp);
inline control::action<finalize, cp::finalize> finalize_action;

} // namespace muscl::action

#endif // MUSCL_FINALIZE_HH
