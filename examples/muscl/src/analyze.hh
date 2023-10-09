#ifndef MUSCL_ANALYZE_HH
#define MUSCL_ANALYZE_HH

#include "control.hh"

namespace muscl::action {

void analyze(control_policy & cp);
inline control::action<analyze, cp::analyze> analyze_action;

} // namespace muscl::action

#endif // MUSCL_ANALYZE_HH
