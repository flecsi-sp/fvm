#ifndef MUSCL_INITIALIZE_HH
#define MUSCL_INITIALIZE_HH

#include "control.hh"

namespace muscl {
namespace action {
void initialize(muscl::control_policy & cp);
inline control::action<initialize, cp::initialize> initialize_action;
} // namespace action
} // namespace muscl

#endif // MUSCL_INITIALIZE_HH
