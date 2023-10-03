#ifndef MUSCL_TYPES_HH
#define MUSCL_TYPES_HH

#include "fvm/mesh.hh"

#include <flecsi/data.hh>

namespace muscl {
inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fvm::mesh;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;

} // namespace muscl
  //
#endif // MUSCL_TYPES_HH
