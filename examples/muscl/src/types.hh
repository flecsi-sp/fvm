#ifndef MUSCL_TYPES_HH
#define MUSCL_TYPES_HH

#include <flecsi/data.hh>
#include <fvm/mesh.hh>

namespace muscl {
inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fvm::mesh;
using index = flecsi::topo::index;
using global = flecsi::topo::global;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;
template<typename T>
using single = field<T, flecsi::data::single>;

struct vec3 {
  double x, y, z;
};

inline std::ostream &
operator<<(std::ostream & s, vec3 const & v) {
  s << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  return s;
} // operator<<

} // namespace muscl

#endif // MUSCL_TYPES_HH
