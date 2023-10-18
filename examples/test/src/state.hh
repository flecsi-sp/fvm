#ifndef TEST_STATE_HH
#define TEST_STATE_HH

#include <flecsi/data.hh>
#include <fvm/mesh.hh>

namespace test {

inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fvm::mesh;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;

inline mesh::slot m;
inline const field<double>::definition<mesh, mesh::cells> f;

} // namespace test

#endif // TEST_STATE_HH
