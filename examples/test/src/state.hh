#ifndef TEST_STATE_HH
#define TEST_STATE_HH

#include <flecsi/data.hh>
#include <fvm/mesh.hh>

namespace test {

inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fvm::mesh;
using global = flecsi::topo::global;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;
template<typename T>
using single = field<T, flecsi::data::single>;

inline mesh::slot m;
// Convenience slot (shorter name) using flecsi's global topology instance.
inline global::slot & gt = flecsi::global_topology;
inline const field<double>::definition<mesh, mesh::cells> f;
inline const single<mesh::bmap>::definition<flecsi::topo::global> bmap;

} // namespace test

#endif // TEST_STATE_HH
