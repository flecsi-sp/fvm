#ifndef FVM_MESH_HH
#define FVM_MESH_HH

#include <flecsi/data.hh>
#include <flecsi/execution.hh>
#include <flecsi/flog.hh>
#include <flecsi/topo/narray/coloring_utils.hh>

namespace fvm {

struct mesh : flecsi::topo::specialization<flecsi::topo::narray, mesh> {

  /*--------------------------------------------------------------------------*
    Policy Information.
   *--------------------------------------------------------------------------*/

  enum index_space { cells };
  using index_spaces = has<cells>;
  enum domain { quantities, predictor, corrector, all, global };
  enum axis { x_axis, y_axis, z_axis };
  using axes = has<x_axis, y_axis, z_axis>;
  enum boundary { low, high };
  enum boundary_type { inflow, outflow, reflecting, periodic };

  using coord = base::coord;
  using gcoord = base::gcoord;
  using bmap = std::array<std::array<boundary_type, 2>, 3>;
  using grect = std::array<std::array<double, 3>, 3>;
  using colors = base::colors;
  using axis_definition = base::axis_definition;
  using index_definition = base::index_definition;

  struct meta_data {
    double xdelta;
    double ydelta;
    double zdelta;
  };

  static constexpr std::size_t dimension = 3;

  template<auto>
  static constexpr std::size_t privilege_count = 2;

  /*--------------------------------------------------------------------------*
    Interface.
   *--------------------------------------------------------------------------*/

  template<class B>
  struct interface : B {

    template<axis A, domain DM = quantities>
    std::size_t size() {
      if constexpr(DM == quantities) {
        return B::template size<mesh::cells, A, base::domain::logical>();
      }
      else if constexpr(DM == predictor) {
        return B::template size<mesh::cells, A, base::domain::logical>() + 2;
      }
      else if constexpr(DM == corrector) {
        return B::template size<mesh::cells, A, base::domain::logical>() + 1;
      }
      else if constexpr(DM == all) {
        return B::template size<mesh::cells, A, base::domain::all>();
      }
      else if constexpr(DM == global) {
        return B::template size<mesh::cells, A, base::domain::global>();
      } // if
    } // size

    template<axis A, domain DM = quantities>
    FLECSI_INLINE_TARGET auto cells() const {
      if constexpr(DM == quantities) {
        return B::template range<mesh::cells, A, base::domain::logical>();
      }
      else if constexpr(DM == predictor) {
        return flecsi::topo::make_ids<mesh::cells>(
          flecsi::util::iota_view<flecsi::util::id>(
            1, B::template size<mesh::cells, A, base::domain::all>() - 1));
      }
      else if constexpr(DM == corrector) {
        return flecsi::topo::make_ids<mesh::cells>(
          flecsi::util::iota_view<flecsi::util::id>(
            2, B::template size<mesh::cells, A, base::domain::all>() - 1));
      }
      else if constexpr(DM == all) {
        return B::template range<mesh::cells, A, base::domain::all>();
      } // if
    } // cells

    void set_geometry(double x, double y, double z) { // available if writable
      this->policy_meta() = {x, y, z};
    } // set_geometry

    template<axis A>
    FLECSI_INLINE_TARGET std::size_t global_id(std::size_t i) const {
      return B::template global_id<mesh::cells, A>(i);
    } // global_id

    template<axis A>
    FLECSI_INLINE_TARGET double delta() const {
      if constexpr(A == x_axis) {
        return this->policy_meta().xdelta;
      }
      else if constexpr(A == y_axis) {
        return this->policy_meta().ydelta;
      }
      else if constexpr(A == z_axis) {
        return this->policy_meta().zdelta;
      } // if
    } // delta

    /*!
      Return the cell head for the given axis and id. The head is the trailing
      interface of the cell.
     */
    template<axis A>
    FLECSI_INLINE_TARGET double head(std::size_t i) const {
      return center<A>(i) - 0.5 * delta<A>();
    } // center

    /*!
      Return the cell center for the given axis and id.
     */
    template<axis A>
    FLECSI_INLINE_TARGET double center(std::size_t i) const {
      return delta<A>() * global_id<A>(i) + 0.5 * delta<A>();
    } // center

    /*!
      Return the cell tail for the given axis and id. The tail is the leading
      interface of the cell.
     */
    template<axis A>
    FLECSI_INLINE_TARGET double tail(std::size_t i) const {
      return center<A>(i) + 0.5 * delta<A>();
    } // center

  }; // interface

  /*--------------------------------------------------------------------------*
    Color Task.
   *--------------------------------------------------------------------------*/

  static coloring color(std::size_t num_colors, gcoord axis_extents) {
    index_definition idef;
    idef.axes = flecsi::topo::narray_utils::make_axes(num_colors, axis_extents);
    std::size_t ai{0};

    for(auto & a : idef.axes) {
      a.hdepth = 2;
      a.bdepth = 2;
    } // for

    return {{idef}};
  } // color

  /*--------------------------------------------------------------------------*
    Initialization.
   *--------------------------------------------------------------------------*/

  // FIXME: colors != processes
  static void set_geometry(mesh::accessor<flecsi::rw> sm, grect const & g) {
    sm.set_geometry((g[0][1] - g[0][0]) / sm.size<x_axis, global>(),
      (g[1][1] - g[1][0]) / sm.size<y_axis, global>(),
      (g[2][1] - g[2][0]) / sm.size<z_axis, global>());
  } // set_geometry

  static void initialize(flecsi::data::topology_slot<mesh> & s,
    coloring const &,
    grect const & geometry) {
    flecsi::execute<set_geometry, flecsi::mpi>(s, geometry);
  } // initialize

}; // struct mesh

} // namespace fvm

#endif // FVM_MESH_HH
