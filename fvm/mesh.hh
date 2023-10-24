#ifndef FVM_MESH_HH
#define FVM_MESH_HH

#include <flecsi/data.hh>
#include <flecsi/execution.hh>
#include <flecsi/flog.hh>
#include <flecsi/topo/narray/coloring_utils.hh>

#include <ranges>

namespace fvm {

/*!
  FIXME
 */
struct mesh : flecsi::topo::specialization<flecsi::topo::narray, mesh> {

  /*--------------------------------------------------------------------------*
    Policy Information.
   *--------------------------------------------------------------------------*/

  enum index_space { cells };
  using index_spaces = has<cells>;

  /// Mesh domains.
  /// The domain identifies the supported iteration spaces on the mesh.
  enum domain {
    /// This domain includes the locations of the unknowns of the problem state.
    quantities,
    /// This domain includes all mesh locations where it is possible to compute
    /// a 2nd-order slope approximation.
    predictor,
    /// This domain includes all mesh locations where it is possible to compute
    /// a flux across an interface from reconstructed predictor quantities.
    corrector,
    /// This domain includes all local cells, including halo and boundary cells.
    all,
    /// This domain includes all global cells.
    global
  };

  /// Mesh axes.
  /// The axis identifies a Cartesian coordinate axis on the mesh.
  enum axis {
    /// X-coordinate axis.
    x_axis,
    /// Y-coordinate axis.
    y_axis,
    /// Z-coordinate axis.
    z_axis
  };

  using axes = has<x_axis, y_axis, z_axis>;

  /// Boundary.
  /// Identifies the low or high boundary along a particular axis.
  enum boundary {
    /// Low axis boundary.
    low,
    /// High axis boundary.
    high
  };

  /// Boundary type.
  /// The supported boundary types are derived from Leveque's
  /// [Finite Volume Methods for Hyperbolc
  /// Problems](https://www.amazon.com/Methods-Hyperbolic-Problems-Cambridge-Mathematics/dp/0521009243).
  /// All boundary conditions are implemented using the \em ghost-cell approach
  /// outlined in the text.
  enum boundary_type {
    /// Zero-order extrapolation from the interior solution.
    flow,
    /// Cauchy problem solution for solid wall boundary.
    reflecting,
    /// Substitution with opposite interior solution.
    periodic
  };

  using periodic_axes = std::array<bool, 3>;
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

  template<auto S, typename T>
  static auto make_ids(T && t) {
    return std::ranges::transform_view(std::forward<T>(t),
      [](auto const & i) { return flecsi::topo::id<S>(i); });
  } // make_ids

  /*--------------------------------------------------------------------------*
    Interface.
   *--------------------------------------------------------------------------*/

  /// Mesh Interface.
  template<class B>
  struct interface : B {

    /// Return the size for the given axis and domain.
    /// @tparam A  The mesh axis.
    /// @tparam DM The mesh domain.
    template<axis A, domain DM = quantities>
    auto size() {
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

    /// Return a range over the given axis and domain.
    /// @tparam A  The mesh axis.
    /// @tparam DM The mesh domain (excluding \em global).
    ///
    /// The range can be used to iterate over the cells in the given domain,
    /// e.g.:
    /// @code
    ///   for(auto k: m.cells<mesh::z_axis, mesh::quantities>()) {
    ///     for(auto j: m.cells<mesh::y_axis, mesh::quantities>()) {
    ///       for(auto i: m.cells<mesh::x_axis, mesh::quantities>()) {
    ///         u[k][j][i] = 1.0;
    ///       } // for
    ///     } // for
    ///   } // for
    /// @endcode
    template<axis A, domain DM = quantities, bool R = false>
    FLECSI_INLINE_TARGET auto cells() const {
      flecsi::util::id b, e;

      if constexpr(DM == quantities) {
        b = 2;
        e = B::template size<mesh::cells, A, base::domain::all>() - 2;
      }
      else if constexpr(DM == predictor) {
        b = 1;
        e = B::template size<mesh::cells, A, base::domain::all>() - 1;
      }
      else if constexpr(DM == corrector) {
        b = 2;
        e = B::template size<mesh::cells, A, base::domain::all>() - 1;
      }
      else if constexpr(DM == all) {
        b = 0;
        e = B::template size<mesh::cells, A, base::domain::all>();
      }
      else if(DM == global) {
        flog_fatal(
          "illegal domain: you cannot iterate over the global domain.");
      } // if

      if constexpr(R) {
        return make_ids<mesh::cells>(
          std::ranges::iota_view{b, e} | std::views::reverse);
      }
      else {
        return make_ids<mesh::cells>(std::ranges::iota_view{b, e});
      } // if
    } // cells

    void set_geometry(double x, double y, double z) { // available if writable
      this->policy_meta() = {x, y, z};
    } // set_geometry

    /// Return the global id of \em i for the given axis.
    /// @tparam A The coordinate axis.
    template<axis A>
    FLECSI_INLINE_TARGET std::size_t global_id(std::size_t i) const {
      return B::template global_id<mesh::cells, A>(i);
    } // global_id

    /// Return the mesh spacing for the given axis.
    /// @tparam A  The coordinate axis.
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

    /// Return the cell head for the given axis and id. The head is the trailing
    /// interface of the cell.
    /// @tparam A  The coordinate axis.
    template<axis A>
    FLECSI_INLINE_TARGET double head(std::size_t i) const {
      return center<A>(i) - 0.5 * delta<A>();
    } // center

    /// Return the cell center for the given axis and id.
    /// @tparam A  The coordinate axis.
    template<axis A>
    FLECSI_INLINE_TARGET double center(std::size_t i) const {
      return delta<A>() * global_id<A>(i) + 0.5 * delta<A>();
    } // center

    /// Return the cell tail for the given axis and id. The tail is the leading
    /// interface of the cell.
    /// @tparam A  The coordinate axis.
    template<axis A>
    FLECSI_INLINE_TARGET double tail(std::size_t i) const {
      return center<A>(i) + 0.5 * delta<A>();
    } // center

    /// Return true if the current color is a low-edge partition for the given
    /// axis.
    /// @tparam A  The coordinate axis.
    template<axis A>
    bool is_low() {
      return B::template is_low<mesh::cells, A>();
    } // is_low

    /// Return true if the current color is a high-edge partition for the given
    /// axis.
    /// @tparam A  The coordinate axis.
    template<axis A>
    bool is_high() {
      return B::template is_high<mesh::cells, A>();
    } // is_high

  }; // interface

  /*--------------------------------------------------------------------------*
    Color Task.
   *--------------------------------------------------------------------------*/

  static coloring
  color(std::size_t num_colors, gcoord axis_extents, periodic_axes p) {
    index_definition idef;
    idef.axes = flecsi::topo::narray_utils::make_axes(num_colors, axis_extents);
    std::size_t ai{0};

    for(auto & a : idef.axes) {
      a.hdepth = 2;
      a.bdepth = 2;
      a.periodic = p[ai];
    } // for

    return {{idef}};
  } // color

  /*--------------------------------------------------------------------------*
    Initialization.
   *--------------------------------------------------------------------------*/

  static void set_geometry(mesh::accessor<flecsi::rw> m, grect const & g) {
    m.set_geometry((g[0][1] - g[0][0]) / m.size<x_axis, global>(),
      (g[1][1] - g[1][0]) / m.size<y_axis, global>(),
      (g[2][1] - g[2][0]) / m.size<z_axis, global>());
  } // set_geometry

  static void initialize(flecsi::data::topology_slot<mesh> & s,
    coloring const &,
    grect const & geometry) {
    flecsi::execute<set_geometry>(s, geometry);
  } // initialize

}; // struct mesh

} // namespace fvm

#endif // FVM_MESH_HH
