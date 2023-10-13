#ifndef MUSCL_BOUNDARY_HH
#define MUSCL_BOUNDARY_HH

#include "state.hh"

namespace muscl {

template<mesh::axis A, mesh::boundary B>
void
flow(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto ru = m.mdspan<mesh::cells>(ru_a);
  auto rE = m.mdspan<mesh::cells>(rE_a);

  if constexpr(A == mesh::axis::x_axis) {
    if(B == mesh::boundary::low) {
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
          r[k][j][0] = r[k][j][2];
          r[k][j][1] = r[k][j][2];
          ru[k][j][0].x = ru[k][j][2].x;
          ru[k][j][1].x = ru[k][j][2].x;
          ru[k][j][0].y = ru[k][j][2].y;
          ru[k][j][1].y = ru[k][j][2].y;
          ru[k][j][0].z = ru[k][j][2].z;
          ru[k][j][1].z = ru[k][j][2].z;
          rE[k][j][0] = rE[k][j][2];
          rE[k][j][1] = rE[k][j][2];
        } // for
      } // for
    }
    else {
      const std::size_t i = m.size<mesh::x_axis, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
          r[k][j][i - 1] = r[k][j][i - 3];
          r[k][j][i - 2] = r[k][j][i - 3];
          ru[k][j][i - 1].x = ru[k][j][i - 3].x;
          ru[k][j][i - 2].x = ru[k][j][i - 3].x;
          ru[k][j][i - 1].y = ru[k][j][i - 3].y;
          ru[k][j][i - 2].y = ru[k][j][i - 3].y;
          ru[k][j][i - 1].z = ru[k][j][i - 3].z;
          ru[k][j][i - 2].z = ru[k][j][i - 3].z;
          rE[k][j][i - 1] = rE[k][j][i - 3];
          rE[k][j][i - 2] = rE[k][j][i - 3];
        } // for
      } // for
    } // if
  }
  else if constexpr(A == mesh::axis::y_axis) {
    if(B == mesh::boundary::low) {
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[k][0][i] = r[k][2][i];
          r[k][1][i] = r[k][2][i];
          ru[k][0][i].x = ru[k][2][i].x;
          ru[k][1][i].x = ru[k][2][i].x;
          ru[k][0][i].y = ru[k][2][i].y;
          ru[k][1][i].y = ru[k][2][i].y;
          ru[k][0][i].z = ru[k][2][i].z;
          ru[k][1][i].z = ru[k][2][i].z;
          rE[k][0][i] = rE[k][2][i];
          rE[k][1][i] = rE[k][2][i];
        } // for
      } // for
    }
    else {
      const std::size_t j = m.size<mesh::x_axis, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[k][j - 1][i] = r[k][j - 3][i];
          r[k][j - 2][i] = r[k][j - 3][i];
          ru[k][j - 1][i].x = ru[k][j - 3][i].x;
          ru[k][j - 2][i].x = ru[k][j - 3][i].x;
          ru[k][j - 1][i].y = ru[k][j - 3][i].y;
          ru[k][j - 2][i].y = ru[k][j - 3][i].y;
          ru[k][j - 1][i].z = ru[k][j - 3][i].z;
          ru[k][j - 2][i].z = ru[k][j - 3][i].z;
          rE[k][j - 1][i] = rE[k][j - 3][i];
          rE[k][j - 2][i] = rE[k][j - 3][i];
        } // for
      } // for
    } // if
  }
  else {
    if(B == mesh::boundary::low) {
      for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[0][j][i] = r[2][j][i];
          r[1][j][i] = r[2][j][i];
          ru[0][j][i].x = ru[2][j][i].x;
          ru[1][j][i].x = ru[2][j][i].x;
          ru[0][j][i].y = ru[2][j][i].y;
          ru[1][j][i].y = ru[2][j][i].y;
          ru[0][j][i].z = ru[2][j][i].z;
          ru[1][j][i].z = ru[2][j][i].z;
          rE[0][j][i] = rE[2][j][i];
          rE[1][j][i] = rE[2][j][i];
        } // for
      } // for
    }
    else {
      const std::size_t k = m.size<mesh::z_axis, mesh::all>();
      for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[k - 1][j][i] = r[k - 3][j][i];
          r[k - 2][j][i] = r[k - 3][j][i];
          ru[k - 1][j][i].x = ru[k - 3][j][i].x;
          ru[k - 2][j][i].x = ru[k - 3][j][i].x;
          ru[k - 1][j][i].y = ru[k - 3][j][i].y;
          ru[k - 2][j][i].y = ru[k - 3][j][i].y;
          ru[k - 1][j][i].z = ru[k - 3][j][i].z;
          ru[k - 2][j][i].z = ru[k - 3][j][i].z;
          rE[k - 1][j][i] = rE[k - 3][j][i];
          rE[k - 2][j][i] = rE[k - 3][j][i];
        } // for
      } // for
    } // if
  } // if
} // flow

template<mesh::axis A, mesh::boundary B>
void
reflecting(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  field<vec3>::accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a) {
  auto r = m.mdspan<mesh::cells>(r_a);
  auto ru = m.mdspan<mesh::cells>(ru_a);
  auto rE = m.mdspan<mesh::cells>(rE_a);

  if constexpr(A == mesh::axis::x_axis) {
    if(B == mesh::boundary::low) {
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
          r[k][j][0] = r[k][j][3];
          r[k][j][1] = r[k][j][2];
          ru[k][j][0].x = ru[k][j][3].x;
          ru[k][j][1].x = ru[k][j][2].x;
          ru[k][j][0].y = ru[k][j][3].y;
          ru[k][j][1].y = ru[k][j][2].y;
          ru[k][j][0].z = ru[k][j][3].z;
          ru[k][j][1].z = ru[k][j][2].z;
          rE[k][j][0] = rE[k][j][3];
          rE[k][j][1] = rE[k][j][2];
        } // for
      } // for
    }
    else {
      const std::size_t i = m.size<mesh::x_axis, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
          r[k][j][i - 1] = r[k][j][i - 4];
          r[k][j][i - 2] = r[k][j][i - 3];
          ru[k][j][i - 1].x = ru[k][j][i - 4].x;
          ru[k][j][i - 2].x = ru[k][j][i - 3].x;
          ru[k][j][i - 1].y = ru[k][j][i - 4].y;
          ru[k][j][i - 2].y = ru[k][j][i - 3].y;
          ru[k][j][i - 1].z = ru[k][j][i - 4].z;
          ru[k][j][i - 2].z = ru[k][j][i - 3].z;
          rE[k][j][i - 1] = rE[k][j][i - 4];
          rE[k][j][i - 2] = rE[k][j][i - 3];
        } // for
      } // for
    } // if
  }
  else if constexpr(A == mesh::axis::y_axis) {
    if(B == mesh::boundary::low) {
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[k][0][i] = r[k][3][i];
          r[k][1][i] = r[k][2][i];
          ru[k][0][i].x = ru[k][3][i].x;
          ru[k][1][i].x = ru[k][2][i].x;
          ru[k][0][i].y = ru[k][3][i].y;
          ru[k][1][i].y = ru[k][2][i].y;
          ru[k][0][i].z = ru[k][3][i].z;
          ru[k][1][i].z = ru[k][2][i].z;
          rE[k][0][i] = rE[k][3][i];
          rE[k][1][i] = rE[k][2][i];
        } // for
      } // for
    }
    else {
      const std::size_t j = m.size<mesh::x_axis, mesh::all>();
      for(auto k : m.cells<mesh::z_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[k][j - 1][i] = r[k][j - 4][i];
          r[k][j - 2][i] = r[k][j - 3][i];
          ru[k][j - 1][i].x = ru[k][j - 4][i].x;
          ru[k][j - 2][i].x = ru[k][j - 3][i].x;
          ru[k][j - 1][i].y = ru[k][j - 4][i].y;
          ru[k][j - 2][i].y = ru[k][j - 3][i].y;
          ru[k][j - 1][i].z = ru[k][j - 4][i].z;
          ru[k][j - 2][i].z = ru[k][j - 3][i].z;
          rE[k][j - 1][i] = rE[k][j - 4][i];
          rE[k][j - 2][i] = rE[k][j - 3][i];
        } // for
      } // for
    } // if
  }
  else {
    if(B == mesh::boundary::low) {
      for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[0][j][i] = r[3][j][i];
          r[1][j][i] = r[2][j][i];
          ru[0][j][i].x = ru[3][j][i].x;
          ru[1][j][i].x = ru[2][j][i].x;
          ru[0][j][i].y = ru[3][j][i].y;
          ru[1][j][i].y = ru[2][j][i].y;
          ru[0][j][i].z = ru[3][j][i].z;
          ru[1][j][i].z = ru[2][j][i].z;
          rE[0][j][i] = rE[3][j][i];
          rE[1][j][i] = rE[2][j][i];
        } // for
      } // for
    }
    else {
      const std::size_t k = m.size<mesh::z_axis, mesh::all>();
      for(auto j : m.cells<mesh::y_axis, mesh::all>()) {
        for(auto i : m.cells<mesh::x_axis, mesh::all>()) {
          r[k - 1][j][i] = r[k - 4][j][i];
          r[k - 2][j][i] = r[k - 3][j][i];
          ru[k - 1][j][i].x = ru[k - 4][j][i].x;
          ru[k - 2][j][i].x = ru[k - 3][j][i].x;
          ru[k - 1][j][i].y = ru[k - 4][j][i].y;
          ru[k - 2][j][i].y = ru[k - 3][j][i].y;
          ru[k - 1][j][i].z = ru[k - 4][j][i].z;
          ru[k - 2][j][i].z = ru[k - 3][j][i].z;
          rE[k - 1][j][i] = rE[k - 4][j][i];
          rE[k - 2][j][i] = rE[k - 3][j][i];
        } // for
      } // for
    } // if
  } // if
} // reflecting

} // namespace muscl

#endif // MUSCL_BOUNDARY_HH
