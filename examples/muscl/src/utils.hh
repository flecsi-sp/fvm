#ifndef MUSCL_UTILS_HH
#define MUSCL_UTILS_HH

#include <fvm/mesh.hh>

namespace muscl::utils {

template<typename T>
T
sqr(T t) {
  return t * t;
} // sq

inline fvm::mesh::boundary_type
mesh_boundary(std::string const & b) {
  if(b == "flow") {
    return fvm::mesh::boundary_type::flow;
  }
  if(b == "reflecting") {
    return fvm::mesh::boundary_type::reflecting;
  }
  if(b == "periodic") {
    return fvm::mesh::boundary_type::periodic;
  }
  else {
    flog_fatal("invalid boundary type(" << b << ")");
  } // if
} // mesh_boundary

} // namespace muscl::utils

#endif // MUSCL_UTILS_HH
