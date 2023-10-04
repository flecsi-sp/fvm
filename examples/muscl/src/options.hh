#ifndef MUSCL_OPTIONS_HH
#define MUSCL_OPTIONS_HH

#include <flecsi/execution.hh>

namespace muscl::opt {
inline flecsi::program_option<std::size_t>
  x_extents("x-extents", "The x extents of the mesh.", 1);
inline flecsi::program_option<std::size_t>
  y_extents("y-extents", "The y extents of the mesh.", 1);
inline flecsi::program_option<std::size_t>
  z_extents("z-extents", "The z extents of the mesh.", 1);

inline flecsi::program_option<int> colors("MUSCL Options",
  "colors,c",
  "Specify the number of colors (default: num processes).",
  {{flecsi::option_default, -1}});
} // namespace muscl::opt

#endif // MUSCL_OPTIONS_HH
