#ifndef TEST_OPTIONS_HH
#define TEST_OPTIONS_HH

#include <flecsi/execution.hh>

namespace test::opt {

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

} // namespace test::opt

#endif // TEST_OPTIONS_HH
