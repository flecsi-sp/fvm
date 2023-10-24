#ifndef TEST_OPTIONS_HH
#define TEST_OPTIONS_HH

#include <flecsi/execution.hh>

namespace test::opt {
inline flecsi::program_option<std::string> config("yaml file",
  "The yaml config file.",
  1,
  [](flecsi::any const & v, std::stringstream & ss) {
    const std::string value = flecsi::option_value<std::string>(v);
    return value.find(".yaml") != std::string::npos
             ? true
             : (ss << "file(" << value << ") has invalid suffix") && false;
  });

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
