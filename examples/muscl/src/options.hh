#ifndef MUSCL_OPTIONS_HH
#define MUSCL_OPTIONS_HH

#include <flecsi/execution.hh>

namespace muscl::opt {

inline flecsi::program_option<std::string> config("yaml file",
  "The yaml config file.",
  1,
  [](std::string const & v, std::stringstream & ss) {
    return v.find(".yaml") != std::string::npos
             ? true
             : (ss << "file(" << v << ") has invalid suffix") && false;
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

inline flecsi::program_option<std::string> flog_tags("FLOG Options",
  "tags,t",
  "Specify the flog tags to enable.",
  {{flecsi::option_default, "all"}});

inline flecsi::program_option<int> flog_verbose("FLOG Options",
  "verbose,v",
  "Enable verbose output. Passing '-1' will strip any additional"
  " decorations added by flog and will only output the user's message.",
  {{flecsi::option_default, 0}});

} // namespace muscl::opt

#endif // MUSCL_OPTIONS_HH
