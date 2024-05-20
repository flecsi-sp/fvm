from spack.package import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI.
    """
    version("2.3-beta", commit="45050e422eb85b0d6bfd81cb2c40f9b1d3f9ca85")
    patch("get_axis.patch", when="@2.3-beta")
