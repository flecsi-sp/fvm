from spack.package import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI.
    """
    version("2.3-beta", commit="4c22b865904686647205b28518dd649005358d4d")
    patch("get_axis.patch", when="@2.3-beta")
