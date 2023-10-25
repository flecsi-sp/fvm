from spack.package import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    git = "https://github.com/tuxfan/flecsi.git"
    version("2.3-fvm", commit="c919855bdaa2add3144b281541b99c1612432c75")

    # since we compile with -Wextra -Werror for development
    # and https://gitlab.kitware.com/cmake/cmake/-/issues/23141
    depends_on('cmake@3.19:3.21,3.22.3:')

    depends_on("legion@cr-16:cr-99", when="backend=legion")
