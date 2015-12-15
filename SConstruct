import os
import SCons.Conftest

exports = {}

import sys
print ("Python version is %s" % str (sys.version))

exports ["include_mpi"] = False
exports ["include_vt"] = False
exports ["include_mp"] = True
exports ["dimensions"] = 2
exports ["quiet"] = -1
exports ["vertical_grid"] = 'chebyshev'
exports ["logger"] = None
exports ["defaults"] = "src/defaults.yaml"

exports ["cppcompiler"] = None
exports ["ccompiler"] = None
exports ["compiler_flags"] = []
# exports ["compiler_flags"] = ["-mt_mpi"]
exports ["cxxtestgen_location"] = "/usr/local/Cellar/cxxtest/4.4/libexec/bin/cxxtestgen"
exports ["cxxtest"] = False
exports ["include_dirs"] = ["/usr/local/include", "/usr/local/opt/lapack/include"]
exports ["library_dirs"] = ["/usr/local/lib", "/usr/local/opt/lapack/lib"]
exports ["env_vars"] = {}

exports ["root"] = os.path.abspath ("..") + "/"
exports ["defaults"] = os.path.join (exports ["root"], exports ["defaults"])

exports ["env"] = Environment (tools = ["default"], **os.environ)
exports ["env"].Tool ("cxxtest", CXXTEST = exports ["cxxtestgen_location"])

exports ["env"].SConscript ("src/SConscript", variant_dir = "build/", duplicate = 0, exports = exports)
