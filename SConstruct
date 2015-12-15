import os
import sys
import SCons.Conftest
import yaml

exports = {}

with open("build_defaults.yaml", 'r') as stream:
    exports = yaml.load(stream)

try:
	with open("build_custom.yaml", 'r') as stream:
		new_exports = yaml.load(stream)
		for key in new_exports:
		    exports [key] = new_exports [key]
except IOError:
	pass

exports ["root"] = os.path.abspath (".") + "/"
exports ["defaults"] = os.path.join (exports ["root"], exports ["defaults"])

exports ["env"] = Environment (tools = ["default"], **os.environ)

exports ["cxxtest"] = False
if "check" in sys.argv:
	exports ["cxxtest"] = True
	exports ["cxxtestgen_location"] = "/usr/local/Cellar/cxxtest/4.4/libexec/bin/cxxtestgen"
	exports ["env"].Tool ("cxxtest", CXXTEST = exports ["cxxtestgen_location"])

exports ["env"].SConscript ("src/SConscript", variant_dir = "build/", duplicate = 0, exports = exports)
