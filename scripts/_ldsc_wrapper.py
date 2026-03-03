# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "numpy==1.21.5",
#   "pandas==1.3.3",
#   "scipy==1.7.3",
#   "bitarray==2",
#   "python-dateutil==2.8.2",
#   "pytz==2022.4",
#   "six==1.16.0",
# ]
# ///
import os
import runpy
import sys

script = os.environ["LDSC_SCRIPT"]
project_root = os.path.dirname(script)
if project_root not in sys.path:
    sys.path.insert(0, project_root)

sys.argv = [script] + sys.argv[1:]
runpy.run_path(script, run_name="__main__")
