"""Compatibility entry point for the adopted three-panel A549 / DU145 atlases."""
from pathlib import Path
import subprocess,sys
if len(sys.argv)<2 or sys.argv[1] not in ['A549','DU145']:
    raise SystemExit('Usage: python make_cell_line_lpwgs_atlases.py A549|DU145 [--output-dir DIR]')
subprocess.run([sys.executable,str(Path(__file__).with_name('make_regional_y_atlas.py')),*sys.argv[1:]],check=True)
