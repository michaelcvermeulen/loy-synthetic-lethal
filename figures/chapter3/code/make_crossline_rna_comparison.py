"""Current cross-line figure entry point.

The former five-panel implementation is archived at
archive/crossline_condense_20261008_175215/figures/chapter3/code/make_crossline_rna_comparison.py.
"""
from pathlib import Path
import runpy

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).with_name('make_crossline_concise.py')), run_name='__main__')
