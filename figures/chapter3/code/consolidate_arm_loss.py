"""Current chromosome-loss figure entry point.

The current four-panel supplementary figure uses the compact 30 October layout.
Previous exports are preserved in archive/arm_loss_whitespace_before_20261030_153759.
"""
from pathlib import Path
import runpy
if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).with_name('make_arm_loss_compact.py')), run_name='__main__')
