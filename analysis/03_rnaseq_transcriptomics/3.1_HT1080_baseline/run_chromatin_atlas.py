"""Rebuild the measured reference atlas from frozen local inputs."""
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
R='C:/Program Files/R/R-4.5.2/bin/Rscript.exe'
LINUX_PY='/home/micha/.local/share/atacseq-sept08-2026-venv/bin/python'


def run(args):
    print('Running:', ' '.join(map(str,args)),flush=True)
    subprocess.run(list(map(str,args)),cwd=ROOT,check=True)


def main():
    py=[sys.executable,'-X','utf8','-u']
    run(py+[ROOT/'prepare_chromatin_atlas.py'])
    linux=str(ROOT/'summarize_chromatin_signal.py').replace('\\','/')
    linux='/mnt/'+linux[0].lower()+linux[2:]
    wsl=['wsl.exe','-d','Ubuntu-24.04','--',LINUX_PY,linux]
    run(wsl+['--metadata-only'])
    run([R,'--vanilla',ROOT/'build_chromatin_peak_features.R'])
    run(wsl)
    run(py+[ROOT/'finalize_chromatin_atlas.py'])
    run(py+[ROOT/'validate_chromatin_atlas.py'])
    run(wsl+['--validate-only'])
    run(py+[ROOT/'plot_chromatin_atlas.py'])
    run(py+[ROOT/'complete_chromatin_atlas.py'])


if __name__=='__main__':main()
