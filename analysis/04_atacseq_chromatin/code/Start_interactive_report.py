"""Start the localhost-only report server without a visible background console."""
from pathlib import Path
import json
import os
import subprocess
import sys
import urllib.request
import time
import webbrowser

ROOT=Path(__file__).resolve().parents[1]
URL='http://127.0.0.1:8765/report/interactive_ATAC_RNA.html'
def open_browser():
    print(f'Opening your browser: {URL}', flush=True)
    try:
        if os.name == 'nt':
            os.startfile(URL)
        elif not webbrowser.open(URL, new=2):
            raise RuntimeError('No default browser was found.')
    except Exception as error:
        print(f'Could not open the browser: {error}\nOpen this address manually: {URL}', file=sys.stderr)
        sys.exit(1)

try:
    with urllib.request.urlopen('http://127.0.0.1:8765/api/health',timeout=2) as r:
        if json.load(r).get('application')=='HT1080_ATAC_RNA':
            print(URL)
            if '--open' in sys.argv:open_browser()
            sys.exit(0)
except (OSError,ValueError):
    pass
with (ROOT/'logs/interactive_viewer_stdout.log').open('w') as stdout,(ROOT/'logs/interactive_viewer_stderr.log').open('w') as stderr:
    process=subprocess.Popen([sys.executable,'-u',str(ROOT/'code/51_interactive_server.py')],cwd=ROOT,
        stdout=stdout,stderr=stderr,stdin=subprocess.DEVNULL,
        creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0,env=dict(os.environ))
(ROOT/'integration/interactive_atac_rna/server.pid').write_text(str(process.pid))
print(f'Started local viewer (PID {process.pid}): {URL}')
if '--open' in sys.argv:
    for _ in range(40):
        try:
            with urllib.request.urlopen('http://127.0.0.1:8765/api/health',timeout=1) as r:
                status = json.load(r)
                if status.get('status')=='ready' and status.get('application')=='HT1080_ATAC_RNA':
                    open_browser()
                    break
        except (OSError,ValueError):time.sleep(.5)
        if process.poll() is not None:
            print('The report server stopped during startup. Error log:', file=sys.stderr)
            print((ROOT/'logs/interactive_viewer_stderr.log').read_text(encoding='utf-8', errors='replace')[-4000:], file=sys.stderr)
            sys.exit(1)
    else:
        print(f'The server did not become ready. See {ROOT / "logs/interactive_viewer_stderr.log"}', file=sys.stderr)
        sys.exit(1)
