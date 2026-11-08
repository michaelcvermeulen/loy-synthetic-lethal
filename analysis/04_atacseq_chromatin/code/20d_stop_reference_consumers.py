import os,signal,subprocess
root='/mnt/a/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS/code/'
lines=subprocess.check_output(['ps','-eo','pid=,ppid=,args='],text=True).splitlines()
procs={}
for line in lines:
 p,p_parent,args=line.strip().split(None,2);procs[int(p)]=(int(p_parent),args)
targets={p for p,(_,args) in procs.items() if root+'25_run_homer.py' in args or root+'27_continue_tobias.py' in args}
targets.discard(os.getpid())
while True:
 children={p for p,(parent,_) in procs.items() if parent in targets}
 if children.issubset(targets):break
 targets|=children
for p in sorted(targets,reverse=True):
 print('Stopping owned reference consumer',p,procs[p][1],flush=True)
 try:os.kill(p,signal.SIGTERM)
 except ProcessLookupError:pass
