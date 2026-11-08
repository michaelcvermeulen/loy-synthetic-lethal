"""Fix the three source notebooks; preserve originals beside this rerun."""
from pathlib import Path
import re
root = Path(__file__).resolve().parents[4]
out = Path(__file__).resolve().parent
for cohort in ('HT1080', 'A549', 'LNCaP'):
    path = root / 'WGS/code' / f'coverage_report_{cohort}_QDNAseq.Rmd'
    old = path.read_bytes()
    backup = out / 'inputs' / path.name
    if not backup.exists():
        backup.write_bytes(old)
    text = old.decode('utf-8-sig')
    pattern = r"# Disable all (?:of QDNAseq's|QDNAseq) per-bin filters.*?chromosomes = NA\)"
    replacement = '''# Preserve the original T2T eligibility mask (bases >= 60%,
# mappability >= 55%, blacklist < 20%). An all-disabled applyFilters()
# call resets use flags rather than preserving them. Keep X and Y.
source("preserve_t2t_bin_mask.R")
set.seed(20261027) # reproducible CBS segmentation and subsequent normalization
rc_f <- preserve_t2t_bin_mask(rc, bins)
stopifnot(sum(fData(rc_f)$use) == 22682,
          sum(fData(rc_f)$use & fData(rc_f)$chromosome == "Y") == 58)'''
    text, n = re.subn(pattern, lambda _: replacement, text, count=1, flags=re.S)
    if n != 1:
        if 'rc_f <- preserve_t2t_bin_mask(rc, bins)' in text:
            print('Already patched:', path.name)
            continue
        raise RuntimeError(f'Expected one filtering block in {path}')
    text = text.replace('dependson="run_qdnaseq"}',
                        'dependson="run_qdnaseq", cache.extra=tools::md5sum("preserve_t2t_bin_mask.R")}')
    text = text.replace('cn <- if (calls_ok) cn_called else cn_seg',
                        'cn <- if (calls_ok) cn_called else cn_seg\nstopifnot(identical(fData(cn)$use, fData(rc_f)$use))')
    path.write_text(text, encoding='utf-8', newline='')
    print('Patched:', path.name)
