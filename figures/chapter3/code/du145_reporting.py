"""Author-approved DU145 reporting subset, 26 September 2026; frozen model unchanged."""
def annotate(de):
    de=de.copy()
    de['reporting_subset']=(de.FDR<.05)&(de.logFC.abs()>=.5)&(de.pairs_agreeing_with_logFC==3)
    return de
LABEL='FDR <0.05, |log2FC| \u22650.5,\nconsistent direction'
