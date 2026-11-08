"""Generate selected DU145 F04. The former S04/S03c regulator panel was retired on 9 October 2026.

The existing --main-only invocation remains accepted for compatibility.
Regenerate S03 separately with make_du145_rna_gene_set_themes.py.
"""
from du145_programmes_matched_main import render_main

if __name__ == '__main__':
    render_main()
