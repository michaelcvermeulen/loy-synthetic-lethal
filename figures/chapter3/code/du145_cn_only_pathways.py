"""Presentation of the saved CN-only GSEA sensitivity; frozen source flags preserved."""
import pandas as pd
import numpy as np
from du145_rna_common import T, THEMES, THEME_LABEL as ORIGINAL_THEME_LABEL
STATUS_ORDER = ['retained', 'sensitive']
STATUS_COLOUR = {'retained': '#E64B35', 'sensitive': '#F0B27A'}
STATUS_LABEL = {'retained': 'Retained after CN exclusion', 'sensitive': 'Not retained after CN exclusion'}
THEME_LABEL = dict(ORIGINAL_THEME_LABEL, positional='Chromosomal sets', GRHL2_epithelial='Epithelial')
def cn_only_themes():
    th = pd.read_csv(THEMES).rename(columns={'status': 'original_screen_status'})
    fg = pd.read_csv(T / 'fgsea_all_collections.csv')
    th = th.merge(fg[['pathway','survives_CN_screen']], on='pathway', validate='one_to_one')
    th['status'] = np.where(th.survives_CN_screen, 'retained', 'sensitive')
    assert len(th)==154 and th.survives_CN_screen.sum()==139
    return th
