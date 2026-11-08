"""Deterministic display selection from saved DU145 leading edges and expression effects."""
from collections import Counter
import numpy as np
import pandas as pd
from du145_rna_common import read

BIO_THEMES = ['interferon_inflammatory', 'GRHL2_epithelial', 'mitotic',
              'ribosomal/translation', 'histone cluster (6p22)']


def select_leading_edge_genes(themes, per_theme=6):
    # Use the exact frozen gene-ID/symbol mapping used for the GSEA ranking.
    ranking = read('pathway_TF_effect_ranking')
    assert ranking.gene_name.is_unique and ranking.gene_id.is_unique
    effects = read('DU145_LOY_vs_sgCTRL_primary_effects')
    de = read('DE_paired_voom_primary')
    cols = ['gene_id', 'gene_name', 'chrom', 'CTRL_mean_CPM', 'LOY_mean_CPM',
            'mean_difference', 'pair_P12', 'pair_P13', 'pair_P14']
    mapped = ranking[['gene_id','gene_name']].merge(effects[cols], on=['gene_id','gene_name'],
                                                  validate='one_to_one')
    assert len(mapped) == len(ranking)
    mapped = mapped.merge(de[['gene_id','logFC','FDR']], on='gene_id', validate='one_to_one')
    assert np.allclose(mapped[['pair_P12','pair_P13','pair_P14']].mean(axis=1), mapped.mean_difference)
    candidates, selected, memberships, summaries = [], [], [], []
    for theme in BIO_THEMES:
        sets = themes[themes.theme.eq(theme) & themes.status.eq('retained')]
        assert len(sets) > 0
        counts = Counter()
        for r in sets.itertuples():
            genes = sorted(set(r.leadingEdge.split(';')))
            counts.update(genes)
            memberships.extend(dict(theme=theme, pathway=r.pathway, gene_name=g) for g in genes)
        freq = pd.DataFrame(counts.items(), columns=['gene_name','leading_edge_sets'])
        freq = freq.merge(mapped, on='gene_name', validate='one_to_one')
        assert len(freq) == len(counts), 'Do not silently drop unmapped leading-edge genes.'
        freq['theme'] = theme
        freq['retained_sets'] = len(sets)
        freq['leading_edge_fraction'] = freq.leading_edge_sets / len(sets)
        freq['abs_mean_difference'] = freq.mean_difference.abs()
        freq = freq.sort_values(['leading_edge_sets','abs_mean_difference','gene_name'],
                                ascending=[False,False,True], kind='stable').reset_index(drop=True)
        freq['frequency_rank'] = np.arange(1, len(freq)+1)
        freq['selected_for_display'] = freq.frequency_rank <= per_theme
        freq['gene_FDR_below_005'] = freq.FDR < .05
        candidates.append(freq)
        chosen = freq.head(per_theme).copy()
        assert len(chosen) == per_theme
        selected.append(chosen)
        summaries.append(dict(theme=theme, retained_sets=len(sets), distinct_leading_edge_genes=len(freq),
                              displayed_genes=per_theme, gene_FDR_significant=int(chosen.gene_FDR_below_005.sum())))
    return (pd.concat(selected, ignore_index=True), pd.concat(candidates, ignore_index=True),
            pd.DataFrame(memberships), pd.DataFrame(summaries))
