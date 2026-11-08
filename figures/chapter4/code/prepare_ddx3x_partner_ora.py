"""ORA of the 31 DDX3X partners against the eligible partner universe.

Local frozen human MSigDB 2026.1: GO:BP, GO:CC and Reactome.
One BH family across all size-eligible sets, including zero-overlap sets.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact, false_discovery_control
from prepare_xy_ranked_partners_by_y_status import sha256

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS/'figures/chapter4'
DATA = OUT/'source_data/ddx3x_partner_ora'
BASE = OUT/'source_data/ddx3x_genomewide_y_pattern'
RES = THESIS.parent/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/resources'
COLLECTIONS = {'GO:BP':'GO biological process','GO:CC':'GO cellular component','CP:REACTOME':'Reactome'}
MIN_SIZE,MAX_SIZE = 10,500


def bh(p):
    """Independent simple BH implementation, checked against scipy."""
    p = np.asarray(p,float)
    order = np.argsort(p,kind='stable')
    ranked = p[order]*len(p)/np.arange(1,len(p)+1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(p));out[order] = np.minimum(ranked,1)
    np.testing.assert_allclose(out,false_discovery_control(p,method='bh'),rtol=1e-12,atol=1e-15)
    return out


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    scan = pd.read_csv(BASE/'all_gene_correlations.csv')
    picked = pd.read_csv(BASE/'selected_partners.csv')
    background = set(scan.loc[scan.eligible,'gene'])
    query = set(picked.gene)
    assert len(background)==17786 and len(query)==31 and query<=background
    assert query==set(scan.loc[scan.selected,'gene']) and 'DDX3X' not in background
    (DATA/'query_31.txt').write_text('\n'.join(sorted(query))+'\n')
    (DATA/'background_17786.txt').write_text('\n'.join(sorted(background))+'\n')
    metadata = pd.read_csv(RES/'msigdb_all_set_metadata.csv').fillna('')
    metadata = metadata.loc[metadata.gs_subcollection.isin(COLLECTIONS)].copy()
    assert metadata.gs_name.is_unique and set(metadata.db_version)=={'2026.1.Hs'}
    membership = pd.read_csv(RES/'msigdb_all_membership.csv.gz')
    membership = membership[membership.pathway.isin(set(metadata.gs_name))].drop_duplicates(['pathway','symbol'])
    assert membership.symbol.notna().all()
    rawsets = membership.groupby('pathway').symbol.agg(set).to_dict()
    annotated = {}
    for collection in COLLECTIONS:
        names = set(metadata.loc[metadata.gs_subcollection.eq(collection),'gs_name'])
        annotated[collection] = set().union(*(rawsets[n] for n in names))
    annotation_rows = []
    for gene in sorted(background):
        row = {'gene':gene,'in_query':gene in query}
        row.update({c:gene in annotated[c] for c in COLLECTIONS})
        annotation_rows.append(row)
    pd.DataFrame(annotation_rows).to_csv(DATA/'annotation_coverage.csv',index=False)
    rows,excluded = [],[]
    for term in metadata.itertuples():
        genes = rawsets[term.gs_name]&background
        k = len(genes&query)
        K,N,n = len(genes),len(background),len(query)
        if not MIN_SIZE<=K<=MAX_SIZE:
            excluded.append({'pathway':term.gs_name,'collection':term.gs_subcollection,'background_genes':K,
                             'reason':'Below minimum size' if K<MIN_SIZE else 'Above maximum size'})
            continue
        p = float(hypergeom.sf(k-1,N,K,n))
        annotated_background = background&annotated[term.gs_subcollection]
        annotated_query = query&annotated_background
        p_annotated = float(hypergeom.sf(k-1,len(annotated_background),K,len(annotated_query)))
        rows.append(dict(pathway=term.gs_name,collection=term.gs_subcollection,source=COLLECTIONS[term.gs_subcollection],
                         msigdb_id=term.gs_id,source_url=term.gs_url,db_version=term.db_version,
                         background_n=N,query_n=n,term_background_n=K,overlap_n=k,
                         expected_overlap=n*K/N,fold_enrichment=(k/n)/(K/N),p_value=p,
                         overlap_genes=';'.join(sorted(query&genes)),
                         term_background_genes=';'.join(sorted(genes)),
                         annotated_background_n=len(annotated_background),annotated_query_n=len(annotated_query),
                         p_annotated_background=p_annotated))
    result = pd.DataFrame(rows)
    result['fdr_bh_global'] = bh(result.p_value)
    result['fdr_bh_annotated_background'] = bh(result.p_annotated_background)
    # BY is reported as a dependence-conservative robustness statistic, not the primary criterion.
    harmonic = float(np.sum(1/np.arange(1,len(result)+1)))
    result['fdr_by_global'] = np.minimum(result.fdr_bh_global*harmonic,1)
    result = result.sort_values(['fdr_bh_global','p_value','pathway']).reset_index(drop=True)
    result.to_csv(DATA/'all_tested_terms.csv',index=False)
    pd.DataFrame(excluded).to_csv(DATA/'excluded_term_sizes.csv',index=False)
    significant = result.loc[result.fdr_bh_global.lt(.05)].copy()
    significant.to_csv(DATA/'significant_terms_fdr005.csv',index=False)
    # Selection only affects the summary display, never the test family or FDRs.
    # Rank all significant terms; keep a representative unless >=75% of its
    # query overlap is redundant with a previously retained overlap (Jaccard).
    representatives,assignment = [],[]
    for row in significant.itertuples():
        genes = set(row.overlap_genes.split(';'))
        redundant = next((r for r in representatives
                          if len(genes&r['genes'])/len(genes|r['genes'])>=.75),None)
        if redundant is None:
            representatives.append({'pathway':row.pathway,'genes':genes})
            rep = row.pathway
        else:
            rep = redundant['pathway']
        assignment.append({'pathway':row.pathway,'representative':rep,'retained_representative':redundant is None})
    assignment = pd.DataFrame(assignment)
    assignment.to_csv(DATA/'term_display_redundancy.csv',index=False)
    names = [r['pathway'] for r in representatives][:12]
    display = significant.set_index('pathway').loc[names].reset_index()
    display.insert(0,'display_rank',np.arange(1,len(display)+1))
    display.to_csv(DATA/'display_terms.csv',index=False)
    checks = []
    for row in significant.itertuples():
        table = [[row.overlap_n,row.query_n-row.overlap_n],
                 [row.term_background_n-row.overlap_n,
                  row.background_n-row.query_n-row.term_background_n+row.overlap_n]]
        assert min(sum(table,[]))>=0
        independent = fisher_exact(table,alternative='greater').pvalue
        checks.append(abs(independent-row.p_value))
    assert max(checks,default=0)<1e-12
    sources = [BASE/'all_gene_correlations.csv',BASE/'selected_partners.csv',BASE/'manifest.json',
               RES/'msigdb_all_membership.csv.gz',RES/'msigdb_all_set_metadata.csv',
               RES/'msigdb_all_collection_catalogue.csv',Path(__file__)]
    manifest = dict(
        query_genes=sorted(query),query_size=len(query),background_size=len(background),
        background_definition='All 17,786 non-self genes eligible for the preceding DDX3X co-dependency pattern scan; no genome-wide default universe.',
        database='Local frozen human MSigDB 2026.1.Hs from the Chapter 3 HT1080 analysis resources',
        collections=COLLECTIONS,
        mapping='Exact human gene symbols, without uppercase conversion, aliases or expansion. Coverage is exported; input symbols and gene sets are deduplicated.',
        coverage={c:{'query_annotated':len(query&g),'background_annotated':len(background&g),
                     'query_without_annotation':sorted(query-g)} for c,g in annotated.items()},
        gene_set_size_rule={'minimum':MIN_SIZE,'maximum':MAX_SIZE,'defined_after_intersection_with_background':True},
        statistical_test='One-sided hypergeometric upper tail P(X >= observed overlap), equivalent to one-sided Fisher exact over-representation.',
        multiple_testing='Benjamini-Hochberg across ALL size-eligible terms in the three collections jointly, including zero-overlap terms. Primary threshold FDR <0.05.',
        tested_terms=len(result),tested_by_collection=result.groupby('collection').size().to_dict(),
        significant_terms=len(significant),significant_by_collection=significant.groupby('collection').size().to_dict(),
        significant_with_annotated_background=int(result.fdr_bh_annotated_background.lt(.05).sum()),
        significant_with_by=int(result.fdr_by_global.lt(.05).sum()),
        sensitivities='Also report separately recalculated global BH after restricting each collection background to its annotated eligible genes, and global BY under arbitrary dependence. Neither defines primary calls.',
        display_selection='Up to 12 terms, ordered by global BH then raw P, with greedy overlap-gene Jaccard >=0.75 redundancy removal across collections. All tested terms and assignments retained.',
        display_terms=len(display),
        limitations='Exploratory functional annotation of a correlation-selected gene set. ORA assumes exchangeability of selected genes under its null and does not model co-dependency among genes, unequal annotation, essentiality or selection bias. Enrichment does not establish altered pathway activity or LOY-specific dependency.',
        validation={'query_is_exact_previous_31':True,'query_within_background':True,'BH_matches_scipy':True,
                    'Fisher_checks_on_all_significant_terms':len(checks),'max_Fisher_p_error':max(checks,default=0)},
        sources={str(p):sha256(p) for p in sources})
    (DATA/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({k:manifest[k] for k in ['coverage','tested_terms','tested_by_collection','significant_terms','significant_by_collection','significant_with_annotated_background','significant_with_by','validation']},indent=2))
    print(display[['source','pathway','overlap_n','term_background_n','fold_enrichment','fdr_bh_global','overlap_genes']].to_string(index=False))


if __name__=='__main__':
    main()
