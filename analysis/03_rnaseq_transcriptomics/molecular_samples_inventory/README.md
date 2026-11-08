# Molecular sample map

[rna_library_map.csv](rna_library_map.csv) lists the 44 unique RNA-seq libraries used in Chapter 3, with original and corrected identifiers, condition, collection information, retained count path and membership in baseline/parental and add-back analyses.

The total is 15 HT1080 baseline/parental + 11 A549 + 6 DU145 + 18 add-back − 6 shared May HT1080 libraries = 44. A library is a measurement, not necessarily an independently derived population. Historical LNCaP identifiers are retained alongside the author-confirmed DU145 identity. Relative HT1080 collection indices are not absolute passage numbers; A549 passage matching and independent derivation remain unestablished.

The compact table in Chapter3_methods.md also describes baseline lpWGS and HT1080 ATAC sampling. Individual specimens across DNA, RNA and ATAC are not assumed to match. The RNA crosswalk does not provide a complete ddPCR, drug-optimisation or DNA-sample inventory.

Reproduce with [build_sample_map.py](../../../research/methods_completion_2026-10-09/build_sample_map.py). [provenance.json](provenance.json) records the four source manifests and upstream count headers. This operation only assembles metadata; it does not modify counts or run statistical analyses.
