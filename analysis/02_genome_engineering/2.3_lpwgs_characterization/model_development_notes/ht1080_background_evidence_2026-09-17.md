# HT1080 background and TP53 annotation — 17 September 2026

The author requested TP53 WT in the profile and a short Results rationale covering HT1080's broad research use, pseudodiploidy and chromosome-engineering history. Both profile variants and their metadata now state **TP53 WT (reported in ATCC-derived cells)** and **pseudodiploid; modal chromosome count 46**. The original versions were archived before revision. The shared plotting heading is “Reported genotype”, because wild type is not an alteration.

## Evidence used

| Claim | Source and scope |
|---|---|
| Established in 1972 from an untreated 35-year-old man; first publication in 1974 | [DSMZ ACC 315](https://www.dsmz.de/collection/catalogue/details/culture/ACC-315) and [Rasheed et al., Cancer 33, 1027–1033 (1974)](https://pubmed.ncbi.nlm.nih.gov/4132053/). Publication record checked; the original full text was not retrieved. |
| Pseudodiploidy, modal count 46, transfection suitability | [ATCC CCL-121](https://www.atcc.org/products/ccl-121). ATCC also describes chromosome abnormalities and a 44–48 range in initial seed stock. |
| Broad research use | [Kim et al., PLOS ONE (2011)](https://doi.org/10.1371/journal.pone.0018556) identifies HT1080 among widely used human cell lines and uses it for genetic-engineering experiments. |
| Artificial-chromosome engineering | [Nakano et al., Developmental Cell (2008)](https://doi.org/10.1016/j.devcel.2008.02.001) constructs the alphoid-tetO HAC in HT1080 and manipulates its kinetochore. |
| Microcell-mediated chromosome transfer | [Suzuki et al., Journal of Biological Chemistry (2006)](https://doi.org/10.1074/jbc.M603053200) constructs HACs in HT1080 and transfers them to mouse embryonic stem cells. |
| TP53 WT and nutlin response in ATCC-derived HT1080 | [Moran and Maki, Molecular Cancer Therapeutics 9, 895–905 (2010)](https://doi.org/10.1158/1535-7163.MCT-09-1220). Methods identify ATCC as the cell source; the paper identifies HT1080 as p53 wild type and reports nutlin-induced cell-cycle arrest. The PMC availability date is 2011, but the publication year is 2010. |
| TP53 status differs between reported HT1080 cultures | [1994 chromosome-transfer study, PMID 7519049](https://pubmed.ncbi.nlm.nih.gov/7519049/) reports mutations in both TP53 copies in its HT1080 culture. This prevents an unqualified statement that every HT1080 stock is TP53 WT. |

## Interpretation and drafting choices

The TP53 annotation is a published background property, not new sequencing or functional validation of the author's stocks. Pseudodiploidy does not establish a normal or uniformly stable genome: the repository describes abnormalities, and the new public WGS profile shows copy-number changes. Therefore the prose does not claim that baseline aneuploidy or copy-number confounding is absent. The existing ~30-h DSMZ doubling-time annotation is retained; the unverified 20–24-h range in the pasted material was not substituted.

The Results now introduces HT1080's history and chromosome-research uses before the existing nucleofection paragraph. The broad statement that HT1080 met every paralog-dependency criterion was removed from that opening: current public scores support gene-specific interpretation, and the later 26Q1 profile does not by itself reconstruct historical selection decisions. The experimental guide-testing chronology and numerical RT-qPCR/ddPCR results are preserved.

Only the requested background claims were added from the supplied material; no new claim of a fully intact Y chromosome, expression of every MSY gene, or genotype verification in the author's parental stock was introduced. Active captions remain in `Chapter3_figure_captions.md`, and both metadata files carry claim-specific source URLs.
