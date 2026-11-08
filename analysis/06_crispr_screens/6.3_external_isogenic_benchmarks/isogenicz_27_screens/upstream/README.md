# IsogenicZ
IsogenicZ is an adapted version of the DrugZ software from the Hart Lab. For paired analysis such as a  
CRISPR drug screen, DrugZ analyses two endpoint samples per replicate as input. IsogenicZ uses two  
endpoint samples (wildtype vs. isogenic mutant) and the two corresponding t0 samples per replicate as  
input. All samples are normalized to an equal total read number, and per replicate (endpoint/t0) fold  
changes for each sgRNA are calculated for wildtype and isogenic mutant, multiplied by the combined  
wildtype and isogenic mutant t0 count, and normalized to equal read numbers again to yield t0-corrected  
endpoint samples. These samples are then analyzed using the standard DrugZ analysis. 

This normalization method was used for analysis of pooled full genome CRISPR screens in multiple  
research papers:

[ELOF1 is a transcription-coupled DNA repair factor that directs RNA polymerase II ubiquitylation](https://pubmed.ncbi.nlm.nih.gov/34108663/)

[MMS22L-TONSL functions in sister chromatid cohesion in a pathway parallel to DSCC1-RFC](https://pubmed.ncbi.nlm.nih.gov/36622344/)

[CRISPR screens in sister chromatid cohesion defective cells reveal PAXIP1-PAGR1 as regulator of chromatin association of cohesin](https://pubmed.ncbi.nlm.nih.gov/37702151/)

Below the DrugZ readme.md edited where appropriate: 
  
IsogenicZ analyzes parallel isogenic CRISPR screens for synthetic lethal interactions.  

```
usage: isogenicz.py [-h] [-i sgRNA_count.txt] [-o isogenicz-output.txt]  
                [-f drugz-foldchange.txt] -c wildtype samples t0 -d wildtype samples endpoint
                -x mutant samples t0 -y mutant samples endpoint
                [-r remove genes] [-p pseudocount] [-I INDEX_COLUMN]  
                [--minobs minObs] [--half_window_size half_window_size] [-q]  
  
-i      	Readcount file, tab-delimited text (input)  
-o      	IsogenicZ results file, tab-delimited text (output)  
-f      	DrugZ Z-transformed fold change file (optional)  
-c      	Wildtype samples t0: comma-delimited list of column headers in readcount file
-d      	Wildtype samples endpoint: comma-delimited list of column headers in readcount file
-x      	Mutant samples t0: comma-delimited list of column headers in readcount file
-y      	Mutant samples endpoint: comma-delimited list of column headers in readcount file
-r      	Comma-delimited list of genes to remove before analysis  
-p      	Pseudocount to add to all readcounts; prevents log(0) problems (default=5) 
-I      	Index column (default=0)  
--minobs   	Ignore genes with fewer observations ( gRNA/gene x replicates) (default=1) 
--half_window_size  Size of the first bin and half the size of the inital sample
    (window) to estimate std (default=500) 
-unpaired Unpaired approach: compares mean(treated samples) to mean(control samples) (default=False)
```
  
The input file should be a tab-delimited file with the following format (order of sample columns or unused columns have no effect):

```
sgRNA	Gene	WTt0_1	WTt0_2	WTt0_3	WTt12_1	WTt12_2	WTt12_3	MUTt0_1	MUTt0_2	MUTt0_3	MUTt12_1	MUTt12_2	MUTt12_3
ACTGGCGCCATCGAGAGCCA	A1BG	150	148	161	130	95	247	154	184	128	124	213	169
CAAGAGAAAGACCACGAGCA	A1BG	352	367	344	355	320	445	520	539	484	679	627	532
GCTCAGCTGGGTCCATCCTG	A1BG	673	839	826	624	885	724	857	841	704	711	982	854
GTCGAGCTGATTCTGAGCGA	A1BG	305	264	360	444	193	379	343	361	295	319	394	288
AGTTATGTTAGGTATACCCG	A1CF	380	326	371	510	257	312	351	363	343	458	467	332
ATGACTCTCATACTCCACGA	A1CF	371	457	445	444	472	508	604	590	547	642	667	545
GGTGCAGCATCCCAACCAGG	A1CF	407	468	401	618	449	306	464	545	492	508	494	501
TGCGCTGGACCAGTGCGCGG	A1CF	213	186	251	272	159	285	167	206	138	190	206	253
AACAGCCACACTCACAGCGA	A2M	475	664	642	547	552	501	632	644	622	758	701	529
(etc)
```
Usage from command line, assuming python 3 with numpy and pandas is installed:
```
python isogenicz.py -i elof1_in.txt -o elof1_out.txt -c WTt0_1,WTt0_2,WTt0_3 -d WTt12_1,WTt12_2,WTt12_3 -x MUTt0_1,MUTt0_2,MUTt0_3 -y MUTt12_1,MUTt12_2,MUTt12_3
```
File elof1_in.txt contains the t0 and t12 counts from the screens of the RPE-1 iCas9 cell line ("WT") and the derived ELOF1 knockout line ("MUT"), see [ELOF1 publication](https://pubmed.ncbi.nlm.nih.gov/34108663/)

Performing analysis using these parameters as input will create elof1_out.txt as an output file, with results as in elof1_out_example.txt
```

