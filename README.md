# Cassava Genomics project
Collection of scripts used in the Cassava Genomics Project

| Script                         | Version | Source                                                   | Cite                                                     | Description                                                                          |
|---|---|---|---|---|
| clean_genomic_fasta.py          | 0.15    | https://github.com/bpucker/GenomeAssembly/                | https://doi.org/10.1101/2023.06.27.546741                | Clean contig identifiers to avoid incompatibility issues                             |
| contig_stats.py                 | 1.31    | https://github.com/bpucker/script_collection/             | http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164321 | Calculate contig statistics                                                          |
| genetic_map_to_fasta.py         | 0.2     | -                                                         | TODO                                                     | Create input file for ALLMAPS merge command by mapping genetic markers to assembly contigs |
| cov_plot.py                     | 0.2     | https://github.com/bpucker/At7                            | https://doi.org/10.1371/journal.pone.0164321             | Create assembly coverage plot from coverage file.                                    |
| coverage_te_plot.py  | 0.5     | -                                                         | https://doi.org/10.1371/journal.pone.0164321             | The cov_plot.py script is adjusted to create a coverage plot including density of TE repeats. Script adjusted for coverage plot of *M. esculenta*. |
| RNAseq_cov_analysis.py          | 0.1     | https://github.com/bpucker/GenomeAssembly/                | https://doi.org/10.1101/2023.06.27.546741                | Analyse coverage of predicted polypeptide sequences by RNAseq data                   |
| COL40_TE_repeat_analysis.ipynb  |         |                                                           |                                                          |                                                                                      |

### genetic_map_to_fasta.py
```
python genetic_map_to_fasta.py \
--map <FULL_PATH_TO_GENETIC_MAP_FILE> \
--contigs <FULL_PATH_TO_CONTIGS_FILE>
--output <BASE_PATH_TO_OUTPUT_FILE> \
[--sim <MINIMUM_SIMILARITY_BEST_HIT]
[--score <MINIMUM_SCORE_BEST_HIT]
```
Mapping of genetic markers (--map) to assembly contigs (--contigs). Genetic markers are expected in the format of the composite genetic map from *Manihot esculenta* Crantz by the ICGMC (File S2, https://doi.org/10.1534/g3.114.015008).

**Input:**
- \-\-map: genetic markers
- \-\-fasta: FASTA file containing assembly contigs
- \-\-output: Base path to output files (without extension)

**Output:**
- \<output>.fasta: FASTA file containing marker sequences
- \<output>_blastn_marker_contig_mapping.txt: BLASTN output
- \<output>_mapped_contigs.csv: Mapped markers, compatible with ALLMAPS merge command

### coverage_te_plot.py
```
python coverage_te_plot.py \
--coverage_file <FULL_PATH_TO_COVERAGE_FILE> \
--te_file <FULL_PATH_TO_TE_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE> \
--cov <AVERAGE_COVERAGE> \
[--res <RESOLUTION, WINDOW_SIZE_FOR_COVERAGE_CALCULATION> 1000]
[--sat <SATURATION, CUTOFF_FOR_MAX_COVERAGE_VALUE> 100.0]
[--num_contigs <NUMBER_OF_CONTIGS_TO_PLOT> 18]
[--max_chromosome <MAXIMUM_CHROMOSOME_SIZE_BP 55000000]
```
Creates a coverage plot showing the average coverage in blocks of resolution size (--res). The maximum displayed coverage is defined by the saturation (--sat). Each chromosome is plotted separatly and the average coverage is marked by a red line. For each chromosome, a histogram is created showing the coverage distribution.

**Input:**
- \-\-coverage_file: Coverage file
- \-\-te_file: Repeats TSV created with Circos genomicDensity function
- \-\-cov: Average coverage
- \-\-out: Base path to output files (without extension)

**Output:**
- \<out>.png: Coverage plot
- \<out>_\<contig>.png: Histogram of coverage for each contig/chromosome
- \<out>_coverage_resolution\<res>.tsv: Average coverage per block


## References
TODO
