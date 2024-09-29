# Cassava Genomics project
[![DOI](https://zenodo.org/badge/846598077.svg)](https://doi.org/10.5281/zenodo.13858428)

Collection of scripts used in the Cassava Genomics Project

| Script                                                                                                   | Version | Source                                              | Cite                                                |
|-----------------------------------------------------------------------------------------------------------------------|---------|-----------------------------------------------------|-----------------------------------------------------|
| **clean_genomic_fasta.py**: Clean contig identifiers to avoid incompatibility issues                                   | 0.15    | https://github.com/bpucker/GenomeAssembly/           | [10.1101/2023.06.27.546741](https://doi.org/10.1101/2023.06.27.546741)           |
| **contig_stats.py**: Calculate contig statistics                                                                       | 1.31    | https://github.com/bpucker/script_collection/        | [10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321)     |
| **genetic_map_to_fasta.py**: Create input file for ALLMAPS merge command by mapping genetic markers to assembly contigs | 0.2     | -                                                   | TODO                                                |
| **cov_plot.py**: Create assembly coverage plot from coverage file.                                                     | 0.2     | https://github.com/bpucker/At7                      | [10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321)     |
| **coverage_te_plot.py**: Adjusted to create a coverage plot including density of TE repeats for *M. esculenta*.        | 0.5     | -                                                   | [10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321)     |
| **RNAseq_cov_analysis.py**: Analyse coverage of predicted polypeptide sequences by RNAseq data                         | 0.1     | https://github.com/bpucker/GenomeAssembly/           | [10.1101/2023.06.27.546741](https://doi.org/10.1101/2023.06.27.546741)           |
| **TE_repeat_analysis.R**: Analyse repeat density of EDTA results | 0.1       | -                                                   | [10.1038/s41597-023-02800-0](https://doi.org/10.1038/s41597-023-02800-0)                                                   |



## genetic_map_to_fasta.py
```
python genetic_map_to_fasta.py \
--map <FULL_PATH_TO_GENETIC_MAP_FILE> \
--contigs <FULL_PATH_TO_CONTIGS_FILE>
--output <BASE_PATH_TO_OUTPUT_FILE>
[--sim <MINIMUM_SIMILARITY_BEST_HIT]
[--score <MINIMUM_SCORE_BEST_HIT]
```
Mapping of genetic markers (`--map`) to assembly contigs (`--contigs`). Genetic markers are expected in the format of the composite genetic map from *Manihot esculenta* Crantz by the ICGMC (File S2, https://doi.org/10.1534/g3.114.015008).

**Input:**
- `--map`: genetic markers
- `--fasta`: FASTA file containing assembly contigs
- `--output`: Base path to output files (without extension)

**Output:**
- `<output>.fasta`: FASTA file containing marker sequences
- `<output>_blastn_marker_contig_mapping.txt`: BLASTN output
- `<output>_mapped_contigs.csv`: Mapped markers, compatible with ALLMAPS merge command

## coverage_te_plot.py
```
python coverage_te_plot.py \
--coverage_file <FULL_PATH_TO_COVERAGE_FILE> \
--te_file <FULL_PATH_TO_TE_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE> \
--cov <AVERAGE_COVERAGE>
[--res <RESOLUTION, WINDOW_SIZE_FOR_COVERAGE_CALCULATION> 1000]
[--sat <SATURATION, CUTOFF_FOR_MAX_COVERAGE_VALUE> 100.0]
[--num_contigs <NUMBER_OF_CONTIGS_TO_PLOT> 18]
[--max_cov <MAXIMUM_COVERAGE 600]
[--max_chromosome <MAXIMUM_CHROMOSOME_SIZE_BP 55000000]
```
Creates a coverage plot showing the average coverage in blocks of resolution size (`--res`). The maximum displayed coverage is defined by the saturation (`--sat`). Each chromosome is plotted separatly and the average coverage is marked by a red line. For each chromosome, a histogram is created showing the coverage distribution.

**Input:**
- `--coverage_file`: Coverage file
- `--te_file`: Repeats TSV created with Circos genomicDensity function
- `--cov`: Average coverage
- `--out`: Base path to output files (without extension)

**Output:**
- `<out>.png`: Coverage plot
- `<out>_<contig>.png`: Histogram of coverage for each contig/chromosome
- `<out>_coverage_resolution<res>.tsv`: Average coverage per block

## TE_repeat_analysis.R
```
Rscript TE_repeat_analysis.R \
--repeat_gff3 <FULL_PATH_TO_REPEAT_GFF3_FILE> \
--repeat_gff3_intact <FULL_PATH_TO_INTACT_REPEAT_GFF3_FILE> \
--gene_gff3 <FULL_PATH_TO_GENE_GFF3_FILE> \
--chr_length_A <FULL_PATH_TO_CHR_LENGTH_A_FILE> \
--chr_length_B <FULL_PATH_TO_CHR_LENGTH_B_FILE> \
--output_dir <FULL_PATH_TO_OUTPUT_DIRECTORY>
```
Analyzes EDTA annotation results by classifying repeats. A bar plot is generated to show the distribution of transposable element (TE) families across chromosomes for each haplophase separately. Additionally, a circos density plot is produced, displaying the genomic density of TE repeats, intact TE repeats, and predicted coding sequences. The inner track of the plot includes a rainfall (rainbow) plot that illustrates the minimal distance between neighboring repeats for each TE family.

**Input**:
- `--repeat_gff3`: Path to repeat GFF3 file from EDTA results
- `--repeat_gff3_intact`: Path to intact repeat GFF3 file from EDTA results
- `--gene_gff3`: Path to GFF3 file annotating predicted coding sequences
- `--chr_length_A`: Path to TSV mapping chromosome IDs to chromosome length for haplophase A
- `--chr_length_B`: Path to TSV mapping chromosome IDs to chromosome length for haplophase B
- `--output_dir`: Directory for output files

**Output**:
- `<output_dir>/tables/TEs_barplot.png`: Barplot showing distribution of TE families for each chromosome
- `<output_dir>/tables/density_repeats_A.tsv`: Density of repeats for haplophase A
- `<output_dir>/tables/density_repeats_B.tsv`: Density of repeats for haplophase B
- `<output_dir>/plots/circos_genomic_density_A.png`: Circos plot for haplophase A
- `<output_dir>/plots/circos_genomic_density_B.png`: Circos plot for haplophase B
