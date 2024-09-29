# Cassava Genomics Project at PuckerLab
# https://www.tu-braunschweig.de/en/ifp/pbb
# Function for TE classification and circos plot taken and adjusted from: https://github.com/LandiMi2/GenomeAssemblyTMEB117
# Cassava Genomics Projekt repository: https://github.com/c-thoben/CassavaGenomicsProject

# Version v0.1

library(optparse)
library(tidyverse)
library(circlize)
library(data.table)
library(RColorBrewer)
library(cowplot)
library(ComplexHeatmap)


option_list = list(
  make_option("--repeat_gff3", type="character", help="Path to repeat GFF3 file from EDTA results"),
  make_option("--repeat_gff3_intact", type="character", help="Path to intact repeat GFF3 file from EDTA results"),
  make_option("--gene_gff3", type="character", help="Path to GFF3 file annotating predicted coding sequences"),
  make_option("--chr_length_A", type="character", help="Path to TSV mapping chromosome IDs to chromosome length for haplophase A"),
  make_option("--chr_length_B", type="character", help="Path to TSV mapping chromosome IDs to chromosome length for haplophase B"),
  make_option("--output_dir", type="character", help="Directory for output files")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# print usage in case of missing arguments
print_usage <- function() {
  print_help(opt_parser)
  quit(status = 1)
}
mandatory_args <- c("repeat_gff3", "repeat_gff3_intact", "gene_gff3", "chr_length_A", "chr_length_B", "output_dir")

missing_args <- mandatory_args[!mandatory_args %in% names(opt) || sapply(opt[mandatory_args], is.null)]
if (length(missing_args) > 0) {
  cat(paste("Error: The following arguments are mandatory and missing:", paste(missing_args, collapse = ", ")), "\n\n")
  print_usage()
}


# Create output directories if they don't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}
if (!dir.exists(file.path(opt$output_dir, "plots"))) {
  dir.create(file.path(opt$output_dir, "plots"), recursive = TRUE)
}
if (!dir.exists(file.path(opt$output_dir, "tables"))) {
  dir.create(file.path(opt$output_dir, "tables"), recursive = TRUE)
}



# --- Load TE elements and classify repeats
message("Loading and classifying TE elements...")

TE <- fread(opt$repeat_gff3, header = FALSE)
TE <- TE[, .(V1, V4, V5, V2, V9, V3)]
colnames(TE) <- c("chr", "start", "end", "id", "info", "classification")

# Classify repeats and filter data
REPEATS <- TE %>% mutate(
  tmp_name = str_extract(info, "Name=[^;]*"),
  family = str_extract(tmp_name, "[^=]*$"),
  tmp_subclass = str_extract(info, "Classification=[^;]*"),
  subclass = str_extract(tmp_subclass, "[^=]*$")
) %>%
  select(chr, start, end, id, family, subclass, classification) %>%
  filter(str_detect(classification, "TIR|LTR|helitron|LINE_element|repeat_region")) %>%
  separate(subclass, c("class", "superfamily"), sep = "/") %>%
  mutate(
    order = case_when(
      str_detect(classification, "TIR") & !str_detect(class, "MITE") ~ "TIR",
      str_detect(class, "MITE") ~ "MITE",
      str_detect(classification, "helitron") ~ "Helitron",
      str_detect(classification, "LTR") ~ "LTR-RT",
      str_detect(classification, "LINE") ~ "LINE",
      str_detect(classification, "repeat_region") ~ "Unclassified repeat"
    )
  ) %>%
  mutate(
    class = case_when(
      order == "TIR" | order == "MITE" | order == "Helitron" ~ "DNA",
      order == "LTR-RT" ~ "retrotransposon",
      order == "LINE" ~ "non-LTR retrotransposon",
      order == "Unclassified repeat" ~ "Other"
    ),
    superfamily = if_else(is.na(superfamily), order, superfamily),
    superfamily = case_when(
      str_detect(superfamily, "Gypsy") ~ "Gypsy",
      str_detect(superfamily, "Copia") ~ "Copia",
      str_detect(superfamily, "Helitron") ~ "Helitron",
      str_detect(superfamily, "LINE") ~ "LINE",
      TRUE ~ superfamily
    )
  )

message("Filtering TE elements...")
REPEATS_all <- REPEATS[REPEATS$order != "Unclassified repeat",]

# Split by haplophase A and B
REPEATS_A <- REPEATS_all %>%
  filter(str_detect(chr, "chr\\d+_A"))
REPEATS_B <- REPEATS_all %>%
  filter(str_detect(chr, "chr\\d+_B"))



# --- Load TE elements and classify intact repeats
message("Loading and classifying intact TE elements...")
TE <- fread( opt$repeat_gff3_intact, h=F)
TE <- TE[, .(V1,V4,V5,V2,V9,V3)]
colnames(TE) <- c("chr", "start", "end", "id", "info", "classification")

# Classification of the annotated TEs - using the classification column
REPEATS_intact <- TE %>% mutate(
  # Extract 'Name' and 'Classification' information from the 'info' column
  tmp_name = str_extract(info, "Name=[^;]*"),
  family = str_extract(tmp_name, "[^=]*$"),
  tmp_subclass = str_extract(info, "Classification=[^;]*"),
  subclass = str_extract(tmp_subclass, "[^=]*$") )%>%
   # Select relevant columns and filter rows based on 'classification'
  select(chr, start, end, id, family, subclass, classification) %>%
  filter(str_detect(classification, "TIR|LTR|helitron|LINE_element|repeat_region")) %>%
  # Separate 'subclass' into 'class' and 'superfamily'
  separate(subclass, c("class", "superfamily"), sep = "/") %>%
   # Classify 'order' 
  mutate(
    order = case_when(
      str_detect(classification, "TIR") & !str_detect(class, "MITE") ~ "TIR",
      str_detect(class, "MITE") ~ "MITE",
      str_detect(classification, "helitron") ~ "Helitron",
      str_detect(classification, "LTR") ~ "LTR-RT",
      str_detect(classification, "LINE") ~ "LINE",
      str_detect(classification, "repeat_region") ~ "Unclassified repeat"
    )
  )%>%mutate(  
    class = case_when(
      order == "TIR" | order == "MITE" | order == "Helitron" ~ "DNA",
      order == "LTR-RT" ~ "retrotransposon",
      order == "LINE" ~ "non-LTR retrotranspson",
      order == "Unclassified repeat" ~ "Other"
    ),
    superfamily = if_else(is.na(superfamily), order, superfamily),
    superfamily = case_when(
      str_detect(superfamily, "Gypsy") ~ "Gypsy",
      str_detect(superfamily, "Copia") ~ "Copia",
      str_detect(superfamily, "Helitron") ~ "Helitron",
      str_detect(superfamily, "LINE") ~ "LINE",
      str_detect(superfamily, "unknown") ~ "unknown",
      str_detect(superfamily, "DTM|Mutator") ~ "MUDR-Mutator",
      str_detect(superfamily, "DTH|Harbinger") ~ "PIF-Harbinger",
      str_detect(superfamily, "DTA|hAT") ~ "hAT",
      str_detect(superfamily, "DTC|CACTA") ~ "CACTA",
      str_detect(superfamily, "DTT|Mariner") ~ "Tc1-Mariner"
    ),
    superfamily = if_else(is.na(superfamily), "-", superfamily)
  )

# Filter: chromosomes haplophase A and only classified repeats
REPEATS_all <- REPEATS_intact[REPEATS_intact$order!="Unclassified repeat",]
REPEATS_A_intact <- REPEATS_all %>%
  filter(str_detect(chr, "chr\\d+_A"))
REPEATS_B_intact <- REPEATS_all %>%
  filter(str_detect(chr, "chr\\d+_B"))



# --- TE Barplot for Hap A and Hap B
message("Creating bar plot...")
colors <- c("#FDB172", "#AAD9A7", "#97809A", "#A0BFDB")

bar_A <- REPEATS_A %>%
  group_by(chr, order) %>%
  summarise(num = n()) %>%
  arrange(chr, order)

bar_B <- REPEATS_B %>%
  group_by(chr, order) %>%
  summarise(num = n()) %>%
  arrange(chr, order)

# Generate bar plots
plota <- ggplot(bar_A, aes(x=fct_reorder(chr, parse_number(chr), .desc = TRUE), y=num, fill = order)) +
  geom_bar(stat="identity") +
  coord_flip() + 
  ggtitle("Hap A") +
  xlab("Chromosomes") + ylab("TE_count") +
  theme_bw() + 
  scale_fill_manual(values = colors)

plotb <- ggplot(bar_B, aes(x=fct_reorder(chr, parse_number(chr), .desc = TRUE), y=num, fill = order)) +
  geom_bar(stat="identity") +
  coord_flip() + 
  ggtitle("Hap B") +
  xlab("Chromosomes") + ylab("TE_count") +
  theme_bw() + 
  scale_fill_manual(values = colors)

cowplot::plot_grid(plota, plotb, labels = "AUTO", greedy = TRUE, rel_widths = c(1, 1.17))
ggsave(paste0(opt$output_dir, "/plots/TEs_barplot.png"), device="png", width=4000, height=3000, units="px")



# --- Calculate genomic densities for repeats (HapA and HapB)
message(" Calculating repeat densities...")

den_a  = genomicDensity(REPEATS_A, window.size = 1e6)
write.table(den_a, paste0(opt$output_dir, "/tables/density_repeats_A.tsv"), 
row.names = FALSE, 
sep="\t")

den_b  = genomicDensity(REPEATS_B, window.size = 1e6)
write.table(den_b, paste0(opt$output_dir, "/tables/density_repeats_B.tsv"), 
row.names = FALSE, 
sep="\t")



# --- Circos plot for HapA
message("Creating Circos plot for HapA...")
cols <- c("#AAD9A7", "#A0BFDB", "#FDB172", "#97809A")

# Prepare chromosome length data
chrlenA <- fread(opt$chr_length_A)
colnames(chrlenA) <- c("chr", "len")
chrlenA <- chrlenA %>%
  mutate(c = as.integer(gsub(".*chr([0-9]+)_.*", "\\1", chr))) %>%
  arrange(c) %>%
  select(chr, len)

chrlenB <- fread(opt$chr_length_B)
colnames(chrlenB) <- c("chr", "len")
chrlenB <- chrlenB %>%
  mutate(c = as.integer(gsub(".*chr([0-9]+)_.*", "\\1", chr))) %>%
  arrange(c) %>%
  select(chr, len)

# Prepare gene data
gff1 <- fread(opt$gene_gff3, h = FALSE)
gff1 <- gff1[, .(V1, V4, V5, V3)]
colnames(gff1) <- c("chr", "start", "end", "name")

genesA <- gff1 %>%
  filter(name == "gene") %>%
  filter(str_detect(chr, "chr\\d+_A")) %>%
  mutate(chr = gsub("^col40_", "", chr))

genesB <- gff1 %>%
  filter(name == "gene") %>%
  filter(str_detect(chr, "chr\\d+_B")) %>%
  mutate(chr = gsub("^col40_", "", chr))


# Initialize Circos plot
circos.clear()
circos.par("track.height" = 0.8, gap.degree = 1, cell.padding = c(0, 0, 0, 0))

png(paste0(opt$output_dir, "/plots/circos_genomic_density_A.png"), width = 3000, height = 3000, res = 300)
circos.initialize(
  factors = c(
    "chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A",
    "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A",
    "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A",
    "chr16_A", "chr17_A", "chr18_A"
  ),
  xlim = matrix(c(rep(0, 18), chrlenA$len), ncol = 2)
)

# chromosomes
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr,
              cex = 1, col = "black",
              facing = "outside", niceFacing = TRUE
  )
}, track.height = 0.06, bg.col = "#f3f3f3") # bg.col = "grey90",  bg.border = F, 

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h = "top", labels.cex = 0.5, col = "#000000")
})

#plot
circos.genomicDensity(REPEATS_A,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1e6,
                      col = cols[1],
                      track.height = 0.12,
                      bg.col = "#fdfdfd"
) 

circos.genomicDensity(REPEATS_A_intact,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1e6,
                      col = cols[2],
                      track.height = 0.12,
                      bg.col = "#fdfdfd"
) 

circos.genomicDensity(genesA,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1e6,
                      col = cols[3],
                      track.height = 0.12,
                      bg.col = "#fdfdfd"
) 

# TE order density
repA = REPEATS_A %>%
  split(f = as.factor(.$order))

bed_list <- list(repA$"LTR-RT", repA$TIR, repA$MITE, repA$Helitron)
bed_col <-  c("#FCBA04", "#F87575", "#6DD3CE", "#b388eb")
circos.genomicRainfall(bed_list, pch = 16, cex = 0.1, col = bed_col, track.height = 0.12)

# legend
lgd_lines = Legend(at = c("TE", "TE intact", "predicted coding sequences"),  legend_gp = gpar(fill= c(cols[1],cols[2],cols[3]), lwd = 3, fontsize=30, title_gp = gpar(fontsize=30)), 
  title_position = "topleft", title = "Densities",
  row_gap = unit(1, "mm"))
lgd_list_vertical = packLegend(lgd_lines)
draw(lgd_list_vertical,  x = unit(0.6, "npc") , y = unit(1, "npc") - unit(11.7, "cm"), just = c("right", "top"))
dev.off()



# --- Circos plot for HapB
message("Creating Circos plot for HapB...")

circos.clear()
circos.par("track.height" = 0.8, gap.degree = 1, cell.padding = c(0, 0, 0, 0))

png(paste0(opt$output_dir, "/plots/circos_genomic_density_B.png"), width = 3000, height = 3000, res = 300)
circos.initialize(
  factors = c(
    "chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B",
    "chr6_B", "chr7_B", "chr8_B", "chr9_B", "chr10_B",
    "chr11_B", "chr12_B", "chr13_B", "chr14_B", "chr15_B",
    "chr16_B", "chr17_B", "chr18_B"
  ),
  xlim = matrix(c(rep(0, 18), chrlenB$len), ncol = 2)
)

# chromosomes
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr,
              cex = 1, col = "black",
              facing = "outside", niceFacing = TRUE
  )
}, track.height = 0.06, bg.col = "#f3f3f3") # bg.col = "grey90",  bg.border = F, 

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h = "top", labels.cex = 0.5, col = "#000000")
})

#plot
circos.genomicDensity(REPEATS_B,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1e6,
                      col =  cols[1],
                      track.height = 0.12,
                      bg.col = "#fdfdfd"
) 

circos.genomicDensity(REPEATS_B_intact,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1e6,
                      col = cols[2],
                      track.height = 0.12,
                      bg.col = "#fdfdfd"
) 

circos.genomicDensity(genesB,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1e6,
                      col = cols[3],
                      track.height = 0.12,
                      bg.col = "#fdfdfd"
)

# TE order densities
repB = REPEATS_B %>%
  split(f = as.factor(.$order))

bed_list <- list(repB$"LTR-RT", repB$TIR, repB$MITE, repB$Helitron)
bed_col <-  c("#FCBA04", "#F87575", "#6DD3CE", "#b388eb")
circos.genomicRainfall(bed_list, pch = 16, cex = 0.1, col = bed_col, track.height = 0.12)

#  legend
lgd_lines = Legend(at = c("LTR-RT", "TIR", "MITE", "Helitron"),  legend_gp = gpar(fill= bed_col, lwd = 3, fontsize=30, title_gp = gpar(fontsize=23)), 
  title_position = "topleft", title = "TE order densities",
  row_gap = unit(1, "mm"))
lgd_list_vertical = packLegend(lgd_lines)
draw(lgd_list_vertical, x = unit(0.57, "npc") , y = unit(1, "npc") - unit(14, "cm"), just = c("right", "bottom"))

dev.off()
