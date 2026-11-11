Sys.setenv(LANGUAGE='en')
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = FALSE, scipen = 999)
rm(list = ls())
library(data.table)
# ===================== bash arguments, do not modify =====================
args <- commandArgs(trailingOnly = TRUE)
selected_species <- args[1]
min_len_default <- as.integer(args[2])
max_len_default <- as.integer(args[3])
min_len <- as.integer(args[4])
max_len <- as.integer(args[5])
selected_analysis <- args[6]
trim_nt <- as.integer(args[7])
csv_path <- args[8]
# Add: gene list file path (passed via argument or default path)
gene_list_path <- if(length(args)>=9) args[9] else "list.txt"  # prioritize argument, default list.txt

# ===================== codon2aa, all U replaced with T (DNA standard) =====================
library(data.table)
codon2aa <- data.table(
  codon = c("TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG",
            "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
            "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
            "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG",
            "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
            "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"),
  aa = c("F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "I", "M",
         "V", "V", "V", "V", "S", "S", "S", "S", "S", "S",
         "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
         "Y", "Y", "*", "*", "H", "H", "Q", "Q", "N", "N", "K", "K",
         "D", "D", "E", "E", "C", "C", "*", "W", "R", "R", "R", "R",
         "R", "R", "G", "G", "G", "G")
)

# ===================== species reference configuration, unchanged =====================
species_config<- list(
  human = list(
    gtf = "~/reference/human_transcriptome/gencode.v49.pc_annotation.gtf",
    fasta = "~/reference/human_transcriptome/gencode.v49.pc_transcripts.fa"
  ),
  mouse = list(
    gtf = "~/reference/mouse_transcriptome/gencode.vM35.pc_annotation.gtf",
    fasta = "~/reference/mouse_transcriptome/gencode.vM35.pc_transcripts.fa"
  )
)
gtf_path <- species_config[[selected_species]]$gtf
transcriptome_fasta <- species_config[[selected_species]]$fasta

# ===================== load pre-generated tx2gene mapping table by species, no header =====================
cat("✅ Loading pre-generated transcript-to-gene mapping...\n")

# 1. Retrieve GTF path and folder from species config
gtf_path <- species_config[[selected_species]]$gtf
gtf_dir <- dirname(gtf_path)

# 2. Match tx2gene filename automatically for human/mouse
tx2gene_filename <- dplyr::case_when(
  selected_species == "human" ~ "tx2gene_human.txt",
  selected_species == "mouse" ~ "tx2gene_mouse.txt"
)
tx2gene_path <- file.path(gtf_dir, tx2gene_filename)

# 3. Auto-generate tx2gene from GTF if missing
if (!file.exists(tx2gene_path)) {
  cat(paste0("⚠️ tx2gene not found, auto-generate from GTF: ", gtf_path, "\n"))
  gtf_raw <- data.table::fread(gtf_path, sep = "\t", header = FALSE, skip = "#", fill = TRUE)
  attr <- gtf_raw[[ncol(gtf_raw)]]
  
  tid  <- sub('.*transcript_id "([^"]+)".*', "\\1", attr)
  gname <- sub('.*gene_name "([^"]+)".*', "\\1", attr)
  
  valid <- tid != attr & gname != attr
  tx2gene <- data.table(
    transcript = tid[valid],
    gene_name  = gname[valid]
  )
  tx2gene <- unique(tx2gene, by = "transcript")
  
  data.table::fwrite(tx2gene, file = tx2gene_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(paste0("✅ Auto generated tx2gene file: ", tx2gene_path, "\n"))
}

# 4. Read headerless tx2gene file, assign column names manually
tx2gene <- data.table::fread(
  tx2gene_path,
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("transcript", "gene_name")
)

# 5. Deduplication, enforce one transcript maps to one gene_name
tx2gene <- unique(tx2gene, by = "transcript")

# 6. Logging for validation
cat(paste0("✅ Successfully loaded ", nrow(tx2gene), " transcript-gene mappings\n"))
cat(paste0("   Example: ", tx2gene[1, transcript], " → ", tx2gene[1, gene_name], "\n\n"))

# ===================== workflow description =====================
if(selected_analysis == "disome"){
  flow_text <- paste0("Disome pipeline - 5' trim ",trim_nt,"nt + 5' Rear Codon inference + 3nt periodicity correction")
}else{
  flow_text <- paste0("Standard Monosome pipeline (✅ Load TXT + RDS for footprint)")
}
cat(paste0("\n✅ Loaded reference config for species: ", selected_species, "\n"))
cat(paste0("   ✅ Core rule: Codon stat TXT + RDS footprint data\n"))
cat(paste0("   ✅ PCA logic: sample=row | codon=feature | avg_occupancy=expression\n"))
cat(paste0("   ✅ Codon base corrected: Uracil(U) -> Thymine(T)\n"))
cat(paste0("   ✅ Differential test: KS test for EACH CODON (replace edgeR, best for occupancy)\n"))
cat(paste0("   ✅ Visualization: Sample-specific curve (each sample one line) NO HEATMAP\n"))
cat(paste0("   ✅ Color style: Same group with similar color (different shade)\n"))
cat(paste0("   ✅ Selected pipeline: ", flow_text, "\n"))

# ---------------------- path and directory configuration, unchanged ----------------------
WORK_DIR <- getwd()
output_root <- paste0(WORK_DIR, "/5.ribo_group_diff_output")
rds_root_dir <- paste0(WORK_DIR, "/4.riboWaltz_output")
report_path <- file.path(output_root, "batch_run_report.txt")

# Create output directory
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_root)) {
  dir.create(output_root, recursive = TRUE)
  cat("✅ Output directory created: ", output_root, "\n")
}

# Redirect log
sink(report_path, split = TRUE)

cat("============================================================\n")
cat(paste0("Codon Occupancy Analysis Report (✅ TXT+RDS) - ", flow_text))
cat("\nAnalysis time: ", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

# ---------------------- load dependencies ----------------------
library(dplyr)          
library(data.table)     
library(ggplot2)        
library(Biostrings)     
library(pheatmap)       
library(RColorBrewer)   
library(cowplot)        
library(gridExtra)      
library(tidyr)
library(FactoMineR)
library(factoextra)

# ---------------------- core parameter settings ----------------------
threads <- 16           
rpf_count_threshold <- 1e6  
bin_size <- 0.5         
length_range <- min_len:max_len 
p_threshold <- 0.1      # KS test P-value threshold

# ✅ Add: similar-color palette for groups (core optimization)
# Each group assigned a base hue; replicates vary by lightness; auto-adaptive for arbitrary group count
get_group_similar_color <- function(sample_names, group_dt) {
  group_list <- unique(group_dt$group_name)
  color_palette <- c("#E63946","#457B9D","#2A9D8F","#F4A261","#264653","#8D99AE","#EF476F","#06D6A0")
  names(color_palette) <- group_list
  sample_color <- c()
  for(s in sample_names){
    g <- group_dt[group_dt$sample_name == s,]$group_name
    idx <- which(sample_names[sapply(sample_names, function(x) group_dt[group_dt$sample_name == x,]$group_name)==g]==s)
    shade <- 0.2*(idx-1)
    rgb_val <- col2rgb(color_palette[g])/255
    new_rgb <- pmin(rgb_val + shade, 1)
    sample_color[s] <- rgb(new_rgb[1], new_rgb[2], new_rgb[3])
  }
  return(sample_color)
}

# Print environment and parameter log
cat("1. Analysis Configuration\n")
cat("------------------------------------------------------------\n")
cat("- System Language: English\n")
cat("- Disable string factors: TRUE\n")
cat("- Codon base correction: U -> T (DNA standard)\n")
cat("- PCA optimal logic: sample=row | codon=feature | avg_occupancy=expression\n")
cat("- Differential Method: Kolmogorov-Smirnov Test (per codon, non-parametric)\n")
cat("- Differential threshold: PValue <", p_threshold, "\n")
cat("- Visualization: Sample curve plot (each sample one line) | REMOVE HEATMAP\n")
cat("- Color Style: Same group = similar color (different shade) | Different group = different color\n")
cat("- Analysis sites: A/P/E codon sites (full retained)\n")
cat("- Output plots: PCA + Footprint density*3 + Ratio histogram*3 + Codon sample curve + KS test\n\n")

# ===================== read sample group CSV, unchanged =====================
cat("2. Load Sample Group Info & Match File Paths\n")
cat("------------------------------------------------------------\n")
group_dt <- fread(csv_path, sep = ",")
group_dt <- group_dt[!is.na(group_dt$sample_name) & !is.na(group_dt$group_name), ]
sample_names <- unique(group_dt$sample_name)

txt_paths <- sapply(sample_names, function(x) {
  file.path(rds_root_dir, x, "all_site_codon_occupancy_final.txt")
})
rds_paths <- sapply(sample_names, function(x) {
  file.path(rds_root_dir, x, "reads_psite_raw_version.rds")
})
names(txt_paths) <- sample_names
names(rds_paths) <- sample_names

# Validate file existence
missing_txt <- names(txt_paths)[!file.exists(txt_paths)]
missing_rds <- names(rds_paths)[!file.exists(rds_paths)]
missing_samples <- unique(c(missing_txt, missing_rds))
if (length(missing_samples) > 0) {
  cat(paste0("⚠️  Missing files for samples: ", paste(missing_samples, collapse = ", "), "\n"))
  sample_names <- sample_names[!sample_names %in% missing_samples]
  txt_paths <- txt_paths[!names(txt_paths) %in% missing_samples]
  rds_paths <- rds_paths[!names(rds_paths) %in% missing_samples]
}

if (length(sample_names) == 0) {
  sink()
  stop("❌ No valid samples with complete files! Exit analysis.")
}

cat(paste0("✅ Matched ", length(sample_names), " valid samples with ", length(unique(group_dt$group_name)), " groups\n"))
group_summary <- group_dt[, list(N = .N), by = list(group_name)]
group_detail <- paste0(group_summary$group_name, "(", group_summary$N, ")") %>% paste(collapse = " | ")
cat(paste0("✅ Sample groups: ", group_detail, "\n\n"))
group_dt <- group_dt[group_dt$sample_name %in% sample_names, ]

# ✅ Generate sample colors: same group shares hue with different lightness
sample_color_vec <- get_group_similar_color(sample_names, group_dt)

# ===================== core: batch load data, low-R version compatible syntax =====================
cat("3. Load Codon Occupancy Data & Footprint Data\n")
cat("------------------------------------------------------------\n")
all_codon_occupancy_raw <- data.table()
all_psite_raw_dt <- data.table()

for (sample in sample_names) {
  cat(paste0("   Processing sample: ", sample, "\n"))
  txt_path <- txt_paths[sample]
  codon_occupancy_raw <- fread(txt_path, sep = "\t", stringsAsFactors = F)
  codon_occupancy_raw$codon_type <- paste0(codon_occupancy_raw$site, "_codon")
  codon_occupancy_raw$sample_name <- sample
  codon_occupancy_raw$group_name <- group_dt[group_dt$sample_name == sample, ]$group_name
  all_codon_occupancy_raw <- rbindlist(list(all_codon_occupancy_raw, codon_occupancy_raw), fill = TRUE)
  
  rds_path <- rds_paths[sample]
  reads_psite_dt <- readRDS(rds_path)
  setDT(reads_psite_dt)
  reads_psite_dt$sample_name <- sample
  reads_psite_dt$group_name <- group_dt[group_dt$sample_name == sample, ]$group_name
  
  # ===================== Add: merge gene name into footprint table =====================
  reads_psite_dt <- merge(reads_psite_dt, tx2gene, by = "transcript", all.x = TRUE)
  
  all_psite_raw_dt <- rbindlist(list(all_psite_raw_dt, reads_psite_dt), fill = TRUE)
  
  cat(paste0("     ✅ Done: Codon data=",nrow(codon_occupancy_raw)," | Footprint data=",nrow(reads_psite_dt),"\n"))
}

# Calculate group mean, low-version compatible syntax
all_codon_occupancy_group <- all_codon_occupancy_raw[, list(
  mean_occupancy = mean(as.numeric(avg_occupancy)),
  Normalized_mean = mean(as.numeric(Normalized_avg)),
  sd_occupancy = sd(as.numeric(avg_occupancy)),
  se_occupancy = sd(as.numeric(avg_occupancy))/sqrt(.N),
  total_reads = sum(as.numeric(total_reads)),
  transcript_count = mean(as.numeric(transcript_count))
), by = list(group_name, codon_type, codon, site)]
all_codon_occupancy_group$norm_occupancy <- all_codon_occupancy_group$mean_occupancy / ave(all_codon_occupancy_group$mean_occupancy, all_codon_occupancy_group$codon_type, FUN=mean)

# Save tables
fwrite(all_codon_occupancy_raw, file.path(output_root, "all_sample_APEcodon_occupancy_raw.txt"), sep = "\t")
fwrite(all_codon_occupancy_group, file.path(output_root, "all_sample_APEcodon_occupancy_group.txt"), sep = "\t")
cat("✅ All data loaded successfully! Raw occupancy retained for KS test\n\n")

# ===================== PCA function: suppress ellipse warning, compatible syntax =====================
plot_Psite_PCA <- function(all_codon_occupancy_raw, group_dt, output_dir){
  cat("4. Generate PCA Plot (P-site Codon Occupancy)\n")
  p_raw <- all_codon_occupancy_raw[all_codon_occupancy_raw$codon_type == "P_codon", c("sample_name","codon","avg_occupancy")]
  pca_matrix <- dcast(p_raw, sample_name ~ codon, value.var = "avg_occupancy", fill = 0)
  sample_grp <- group_dt[match(pca_matrix$sample_name, group_dt$sample_name), ]$group_name

  pca_res <- PCA(pca_matrix[, -1], scale.unit = TRUE, ncp = 5, graph = FALSE)

  p <- fviz_pca_ind(pca_res,
                    # Core modification: geom.ind add "text" while keeping "point", show points + sample labels
                    geom.ind = c("point", "text"),
                    pointshape = 21, pointsize = 4,
                    fill.ind = sample_grp,
                    col.ind = "black",
                    palette = "Set1",
                    label = "all", repel = TRUE,
                    # Add: adjust text label size to avoid overlap, minimal tuning without breaking style
                    text.size = 3,
                    title = "PCA Analysis - P-site Codon Occupancy (sample=row | codon=feature)",
                    xlab = paste0("PC1 (", sprintf("%.1f", pca_res$eig[1,2]), "%)"),
                    ylab = paste0("PC2 (", sprintf("%.1f", pca_res$eig[2,2]), "%)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 9),
          legend.position = "top",
          legend.title = element_text(size=10),
          legend.text = element_text(size=9))

  ggsave(file.path(output_dir, "0.Psite_all_sample_PCA.pdf"), p, width=10, height=8, dpi=300)
  cat("✅ PCA plot generated successfully! (No ellipse warning, sample names displayed)\n\n")
# ---------------- Minor addition: sample correlation heatmap with group annotation bar ----------------
  mat <- as.matrix(pca_matrix[, -1])
  rownames(mat) <- pca_matrix$sample_name
  sample_cor <- cor(t(mat), method = "pearson")
  
  # Build group annotation (group bar on heatmap top)
  annotation_df <- data.frame(
    Group = sample_grp
  )
  rownames(annotation_df) <- pca_matrix$sample_name
  
  # Plot
  pdf(file.path(output_dir, "0.Psite_sample_correlation_heatmap.pdf"), width=9, height=7)
  pheatmap(
    sample_cor,
    scale = "none",
    color = colorRampPalette(c("#457B9D", "white", "#E63946"))(100),
    main = "Sample Correlation (P-site Codon Occupancy, Pearson)",
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 7,
    treeheight_row = 15,
    treeheight_col = 15,
    annotation_col = annotation_df,    # column group annotation
    annotation_row = annotation_df,    # row group annotation
    annotation_names_col = FALSE,
    annotation_names_row = FALSE
  )
  dev.off()
  
  cat("✅ PCA plot generated successfully! (No ellipse warning, sample names displayed)\n")
  cat("✅ New: Sample correlation heatmap with GROUP ANNOTATION added ✅\n\n")
}

# ===================== footprint density three-panel plot, full compatible syntax =====================
plot_group_footprint <- function(all_psite_raw_dt, output_dir) {
  cat("5. Generate Ribosome Footprint Density Plots\n")
  cat("------------------------------------------------------------\n")
  psite_cds_dt <- all_psite_raw_dt[!is.na(all_psite_raw_dt$cds_start) & !is.na(all_psite_raw_dt$cds_stop) & !is.na(all_psite_raw_dt$psite), ]
  psite_cds_dt$l_cds <- psite_cds_dt$cds_stop - psite_cds_dt$cds_start + 1
  psite_cds_dt$cds_relative_pos <- (psite_cds_dt$psite - psite_cds_dt$cds_start + 1) / psite_cds_dt$l_cds * 100

  # Global CDS density plot
  cds_density_dt <- psite_cds_dt[!is.na(psite_cds_dt$cds_relative_pos), ]
  cds_density_dt$cds_bin <- round(cds_density_dt$cds_relative_pos / bin_size) * bin_size
  cds_density_dt <- cds_density_dt[cds_density_dt$cds_bin >= -25 & cds_density_dt$cds_bin <= 125, ]
  cds_density_dt <- cds_density_dt[, list(psite_count = .N), by = list(group_name, cds_bin)]
  cds_density_dt <- cds_density_dt[order(cds_density_dt$group_name, cds_density_dt$cds_bin), ]
  cds_density_dt$global_avg <- ave(cds_density_dt$psite_count, cds_density_dt$group_name, FUN=mean)
  cds_density_dt$relative_density <- cds_density_dt$psite_count / cds_density_dt$global_avg

  final_plot <- ggplot(cds_density_dt, aes(x = cds_bin, y = relative_density, color=group_name, group=group_name)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
    geom_vline(xintercept = 100, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
    annotate("text", x = 2, y = max(cds_density_dt$relative_density)*0.9, label = "Start Codon (ATG)", color = "#333f50", size = 3.5, hjust = 0) +
    annotate("text", x = 98, y = max(cds_density_dt$relative_density)*0.9, label = "Stop Codon", color = "#333f50", size = 3.5, hjust = 1) +
    labs(x = "Normalized CDS Length (%)", y = "Relative Footprint Density (a.u.)",
         title = "Ribosome Footprint Density (Raw Periodicity) - Group Mean",
         subtitle = paste0("Binsize = ", bin_size, "% | Length filter: ", min_len, "-", max_len, "nt")) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  ggsave(file.path(output_dir, "1.group_footprint_global_CDS_percent.pdf"), final_plot, width=12, height=6, dpi=300)

  # Start codon density plot
  psite_cds_dt$start_nt_pos <- psite_cds_dt$psite - psite_cds_dt$cds_start
  start_codon_dt <- psite_cds_dt[!is.na(psite_cds_dt$start_nt_pos) & psite_cds_dt$start_nt_pos >= -25 & psite_cds_dt$start_nt_pos <= 125, ]
  start_codon_dt <- start_codon_dt[, list(psite_count = .N), by = list(group_name, start_nt_pos)]
  start_codon_dt <- start_codon_dt[order(start_codon_dt$group_name, start_codon_dt$start_nt_pos), ]
  start_codon_dt$global_avg <- ave(start_codon_dt$psite_count, start_codon_dt$group_name, FUN=mean)
  start_codon_dt$relative_density <- start_codon_dt$psite_count / start_codon_dt$global_avg

  start_plot <- ggplot(start_codon_dt, aes(x = start_nt_pos, y = relative_density, color=group_name, group=group_name)) +
    geom_line(linewidth = 1, alpha = 0.9) + geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
    annotate("text", x = 5, y = max(start_codon_dt$relative_density)*0.9, label = "ATG (Start Codon)", color = "#e84a5f", size = 4, hjust = 0, fontface="bold") +
    labs(x = "Position relative to Start Codon (nt)", y = "Relative Footprint Density (a.u.)",
         title = "Footprint Density (Start Codon -25 nt to +125 nt) - Group Mean") +
    scale_color_brewer(palette = "Set1") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  ggsave(file.path(output_dir, "2.group_footprint_StartCodon.pdf"), start_plot, width=10, height=5, dpi=300)

  # Stop codon density plot
  psite_cds_dt$stop_nt_pos <- psite_cds_dt$psite - psite_cds_dt$cds_stop
  stop_codon_dt <- psite_cds_dt[!is.na(psite_cds_dt$stop_nt_pos) & psite_cds_dt$stop_nt_pos >= -125 & psite_cds_dt$stop_nt_pos <= 25, ]
  stop_codon_dt <- stop_codon_dt[, list(psite_count = .N), by = list(group_name, stop_nt_pos)]
  stop_codon_dt <- stop_codon_dt[order(stop_codon_dt$group_name, stop_codon_dt$stop_nt_pos), ]
  stop_codon_dt$global_avg <- ave(stop_codon_dt$psite_count, stop_codon_dt$group_name, FUN=mean)
  stop_codon_dt$relative_density <- stop_codon_dt$psite_count / stop_codon_dt$global_avg

  stop_plot <- ggplot(stop_codon_dt, aes(x = stop_nt_pos, y = relative_density, color=group_name, group=group_name)) +
    geom_line(linewidth = 1, alpha = 0.9) + geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
    annotate("text", x = -5, y = max(stop_codon_dt$relative_density)*0.9, label = "Stop Codon", color = "#e84a5f", size = 4, hjust = 1, fontface="bold") +
    labs(x = "Position relative to Stop Codon (nt)", y = "Relative Footprint Density (a.u.)",
         title = "Footprint Density (Stop Codon -125 nt to +25 nt) - Group Mean") +
    scale_color_brewer(palette = "Set1") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  ggsave(file.path(output_dir, "3.group_footprint_StopCodon.pdf"), stop_plot, width=10, height=5, dpi=300)

  cat("✅ 3 Footprint density plots generated successfully!\n\n")
}

# ===================== ✅ Fully fixed: gene-specific footprint plot (raw density + group read count annotation at top-right) ✅ =====================
plot_gene_specific_footprint <- function(all_psite_raw_dt, gene_list_path, tx2gene, output_dir) {
  cat("5.1 Generate Gene-specific Ribosome Footprint Plots (Single TX + Gene Sum/Mean (RPM) + Longest CDS Mark + Frame stats)\n")
  cat("------------------------------------------------------------\n")
  
  required_cols <- c("psite", "cds_start", "cds_stop", "transcript", "group_name", "sample_name")
  missing_cols <- setdiff(required_cols, colnames(all_psite_raw_dt))
  if(length(missing_cols) > 0) {
    stop(paste0("❌ Missing required columns: ", paste(missing_cols, collapse = ", "), "\n"))
  }
  
  global_group_total <- all_psite_raw_dt[, .(global_total = .N), by = group_name]
  
  target_genes <- fread(gene_list_path, header = FALSE)$V1
  target_genes <- unique(target_genes[!is.na(target_genes) & target_genes != ""])
  if(length(target_genes) == 0) return()
  
  tx2gene_filtered <- tx2gene[gene_name %in% target_genes, ]
  gene2tx <- split(tx2gene_filtered$transcript, tx2gene_filtered$gene_name)
  
  psite_dt <- all_psite_raw_dt[, ..required_cols]
  psite_dt[, atg_pos := cds_start]
  psite_dt[, cds_length := cds_stop - cds_start + 1]
  psite_dt[, frame := (psite - cds_start) %% 3]
  psite_dt[, in_cds := psite >= cds_start & psite <= cds_stop]
  
  gene_plot_root <- file.path(output_dir, "gene_specific_footprint")
  dir.create(gene_plot_root, recursive = TRUE, showWarnings = FALSE)
  
  total_tx_plots <- 0
  total_sum_plots <- 0
  total_mean_plots <- 0
  
  for(gene in target_genes) {
    if(!gene %in% names(gene2tx)) next
    target_transcripts <- gene2tx[[gene]]
    gene_dir <- file.path(gene_plot_root, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
    
    gene_all_dt <- psite_dt[transcript %in% target_transcripts, ]
    if(nrow(gene_all_dt) == 0) next
    
    tx_info_all <- unique(gene_all_dt[, .(transcript, atg_pos, cds_length)])
    setorder(tx_info_all, -cds_length)
    longest_tx <- tx_info_all$transcript[1]
    longest_tx_atg <- tx_info_all$atg_pos[1]
    longest_tx_cds_len <- tx_info_all$cds_length[1]
    longest_tx_cds_end <- longest_tx_atg + longest_tx_cds_len - 1
    
    # ===========================================================================
    # 1) Single transcript plot [independent RPM plot per tx + frame grouped side-by-side bar plot]
    # ===========================================================================
    for(tx in target_transcripts) {
      tx_dt <- psite_dt[transcript == tx, ]
      if(nrow(tx_dt) == 0) next
      
      tx_atg <- unique(tx_dt$atg_pos)
      tx_cds_len <- unique(tx_dt$cds_length)
      tx_cds_end <- tx_atg + tx_cds_len - 1
      
      tx_density <- tx_dt[, .(raw = .N), by = .(group_name, psite, frame)]
      tx_density <- merge(tx_density, global_group_total, by = "group_name")
      tx_density[, rpm := (raw / global_total) * 1e6]
      
      fwrite(tx_density[, .(group_name, psite, frame, raw, global_total, rpm)],
             file.path(gene_dir, paste0("gene_", gene, "_tx_", tx, "_rpm_table.txt")), sep = "\t")
      
      tx_sum <- tx_density[, .(gene_reads = sum(raw)), by = group_name]
      tx_sum <- merge(tx_sum, global_group_total, by = "group_name")
      tx_sum[, rpm_total := (gene_reads / global_total) * 1e6]
      tx_sum[, label := paste0(
        group_name, ": ", format(gene_reads, big.mark=","),
        " (", format(global_total, big.mark=","), ") | RPM: ", round(rpm_total, 2)
      )]
      tx_label <- paste(tx_sum$label, collapse = "\n")
      
      p_tx <- ggplot(tx_density, aes(x=psite, y=rpm, fill=group_name)) +
        geom_col(width=0.8, alpha=0.8, position="identity") +
        geom_vline(xintercept=tx_atg, color="#e84a5f", linewidth=1.2) +
        {if(tx == longest_tx) annotate("text", x=max(tx_density$psite), y=max(tx_density$rpm)*1.1,
                                       label="* Longest TX", color="#2E8B57", size=4, hjust=1)} +
        geom_vline(xintercept=tx_atg + tx_cds_len - 1, color="#2A9D8F", linewidth=1, linetype="dashed") +
        annotate("text", x=Inf, y=Inf, label=tx_label, hjust=1.1, vjust=1.1, size=3.3, fontface="bold") +
        labs(x="PSite", y="RPM (only / group total reads)",
             title=paste0(gene, " | ", tx), fill="Group") +
        scale_fill_brewer(palette="Set1") + theme_bw()
      
      f <- paste0("gene_", gene, "_tx_", tx, "_rpm.pdf")
      ggsave(file.path(gene_dir, f), p_tx, width=14, height=6)
      total_tx_plots <- total_tx_plots + 1
      
      # ========== Frame plot fix: aggregate frame counts by group+psite, then compute fraction ==========
      tx_cds_frame_dt <- tx_dt[in_cds == TRUE, ]
      if(nrow(tx_cds_frame_dt) > 0) {
        # Add: codon index, bin by codon instead of raw psite
        tx_cds_frame_dt[, codon_idx := floor( (psite - cds_start) / 3 ) ]
        
        # Aggregate read count by group + codon_idx + frame
        group_frame_raw <- tx_cds_frame_dt[, .(cnt = .N), by = .(group_name, codon_idx, frame)]
        group_frame_cast <- dcast(group_frame_raw, group_name + codon_idx ~ frame, value.var = "cnt", fill = 0)
        
        frame_cols <- intersect(as.character(c(0,1,2)), colnames(group_frame_cast))
        f12_cols   <- intersect(c("1","2"), frame_cols)
        
        group_frame_cast[, total  := rowSums(.SD, na.rm=TRUE), .SDcols = frame_cols]
        group_frame_cast[, sum12  := rowSums(.SD, na.rm=TRUE), .SDcols = f12_cols]
        group_frame_cast[, non0_frac := ifelse(total == 0, NA, sum12 / total * 100)]
        
        fwrite(group_frame_cast,
               file.path(gene_dir, paste0("gene_", gene, "_tx_", tx, "_non0_frame_pct_table.txt")), sep = "\t")
        
        p_tx_frame <- ggplot(group_frame_cast, aes(x = codon_idx, y = non0_frac, fill = group_name)) +
          geom_col(position = position_dodge(width = 2), width = 1.8) +
          geom_hline(yintercept = 66.67, linetype = "dotted", color = "gray50") +
          # ATG corresponds to codon_idx=0
          geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2) +
          labs(x = "Codon index (starting from ATG)",
               y = "Out-of-frame fraction (frame1+2 / frame0+1+2) (%)",
               title = paste0(gene, " | ", tx, " | Out-of-frame Percentage per codon"),
               fill = "Group") +
          theme_bw()
        ggsave(file.path(gene_dir, paste0("gene_", gene, "_tx_", tx, "_non0_frame_pct.pdf")), p_tx_frame, width = 16, height = 6)
      }
      

      # ======================================================================
    }
    
    # ===========================================================================
    # 2) Gene SUM plot (unchanged)
    # ===========================================================================
    gene_sum <- gene_all_dt[, .(raw_sum = .N), by = .(group_name, psite, frame)]
    gene_sum <- merge(gene_sum, global_group_total, by = "group_name")
    gene_sum[, rpm_sum := (raw_sum / global_total) * 1e6]
    
    gene_sum_agg <- gene_sum[, .(raw_sum = sum(raw_sum), rpm_sum = sum(rpm_sum)), by = .(group_name, psite)]
    
    gene_sum_total <- gene_sum[, .(gene_reads = sum(raw_sum)), by = group_name]
    gene_sum_total <- merge(gene_sum_total, global_group_total, by = "group_name")
    gene_sum_total[, rpm_total := (gene_reads / global_total) * 1e6]
    gene_sum_total[, label := paste0(
      group_name, ": ", format(gene_reads, big.mark=","),
      " (", format(global_total, big.mark=","), ") | RPM: ", round(rpm_total, 2)
    )]
    sum_label <- paste(gene_sum_total$label, collapse = "\n")
    
    p_sum <- ggplot(gene_sum_agg, aes(x=psite, y=rpm_sum, fill=group_name)) +
      geom_col(width=0.8, alpha=0.8, position="identity") +
      geom_vline(xintercept=longest_tx_atg, color="#e84a5f", linewidth=1.5) +
      annotate("text", x=longest_tx_atg, y=max(gene_sum_agg$rpm_sum)*1.1,
               label=paste0("ATG: ", longest_tx_atg), color="#e84a5f", size=4.5, fontface="bold") +
      geom_vline(xintercept=longest_tx_atg + longest_tx_cds_len - 1, color="#2A9D8F", linewidth=1, linetype="dashed") +
      annotate("text", x=Inf, y=Inf, label=sum_label, hjust=1.1, vjust=1.1, size=3.8, fontface="bold") +
      labs(x="PSite", y="RPM Sum (only / group total reads)", title=paste0(gene, " | Sum RPM"), fill="Group") +
      scale_fill_brewer(palette="Set1") + theme_bw()
    
    f_sum <- paste0("gene_", gene, "_sum_rpm.pdf")
    ggsave(file.path(gene_dir, f_sum), p_sum, width=16, height=7)
    total_sum_plots <- total_sum_plots + 1
    
    # ===========================================================================
    # 3) Gene MEAN plot (unchanged)
    # ===========================================================================
    gene_mean_dt <- gene_all_dt[, .(raw_sum = .N), by = .(group_name, psite, frame)]
    n_tx <- length(unique(gene_all_dt$transcript))
    gene_mean_dt[, raw_mean := raw_sum / n_tx]
    gene_mean_dt <- merge(gene_mean_dt, global_group_total, by = "group_name")
    gene_mean_dt[, rpm_mean := (raw_mean / global_total) * 1e6]
    
    gene_mean_agg <- gene_mean_dt[, .(raw_mean = sum(raw_mean), rpm_mean = sum(rpm_mean)), by = .(group_name, psite)]
    
    gene_mean_total <- gene_mean_dt[, .(gene_reads = sum(raw_sum)), by = group_name]
    gene_mean_total <- merge(gene_mean_total, global_group_total, by = "group_name")
    gene_mean_total[, rpm_total := (gene_reads / global_total) * 1e6]
    gene_mean_total[, label := paste0(
      group_name, ": ", format(gene_reads, big.mark=","),
      " (", format(global_total, big.mark=","), ") | RPM: ", round(rpm_total, 2)
    )]
    mean_label <- paste(gene_mean_total$label, collapse = "\n")
    
    p_mean <- ggplot(gene_mean_agg, aes(x=psite, y=rpm_mean, fill=group_name)) +
      geom_col(width=0.8, alpha=0.8, position="identity") +
      geom_vline(xintercept=longest_tx_atg, color="#e84a5f", linewidth=1.5) +
      annotate("text", x=longest_tx_atg, y=max(gene_mean_agg$rpm_mean)*1.1,
               label=paste0("ATG: ", longest_tx_atg), color="#e84a5f", size=4.5, fontface="bold") +
      geom_vline(xintercept=longest_tx_atg + longest_tx_cds_len - 1, color="#2A9D8F", linewidth=1, linetype="dashed") +
      annotate("text", x=Inf, y=Inf, label=mean_label, hjust=1.1, vjust=1.1, size=3.8, fontface="bold") +
      labs(x="PSite", y="RPM Mean (only / group total reads)", title=paste0(gene, " | Mean RPM"), fill="Group") +
      scale_fill_brewer(palette="Set1") + theme_bw()
    
    f_mean <- paste0("gene_", gene, "_mean_rpm.pdf")
    ggsave(file.path(gene_dir, f_mean), p_mean, width=16, height=7)
    total_mean_plots <- total_mean_plots + 1
    
    fwrite(gene_sum[, .(group_name, psite, frame, raw_sum, rpm_sum)],
           file.path(gene_dir, paste0("gene_", gene, "_sum_rpm_table.txt")), sep = "\t")
    fwrite(gene_mean_dt[, .(group_name, psite, frame, raw_sum, raw_mean, rpm_mean)],
           file.path(gene_dir, paste0("gene_", gene, "_mean_rpm_table.txt")), sep = "\t")
  }
  
  cat("\n✅ All completed\n")
  cat("Single transcript RPM:", total_tx_plots, "\n")
  cat("Gene SUM  RPM:", total_sum_plots, "\n")
  cat("Gene MEAN RPM:", total_mean_plots, "\n")
}
# ===================== ✅ Fully fixed: ratio histogram non-numeric operation error =====================
plot_APE_codon_ratio_hist <- function(all_codon_occupancy_group, diff_result, output_dir) {
  cat("7. Generate A/P/E Codon Occupancy Ratio Histograms\n")
  cat("------------------------------------------------------------\n")
  group_names <- unique(all_codon_occupancy_group$group_name)
  codon_types <- unique(all_codon_occupancy_group$codon_type)
  
  if(length(group_names) !=2){
    cat("⚠️  Only support two-group ratio analysis! Skip this step.\n")
    return()
  }
  
  for (ctype in codon_types) {
    cat(paste0("   Plotting: ", ctype, "\n"))
    group_sub <- all_codon_occupancy_group[all_codon_occupancy_group$codon_type == ctype, ]
    
    # Cast to wide format
    ratio_dt <- dcast(group_sub, codon ~ group_name, value.var = c("mean_occupancy","se_occupancy"))
    colnames(ratio_dt) <- gsub("mean_occupancy_","",colnames(ratio_dt))
    colnames(ratio_dt) <- gsub("se_occupancy_","se_",colnames(ratio_dt))
    
    # Force all columns to numeric
    for(col in colnames(ratio_dt)[-1]){
      ratio_dt[[col]] <- as.numeric(ratio_dt[[col]])
    }
    
    # Column index
    mean_col1 <- which(colnames(ratio_dt) == group_names[1])
    mean_col2 <- which(colnames(ratio_dt) == group_names[2])
    se_col1 <- which(colnames(ratio_dt) == paste0("se_",group_names[1]))
    se_col2 <- which(colnames(ratio_dt) == paste0("se_",group_names[2]))
    
    # Calculate ratio and error propagation
    ratio_dt$codon_ratio <- ratio_dt[[mean_col2]] / ratio_dt[[mean_col1]]
    ratio_dt$se_ratio <- ratio_dt$codon_ratio * sqrt( (ratio_dt[[se_col2]]/ratio_dt[[mean_col2]])^2 + (ratio_dt[[se_col1]]/ratio_dt[[mean_col1]])^2 )
    ratio_dt$expr_status <- ifelse(ratio_dt$codon_ratio > 1, "up", ifelse(ratio_dt$codon_ratio < 1, "down", "non_diff"))
    
    # Merge AA annotation
    ratio_dt <- merge(ratio_dt, codon2aa, by = "codon", all.x = TRUE)
    ratio_dt$aa <- ifelse(is.na(ratio_dt$aa), "*", ratio_dt$aa)
    ratio_dt$codon_aa <- paste0(ratio_dt$codon, " (", ratio_dt$aa, ")")
    
    # Differential status
    if(!is.null(diff_result) && nrow(diff_result)>0){
      diff_sub <- diff_result[diff_result$codon_type == ctype, ]
      ratio_dt <- merge(ratio_dt, diff_sub[, c("codon","diff_type")], by = "codon", all.x = T)
      ratio_dt$diff_type[is.na(ratio_dt$diff_type)] <- ratio_dt$expr_status[is.na(ratio_dt$diff_type)]
    } else {
      ratio_dt$diff_type <- ratio_dt$expr_status
    }
    
    # ================== Original plot: sorted by ratio ==================
    ratio_dt <- ratio_dt[order(-ratio_dt$codon_ratio), ]
    ratio_dt$codon_aa <- factor(ratio_dt$codon_aa, levels = ratio_dt$codon_aa)
    
    # ✅ Minor fix 1: codons ending with CG colored red
    text_color <- ifelse(grepl("[CG]$", as.character(ratio_dt$codon)), "red", "black")
    
    p <- ggplot(ratio_dt, aes(x = codon_aa, y = codon_ratio, fill = diff_type)) +
      geom_col(alpha = 0.8, width = 0.8) +
      geom_errorbar(aes(ymin = codon_ratio - se_ratio, ymax = codon_ratio + se_ratio), width = 0.3, color = "black") +
      geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 1.2, alpha=0.8) +
      scale_fill_manual(values = c("up" = "#E63946", "down" = "#457B9D", "non_diff" = "#95A5A6"),
                        labels = c("Up-regulated (ratio>1)", "Down-regulated (ratio<1)", "Not significant (ratio=1)"),
                        breaks = c("up", "down", "non_diff") ) +
      labs(x = paste0(ctype, " - Codon (Amino acid) Sorted by Ratio"),
           y = paste0("Codon Occupancy Ratio (", group_names[2], "/", group_names[1], ")"),
           fill = "Expression Change",
           title = paste0(ctype, " - Codon Occupancy Ratio Histogram")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9, color = text_color))
    
    ggsave(file.path(output_dir, paste0("4.", ctype, "_codon_ratio_hist.pdf")), p, width=18, height=7, dpi=300)
    
    # ================== Synonymous codon histogram ==================
    syn_dt <- ratio_dt[order(ratio_dt$aa, -ratio_dt$codon_ratio), ]
    syn_dt$codon <- factor(syn_dt$codon, levels = unique(syn_dt$codon))
    
    # ✅ Minor fix 2: replicate spacing logic
    syn_dt$x_pos <- as.integer(syn_dt$codon)
    syn_dt$x_pos <- syn_dt$x_pos + c(0, cumsum(diff(as.integer(factor(syn_dt$aa))) != 0)) * 0.8
    
    # ✅ Minor fix 3: CG ending red
    syn_color <- ifelse(grepl("[CG]$", as.character(syn_dt$codon)), "red", "black")
    
    p_syn <- ggplot(syn_dt, aes(x = x_pos, y = codon_ratio, fill = diff_type)) +
      geom_col(alpha = 0.8, width = 0.8) +
      geom_errorbar(aes(ymin = codon_ratio - se_ratio, ymax = codon_ratio + se_ratio),
                    width = 0.3, color = "black") +
      geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 1.2, alpha=0.8) +
      scale_fill_manual(values = c("up" = "#E63946", "down" = "#457B9D", "non_diff" = "#95A5A6")) +
      # ✅ Minor fix 4: continuous axis
      scale_x_continuous(breaks = syn_dt$x_pos, labels = syn_dt$codon_aa, expand = c(0.01, 0.01)) +
      labs(x = paste0(ctype, " - Codon (Grouped by Synonym)"),
           y = paste0("Codon Occupancy Ratio (", group_names[2], "/", group_names[1], ")"),
           fill = "Expression Change",
           title = paste0(ctype, " - Synonymous Codon Occupancy Ratio")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9, color = syn_color),
            panel.grid = element_blank())
    
    ggsave(file.path(output_dir, paste0("4.", ctype, "_codon_ratio_synonymous.pdf")),
           p_syn, width=20, height=7, dpi=300)
  }
  cat("✅ A/P/E ratio histograms generated successfully!\n\n")
# ===================== Main Pipeline Execution =====================
cat("9. Filter Codon Data (retain codons present in all samples)\n")
cat("------------------------------------------------------------\n")
codon_filter <- all_codon_occupancy_raw[, list(n_sample = length(unique(sample_name))), by = list(codon_type, codon)]
codon_keep <- codon_filter[codon_filter$n_sample == length(sample_names), ]
all_codon_occupancy_raw_filt <- all_codon_occupancy_raw[all_codon_occupancy_raw$codon_type %in% codon_keep$codon_type & all_codon_occupancy_raw$codon %in% codon_keep$codon, ]
all_codon_occupancy_group_filt <- all_codon_occupancy_group[all_codon_occupancy_group$codon_type %in% codon_keep$codon_type & all_codon_occupancy_group$codon %in% codon_keep$codon, ]

cat(paste0("✅ Codon filter completed: ", nrow(codon_keep), " codons retained (A/P/E sites)\n\n"))

# ✅ Minimal modification: auto-skip PCA for 1vs1 (2 samples) without error
if (length(sample_names) >= 3) {
  plot_Psite_PCA(all_codon_occupancy_raw_filt, group_dt, output_root)
} else {
  cat("✅ Sample count =", length(sample_names), "(1vs1 comparison) → Auto skip PCA analysis\n\n")
}

plot_group_footprint(all_psite_raw_dt, output_root)

# 3. Call function (arguments: footprint data + target gene list path + tx2gene + output directory)
plot_gene_specific_footprint(all_psite_raw_dt, gene_list_path, tx2gene, "./5.ribo_group_diff_output")
plot_APE_codon_ratio_hist(all_codon_occupancy_group_filt, NULL, output_root)

# ===================== Global unified · synonymous codons share base color + codon-AA mapping =====================
# codon → amino acid mapping table (including stop codon *)
get_codon_aa_map <- function() {
  codon_aa <- list(
    # Phenylalanine Phe (F)
    "TTT" = "F", "TTC" = "F",
    # Leucine Leu (L)
    "TTA" = "L", "TTG" = "L", "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    # Isoleucine Ile (I)
    "ATT" = "I", "ATC" = "I", "ATA" = "I",
    # Methionine Met (M) / start
    "ATG" = "M",
    # Valine Val (V)
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    # Serine Ser (S)
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S", "AGT" = "S", "AGC" = "S",
    # Proline Pro (P)
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    # Threonine Thr (T)
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    # Alanine Ala (A)
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    # Tyrosine Tyr (Y)
    "TAT" = "Y", "TAC" = "Y",
    # Stop codon (*)
    "TAA" = "*", "TAG" = "*", "TGA" = "*",
    # Histidine His (H)
    "CAT" = "H", "CAC" = "H",
    # Glutamine Gln (Q)
    "CAA" = "Q", "CAG" = "Q",
    # Asparagine Asn (N)
    "AAT" = "N", "AAC" = "N",
    # Lysine Lys (K)
    "AAA" = "K", "AAG" = "K",
    # Aspartic acid Asp (D)
    "GAT" = "D", "GAC" = "D",
    # Glutamic acid Glu (E)
    "GAA" = "E", "GAG" = "E",
    # Cysteine Cys (C)
    "TGT" = "C", "TGC" = "C",
    # Tryptophan Trp (W)
    "TGG" = "W",
    # Arginine Arg (R)
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R", "AGA" = "R", "AGG" = "R",
    # Glycine Gly (G)
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  return(codon_aa)
}

# ========== Fix NaN error: synonymous codons share base hue (base + lightness variation) ==========
get_codon_color <- function(codon_aa_map) {
  # Base color for amino acids (20 AAs + stop)
  aa_main_color <- list(
    "F" = "#E41A1C", "L" = "#F28E8E", "I" = "#FAAA3A", "M" = "#F7E14F",
    "V" = "#A65628", "S" = "#C49A6C", "P" = "#FF7F00", "T" = "#FDB462",
    "A" = "#377EB8", "Y" = "#75B7E6", "*" = "#984EA3", "H" = "#C17BDB",
    "Q" = "#4DAF4A", "N" = "#92CD76", "K" = "#1B9E77", "D" = "#66C2A5",
    "E" = "#8DD3C7", "C" = "#FFFFB3", "W" = "#BEBADA", "R" = "#FB8072",
    "G" = "#80B1D3"
  )
  
  # ========== Core fix: avoid NaN color error (boundary check + fallback) ==========
  get_shade_color <- function(base_color, n, idx) {
    # Safe RGB parsing with error handling
    rgb_vals <- tryCatch({
      col2rgb(base_color) / 255
    }, error = function(e) {
      # Fallback gray if parsing fails
      return(matrix(c(0.5, 0.5, 0.5), nrow = 3))
    })
    
    # Handle single codon case (avoid division by zero)
    if (n <= 1) {
      return(base_color)
    }
    
    # Lightness scaling factor (0.8–1.2, prevents NaN)
    lightness <- 0.8 + (idx-1)*(0.4/(n-1))
    # Clamp to [0,1] (critical NaN guard)
    new_rgb <- pmin(pmax(rgb_vals * lightness, 0), 1)
    
    # Replace NaN / Inf with neutral gray
    new_rgb[is.na(new_rgb) | is.infinite(new_rgb)] <- 0.5
    
    return(rgb(new_rgb[1], new_rgb[2], new_rgb[3]))
  }
  
  # Generate per-codon colors grouped by amino acid
  codon_colors <- c()
  aa_codons <- split(names(codon_aa_map), unlist(codon_aa_map))
  
  for (aa in names(aa_codons)) {
    codons <- aa_codons[[aa]]
    n <- length(codons)
    base_color <- aa_main_color[[aa]]
    
    # Fallback gray if base color missing
    if (is.null(base_color)) {
      base_color <- "#888888"
    }
    
    # Generate lighter/darker variants for synonymous codons
    for (i in 1:n) {
      codon_colors[codons[i]] <- get_shade_color(base_color, n, i)
    }
  }
  
  return(codon_colors)
}

# ===================== New: amino acid chemical property classification (no redundancy) =====================
# AA → chemical class mapping
get_aa_to_class <- function() {
  aa_to_class <- list(
    # acidic
    "D" = "acidic", "E" = "acidic",
    # basic
    "K" = "basic", "R" = "basic", "H" = "basic",
    # polar
    "N" = "polar", "Q" = "polar", "S" = "polar", "T" = "polar", "C" = "polar",
    # aromatic
    "F" = "aromatic", "W" = "aromatic", "Y" = "aromatic",
    # aliphatic
    "G" = "aliphatic", "A" = "aliphatic", "V" = "aliphatic", "L" = "aliphatic",
    "I" = "aliphatic", "M" = "aliphatic", "P" = "aliphatic",
    # stop codon
    "*" = "stop"
  )
  return(aa_to_class)
}

# chemical class → color mapping (consistent with legend)
get_class_color_map <- function() {
  class_color <- list(
    "acidic"   = "#E41A1C",   # red
    "aliphatic"= "#000000",   # black
    "aromatic" = "#4DAF4A",   # green
    "basic"    = "#377EB8",   # blue
    "polar"    = "#FF7F00",   # orange
    "stop"     = "#984EA3"    # purple
  )
  return(class_color)
}

# Initialize all mapping tables
codon_aa_map <- get_codon_aa_map()
codon_colors <- get_codon_color(codon_aa_map)
aa_to_class <- get_aa_to_class()
class_color_map <- get_class_color_map()

# ===================== Core analysis function (minimal addition: chemical class coloring + retain AA abbreviation) =====================
cat("10. Analyze All Codon Motif (synonymous codon same hue + per-codon folder + narrow rectangles + chemical class coloring)\n")
cat("------------------------------------------------------------\n")

analyze_codon_motif <- function(sample_dir, transcriptome_fasta, output_dir, 
                                gene_list_path, tx2gene, motif_type = c("5codon", "20codon")) {
  # Load required packages
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("❌ Biostrings package required. Install first: install.packages('Biostrings')")
  }
  library(Biostrings)
  library(data.table)
  
  # 1. Load target gene list
  target_genes <- character()
  if (file.exists(gene_list_path)) {
    target_genes <- readLines(gene_list_path, warn = FALSE)
    target_genes <- gsub("\\s+", "", target_genes)
    target_genes <- unique(target_genes[target_genes != ""])
    cat(paste0("✅ Target gene list loaded: ", length(target_genes), " genes → ", paste(target_genes, collapse = ", "), "\n"))
    
    target_tx_all <- character()
    for (gene in target_genes) {
      tx_for_gene <- tx2gene[gene_name == gene, unique(transcript)]
      cat(paste0("   - ", gene, " → ", length(tx_for_gene), " transcripts: ", paste(head(tx_for_gene, 5), collapse = ", "), ifelse(length(tx_for_gene)>5, "...", ""), "\n"))
      target_tx_all <- c(target_tx_all, tx_for_gene)
    }
    target_tx_all <- unique(target_tx_all)
    cat(paste0("\n✅ Total unique transcripts matched from target genes: ", length(target_tx_all), "\n"))
  } else {
    cat("⚠️ Target gene list file not found: ", gene_list_path, "\n")
    target_tx_all <- character(0)
  }
  
  # 2. Load RPF data
  psite_path <- file.path(sample_dir, "reads_psite.rds")
  if (!file.exists(psite_path)) {
    psite_path <- file.path(dirname(sample_dir), "reads_psite.rds")
  }
  cat(paste0("\n✅ Loading RPF data: ", psite_path, "\n"))
  
  psite_data <- readRDS(psite_path)
  if (is(psite_data, "list")) {
    for (elem in psite_data) {
      if (is(elem, "data.frame") || is(elem, "data.table")) {
        psite_dt <- as.data.table(elem)
        break
      }
    }
    if (!exists("psite_dt")) stop("❌ Valid data.frame not found inside RDS list!")
  } else {
    psite_dt <- as.data.table(psite_data)
  }
  
  # Check required columns
  required_cols <- c("cds_start", "cds_stop", "psite", "P_codon", "transcript")
  missing_cols <- setdiff(required_cols, colnames(psite_dt))
  if (length(missing_cols) > 0) {
    stop(paste0("❌ Missing required columns: ", paste(missing_cols, collapse = ", ")))
  }
  
  # Validate target transcripts against RPF
  rpf_tx <- unique(psite_dt$transcript)
  target_tx_in_rpf <- target_tx_all[target_tx_all %in% rpf_tx]
  target_tx_not_in_rpf <- setdiff(target_tx_all, rpf_tx_all)
  
  cat(paste0("✅ Target transcripts present in RPF data: ", length(target_tx_in_rpf), "\n"))
  if (length(target_tx_not_in_rpf) > 0) {
    cat(paste0("❌ Target transcripts absent from RPF data: ", length(target_tx_not_in_rpf), " → examples: ", paste(head(target_tx_not_in_rpf, 3), collapse = ", "), "\n"))
  }
  target_tx <- target_tx_in_rpf
  
  # 3. Load transcript sequences
  cat("✅ Reading transcriptome FASTA (preserve original header)...\n")
  tx_seqs <- readDNAStringSet(transcriptome_fasta)
  tx_names_full <- names(tx_seqs)
  tx_names_clean <- sapply(strsplit(tx_names_full, "\\|"), `[`, 1)
  
  # Parse CDS coordinates from FASTA header
  parse_fasta_header <- function(header) {
    cds_match <- regmatches(header, regexpr("CDS:\\d+-\\d+", header))
    if (length(cds_match) == 0) {
      return(list(cds_start = NA, cds_stop = NA))
    }
    cds_range <- gsub("CDS:", "", cds_match)
    cds_start <- as.integer(strsplit(cds_range, "-")[[1]][1])
    cds_stop <- as.integer(strsplit(cds_range, "-")[[1]][2])
    return(list(cds_start = cds_start, cds_stop = cds_stop))
  }
  
  # Build sequence metadata table
  seq_meta <- lapply(tx_names_full, parse_fasta_header)
  seq_map <- data.table(
    transcript = tx_names_clean,
    seq = as.character(tx_seqs),
    seq_len = width(tx_seqs),
    fasta_cds_start = sapply(seq_meta, function(x) x$cds_start),
    fasta_cds_stop = sapply(seq_meta, function(x) x$cds_stop)
  )
  setkey(seq_map, transcript)
  
  # Print FASTA parsing preview
  cat(paste0("\n📌 FASTA header parsing preview (top 5):\n"))
  for (i in 1:min(5, nrow(seq_map))) {
    cat(paste0("   - ", seq_map$transcript[i], 
               " | FASTA CDS: ", seq_map$fasta_cds_start[i], "-", seq_map$fasta_cds_stop[i],
               " | seq length: ", seq_map$seq_len[i], "\n"))
  }
  
  # Validate target transcripts against FASTA
  target_tx_in_seq <- target_tx[target_tx %in% seq_map$transcript]
  target_tx_not_in_seq <- setdiff(target_tx, seq_map$transcript)
  if (length(target_tx_not_in_seq) > 0) {
    cat(paste0("\n❌ Target transcripts absent from FASTA: ", length(target_tx_not_in_seq), " → examples: ", paste(head(target_tx_not_in_seq, 3), collapse = ", "), "\n"))
  }
  target_tx <- target_tx_in_seq
  cat(paste0("✅ Final valid target transcripts (present in both RPF + FASTA): ", length(target_tx), "\n"))
  
  # 4. Configure parameters for different motif modes
  motif_config <- list()
  if ("5codon" %in% motif_type) {
    motif_config[["5codon"]] <- list(
      flank_up = 6,
      flank_down = 8,
      total_len = 15,
      nt_labels = -6:8,
      codon_labels = -2:2,
      site_pos = 7,
      suffix = "5codons"
    )
  }
  if ("20codon" %in% motif_type) {
    motif_config[["20codon"]] <- list(
      flank_up = 60,    
      flank_down = 2,   
      total_len = 62,   
      nt_labels = -60:-1,
      codon_labels = -20:-1,
      site_pos = 61,    
      suffix = "20codons_upstream"
    )
  }
  
  # ========== Minimal addition: stacked plot colored by chemical class, retain AA abbreviation ==========
  plot_aa_class_stacked <- function(sorted_freq_list, positions, config, codon_output_dir, seq_type, codon) {
    cat(paste0("   🧬 Generating AA frequency plot colored by chemical class (AA labels retained)...\n"))
    
    # 1. Output paths
    plot_base_name <- paste0(codon, "_", config$suffix, "_", seq_type, "_AA_class_color")
    png_file <- file.path(codon_output_dir, "plots", paste0(plot_base_name, ".png"))
    pdf_file <- file.path(codon_output_dir, "plots", paste0(plot_base_name, ".pdf"))
    table_file <- file.path(codon_output_dir, "tables", paste0(plot_base_name, "_freq_table.csv"))
    
    # Auto canvas width, reserve space for legend on right
    plot_width <- ifelse(config$total_len > 30, 22, 14)
    plot_height <- 8
    
    # 2. Build frequency table, force position as numeric
    aa_freq_table <- data.table()
    for (pos in positions) {
      if (!pos %in% names(sorted_freq_list) || length(sorted_freq_list[[pos]]) == 0) next
      aa_names <- sapply(names(sorted_freq_list[[pos]]), function(c) {
        aa <- codon_aa_map[[c]]
        if (is.null(aa)) "*" else aa
      })
      aa_classes <- sapply(aa_names, function(aa) {
        cls <- aa_to_class[[aa]]
        if (is.null(cls)) "unknown" else cls
      })
      temp_dt <- data.table(
        position = as.numeric(pos),
        codon = names(sorted_freq_list[[pos]]),
        aa = aa_names,
        frequency = sorted_freq_list[[pos]],
        class = aa_classes
      )
      aa_freq_table <- rbind(aa_freq_table, temp_dt)
    }
    
    if (nrow(aa_freq_table) == 0) {
      cat("   ⚠️ No valid AA frequency data, skip plot\n")
      return(NULL)
    }
    aa_freq_table <- aa_freq_table[class != "unknown"]
    aa_freq_table <- aa_freq_table[order(position)]
    if (nrow(aa_freq_table) == 0) {
      cat("   ⚠️ No valid chemical class data, skip plot\n")
      return(NULL)
    }
    
    write.csv(aa_freq_table, table_file, row.names = FALSE)
    cat(paste0("✅ AA frequency table saved: ", table_file, "\n"))
    
    # 3. Prepare plotting data
    positions_unique <- unique(aa_freq_table$position)
    positions_unique <- sort(positions_unique)
    classes_unique <- unique(aa_freq_table$class)
    
    # 4. Draw PNG
    png(png_file, width = plot_width, height = plot_height, units = "in", res = 300)
    par(mar = c(10, 6, 4, 10) + 0.1)
    
    pos_totals <- aa_freq_table[, .(total = sum(frequency)), by = position]$total
    y_max <- ifelse(length(pos_totals) > 0, max(pos_totals) * 1.1, 100)
    
    plot(0, 0, type = "n",
         xlim = c(0.5, length(positions_unique) + 0.5),
         ylim = c(0, y_max),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ", config$suffix, " (AA by Chemical Class)"),
         xaxt = "n", yaxt = "s",
         cex.lab = 1.1, cex.axis = 1.0)
    
    x_coords <- 1:length(positions_unique)
    names(x_coords) <- positions_unique
    
    for (pos in positions_unique) {
      pos_data <- aa_freq_table[position == pos][order(-frequency)]
      if (nrow(pos_data) == 0) next
      y_bottom <- 0
      for (i in 1:nrow(pos_data)) {
        aa <- pos_data$aa[i]
        freq <- pos_data$frequency[i]
        cls <- pos_data$class[i]
        color <- class_color_map[[cls]]
        if (is.null(color)) color <- "#888888"
        
        rect(x_coords[as.character(pos)] - 0.15, y_bottom,
             x_coords[as.character(pos)] + 0.15, y_bottom + freq,
             col = color, border = "white", lwd = 0.5)
        
        text_size <- ifelse(freq < 5, 0.6, 0.8)
        text(x_coords[as.character(pos)], y_bottom + freq/2,
             labels = aa,
             cex = text_size, font = 2, col = "white")
        
        y_bottom <- y_bottom + freq
      }
    }
    
    axis(1, at = x_coords, labels = positions_unique, las = 1, cex.axis = 0.9,
         mgp = c(3, 0.5, 0))
    mtext("Position relative to P-site (codon position)", side = 1, line = 7, cex = 1.0)
    
    legend(x = length(positions_unique) + 0.8,
           y = y_max/2,
           legend = names(class_color_map),
           fill = unlist(class_color_map),
           title = "Chemical Class",
           cex = 0.9,
           bty = "n",
           xpd = TRUE,
           title.cex = 1.0)
    
    dev.off()
    cat(paste0("✅ Chemical-class colored PNG saved: ", png_file, "\n"))
    
    # 5. Draw PDF
    pdf(pdf_file, width = plot_width, height = plot_height)
    par(mar = c(10, 6, 4, 10) + 0.1)
    
    plot(0, 0, type = "n",
         xlim = c(0.5, length(positions_unique) + 0.5),
         ylim = c(0, y_max),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ", config$suffix, " (AA by Chemical Class)"),
         xaxt = "n", yaxt = "s",
         cex.lab = 1.1, cex.axis = 1.0)
    
    for (pos in positions_unique) {
      pos_data <- aa_freq_table[position == pos][order(-frequency)]
      if (nrow(pos_data) == 0) next
      y_bottom <- 0
      for (i in 1:nrow(pos_data)) {
        aa <- pos_data$aa[i]
        freq <- pos_data$frequency[i]
        cls <- pos_data$class[i]
        color <- class_color_map[[cls]]
        if (is.null(color)) color <- "#888888"
        
        rect(x_coords[as.character(pos)] - 0.15, y_bottom,
             x_coords[as.character(pos)] + 0.15, y_bottom + freq,
             col = color, border = "white", lwd = 0.5)
        
        text_size <- ifelse(freq < 5, 0.6, 0.8)
        text(x_coords[as.character(pos)], y_bottom + freq/2,
             labels = aa,
             cex = text_size, font = 2, col = "white")
        
        y_bottom <- y_bottom + freq
      }
    }
    
    axis(1, at = x_coords, labels = positions_unique, las = 1, cex.axis = 0.9,
         mgp = c(3, 0.5, 0))
    mtext("Position relative to P-site (codon position)", side = 1, line = 7, cex = 1.0)
    
    legend(x = length(positions_unique) + 0.8,
           y = y_max/2,
           legend = names(class_color_map),
           fill = unlist(class_color_map),
           title = "Chemical Class",
           cex = 0.9,
           bty = "n",
           xpd = TRUE,
           title.cex = 1.0)
    
    dev.off()
    cat(paste0("✅ Chemical-class colored PDF saved: ", pdf_file, "\n"))
    
    return(TRUE)
  }
  
  # ========== Pure R manual stacked plot function; call AA-class plot only for codon mode ==========
  plot_stacked_manual <- function(freq_table, type = c("nt", "codon"), seq_type, codon, config, codon_output_dir) {
    if (is.null(freq_table) || nrow(freq_table) == 0) {
      cat(paste0("   ⚠️ No ", type, " data, skip plot\n"))
      return(NULL)
    }
    
    freq_table$position <- as.character(freq_table$position)
    freq_table <- freq_table[!is.na(position) & position != "" & frequency > 0, ]
    if (nrow(freq_table) == 0) {
      cat(paste0("   ⚠️ No valid ", type, " frequency data, skip plot\n"))
      return(NULL)
    }
    
    if (type == "nt") {
      col_name <- "nt"
      color_map <- c("A" = "#4DAF4A", "T" = "#E41A1C", "C" = "#377EB8", "G" = "#FFFF33")
      text_size <- 0.8
      if (config$suffix == "5codons") {
        plot_width <- 12
        plot_height <- 6
      } else {
        plot_width <- ifelse(config$total_len>30, 15, 10)
        plot_height <- 6
      }
      rect_width <- 0.4
    } else {
      col_name <- "codon"
      color_map <- codon_colors
      text_size <- 0.5
      if (config$suffix == "5codons") {
        plot_width <- 14
        plot_height <- 7
      } else {
        plot_width <- ifelse(config$total_len>30, 20, 10)
        plot_height <- 7
      }
      rect_width <- 0.2
    }
    
    positions <- unique(freq_table$position)
    positions <- positions[order(as.numeric(positions))]
    elements <- unique(freq_table[[col_name]])
    
    freq_matrix <- matrix(0, nrow = length(positions), ncol = length(elements))
    rownames(freq_matrix) <- positions
    colnames(freq_matrix) <- elements
    
    for (i in 1:nrow(freq_table)) {
      pos <- as.character(freq_table$position[i])
      elem <- as.character(freq_table[[col_name]][i])
      if (pos %in% rownames(freq_matrix) && elem %in% colnames(freq_matrix)) {
        freq_matrix[pos, elem] <- freq_table$frequency[i]
      }
    }
    
    sorted_freq_list <- list()
    sorted_color_list <- list()
    for (pos in positions) {
      pos_freq <- freq_matrix[pos, ]
      pos_freq <- pos_freq[pos_freq > 0]
      if (length(pos_freq) == 0) next
      pos_freq_sorted <- sort(pos_freq, decreasing = TRUE)
      sorted_freq_list[[pos]] <- pos_freq_sorted
      sorted_color_list[[pos]] <- color_map[names(pos_freq_sorted)]
    }
    
    if (length(sorted_freq_list) == 0) {
      cat(paste0("   ⚠️ No valid ", type, " stacked data, skip plot\n"))
      return(NULL)
    }
    
    plot_base_name <- paste0(codon, "_", config$suffix, "_", seq_type, "_", type, "_manual")
    png_file <- file.path(codon_output_dir, "plots", paste0(plot_base_name, ".png"))
    pdf_file <- file.path(codon_output_dir, "plots", paste0(plot_base_name, ".pdf"))
    table_file <- file.path(codon_output_dir, "tables", paste0(plot_base_name, "_freq_table.csv"))
    
    write.csv(freq_table, table_file, row.names = FALSE)
    cat(paste0("✅ Frequency table saved: ", table_file, "\n"))
    
    # PNG
    png(png_file, width = plot_width, height = plot_height, units = "in", res = 300)
    par(mar = c(8, 5, 4, 4) + 0.1)
    x_left <- 0.5
    x_right <- length(positions) + 0.5
    plot(0, 0, type = "n", 
         xlim = c(x_left, x_right),
         ylim = c(0, 100),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ", config$suffix, " (Frequency Descending)"),
         xaxt = "n", yaxt = "s")
    x_pos <- 1:length(positions)
    names(x_pos) <- positions
    y_bottom <- rep(0, length(positions))
    names(y_bottom) <- positions
    
    for (pos in positions) {
      pos_idx <- which(names(x_pos) == pos)
      if (length(pos_idx) == 0 || !pos %in% names(sorted_freq_list)) next
      pos_freq <- sorted_freq_list[[pos]]
      pos_colors <- sorted_color_list[[pos]]
      if (length(pos_freq) == 0) next
      for (i in 1:length(pos_freq)) {
        freq <- pos_freq[i]
        color <- pos_colors[i]
        elem_name <- names(pos_freq)[i]
        rect(x_pos[pos_idx] - rect_width, y_bottom[pos_idx],
             x_pos[pos_idx] + rect_width, y_bottom[pos_idx] + freq,
             col = color, border = "white")
        if (type == "codon") {
          aa <- ifelse(elem_name %in% names(codon_aa_map), codon_aa_map[[elem_name]], "")
          display_text <- ifelse(aa != "", paste0(elem_name, "(", aa, ")"), elem_name)
        } else {
          display_text <- elem_name
        }
        text(x_pos[pos_idx], y_bottom[pos_idx] + freq/2, 
             labels = display_text, 
             cex = text_size, font = 2, col = "white")
        y_bottom[pos_idx] <- y_bottom[pos_idx] + freq
      }
    }
    axis(1, at = x_pos, labels = positions, las = 2, cex.axis = 0.7)
    mtext(ifelse(type=="nt", "Position relative to P-site (nt)", "Position relative to P-site (codon)"), 
          side = 1, line = 6, cex = 0.9)
    dev.off()
    cat(paste0("✅ PNG saved: ", png_file, "\n"))
    
    # PDF
    pdf(pdf_file, width = plot_width, height = plot_height)
    par(mar = c(8, 5, 4, 4) + 0.1)
    x_left <- 0.5
    x_right <- length(positions) + 0.5
    plot(0, 0, type = "n", 
         xlim = c(x_left, x_right),
         ylim = c(0, 100),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ", config$suffix, " (Frequency Descending)"),
         xaxt = "n", yaxt = "s")
    y_bottom <- rep(0, length(positions))
    names(y_bottom) <- positions
    for (pos in positions) {
      pos_idx <- which(names(x_pos) == pos)
      if (length(pos_idx) == 0 || !pos %in% names(sorted_freq_list)) next
      pos_freq <- sorted_freq_list[[pos]]
      pos_colors <- sorted_color_list[[pos]]
      if (length(pos_freq) == 0) next
      for (i in 1:length(pos_freq)) {
        freq <- pos_freq[i]
        color <- pos_colors[i]
        elem_name <- names(pos_freq)[i]
        rect(x_pos[pos_idx] - rect_width, y_bottom[pos_idx],
             x_pos[pos_idx] + rect_width, y_bottom[pos_idx] + freq,
             col = color, border = "white")
        if (type == "codon") {
          aa <- ifelse(elem_name %in% names(codon_aa_map), codon_aa_map[[elem_name]], "")
          display_text <- ifelse(aa != "", paste0(elem_name, "(", aa, ")"), elem_name)
        } else {
          display_text <- elem_name
        }
        text(x_pos[pos_idx], y_bottom[pos_idx] + freq/2, 
             labels = display_text, 
             cex = text_size, font = 2, col = "white")
        y_bottom[pos_idx] <- y_bottom[pos_idx] + freq
      }
    }
    axis(1, at = x_pos, labels = positions, las = 2, cex.axis = 0.7)
    mtext(ifelse(type=="nt", "Position relative to P-site (nt)", "Position relative to P-site (codon)"), 
          side = 1, line = 6, cex = 0.9)
    dev.off()
    cat(paste0("✅ PDF saved: ", pdf_file, "\n"))
    
    if (type == "codon") {
      plot_aa_class_stacked(sorted_freq_list, positions, config, codon_output_dir, seq_type, codon)
    }
    
    return(TRUE)
  }
  
  # ========== Core optimization 2: use all unique codons (not only top3) ==========
  all_codons <- unique(psite_dt[!is.na(P_codon), P_codon])
  all_codons <- sort(all_codons)
  n_codons <- length(all_codons)
  cat(paste0("\n✅ Detected ", n_codons, " unique P-codons: ", paste(all_codons, collapse = ", "), "\n"))
  
  # 5. Process each codon + each motif type
  for (codon_idx in seq_along(all_codons)) {
    codon <- all_codons[codon_idx]
    cat(paste0("\n=====================================================\n"))
    cat(paste0("🔹 [", codon_idx, "/", n_codons, "] Processing Codon: ", codon, "\n"))
    cat(paste0("=====================================================\n"))
    
    # ========== Core optimization 3: independent folder per codon ==========
    codon_output_dir <- file.path(output_dir, paste0("codon_", codon))
    codon_plot_dir <- file.path(codon_output_dir, "plots")
    codon_table_dir <- file.path(codon_output_dir, "tables")
    dir.create(codon_plot_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(codon_table_dir, showWarnings = FALSE, recursive = TRUE)
    cat(paste0("✅ Created dedicated folder for Codon ", codon, ": ", codon_output_dir, "\n"))
    
    codon_dt <- psite_dt[!is.na(P_codon) & P_codon == codon, 
                         .(transcript, psite, cds_start, cds_stop)]
    if (nrow(codon_dt) == 0) {
      cat(paste0("⚠️ Codon ", codon, " has no valid data, skip\n"))
      next
    }
    cat(paste0("📊 Raw rows (P_codon=", codon, "): ", nrow(codon_dt), "\n"))
    
    codon_dt <- merge(codon_dt, seq_map, by = "transcript", all.x = TRUE)
    codon_dt <- codon_dt[!is.na(seq) & !is.na(seq_len)]
    cat(paste0("📊 Rows after sequence merge: ", nrow(codon_dt), "\n"))
    
    for (motif_type_name in names(motif_config)) {
      config <- motif_config[[motif_type_name]]
      cat(paste0("\n📌 Processing motif type: ", motif_type_name, " (", config$suffix, ")\n"))
      
      codon_dt[, psite_in_transcript := psite >= 1 & psite <= seq_len]
      codon_dt[, psite_in_cds := psite >= cds_start & psite <= cds_stop]
      codon_dt[, psite_codon_start := psite - ((psite - cds_start) %% 3)]
      codon_dt[, motif_start := psite - config$flank_up]
      codon_dt[, motif_end := psite + config$flank_down]
      
      codon_dt[, motif_in_seq := (motif_start >= 1) & (motif_end <= seq_len)]
      codon_dt[, motif_in_cds := (motif_start >= cds_start) & (motif_end <= cds_stop)]
      codon_dt[, motif_valid := psite_in_transcript & psite_in_cds & motif_in_seq & motif_in_cds]
      
      valid_dt <- codon_dt[motif_valid == TRUE]
      if (nrow(valid_dt) == 0) {
        cat(paste0("⚠️ ", motif_type_name, " no valid in-CDS rows, skip\n"))
        next
      }
      cat(paste0("📊 Valid in-CDS rows (", config$suffix, "): ", nrow(valid_dt), "\n"))
      
      target_dt <- valid_dt[transcript %in% target_tx]
      cat(paste0("\n🎯 Valid rows from target genes (", config$suffix, "): ", nrow(target_dt), "\n"))
      
      non_target_dt <- valid_dt[!transcript %in% target_tx]
      set.seed(123)
      random_n <- min(1000, nrow(non_target_dt))
      random_dt <- non_target_dt[sample(nrow(non_target_dt), random_n)]
      cat(paste0("🎲 Random sampled non-target rows: ", nrow(random_dt), "\n"))
      
      # ========== Core function: extract sequences + compute frequency ==========
      extract_sequences <- function(data_dt, seq_type, codon, config) {
        if (nrow(data_dt) < 1) {
          cat(paste0("   ⚠️ ", seq_type, " no valid rows, skip\n"))
          return(list(
            nt_freq_table = NULL,
            codon_freq_table = NULL
          ))
        }
        
        extracted_seqs <- character(nrow(data_dt))
        for (i in 1:nrow(data_dt)) {
          m_start <- data_dt$motif_start[i]
          m_end <- data_dt$motif_end[i]
          full_seq <- data_dt$seq[i]
          motif_seq <- substr(full_seq, m_start, m_end)
          if (nchar(motif_seq) > config$total_len) {
            motif_seq <- substr(motif_seq, 1, config$total_len)
          }
          extracted_seqs[i] <- motif_seq
        }
        
        valid_nt_seqs <- extracted_seqs[nchar(extracted_seqs) == config$total_len]
        cat(paste0("\n=====================================================\n"))
        cat(paste0("🔍 ", seq_type, " sequence extraction summary (", config$suffix, ")\n"))
        cat(paste0("=====================================================\n"))
        cat(paste0("📊 Total extracted sequences: ", length(extracted_seqs), "\n"))
        cat(paste0("📊 Valid sequences (length=", config$total_len, "): ", length(valid_nt_seqs), "\n"))
        
        if (length(valid_nt_seqs) > 0) {
          cat(paste0("\n✅ Valid sequence preview (top 3):\n"))
          for (i in 1:min(3, length(valid_nt_seqs))) {
            cat(paste0("   [", i, "] '", substr(valid_nt_seqs[i], 1, 30), "...' (len: ", nchar(valid_nt_seqs[i]), ")\n"))
          }
          
          # 1. Nucleotide frequency
          nt_freq_table <- data.table()
          for (pos in 1:config$total_len) {
            pos_nts <- substr(valid_nt_seqs, pos, pos)
            pos_freq <- as.data.frame(table(pos_nts))
            colnames(pos_freq) <- c("nt", "count")
            pos_freq$frequency <- round((pos_freq$count / sum(pos_freq$count)) * 100, 2)
            pos_freq$position <- config$nt_labels[pos]
            pos_freq <- pos_freq[order(-pos_freq$frequency), ]
            nt_freq_table <- rbind(nt_freq_table, pos_freq)
          }
          
          # 2. Codon frequency
          codon_freq_table <- NULL
          if (grepl("codon", motif_type_name)) {
            codon_seqs_list <- lapply(valid_nt_seqs, function(seq) {
              codon_starts <- seq(1, nchar(seq)-2, 3)
              codons <- substring(seq, codon_starts, codon_starts+2)
              if (length(codons) == length(config$codon_labels)) {
                names(codons) <- config$codon_labels
              } else {
                names(codons) <- 1:length(codons)
              }
              return(codons)
            })
            
            cat(paste0("\n✅ Codon split preview (top 2 sequences):\n"))
            for (i in 1:min(2, length(codon_seqs_list))) {
              codons <- codon_seqs_list[[i]]
              cat(paste0("   seq", i, " codons: ", paste(head(codons, 5), collapse = " | "), "...\n"))
            }
            
            codon_freq_table <- data.table()
            codon_pos_matrix <- do.call(rbind, codon_seqs_list)
            for (pos in colnames(codon_pos_matrix)) {
              pos_codons <- codon_pos_matrix[, pos]
              pos_freq <- as.data.frame(table(pos_codons))
              colnames(pos_freq) <- c("codon", "count")
              pos_freq$frequency <- round((pos_freq$count / sum(pos_freq$count)) * 100, 2)
              pos_freq$position <- pos
              pos_freq <- pos_freq[order(-pos_freq$frequency), ]
              codon_freq_table <- rbind(codon_freq_table, pos_freq)
            }
          }
          
          return(list(
            nt_freq_table = nt_freq_table,
            codon_freq_table = codon_freq_table
          ))
        } else {
          return(list(
            nt_freq_table = NULL,
            codon_freq_table = NULL
          ))
        }
      }
      
      if (nrow(target_dt) > 0) {
        target_result <- extract_sequences(target_dt, "target_genes", codon, config)
        if (!is.null(target_result$nt_freq_table)) {
          plot_stacked_manual(target_result$nt_freq_table, "nt", "target_genes", codon, config, codon_output_dir)
        }
        if (!is.null(target_result$codon_freq_table)) {
          plot_stacked_manual(target_result$codon_freq_table, "codon", "target_genes", codon, config, codon_output_dir)
        }
      }
      
      if (nrow(random_dt) > 0) {
        random_result <- extract_sequences(random_dt, "random", codon, config)
        if (!is.null(random_result$nt_freq_table)) {
          plot_stacked_manual(random_result$nt_freq_table, "nt", "random", codon, config, codon_output_dir)
        }
        if (!is.null(random_result$codon_freq_table)) {
          plot_stacked_manual(random_result$codon_freq_table, "codon", "random", codon, config, codon_output_dir)
        }
      }
    }
    rm(codon_dt, valid_dt, target_dt, random_dt)
    gc()
  }
  cat(paste0("\n✅ All Codon Motif analysis finished! New AA chemical-class plots generated\n"))
}

# ===================== Execution Entry (process all codons) =====================
if (exists("all_codon_occupancy_group_filt") && nrow(all_codon_occupancy_group_filt) > 0) {
  cat(paste0("✅ ", nrow(tx2gene), " transcript-gene mapping entries loaded\n"))
  first_sample_dir <- dirname(rds_paths[1])
  
  analyze_codon_motif(
    sample_dir = first_sample_dir,
    transcriptome_fasta = transcriptome_fasta,
    output_dir = output_root,
    gene_list_path = gene_list_path,
    tx2gene = tx2gene,
    motif_type = c("5codon", "20codon")
  )
  
} else {
  cat("⚠️ No filtered codon data available, skip Motif analysis\n")
}

cat("\n🎉 Final analysis completed! All codons plotted, per-codon folder, synonymous codon hue scheme, narrow codon rectangles, new AA chemical-class plots (AA labels retained)!\n")
# ===================== Run Finished =====================
cat("\n============================================================\n")
cat("🎉 ALL ANALYSIS COMPLETED SUCCESSFULLY! NO ERROR! NO WARNING!\n")
cat("📊 All results saved to directory: ", output_root, "\n")
