# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                  amplification characteristics                      # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
####################################################################    library

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggplot2)
library(reshape2)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(cowplot)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

###############################################################################
#################################################################### file 

sample_table_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
ecdna_path              <- "/re_gecip/cancer_lung/CBailey1/pancan_ecdna_analysis/data/summary_objects/ALL_CANCER_TYPES_BED_LIST_2023-01-18.RDS"
copy_number_driver_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copy_number_driver_genes.txt"
output_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/genomic_disruption/"

amp_gene_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/all_hist_amplification_driver_genes.pdf"
loss_gene_path          <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/all_hist_homozygous_deletion_driver_genes.pdf"

###############################################################################
####################################################################      MAIN

# read in the data
sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

amp_data <- readRDS(ecdna_path)
amp_data <- amp_data[amp_data$participant_id %in% sample_table$participant_id, ]
amp_data <- left_join(amp_data, sample_table[, c("participant_id", "histology")])

driver <- read.table(copy_number_driver_path, head = T, sep = "\t")
# only keep amplification driver for this analysis
driver <- driver[driver$type == "amplification", ]

# get gene locations 
# get locations of genes
txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes                 <- genes(txdb, single.strand.genes.only=FALSE)
genes                 <- data.frame(genes)

# extract HUGO symbols
gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)

# merge with gene locations
gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))
gene_names            <- gene_names[!is.na(gene_names$symbol), ]

# add gene locations to driver genes
driver <- left_join(driver, gene_names[, c("symbol", "seqnames", "start", "end")], by = c("gene" = "symbol"))
driver <- driver[grep("alt", driver$seqnames, invert = T), ]

# for each histology individually overlap the amplifications with the copy number driver genes

all_amp_driver <- list()
for(hist in c("ADENOCARCINOMA", "SQUAMOUS_CELL")){
  
  # get the histology specific driver genes
  hist_driver <- driver[driver$histology == hist, ]
  hist_driver_GR <- GRanges(seqnames = hist_driver$seqnames, IRanges(start = hist_driver$start, end = hist_driver$end))
  mcols(hist_driver_GR) <- data.frame(hist_driver[, c("gene", "is_cancer_gene", "histology", "samples_per_bin", "sample_per_bin_prop", "type")])
  
  # get the amplifications
  if(hist %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")){
    hist_amp_data <- amp_data[amp_data$histology == hist, ]
  } else {hist_amp_data <- amp_data}
  
  hist_amp_GR <- GRanges(seqnames = hist_amp_data$chr, IRanges(start = hist_amp_data$start, end = hist_amp_data$end))
  mcols(hist_amp_GR) <- data.frame(hist_amp_data[, c("participant_id", "feature_type", "histology")])
  
  overlap         <-  findOverlaps(hist_amp_GR, hist_driver_GR)
  intersect_df    <-  cbind(hist_amp_data[queryHits(overlap), ], hist_driver[subjectHits(overlap), ])

  all_amp_driver[[hist]] <- intersect_df
}

all_amp_driver <- do.call(rbind, all_amp_driver)

amp_drivers <- unique(all_amp_driver[, c("participant_id", "histology", "feature_type", "gene", "is_cancer_gene", "samples_per_bin")])

amp_driver_count <- list()
for(hist in c("ADENOCARCINOMA", "SQUAMOUS_CELL")){
  
  hist_amp_driver <- amp_drivers[amp_drivers$histology == hist, ]
  hist_amp_driver_count <- data.frame(table(hist_amp_driver[, c("feature_type", "gene")]))
  hist_amp_driver_count$histology <- hist
  
  hist_amp_driver_count <- left_join(hist_amp_driver_count, unique(hist_amp_driver[, c("gene", "samples_per_bin")]))
  
  amp_driver_count[[hist]] <- hist_amp_driver_count
}
amp_driver_count <- do.call(rbind, amp_driver_count)

amp_driver_count$prop <- amp_driver_count$Freq / amp_driver_count$samples_per_bin

# calculate proportion again but based on the numbers we get with amplicon architect
amp_driver_count_sum <- list()

for(gene in unique(amp_driver_count$gene)){
  
  df <- amp_driver_count[amp_driver_count$gene == gene, ]
  
  all_amp <- sum(df$Freq)
  
  df$prop_AA <- df$Freq/all_amp
  
  amp_driver_count_sum[[gene]] <- df
  
}

amp_driver_count_sum <- do.call(rbind, amp_driver_count_sum)
amp_driver_count_sum$prop_AA <- amp_driver_count_sum$prop_AA*100


# for each histology make a plot

all_p <- list()
for(hist in c("ADENOCARCINOMA", "SQUAMOUS_CELL")){
  
  amp_driver_count_hist <- amp_driver_count_sum[amp_driver_count_sum$histology == hist, ]
  amp_driver_count_hist <- amp_driver_count_hist[order(amp_driver_count_hist$prop_AA, decreasing = T), ]
  amp_driver_count_hist$gene <- factor(amp_driver_count_hist$gene, levels = unique(amp_driver_count_hist$gene))
  
  if(hist == "ADENOCARCINOMA"){
    color <- "#67001f"
    hist_name <- "Adenocarcinoma"
  }
  
  if(hist == "SQUAMOUS_CELL"){
    color <- "#053061"
    hist_name <- "Squamous cell"
  }

  amp_driver_count_hist$feature_type <- factor(amp_driver_count_hist$feature_type, levels = c("Linear amplification", "ecDNA", "BFB", "Complex non-cyclic", "unknown"))
  
  # add number of tumours with amplification
  amp_driver_count_hist$gene_name <- paste0(amp_driver_count_hist$gene, " (", amp_driver_count_hist$samples_per_bin, ")")
  amp_driver_count_hist <- amp_driver_count_hist[order(amp_driver_count_hist$samples_per_bin, decreasing = T), ]
  amp_driver_count_hist$gene <- factor(amp_driver_count_hist$gene, levels = unique(amp_driver_count_hist$gene))
  
  # order the genes
  amp_driver_count_hist <- amp_driver_count_hist[order(amp_driver_count_hist$Freq, decreasing = T), ]
  amp_driver_count_hist$gene <- factor(amp_driver_count_hist$gene, levels = unique(amp_driver_count_hist$gene))
  
  
  p <- ggplot(amp_driver_count_hist, aes(gene, feature_type, fill = prop_AA, label = Freq)) +
        geom_tile() +
        geom_text() +
        scale_fill_gradient(low = "white", high = color, name = paste0("proportion of\namplifications in\n", hist_name, "s")) +
        xlab("") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  all_p[[hist]] <- p
}

# make them into one plot

plots_combined <- plot_grid(all_p[[1]], all_p[[2]], ncol = 1, rel_widths = c(1, 1))

pdf(paste0(output_path, "amplification_types_LUAD_LUSC.pdf"), width = 6, height = 6)
plots_combined
dev.off()

