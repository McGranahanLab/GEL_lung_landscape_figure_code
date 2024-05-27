# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                   look for known fusion drivers                     # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) #change to hg19 depending on input data
library(org.Hs.eg.db)
library(grid)
library(gridExtra)
library(gtable)
library(rtracklayer)
library(ggbeeswarm)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cowplot)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

###############################################################################
################################################################### file paths

sample_table_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
SV_path                 <- "/re_gecip/cancer_lung/kthol/SV/input/SV_list.Rdata"
out_path                <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/SV/fusions.txt"
driver_path             <- "/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/2023_12_25/results/tables/driver_mutations/"
SV_driver_mut_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/SV_driver_genes_in_samples.txt"
CN_driver_mut_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copynumber_peak_driver_per_patient.txt"
out_path_plot           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/fusion_genes.pdf"

####fusion_count###########################################################################
###################################################################       MAIN

sample_table <- read.table(sample_table_path, head = TRUE, sep = "\t")

load(SV_path)

all_SV_concatenated <- all_SV_concatenated[all_SV_concatenated$sample %in% sample_table$participant_id, ]

all_SV_concatenated$chr1 <- sub("chr", "", all_SV_concatenated$chr1)
all_SV_concatenated$chr2 <- sub("chr", "", all_SV_concatenated$chr2)

all_SV_concatenated[all_SV_concatenated$chr1 == "X", "chr1"] <- 23
all_SV_concatenated[all_SV_concatenated$chr2 == "X", "chr2"] <- 23
all_SV_concatenated[all_SV_concatenated$chr1 == "Y", "chr1"] <- 24
all_SV_concatenated[all_SV_concatenated$chr2 == "Y", "chr2"] <- 24

all_SV_concatenated$chr1 <- as.numeric(all_SV_concatenated$chr1)
all_SV_concatenated$chr2 <- as.numeric(all_SV_concatenated$chr2)

all_SV_concatenated$pos1 <- as.numeric(as.character(all_SV_concatenated$pos1))
all_SV_concatenated$pos2 <- as.numeric(as.character(all_SV_concatenated$pos2))

# lots of them are duplicated, let's remove

all_SV_concatenated <- all_SV_concatenated %>% rowwise() %>% mutate(ID = paste0(min(chr1, chr2), "_", min(pos1, pos2), "_", max(chr1, chr2), "_", max(pos1, pos2)))
all_SV_concatenated <- data.frame(all_SV_concatenated)

# add sample to the ID column
all_SV_concatenated$ID <- paste0(all_SV_concatenated$ID, "_", all_SV_concatenated$sample)

# kick out duplicates
dups <- duplicated(all_SV_concatenated$ID)

all_SV_concatenated <- all_SV_concatenated[dups == FALSE, ]

# filter them by size a little
all_SV_concatenated$width <- all_SV_concatenated$pos2 - all_SV_concatenated$pos1
all_SV_concatenated$width <- abs(all_SV_concatenated$width)

# remove everything smaller than 50kb
all_SV_concatenated <- all_SV_concatenated[all_SV_concatenated$width >= 50, ]

# get rid of some unnecessary columns
all_SV_concatenated <- all_SV_concatenated[, c("chr1", "pos1", "chr2", "pos2", "classification","sample")]
all_SV_concatenated <- unique(all_SV_concatenated)

# get the names and locations of genes
# load gene data
txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes                 <- genes(txdb, single.strand.genes.only=FALSE)
genes                 <- data.frame(genes)

# extract HUGO symbols
gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)

# merge with gene locations
gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))

# gene_names <- gene_names[gene_names$symbol %in% driver$gene, ]
gene_names <- gene_names[!is.na(gene_names$symbol), ]

genes_interest <- gene_names[gene_names$symbol %in% c("ALK", "ROS1", "RET", "NTRK1", "NTRK2", "NTRK3", "NRG1"), ]

# find SVs in these genes
all_SV_concatenated$chr1 <- paste0("chr", all_SV_concatenated$chr1)
all_SV_concatenated$chr2 <- paste0("chr", all_SV_concatenated$chr2)
all_SV_concatenated$SV_ID <- paste(all_SV_concatenated$chr1, all_SV_concatenated$pos1, all_SV_concatenated$chr2, all_SV_concatenated$pos2, all_SV_concatenated$classification, all_SV_concatenated$sample, sep = "_")

SV_1 <- all_SV_concatenated[, c("sample", "classification", "chr1", "pos1", "SV_ID")]
SV_2 <- all_SV_concatenated[, c("sample", "classification", "chr2", "pos2", "SV_ID")]
colnames(SV_1) <- c("sample", "classification", "chr", "pos", "SV_ID")
colnames(SV_2) <- c("sample", "classification", "chr", "pos", "SV_ID")

SVs <- rbind(SV_1, SV_2)

SV_gr <- GRanges(seqnames = SVs$chr, IRanges(start = SVs$pos, end = SVs$pos))
gene_interest_gr <- GRanges(seqnames = genes_interest$seqnames, IRanges(start = genes_interest$start, end = genes_interest$end))

overlap         <-  findOverlaps(SV_gr, gene_interest_gr)
fusion_interest_df    <-  cbind(SVs[queryHits(overlap), ], genes_interest[subjectHits(overlap), ])

# and also overlap all SVs with all genes so we can see which genes the genes of interest are fused together with
gene_gr <- GRanges(seqnames = gene_names$seqnames, IRanges(start = gene_names$start, end = gene_names$end))

overlap_all         <-  findOverlaps(SV_gr, gene_gr)
SV_gene_df          <-  cbind(SVs[queryHits(overlap_all), ], gene_names[subjectHits(overlap_all), ])

# the SVs which are on one of the genes of interest, where does thei start/stop fall?
fusion_interest <- unique(fusion_interest_df[, c("sample", "SV_ID", "classification", "chr", "pos", "symbol", "seqnames" , "start", "end")])
colnames(fusion_interest) <- c("sample", "SV_ID", "classification", "chr1", "pos1", "symbol1", "gene_chr1", "gene_start1", "gene_end1")

SV_gene         <- unique(SV_gene_df[, c("sample", "SV_ID", "classification", "chr", "pos", "symbol", "seqnames" , "start", "end")])
colnames(SV_gene) <-  c("sample", "SV_ID", "classification", "chr2", "pos2", "symbol2", "gene_chr2", "gene_start2", "gene_end2")

fusion_pairs <- left_join(fusion_interest, SV_gene)
fusion_pairs$chr_pos1 <- paste0(fusion_pairs$chr1, "_", fusion_pairs$pos1)
fusion_pairs$chr_pos2 <- paste0(fusion_pairs$chr2, "_", fusion_pairs$pos2)

# remove all the ones where chr1 and pos1 are equal to chr2 and pos2
fusion_pairs <- fusion_pairs[fusion_pairs$chr_pos1 != fusion_pairs$chr_pos2, ]

# remove duplications and deletions
fusion_pairs_df <- fusion_pairs[!fusion_pairs$classification %in% c("DEL", "DUP"),]
fusion_pairs_df <- fusion_pairs_df[fusion_pairs_df$symbol1 != fusion_pairs_df$symbol2, ]

write.table(fusion_pairs_df, paste0(out_path), sep = "\t", row.names = F, col.names = T, quote = F)

# count
fusion_samples <- unique(fusion_pairs_df[, c("sample", "symbol1", "symbol2")])
fusion_samples <- fusion_samples[fusion_samples$symbol1 != fusion_samples$symbol2, ]
fusion_samples <- fusion_samples[fusion_samples$symbol2 != "NRG1-IT3", ]
fusion_count <- data.frame(table(fusion_samples$symbol1))

# which other driver genes did patients have with these fusions

# read all the driver mutation files for all histologies
all_files <- list.files(driver_path, full.names = T)

SNV_driver <- list()
for(i in all_files){
  
  print(i)
  df <- read.delim(i)
  df <- df[df$var_class != "amp", ]
  df <- df[df$var_class != "hom.del.", ]
  df <- unique(df[, c("participant_id", "gene_name", "gr_id", "key")])
  
  SNV_driver[[i]] <- df
}

SNV_driver <- do.call(rbind, SNV_driver)
rownames(SNV_driver) <- c(1:nrow(SNV_driver))
SNV_driver <- unique(SNV_driver)
SNV_driver$participant_id <- as.character(SNV_driver$participant_id)
SNV_driver$gene_name <- paste0(SNV_driver$gene_name, ":", SNV_driver$gr_id)
SNV_driver <- SNV_driver[, c("participant_id", "gene_name")]
colnames(SNV_driver) <- c("sample", "SNV_driver")

# SVdriver
SV_driver <- read.table(SV_driver_mut_path, head = T, sep = "\t")
SV_driver$gene_gr <- paste0(SV_driver$gene_name, ":", SV_driver$element_type)
SV_driver <- unique(SV_driver[, c("participant_id", "gene_gr")])
colnames(SV_driver) <- c("sample", "SV_driver")
SV_driver$sample <- as.character(SV_driver$sample)

# CN driver
CN_driver <- read.table(CN_driver_mut_path, head = T, sep = "\t")
CN_driver$peak_gene <- paste0(CN_driver$new_peak_name, ":", CN_driver$significant_CN_genes_on_peak)
CN_driver <- unique(CN_driver[, c("patient", "peak_gene")])
colnames(CN_driver) <- c("sample", "CN_driver")
CN_driver$sample <- as.character(CN_driver$sample)

fusion_samples$fusion <- paste0(fusion_samples$symbol1, ":", fusion_samples$symbol2)
fusion_driver_count <- left_join(fusion_samples, unique(SNV_driver))
fusion_driver_count <- left_join(fusion_driver_count, unique(SV_driver))
fusion_driver_count <- left_join(fusion_driver_count, unique(CN_driver))


plot_df <- data.frame(table(fusion_samples[, c("symbol1", "symbol2")]))
plot_df <- plot_df[order(plot_df$Freq, decreasing = T),]
plot_df$symbol1 <- factor(plot_df$symbol1, levels = unique(plot_df$symbol1))
plot_df$symbol2 <- factor(plot_df$symbol2, levels = unique(plot_df$symbol2))

p <- ggplot(plot_df, aes(symbol1, symbol2, size = Freq)) +
        geom_point() +
        xlab("fusion gene 1") +
        ylab("fusion gene 2") +
        theme_bw()

plot2_df <- unique(fusion_driver_count[, c("sample", "symbol1", "SNV_driver")])
plot2_df <- data.frame(table(plot2_df[, c("symbol1", "SNV_driver")]))
plot2_df <- plot2_df[order(plot2_df$Freq, decreasing = T), ]
plot2_df$symbol1 <- factor(plot2_df$symbol1, levels = unique(plot_df$symbol1))
plot2_df$SNV_driver <- factor(plot2_df$SNV_driver, levels = unique(plot2_df$SNV_driver))

p2 <- ggplot(plot2_df, aes(symbol1, SNV_driver, size = Freq)) +
          geom_point() +
          xlab("") +
          ylab("SNV/indel driver") +
          theme_bw() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank())

p_all <- plot_grid(p2, p, ncol = 1, align = "v")

pdf(out_path_plot, width = 5, height = 7)
p_all
dev.off()

