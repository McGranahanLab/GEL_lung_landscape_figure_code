# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # #                         LOH in histologies                          # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
################################################################### libraries
library(dplyr)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggpubr)
library(ggbeeswarm)

################################################################################
################################################################### file path
cn_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/combined_ascatTable_segments_joined.rds"
sample_table_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
output_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/GenomicOverview/"

################################################################################
################################################################### MAIN

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

cn <- readRDS(cn_path)

# indicate LOH and haploid LOH
cn[cn$nMinor == 0 & cn$nMajor > 0, "LOH"] <- TRUE
cn[is.na(cn$LOH), "LOH"] <- FALSE
cn[cn$nMinor == 0 & cn$nMajor == 1, "haploidLOH"] <- TRUE
cn[is.na(cn$haploidLOH), "haploidLOH"] <- FALSE

# calculate the fraction of the genome thats subject to LOH and haploid LOH
chr_sizes                 <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0('chr', c(1:22, "X"))]
chr_sizes                 <- data.frame(chr_sizes)
genome_length             <- sum(chr_sizes$chr_sizes)

# calculate the fraction of the genome that has LOH
cn$seg_length <- cn$endpos - cn$startpos

loh_frac <- cn %>% group_by(patient, LOH, haploidLOH) %>% summarise(loh_sum = sum(seg_length))
loh_frac$LOH_status <- paste0(loh_frac$LOH, "_", loh_frac$haploidLOH  )


# some are missing, fill those in

LOH_mat <- dcast(loh_frac, patient ~ LOH_status, value.var = "loh_sum")
LOH_mat[is.na(LOH_mat)] <- 0
LOH_mat$FALSE_FALSE <- NULL

colnames(LOH_mat) <- c("patient", "LOH", "haploidLOH")
LOH_mat <- melt(LOH_mat, id.vars = "patient")
LOH_mat$genome_proportion <- LOH_mat$value / genome_length
LOH_mat$genome_proportion <- LOH_mat$genome_proportion*100

# add histology 
LOH_mat <- left_join(LOH_mat, sample_table[, c("participant_id", "histology")], by = c("patient" = "participant_id"))
LOH_mat <- LOH_mat[LOH_mat$variable == "LOH", ]

median_LOH <- LOH_mat %>% group_by(histology, variable) %>% summarise(median = median(genome_proportion))
median_LOH <- median_LOH[order(median_LOH$median, decreasing = T), ]

LOH_mat$histology <- factor(LOH_mat$histology, levels = median_LOH$histology)

p <- ggplot(LOH_mat, aes(histology, genome_proportion, fill = histology)) +
  geom_boxplot() +
  geom_quasirandom(alpha = 0.3) +
  stat_compare_means(label.y = 70, label.x = 3) +
  ylab("% genome with LOH") +
  xlab("") +
  scale_fill_manual(values = c("ADENOCARCINOMA" = "#67001f",
                               "MET_ADENOCARCINOMA" = "#67001f80",
                               "SQUAMOUS_CELL" = "#053061",
                               "MET_SQUAMOUS_CELL" = "#05306180",
                               "ADENOSQUAMOUS" = "#66c2a5",
                               "CARCINOID" = "#fc8d62",
                               "LARGE_CELL" = "#54278f",
                               "MESOTHELIOMA" = "#e3309e",
                               "SMALL_CELL" = "#b89704",
                               "MET_SMALL_CELL" = "#b8970480",
                               "NEUROENDOCRINE_CARCINOMA" = "#0b8ca3",
                               "OTHER" = "#7a7979")) + 
  scale_x_discrete(labels = c("Large cell",
                              "Squamous cell",
                              "Small cell",
                              "Adenocarcinoma\nmetastasis",
                              "Adenosquamous",
                              "Squamous cell\nmetastasis",
                              "Other",
                              "Adenocarcinoma",
                              "Neuroendocrine",
                              "Small cell\nmetastasis",
                              "Mesothelioma",
                              "Carcinoid")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pdf(paste0(output_path, "LOH_genome_percentage.pdf"), width = 5, height = 3)
p
dev.off()






