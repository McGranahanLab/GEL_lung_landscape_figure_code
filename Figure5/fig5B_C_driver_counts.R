# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                               driver numbers                        # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
###################################################################### libraries


.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ComplexUpset)
library(survival)
library(survminer)
library(stringr)

################################################################################
##################################################################### file paths

sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
driver_path           <- "/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/2023_12_25/results/tables/driver_mutations/"
SV_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/SV_driver_genes_in_samples.txt"
CN_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copynumber_peak_driver_per_patient.txt"
output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/"
survival_path         <- "/re_gecip/cancer_lung/shared/pancan_survival_release_v16.RDS"
smoking_path          <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/smoking_data.txt"
fusion_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/SV/fusions.txt"
germline_path         <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/germline/germline_second_hit.txt"

################################################################################
######################################################################      MAIN

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

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
SNV_driver$gene_gr <- paste0(SNV_driver$gene_name, ":", SNV_driver$gr_id)


# first let's recreate the bubble plot of number coding / non-coding driver
# how many coding / non-coding driver mutations does each tumour have?
SNV_driver[SNV_driver$gr_id == "CDS", "class"] <- "coding" 
SNV_driver[SNV_driver$gr_id != "CDS", "class"] <- "non_coding" 

SNV_driver[SNV_driver$class == "coding", "coding_SNVdriver"] <- "codingSNV"
SNV_driver[SNV_driver$class == "non_coding", "noncoding_SNVdriver"] <- "noncodingSNV"

driver_count <- data.frame(table(SNV_driver[, c("participant_id", "class")]))

# now i need to combine this for each patient
driver_count_mat <- dcast(participant_id ~ class, data = driver_count)
driver_count_mat$participant_id <- as.character(driver_count_mat$participant_id)

# add all tumours that are missing as 0/0
missing_tumours <- sample_table[-which(sample_table$participant_id %in% driver_count_mat$participant_id), "participant_id"]
missing_tumours_df <- data.frame(participant_id = missing_tumours,
                                 coding = 0, 
                                 non_coding = 0)

driver_count_mat <- rbind(driver_count_mat, missing_tumours_df)

# calculate average number of drivers
sum(driver_count_mat$coding) / nrow(driver_count_mat)
sum(driver_count_mat$non_coding) / nrow(driver_count_mat)
(sum(driver_count_mat$coding) + sum(driver_count_mat$non_coding)) / nrow(driver_count_mat)
nrow(driver_count_mat[driver_count_mat$coding != 0 | driver_count_mat$non_coding != 0,]) / nrow(driver_count_mat)
nrow(driver_count_mat[driver_count_mat$coding != 0 & driver_count_mat$non_coding != 0,]) / nrow(driver_count_mat)
nrow(driver_count_mat[driver_count_mat$coding != 0 & driver_count_mat$non_coding == 0,]) / nrow(driver_count_mat)
nrow(driver_count_mat[driver_count_mat$coding == 0 & driver_count_mat$non_coding != 0,]) / nrow(driver_count_mat)

driver_count_class <- data.frame(table(driver_count_mat[, c("coding", "non_coding")]))

# no make this into a plot

# bin Freq
driver_count_class <- driver_count_class %>% mutate(Freq_bin = cut(Freq, breaks=c(0, 5, 10, 25, 50, 75, 100, 150, 200, 218)))


p_bubble <- ggplot() +
  geom_point(data = driver_count_class[driver_count_class$Freq > 0, ], aes(coding, non_coding, color = Freq_bin, size = Freq)) +
  geom_text(data = driver_count_class[driver_count_class$Freq >= 5, ], aes(coding, non_coding, label = Freq), color = "black", vjust = -1) +
  scale_color_manual(name = "# tumours", values = c("#41ab5d", "#006837", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")) +
  scale_size_continuous(name = "# tumours") +
  xlab("Number of coding cancer gene mutations") +
  ylab("Number of non coding cancer gene mutations") +
  theme_bw()

pdf(paste0(output_path, "number_coding_vs_noncoding_SNV_driver.pdf"), width = 5, height = 4)
p_bubble
dev.off()

# check histologies
sample_table$participant_id <- as.character(sample_table$participant_id)
driver_count_hist<- left_join(driver_count_mat, sample_table[, c("participant_id", "histology")])

# are these enriched for certain histologies

driver_count_hist[driver_count_hist$non_coding > 0 & driver_count_hist$coding == 0, "only_non_coding"] <- TRUE
driver_count_hist[is.na(driver_count_hist$only_non_coding), "only_non_coding"] <- FALSE




# now let's check what happens when we also take SVs and copy number into account and fusions and germline
SV_driver <- read.table(SV_driver_mut_path, head = T, sep = "\t")
SV_driver$SVdriver <- "SV"
SV_driver <- unique(SV_driver[, c("participant_id", "SVdriver")])
SV_driver$participant_id <- as.character(SV_driver$participant_id)

CN_driver <- read.table(CN_driver_mut_path, head = T, sep = "\t")
CN_driver$CNdriver <- "CN"
CN_driver <- unique(CN_driver[, c("patient", "CNdriver")])
colnames(CN_driver) <- c("participant_id", "CNdriver")
CN_driver$participant_id <- as.character(CN_driver$participant_id)

fusions <- read.table(fusion_path, head = T, sep = "\t")
fusions <- unique(fusions[, c("sample", "symbol1", "symbol2")])
fusions <- fusions[fusions$symbol1 != fusions$symbol2, ]
fusions <- fusions[fusions$symbol2 != "NRG1-IT3", ]
fusions <- unique(fusions[, c("sample", "symbol1")])
fusions$fusion <- "fusion"
fusions$symbol1 <- NULL
colnames(fusions) <- c("participant_id", "fusion")
fusions$participant_id <- as.character(fusions$participant_id)

germline <- read.table(germline_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)
germline$patient <- as.character(germline$patient)
second_hit <- left_join(germline, sample_table[, c("participant_id", "histology")], by = c("patient" = "participant_id"))

second_hit <- second_hit[second_hit$germline_driver_present, ]

# where one patient hase two genes mutated lets split that up into two rows
germline_gene_split <- list()
for(i in 1:nrow(second_hit)){
  
  df <- second_hit[i, ]
  
  somatic <- strsplit(df$germline_genes_somatically_mutated, "_")[[1]]
  sv <- strsplit(df$germline_genes_somatic_SV, "_")[[1]]
  germline_mut <- strsplit(df$germline_driver_gene, "_")[[1]]
  germline_loh <- strsplit(df$germline_driver_LOH_genes, "_")[[1]]
  
  trues <- c(length(somatic) > 1, length(sv) > 1, length(germline_mut) > 1, length(germline_loh) > 1)
  
  if(any(trues)){
    df_out <- data.frame(patient = df$patient,
                         germline_gene_somatic_mutation = df$germline_gene_somatic_mutation,
                         germline_genes_somatically_mutated = somatic,
                         germline_genes_somatic_SV = sv,
                         germline_gene_somatic_SV_present = df$germline_gene_somatic_SV_present,
                         germline_driver_present = df$germline_driver_present,
                         germline_driver_gene = germline_mut,
                         germline_driver_LOH_present = df$germline_driver_LOH_present,
                         germline_driver_LOH_genes = germline_loh,
                         histology = df$histology)
    
  } else {df_out <- df}
  
  germline_gene_split[[i]] <- df_out
}
germline_gene_split <- do.call(rbind, germline_gene_split)

# indicate whether there is actually a second hit
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_driver_LOH_genes), "second_hit"] <- "TRUE"
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_genes_somatically_mutated), "second_hit"] <- "TRUE"
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_genes_somatic_SV), "second_hit"] <- "TRUE"
germline_gene_split[is.na(germline_gene_split$second_hit), "second_hit"] <- "FALSE"

germline_gene_split <- unique(germline_gene_split[, c("patient", "second_hit")])
germline_gene_split <- germline_gene_split[germline_gene_split$second_hit == TRUE, ]
colnames(germline_gene_split) <- c("participant_id", "germline")
germline_gene_split$germline <- "germline"


driver_count_mat[driver_count_mat$coding >= 1, "coding"] <- "coding_SNV"
driver_count_mat[driver_count_mat$non_coding >= 1, "non_coding"] <- "noncoding_SNV"
driver_count_mat[driver_count_mat$coding == 0, "coding"] <- NA
driver_count_mat[driver_count_mat$non_coding == 0, "non_coding"] <- NA

driver <- left_join(driver_count_mat, SV_driver)
driver <- left_join(driver, CN_driver)
driver <- left_join(driver, fusions)
driver <- left_join(driver, germline_gene_split)
driver <- unique(driver)

driver$class <- paste0(driver$coding, "_", driver$non_coding, "_", driver$SVdriver, "_", driver$CNdriver, "_", driver$fusion, "_", driver$germline)
driver$class <- sub("_NA_", "_", driver$class)
driver$class <- sub("_NA", "", driver$class)
driver$class <- sub("NA_", "", driver$class)
driver$class <- sub("_NA", "", driver$class)
driver$class <- sub("NA_", "", driver$class)
driver$class <- sub("_NA", "", driver$class)

write.table(driver, paste0(output_path, "all_driver_binary.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

# add histology
sample_table$participant_id <- as.character(sample_table$participant_id)
driver <- left_join(driver, sample_table[, c("participant_id", "histology")])
driver <- unique(driver)

# make a bar plot only for those tumours that don't have a SNV driver
no_driver <- driver[which(is.na(driver$coding) & is.na(driver$non_coding)), ]

no_driver_count <- data.frame(table(no_driver[, c("class", "histology")]))
no_driver_count$class <- as.character(no_driver_count$class)
no_driver_count[no_driver_count$class == "NA", "class"] <- "no_driver"

# indicate if has no driver
no_driver_count[no_driver_count$class == "no_driver", "group"] <- "no driver"
no_driver_count[no_driver_count$class != "no_driver", "group"] <- "driver"

no_driver_count$group <- factor(no_driver_count$group, levels = c("no driver", "driver"))
no_driver_count$class <- factor(no_driver_count$class, levels = c("no_driver", "CN", "SV_CN", "SV", "fusion", "CN_fusion"))

no_driver_count_class_only <- no_driver_count %>% group_by(class) %>% summarise(sum(Freq))
no_driver_count_class_only[no_driver_count_class_only$class == "no_driver", "group"] <- "no driver"
no_driver_count_class_only[no_driver_count_class_only$class != "no_driver", "group"] <- "driver"
no_driver_count_class_only$group <- factor(no_driver_count_class_only$group, levels = c("no driver", "driver"))

no_driver_count$histology <- factor(no_driver_count$histology, levels = c("ADENOCARCINOMA",   
                                                                          "MET_ADENOCARCINOMA",
                                                                          "SQUAMOUS_CELL",
                                                                          "MET_SQUAMOUS_CELL",
                                                                          "ADENOSQUAMOUS",
                                                                          "CARCINOID",
                                                                          "MESOTHELIOMA",
                                                                          "SMALL_CELL",
                                                                          "MET_SMALL_CELL",
                                                                          "NEUROENDOCRINE_CARCINOMA",
                                                                          "OTHER"))

p_no_driver <- ggplot() +
  geom_bar(data = no_driver_count, aes(class, Freq, fill = histology), stat = "identity") +
  geom_text(data = no_driver_count_class_only, aes(class, `sum(Freq)`, label = `sum(Freq)`, vjust = -0.2)) +
  scale_fill_manual(values = c("ADENOCARCINOMA" = "#67001f",   
                               "MET_ADENOCARCINOMA" = "#67001f80",
                               "SQUAMOUS_CELL" = "#053061",
                               "MET_SQUAMOUS_CELL" = "#05306180",
                               "ADENOSQUAMOUS" = "#66c2a5",
                               "CARCINOID" = "#fc8d62",
                               "MESOTHELIOMA" = "#e3309e",
                               "SMALL_CELL" = "#b89704",
                               "MET_SMALL_CELL" = "#b8970480",
                               "NEUROENDOCRINE_CARCINOMA" = "#0b8ca3",
                               "OTHER" = "#7a7979"),
                    labels = c("Adenocarcinoma",
                               "Adenocarcinoma\nmetastasis",
                               "Squamous cell",
                               "Squamous cell\nmetastasis",
                               "Adenosquamous",
                               "Carcinoid",
                               "Mesothelioma",
                               "Small cell",
                               "Small cell\nmetastasis",
                               "Neuroendocrine\ncarcinoma",
                               "Other",
                               "Other\nmetastasis")) +
  theme_bw() +
  xlab("") +
  ylab("# tumours") +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = c("no_driver" = "no driver",
                              "CN" = "CN", 
                              "SV_CN" = "CN+SV", 
                              "SV" = "SV",
                              "fusion" = "fusion",
                              "CN_fusion" = "CN+fusion")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill="white"))


pdf(paste0(output_path, "number_no_SNV_driver.pdf"), width = 5, height = 5)
p_no_driver
dev.off()
