# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                             ONCOPRINT                               # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
##################################################################### libraries
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)

###############################################################################
#####################################################################file paths

sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
SNV_driver_path       <- "/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/2023_12_25/results/tables/driver_mutations/"
fusion_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/SV/fusions.txt"
SV_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/SV_driver_genes_in_samples.txt"
CN_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copynumber_peak_driver_per_patient.txt"
germline_path         <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/germline/germline_second_hit.txt"
output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/"

###############################################################################
#####################################################################     MAIN

sample_table <- read.table(sample_table_path, head = T, sep = "\t")

# read in SNV driver mutations
# read all the driver mutation files for all histologies
all_files <- list.files(SNV_driver_path, full.names = T)

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
SNV_driver[SNV_driver$gr_id == "CDS", "coding"] <- SNV_driver[SNV_driver$gr_id == "CDS", "gene_name"]
SNV_driver[SNV_driver$gr_id != "CDS", "non_coding"] <- SNV_driver[SNV_driver$gr_id != "CDS", "gene_gr"]
SNV_driver <- unique(SNV_driver[, c("participant_id", "coding", "non_coding")])

# load the CN drivers
CN_driver_sample <- read.table(CN_driver_mut_path, head = T, sep = "\t")
CN_driver_sample <- CN_driver_sample[grep("PANCAN", CN_driver_sample$new_peak_name), ]
CN_driver_sample$name <- paste0(CN_driver_sample$new_peak_name, ":", CN_driver_sample$significant_CN_genes_on_peak)
CN_driver_sample$name <- sub("PANCAN_", "", CN_driver_sample$name)
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "TBL1XR1_PIK3CA_SOX2_MAP3K13_BCL6_TP63_MB21D2_MUC4_MECOM_STAG1_PIK3CB_ATR", "gene"] <- "PIK2CA, SOX2"
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "SDHA_TERT_CTNND2_CDH10_DROSHA_IL7R_NIPBL", "gene"] <- "TERT"
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "CSMD1_RPL23AP53_ZNF596_FAM87A_FBXO25_TDRP_ERICH1", "gene"] <- NA
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "SLC16A12_PANK1_FLJ37201_KIF20B_LINC00865_LINC01374", "gene"] <- NA
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "RECQL4_FAM135B_CDH17_PABPC1_UBR5_CSMD3_RAD21_EXT1_MYC", "gene"] <- "MYC"
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "STK11_POLRMT_CD209_ZNF429_EEF2_CREB3L3_MAP2K2_KEAP1", "gene"] <- "STK11, KEAP1"
CN_driver_sample[grep("_", CN_driver_sample$significant_CN_genes_on_peak, invert = T), "gene"] <- CN_driver_sample[grep("_", CN_driver_sample$significant_CN_genes_on_peak, invert = T), "significant_CN_genes_on_peak"]
CN_driver_sample[nchar(CN_driver_sample$significant_CN_genes_on_peak) <= 20, "gene"] <- CN_driver_sample[nchar(CN_driver_sample$significant_CN_genes_on_peak) <= 20, "significant_CN_genes_on_peak"]
CN_driver_sample$short_peak_name <- sub("PANCAN_", "", CN_driver_sample$new_peak_name)
CN_driver_sample$short_peak_name <- sub("_p.*", "", CN_driver_sample$short_peak_name)
CN_driver_sample$short_peak_name <- sub("_q.*", "", CN_driver_sample$short_peak_name)

CN_driver_sample$short_name <- paste0(CN_driver_sample$short_peak_name, ":", CN_driver_sample$gene)
CN_driver_sample$short_name <- sub(":NA", "", CN_driver_sample$short_name)

CN_driver <- unique(CN_driver_sample[, c("patient", "short_name")])
colnames(CN_driver) <- c("participant_id", "CN")

# load the SV drivers
SV_driver <- read.table(SV_driver_mut_path, head = T, sep = "\t")
SV_driver <- unique(SV_driver[, c("participant_id", "gene_name", "element_type", "SV_CLASS")])
colnames(SV_driver) <- c("participant_id", "SV", "SV_gene_element", "SV_class")

# load the fusions
fusions <- read.table(fusion_path, head = T, sep = "\t")
fusions <- unique(fusions[, c("sample", "symbol1", "symbol2")])
fusions <- fusions[fusions$symbol1 != fusions$symbol2, ]
fusions <- fusions[fusions$symbol2 != "NRG1-IT3", ]
fusions$symbol2 <- NULL
colnames(fusions) <- c("participant_id", "fusion")

# load the germline variants
germline <- read.table(germline_path, head = T, sep = "\t")

germline_gene_split <- list()
for(i in 1:nrow(germline)){
  print(i)
  df <- germline[i, ]
  
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
                         germline_driver_LOH_genes = germline_loh)
    
  } else {df_out <- df}
  
  germline_gene_split[[i]] <- df_out
}
germline_gene_split <- do.call(rbind, germline_gene_split)

# indicate whether there is actually a second hit
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_driver_LOH_genes), "second_hit"] <- "TRUE"
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_genes_somatically_mutated), "second_hit"] <- "TRUE"
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_genes_somatic_SV), "second_hit"] <- "TRUE"
germline_gene_split[is.na(germline_gene_split$second_hit), "second_hit"] <- "FALSE"

germline <- unique(germline_gene_split[, c("patient", "germline_driver_gene", "second_hit")])
colnames(germline) <- c("participant_id", "germline", "germline_second_hit")

# make this all into one dataframe
all <- sample_table[, "participant_id", drop = F ]

all <- left_join(all, SNV_driver)
all <- left_join(all, CN_driver)
all <- left_join(all, SV_driver)
all <- left_join(all, fusions)
all <- left_join(all, germline)
all$participant_id <- as.character(all$participant_id)

# make a patient order
# make the patient order based on coding drivers and block-wise per gene
sample_table$participant_id <- as.character(sample_table$participant_id)


all <- left_join(all, sample_table[, c("participant_id", "histology")])
patient_order <- unique(all[, c("participant_id", "coding", "histology")])
patient_order <- na.omit(patient_order)

gene_counts <- data.frame(table(patient_order$coding))
gene_counts <- gene_counts[order(gene_counts$Freq, decreasing = T), "Var1"]

patient_order$coding <- factor(patient_order$coding, levels = gene_counts)
patient_order <- patient_order %>% arrange(coding, histology)
patient_order <- unique(patient_order$participant_id)

# now order the tumours with only non-coding drivers in a similar way
patient_order_nc <- unique(all[, c("participant_id", "non_coding", "histology")])
patient_order_nc <- na.omit(patient_order_nc)
patient_order_nc <- patient_order_nc[!patient_order_nc$participant_id %in% patient_order, ]

gene_counts_nc <- data.frame(table(patient_order_nc$non_coding))
gene_counts_nc <- gene_counts_nc[order(gene_counts_nc$Freq, decreasing = T), "Var1"]

patient_order_nc$non_coding <- factor(patient_order_nc$non_coding, levels = gene_counts_nc)
patient_order_nc <- patient_order_nc %>% arrange(non_coding, histology)
patient_order_nc <- unique(patient_order_nc$participant_id)

patient_order <- c(patient_order, patient_order_nc)

# CN
patient_order_cn <- unique(all[, c("participant_id", "CN", "histology")])
patient_order_cn <- na.omit(patient_order_cn)
patient_order_cn <- patient_order_cn[!patient_order_cn$participant_id %in% patient_order, ]

gene_counts_cn <- data.frame(table(patient_order_cn$CN))
gene_counts_cn <- gene_counts_cn[order(gene_counts_cn$Freq, decreasing = T), "Var1"]

patient_order_cn$CN <- factor(patient_order_cn$CN, levels = gene_counts_cn)
patient_order_cn <- patient_order_cn %>% arrange(CN, histology)
patient_order_cn <- unique(patient_order_cn$participant_id)

patient_order <- c(patient_order, patient_order_cn)

# SV
patient_order_sv <- unique(all[, c("participant_id", "SV", "histology")])
patient_order_sv <- na.omit(patient_order_sv)
patient_order_sv <- patient_order_sv[!patient_order_sv$participant_id %in% patient_order, ]

gene_counts_sv <- data.frame(table(patient_order_sv$SV))
gene_counts_sv <- gene_counts_sv[order(gene_counts_sv$Freq, decreasing = T), "Var1"]

patient_order_sv$SV <- factor(patient_order_sv$SV, levels = gene_counts_sv)
patient_order_sv <- patient_order_sv %>% arrange(SV, histology)
patient_order_sv <- unique(patient_order_sv$participant_id)

patient_order <- c(patient_order, patient_order_sv)


# fusion
patient_order_fusion <- unique(all[, c("participant_id", "fusion", "histology")])
patient_order_fusion <- na.omit(patient_order_fusion)
patient_order_fusion <- patient_order_fusion[!patient_order_fusion$participant_id %in% patient_order, ]

gene_counts_fusion <- data.frame(table(patient_order_fusion$fusion))
gene_counts_fusion <- gene_counts_fusion[order(gene_counts_fusion$Freq, decreasing = T), "Var1"]

patient_order_fusion$fusion <- factor(patient_order_fusion$fusion, levels = gene_counts_fusion)
patient_order_fusion <- patient_order_fusion %>% arrange(fusion, histology)
patient_order_fusion <- unique(patient_order_fusion$participant_id)

patient_order <- c(patient_order, patient_order_fusion)

# germline
patient_order_germline <- unique(all[, c("participant_id", "germline", "histology")])
patient_order_germline <- na.omit(patient_order_germline)
patient_order_germline <- patient_order_germline[!patient_order_germline$participant_id %in% patient_order, ]

gene_counts_germline <- data.frame(table(patient_order_germline$germline))
gene_counts_germline <- gene_counts_germline[order(gene_counts_germline$Freq, decreasing = T), "Var1"]

patient_order_germline$germline <- factor(patient_order_germline$germline, levels = gene_counts_germline)
patient_order_germline <- patient_order_germline %>% arrange(germline, histology)
patient_order_germline <- unique(patient_order_germline$participant_id)

patient_order <- c(patient_order, patient_order_germline)

# now add the last ones that truly dont have a driver
missing <- sample_table[!sample_table$participant_id %in% patient_order, c("participant_id", "histology")]
missing <- missing %>% arrange(histology)
patient_order <- c(patient_order, missing$participant_id)

# plot

# start with coding
coding_df <- unique(all[, c("participant_id", "coding")])
coding_df <- na.omit(coding_df)
coding_df$driver <- "TRUE"

coding_missing <- sample_table[!sample_table$participant_id %in% coding_df$participant_id, "participant_id"]

# add the missing ones
df_add <- data.frame(participant_id = coding_missing,
                     coding = "TP53",
                     driver = "FALSE")

coding_df <- rbind(coding_df, df_add)

# order by how many times a gene is a driver
coding_gene_order <- data.frame(table(coding_df$coding))
coding_gene_order <- coding_gene_order[order(coding_gene_order$Freq, decreasing = T), "Var1"]

coding_df <- dcast(coding_df, participant_id ~ coding)
coding_df[is.na(coding_df)] <- "FALSE"
coding_df <- melt(coding_df, id.vars = "participant_id")
coding_df$participant_id <- factor(coding_df$participant_id, levels = patient_order)
coding_df$value <- factor(coding_df$value, levels = c("TRUE", "FALSE"))
coding_df$variable <- factor(coding_df$variable, levels = rev(coding_gene_order))

p_coding <- ggplot(coding_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  theme_bw() +
  ylab("coding SNVs & indels") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

# now non-coding
noncoding_df <- unique(all[, c("participant_id", "non_coding")])
noncoding_df <- na.omit(noncoding_df)
noncoding_df$driver <- "TRUE"

noncoding_missing <- sample_table[!sample_table$participant_id %in% noncoding_df$participant_id, "participant_id"]


# add the missing ones
df_add <- data.frame(participant_id = noncoding_missing,
                     non_coding = "TP53:ss",
                     driver = "FALSE")

noncoding_df <- rbind(noncoding_df, df_add)

# order by how many times a gene is a driver
noncoding_gene_order <- data.frame(table(noncoding_df$non_coding))
noncoding_gene_order <- noncoding_gene_order[order(noncoding_gene_order$Freq, decreasing = T), "Var1"]

noncoding_df <- dcast(noncoding_df, participant_id ~ non_coding)
noncoding_df[is.na(noncoding_df)] <- "FALSE"
noncoding_df <- melt(noncoding_df, id.vars = "participant_id")
noncoding_df$participant_id <- factor(noncoding_df$participant_id, levels = patient_order)
noncoding_df$value <- factor(noncoding_df$value, levels = c("TRUE", "FALSE"))
noncoding_df$variable <- factor(noncoding_df$variable, levels = rev(noncoding_gene_order))

p_noncoding <- ggplot(noncoding_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  theme_bw() +
  ylab("non-coding SNVs & indels") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank())

# CN
CN_df <- unique(all[, c("participant_id", "CN")])
CN_df <- na.omit(CN_df)
CN_df$driver <- "TRUE"

CN_missing <- sample_table[!sample_table$participant_id %in% CN_df$participant_id, "participant_id"]


# add the missing ones
df_add <- data.frame(participant_id = CN_missing,
                     CN = "homozygous_deletion_chr9:CDKN2A",
                     driver = "FALSE")

CN_df <- rbind(CN_df, df_add)

CN_df$CN <- sub("amplification", "amp", CN_df$CN)
CN_df$CN <- sub("homozygous_deletion", "homdel", CN_df$CN)

# order by how many times a gene is a driver
CN_gene_order <- data.frame(table(CN_df$CN))
CN_gene_order <- CN_gene_order[order(CN_gene_order$Freq, decreasing = T), "Var1"]

CN_df <- dcast(CN_df, participant_id ~ CN)
CN_df[is.na(CN_df)] <- "FALSE"
CN_df <- melt(CN_df, id.vars = "participant_id")
CN_df$participant_id <- factor(CN_df$participant_id, levels = patient_order)
CN_df$value <- factor(CN_df$value, levels = c("TRUE", "FALSE"))
CN_df$variable <- factor(CN_df$variable, levels = rev(CN_gene_order))

p_CN <- ggplot(CN_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  ylab("CN amplifications / homozygous deletions") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank())

# SV
SV_df <- unique(all[, c("participant_id", "SV")])
SV_df <- na.omit(SV_df)
SV_df$driver <- "TRUE"

SV_missing <- sample_table[!sample_table$participant_id %in% SV_df$participant_id, "participant_id"]

# add the missing ones
df_add <- data.frame(participant_id = SV_missing,
                     SV = "CDKN2A",
                     driver = "FALSE")

SV_df <- rbind(SV_df, df_add)

# order by how many times a gene is a driver
SV_gene_order <- data.frame(table(SV_df$SV))
SV_gene_order <- SV_gene_order[order(SV_gene_order$Freq, decreasing = T), "Var1"]

SV_df$SV <- factor(SV_df$SV, levels = rev(SV_gene_order))
SV_df$driver <- factor(SV_df$driver, levels = c("TRUE", "FALSE"))
SV_df$participant_id <- factor(SV_df$participant_id, levels = patient_order)

p_SV <- ggplot(SV_df, aes(participant_id, SV, fill = driver)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        legend.position = "none")

# make a germline plot
germline_df <- unique(all[, c("participant_id", "germline", "germline_second_hit")])
germline_df[!is.na(germline_df$germline), "driver"] <- "TRUE"
germline_df[is.na(germline_df$germline), "driver"] <- "FALSE"

# order by how many times a gene is a driver
germline_gene_order <- data.frame(table(germline_df$germline))
germline_gene_order <- germline_gene_order[order(germline_gene_order$Freq, decreasing = T), "Var1"]

germline_df[is.na(germline_df$germline), "germline"] <- "CHEK2"

germline_df[which(germline_df$driver == "TRUE" & germline_df$germline_second_hit == "FALSE"), "driver"] <- "germline_mutation"
germline_df[which(germline_df$driver == "TRUE" & germline_df$germline_second_hit == "TRUE"), "driver"] <- "germline_mutation_somatic_LOH"

germline_df$participant_id <- factor(germline_df$participant_id, levels = patient_order)
germline_df$driver <- factor(germline_df$driver, levels = c("germline_mutation", "germline_mutation_somatic_LOH", "FALSE"))
germline_df$germline <- factor(germline_df$germline, levels = rev(germline_gene_order))

p_germline <- ggplot(germline_df, aes(participant_id, germline, fill = driver)) +
  geom_tile() +
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3", "white"), labels = c("germline mutation", "germline mutation & somatic LOH", ""), name = "germline driver", na.value = NA) +
  ylab("germline") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank())

# fusions
fusion_df <- unique(all[, c("participant_id", "fusion")])
fusion_df[!is.na(fusion_df$fusion), "driver"] <- "TRUE"
fusion_df[is.na(fusion_df$fusion), "driver"] <- "FALSE"

# order by how many times a gene is a driver
fusion_gene_order <- data.frame(table(fusion_df$fusion))
fusion_gene_order <- fusion_gene_order[order(fusion_gene_order$Freq, decreasing = T), "Var1"]

fusion_df[is.na(fusion_df$fusion), "fusion"] <- "ALK"

fusion_df$participant_id <- factor(fusion_df$participant_id, levels = patient_order)
fusion_df$driver <- factor(fusion_df$driver, levels = c("TRUE", "FALSE"))
fusion_df$fusion <- factor(fusion_df$fusion, levels = rev(fusion_gene_order))

p_fusion <- ggplot(fusion_df, aes(participant_id, fusion, fill = driver)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  ylab("fusion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        legend.position = "none")

# add a plot with histology
hist_df <- unique(all[, c("participant_id", "histology")])
hist_df$participant_id <- factor(hist_df$participant_id, levels = patient_order)

p_hist <- ggplot(hist_df, aes(participant_id, 1, fill = histology)) +
  geom_tile() +
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
                               "OTHER" = "#7a7979"),
                    labels = c("Adenocarcinoma",
                               "Adenocarcinoma\nmetastasis",
                               "Squamous cell",
                               "Squamous cell\nmetastasis",
                               "Adenosquamous",
                               "Carcinoid",
                               "Large cell",
                               "Mesothelioma",
                               "Small cell",
                               "Small cell\nmetastasis",
                               "Neuroendrocrine",
                               "Other")) +
  ylab("histology") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        panel.grid = element_blank(),)

# join plots

p_all <- plot_grid(p_coding ,p_noncoding ,p_CN ,p_SV, p_fusion, p_germline, p_hist, ncol = 1, rel_heights = c(2.5, 3, 2, 1.5, 0.5, 1.5, 0.4), align = "v", axis = "lr")

pdf(paste0(output_path, "oncoplot.pdf"), width = 12, height = 13)
p_all
dev.off()




