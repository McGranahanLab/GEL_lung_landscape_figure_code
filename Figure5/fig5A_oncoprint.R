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
mutburden_path        <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"

###############################################################################
#####################################################################     MAIN

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
mutburden    <- readRDS(mutburden_path)

# read in SNV driver mutations
# read all the driver mutation files for all histologies
all_files <- list.files(SNV_driver_path, full.names = T)

SNV_driver <- list()
for(i in all_files){
  
  print(i)
  df <- read.delim(i)
  df <- df[df$is_driver == TRUE, ]
  df <- df[df$var_class != "amp", ]
  df <- df[df$var_class != "hom.del.", ]
  df <- unique(df[, c("participant_id", "gene_name", "gr_id", "key")])
  
  SNV_driver[[i]] <- df
}

SNV_driver <- do.call(rbind, SNV_driver)
rownames(SNV_driver) <- c(1:nrow(SNV_driver))
SNV_driver <- unique(SNV_driver)
SNV_driver$gene_gr <- paste0(SNV_driver$gene_name, ":", SNV_driver$gr_id)

# need to reformat this
SNV_driver[SNV_driver$gr_id == "CDS", "coding"] <- SNV_driver[SNV_driver$gr_id == "CDS", "gene_name"]
SNV_driver[SNV_driver$gr_id == "3primeUTR", "3primeUTR"] <- SNV_driver[SNV_driver$gr_id == "3primeUTR", "gene_gr"]
SNV_driver[SNV_driver$gr_id == "5primeUTR", "5primeUTR"] <- SNV_driver[SNV_driver$gr_id == "5primeUTR", "gene_gr"]
SNV_driver[SNV_driver$gr_id == "promoter", "promoter"] <- SNV_driver[SNV_driver$gr_id == "promoter", "gene_gr"]
SNV_driver[SNV_driver$gr_id == "ss", "ss"] <- SNV_driver[SNV_driver$gr_id == "ss", "gene_gr"]
SNV_driver[SNV_driver$gr_id == "enhancer", "enhancer"] <- SNV_driver[SNV_driver$gr_id == "enhancer", "gene_gr"]
SNV_driver[SNV_driver$gr_id == "lincRNA_promoter", "lincRNA_promoter"] <- SNV_driver[SNV_driver$gr_id == "lincRNA_promoter", "gene_gr"]
SNV_driver[SNV_driver$gr_id == "lincRNA_ss", "lincRNA_ss"] <- SNV_driver[SNV_driver$gr_id == "lincRNA_ss", "gene_gr"]
SNV_driver <- unique(SNV_driver[, c("participant_id", "coding", "3primeUTR", "5primeUTR", "promoter", "ss", "enhancer", "lincRNA_promoter", "lincRNA_ss")])

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
# 3primeUTR
patient_order_3prime <- unique(all[, c("participant_id", "3primeUTR", "histology")])
patient_order_3prime <- na.omit(patient_order_3prime)
patient_order_3prime <- patient_order_3prime[!patient_order_3prime$participant_id %in% patient_order, ]

gene_counts_3prime <- data.frame(table(patient_order_3prime$`3primeUTR`))
gene_counts_3prime <- gene_counts_3prime[order(gene_counts_3prime$Freq, decreasing = T), "Var1"]

patient_order_3prime$`3primeUTR` <- factor(patient_order_3prime$`3primeUTR`, levels = gene_counts_3prime)
patient_order_3prime <- patient_order_3prime %>% arrange(`3primeUTR`, histology)
patient_order_3prime <- unique(patient_order_3prime$participant_id)

patient_order <- c(patient_order, patient_order_3prime)

# 5primeUTR
patient_order_5prime <- unique(all[, c("participant_id", "5primeUTR", "histology")])
patient_order_5prime <- na.omit(patient_order_5prime)
patient_order_5prime <- patient_order_5prime[!patient_order_5prime$participant_id %in% patient_order, ]

gene_counts_5prime <- data.frame(table(patient_order_5prime$`5primeUTR`))
gene_counts_5prime <- gene_counts_5prime[order(gene_counts_5prime$Freq, decreasing = T), "Var1"]

patient_order_5prime$`5primeUTR` <- factor(patient_order_5prime$`5primeUTR`, levels = gene_counts_5prime)
patient_order_5prime <- patient_order_5prime %>% arrange(`5primeUTR`, histology)
patient_order_5prime <- unique(patient_order_5prime$participant_id)

patient_order <- c(patient_order, patient_order_5prime)

# promoter
patient_order_promoter <- unique(all[, c("participant_id", "promoter", "histology")])
patient_order_promoter <- na.omit(patient_order_promoter)
patient_order_promoter <- patient_order_promoter[!patient_order_promoter$participant_id %in% patient_order, ]

gene_counts_promoter <- data.frame(table(patient_order_promoter$`promoter`))
gene_counts_promoter <- gene_counts_promoter[order(gene_counts_promoter$Freq, decreasing = T), "Var1"]

patient_order_promoter$`promoter` <- factor(patient_order_promoter$`promoter`, levels = gene_counts_promoter)
patient_order_promoter <- patient_order_promoter %>% arrange(`promoter`, histology)
patient_order_promoter <- unique(patient_order_promoter$participant_id)

patient_order <- c(patient_order, patient_order_promoter)

# ss
patient_order_ss <- unique(all[, c("participant_id", "ss", "histology")])
patient_order_ss <- na.omit(patient_order_ss)
patient_order_ss <- patient_order_ss[!patient_order_ss$participant_id %in% patient_order, ]

gene_counts_ss <- data.frame(table(patient_order_ss$`ss`))
gene_counts_ss <- gene_counts_ss[order(gene_counts_ss$Freq, decreasing = T), "Var1"]

patient_order_ss$`ss` <- factor(patient_order_ss$`ss`, levels = gene_counts_ss)
patient_order_ss <- patient_order_ss %>% arrange(`ss`, histology)
patient_order_ss <- unique(patient_order_ss$participant_id)

patient_order <- c(patient_order, patient_order_ss)

# enhancer
patient_order_enhancer <- unique(all[, c("participant_id", "enhancer", "histology")])
patient_order_enhancer <- na.omit(patient_order_enhancer)
patient_order_enhancer <- patient_order_enhancer[!patient_order_enhancer$participant_id %in% patient_order, ]

gene_counts_enhancer <- data.frame(table(patient_order_enhancer$`enhancer`))
gene_counts_enhancer <- gene_counts_enhancer[order(gene_counts_enhancer$Freq, decreasing = T), "Var1"]

patient_order_enhancer$`enhancer` <- factor(patient_order_enhancer$`enhancer`, levels = gene_counts_enhancer)
patient_order_enhancer <- patient_order_enhancer %>% arrange(`enhancer`, histology)
patient_order_enhancer <- unique(patient_order_enhancer$participant_id)

patient_order <- c(patient_order, patient_order_enhancer)

# lincRNA_promoter
patient_order_lincRNA_promoter <- unique(all[, c("participant_id", "lincRNA_promoter", "histology")])
patient_order_lincRNA_promoter <- na.omit(patient_order_lincRNA_promoter)
patient_order_lincRNA_promoter <- patient_order_lincRNA_promoter[!patient_order_lincRNA_promoter$participant_id %in% patient_order, ]

gene_counts_lincRNA_promoter <- data.frame(table(patient_order_lincRNA_promoter$`lincRNA_promoter`))
gene_counts_lincRNA_promoter <- gene_counts_lincRNA_promoter[order(gene_counts_lincRNA_promoter$Freq, decreasing = T), "Var1"]

patient_order_lincRNA_promoter$`lincRNA_promoter` <- factor(patient_order_lincRNA_promoter$`lincRNA_promoter`, levels = gene_counts_lincRNA_promoter)
patient_order_lincRNA_promoter <- patient_order_lincRNA_promoter %>% arrange(`lincRNA_promoter`, histology)
patient_order_lincRNA_promoter <- unique(patient_order_lincRNA_promoter$participant_id)

patient_order <- c(patient_order, patient_order_lincRNA_promoter)

# lincRNA_ss
patient_order_lincRNA_ss <- unique(all[, c("participant_id", "lincRNA_ss", "histology")])
patient_order_lincRNA_ss <- na.omit(patient_order_lincRNA_ss)
patient_order_lincRNA_ss <- patient_order_lincRNA_ss[!patient_order_lincRNA_ss$participant_id %in% patient_order, ]

gene_counts_lincRNA_ss <- data.frame(table(patient_order_lincRNA_ss$`lincRNA_ss`))
gene_counts_lincRNA_ss <- gene_counts_lincRNA_ss[order(gene_counts_lincRNA_ss$Freq, decreasing = T), "Var1"]

patient_order_lincRNA_ss$`lincRNA_ss` <- factor(patient_order_lincRNA_ss$`lincRNA_ss`, levels = gene_counts_lincRNA_ss)
patient_order_lincRNA_ss <- patient_order_lincRNA_ss %>% arrange(`lincRNA_ss`, histology)
patient_order_lincRNA_ss <- unique(patient_order_lincRNA_ss$participant_id)

patient_order <- c(patient_order, patient_order_lincRNA_ss)


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
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("coding SNVs & indels") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.05, 0.1), "cm"))

# now 3primeUTR
primeUTR3_df <- unique(all[, c("participant_id", "3primeUTR")])
primeUTR3_df <- na.omit(primeUTR3_df)
primeUTR3_df$`3primeUTR` <- sub(":.*", "", primeUTR3_df$`3primeUTR`)
primeUTR3_df$driver <- "TRUE"
colnames(primeUTR3_df) <- c("participant_id", "prime3UTR", "driver")

primeUTR3_missing <- sample_table[!sample_table$participant_id %in% primeUTR3_df$participant_id, "participant_id"]


# add the missing ones
df_add <- data.frame(participant_id = primeUTR3_missing,
                     "prime3UTR" = "SLC34A2",
                     driver = "FALSE")

primeUTR3_df <- rbind(primeUTR3_df, df_add)

# order by how many times a gene is a driver
primeUTR3_gene_order <- data.frame(table(primeUTR3_df$prime3UTR))
primeUTR3_gene_order <- primeUTR3_gene_order[order(primeUTR3_gene_order$Freq, decreasing = T), "Var1"]

primeUTR3_df <- dcast(primeUTR3_df, participant_id ~ prime3UTR)
primeUTR3_df[is.na(primeUTR3_df)] <- "FALSE"
primeUTR3_df <- melt(primeUTR3_df, id.vars = "participant_id")
primeUTR3_df$participant_id <- factor(primeUTR3_df$participant_id, levels = patient_order)
primeUTR3_df$value <- factor(primeUTR3_df$value, levels = c("TRUE", "FALSE"))
primeUTR3_df$variable <- factor(primeUTR3_df$variable, levels = rev(primeUTR3_gene_order))

p_primeUTR3 <- ggplot(primeUTR3_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("3'UTR") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

# now 5primeUTR
primeUTR5_df <- unique(all[, c("participant_id", "5primeUTR")])
primeUTR5_df <- na.omit(primeUTR5_df)
primeUTR5_df$`5primeUTR` <- sub(":.*", "", primeUTR5_df$`5primeUTR`)
primeUTR5_df$driver <- "TRUE"
colnames(primeUTR5_df) <- c("participant_id", "prime5UTR", "driver")

primeUTR5_missing <- sample_table[!sample_table$participant_id %in% primeUTR5_df$participant_id, "participant_id"]

# add the missing ones
df_add <- data.frame(participant_id = primeUTR5_missing,
                     "prime5UTR" = "STON2",
                     driver = "FALSE")

primeUTR5_df <- rbind(primeUTR5_df, df_add)

# order by how many times a gene is a driver
primeUTR5_gene_order <- data.frame(table(primeUTR5_df$prime5UTR))
primeUTR5_gene_order <- primeUTR5_gene_order[order(primeUTR5_gene_order$Freq, decreasing = T), "Var1"]

primeUTR5_df <- dcast(primeUTR5_df, participant_id ~ prime5UTR)
primeUTR5_df[is.na(primeUTR5_df)] <- "FALSE"
primeUTR5_df <- melt(primeUTR5_df, id.vars = "participant_id")
primeUTR5_df$participant_id <- factor(primeUTR5_df$participant_id, levels = patient_order)
primeUTR5_df$value <- factor(primeUTR5_df$value, levels = c("TRUE", "FALSE"))
primeUTR5_df$variable <- factor(primeUTR5_df$variable, levels = rev(primeUTR5_gene_order))

p_primeUTR5 <- ggplot(primeUTR5_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("5'UTR") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

# now promoter
promoter_df <- unique(all[, c("participant_id", "promoter")])
promoter_df <- na.omit(promoter_df)
promoter_df$`promoter` <- sub(":.*", "", promoter_df$`promoter`)
promoter_df$driver <- "TRUE"
colnames(promoter_df) <- c("participant_id", "promoter", "driver")

promoter_missing <- sample_table[!sample_table$participant_id %in% promoter_df$participant_id, "participant_id"]

# add the missing ones
df_add <- data.frame(participant_id = promoter_missing,
                     "promoter" = "C11orf58",
                     driver = "FALSE")

promoter_df <- rbind(promoter_df, df_add)

# order by how many times a gene is a driver
promoter_gene_order <- data.frame(table(promoter_df$promoter))
promoter_gene_order <- promoter_gene_order[order(promoter_gene_order$Freq, decreasing = T), "Var1"]

promoter_df <- dcast(promoter_df, participant_id ~ promoter)
promoter_df[is.na(promoter_df)] <- "FALSE"
promoter_df <- melt(promoter_df, id.vars = "participant_id")
promoter_df$participant_id <- factor(promoter_df$participant_id, levels = patient_order)
promoter_df$value <- factor(promoter_df$value, levels = c("TRUE", "FALSE"))
promoter_df$variable <- factor(promoter_df$variable, levels = rev(promoter_gene_order))

p_promoter <- ggplot(promoter_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("promoter") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

# now ss
ss_df <- unique(all[, c("participant_id", "ss")])
ss_df <- na.omit(ss_df)
ss_df$`ss` <- sub(":.*", "", ss_df$`ss`)
ss_df$driver <- "TRUE"
colnames(ss_df) <- c("participant_id", "ss", "driver")

ss_missing <- sample_table[!sample_table$participant_id %in% ss_df$participant_id, "participant_id"]

# add the missing ones
df_add <- data.frame(participant_id = ss_missing,
                     "ss" = "TP53",
                     driver = "FALSE")

ss_df <- rbind(ss_df, df_add)

# order by how many times a gene is a driver
ss_gene_order <- data.frame(table(ss_df$ss))
ss_gene_order <- ss_gene_order[order(ss_gene_order$Freq, decreasing = T), "Var1"]

ss_df <- dcast(ss_df, participant_id ~ ss)
ss_df[is.na(ss_df)] <- "FALSE"
ss_df <- melt(ss_df, id.vars = "participant_id")
ss_df$participant_id <- factor(ss_df$participant_id, levels = patient_order)
ss_df$value <- factor(ss_df$value, levels = c("TRUE", "FALSE"))
ss_df$variable <- factor(ss_df$variable, levels = rev(ss_gene_order))

p_ss <- ggplot(ss_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("ss") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

# now enhancer
enhancer_df <- unique(all[, c("participant_id", "enhancer")])
enhancer_df <- na.omit(enhancer_df)
enhancer_df$`enhancer` <- sub(":.*", "", enhancer_df$`enhancer`)
enhancer_df$driver <- "TRUE"
colnames(enhancer_df) <- c("participant_id", "enhancer", "driver")

enhancer_missing <- sample_table[!sample_table$participant_id %in% enhancer_df$participant_id, "participant_id"]

# add the mienhancering ones
df_add <- data.frame(participant_id = enhancer_mienhancering,
                     "enhancer" = "ZNF578",
                     driver = "FALSE")

enhancer_df <- rbind(enhancer_df, df_add)

# order by how many times a gene is a driver
enhancer_gene_order <- data.frame(table(enhancer_df$enhancer))
enhancer_gene_order <- enhancer_gene_order[order(enhancer_gene_order$Freq, decreasing = T), "Var1"]

enhancer_df <- dcast(enhancer_df, participant_id ~ enhancer)
enhancer_df[is.na(enhancer_df)] <- "FALSE"
enhancer_df <- melt(enhancer_df, id.vars = "participant_id")
enhancer_df$participant_id <- factor(enhancer_df$participant_id, levels = patient_order)
enhancer_df$value <- factor(enhancer_df$value, levels = c("TRUE", "FALSE"))
enhancer_df$variable <- factor(enhancer_df$variable, levels = rev(enhancer_gene_order))

p_enhancer <- ggplot(enhancer_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("enhancer") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

# now lincRNA_promoter
lincRNA_promoter_df <- unique(all[, c("participant_id", "lincRNA_promoter")])
lincRNA_promoter_df <- na.omit(lincRNA_promoter_df)
lincRNA_promoter_df$`lincRNA_promoter` <- sub(":.*", "", lincRNA_promoter_df$`lincRNA_promoter`)
lincRNA_promoter_df$driver <- "TRUE"
colnames(lincRNA_promoter_df) <- c("participant_id", "lincRNA_promoter", "driver")

lincRNA_promoter_missing <- sample_table[!sample_table$participant_id %in% lincRNA_promoter_df$participant_id, "participant_id"]

# add the milincRNA_promotering ones
df_add <- data.frame(participant_id = lincRNA_promoter_missing,
                     "lincRNA_promoter" = "LOC100129917",
                     driver = "FALSE")

lincRNA_promoter_df <- rbind(lincRNA_promoter_df, df_add)

# order by how many times a gene is a driver
lincRNA_promoter_gene_order <- data.frame(table(lincRNA_promoter_df$lincRNA_promoter))
lincRNA_promoter_gene_order <- lincRNA_promoter_gene_order[order(lincRNA_promoter_gene_order$Freq, decreasing = T), "Var1"]

lincRNA_promoter_df <- dcast(lincRNA_promoter_df, participant_id ~ lincRNA_promoter)
lincRNA_promoter_df[is.na(lincRNA_promoter_df)] <- "FALSE"
lincRNA_promoter_df <- melt(lincRNA_promoter_df, id.vars = "participant_id")
lincRNA_promoter_df$participant_id <- factor(lincRNA_promoter_df$participant_id, levels = patient_order)
lincRNA_promoter_df$value <- factor(lincRNA_promoter_df$value, levels = c("TRUE", "FALSE"))
lincRNA_promoter_df$variable <- factor(lincRNA_promoter_df$variable, levels = rev(lincRNA_promoter_gene_order))

p_lincRNA_promoter <- ggplot(lincRNA_promoter_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("lincRNA promoter") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

# now lincRNA_ss
lincRNA_ss_df <- unique(all[, c("participant_id", "lincRNA_ss")])
lincRNA_ss_df <- na.omit(lincRNA_ss_df)
lincRNA_ss_df$`lincRNA_ss` <- sub(":.*", "", lincRNA_ss_df$`lincRNA_ss`)
lincRNA_ss_df$driver <- "TRUE"
colnames(lincRNA_ss_df) <- c("participant_id", "lincRNA_ss", "driver")

lincRNA_ss_missing <- sample_table[!sample_table$participant_id %in% lincRNA_ss_df$participant_id, "participant_id"]

# add the milincRNA_ssing ones
df_add <- data.frame(participant_id = lincRNA_ss_missing,
                     "lincRNA_ss" = "CHL1-AS1",
                     driver = "FALSE")

lincRNA_ss_df <- rbind(lincRNA_ss_df, df_add)

# order by how many times a gene is a driver
lincRNA_ss_gene_order <- data.frame(table(lincRNA_ss_df$lincRNA_ss))
lincRNA_ss_gene_order <- lincRNA_ss_gene_order[order(lincRNA_ss_gene_order$Freq, decreasing = T), "Var1"]

lincRNA_ss_df <- dcast(lincRNA_ss_df, participant_id ~ lincRNA_ss)
lincRNA_ss_df[is.na(lincRNA_ss_df)] <- "FALSE"
lincRNA_ss_df <- melt(lincRNA_ss_df, id.vars = "participant_id")
lincRNA_ss_df$participant_id <- factor(lincRNA_ss_df$participant_id, levels = patient_order)
lincRNA_ss_df$value <- factor(lincRNA_ss_df$value, levels = c("TRUE", "FALSE"))
lincRNA_ss_df$variable <- factor(lincRNA_ss_df$variable, levels = rev(lincRNA_ss_gene_order))

p_lincRNA_ss <- ggplot(lincRNA_ss_df, aes(participant_id, variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  ylab("lincRNA ss") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))


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
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

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
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

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
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

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
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))

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
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))

# add a plot with TMB
burden_df <- mutburden[, c("sample", "SBScount")]
burden_df$SBScount  <- burden_df$SBScount/2800
burden_df$participant_id <- factor(burden_df$sample, levels = patient_order)

p_tmb <- ggplot(burden_df, aes(sample, 1, fill = SBScount)) +
  geom_tile() +
  scale_fill_gradient(name = "SNVs/Mb", low = "white", high = "#ae017e") +
  ylab("SNVs/Mb") +
  theme_bw() +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.1, 0.05, 0.1), "cm"))


# join plots

p_all <- plot_grid(p_coding, p_ss, p_enhancer, p_primeUTR3, p_primeUTR5, p_promoter, p_lincRNA_promoter, p_lincRNA_ss, p_CN, p_SV, p_fusion, p_germline, p_tmb, p_hist,
                   ncol = 1,
                   rel_heights = c(3, 1.5, 0.25, 0.25, 0.7, 1.5, 0.15, 0.15, 2, 1.8, 0.4, 1.7, 0.2, 0.2), align = "v", axis = "lr")

pdf(paste0(output_path, "oncoplot.pdf"), width = 12, height = 13)
p_all
dev.off()




