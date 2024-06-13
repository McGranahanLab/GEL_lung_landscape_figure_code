# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #           barplot / point range number of driver mutations          # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
####################################################################  libraries

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
library(cowplot)

###############################################################################
#################################################################### file paths

sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
SNV_driver_path       <- "/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/2023_12_25/results/tables/driver_mutations/"
SV_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/SV_driver_genes_in_samples.txt"
CN_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copynumber_peak_driver_per_patient.txt"
output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/"
fusion_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/SV/fusions.txt"
germline_path         <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/germline/germline_second_hit.txt"
SNV_driver_count_path <- "/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/figures_for_ms/2023_12_25/f5b_n_coding_vs_noncoding_driver_muts.csv"

driver_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/all_driver_binary.txt"

###############################################################################
####################################################################      MAIN

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

# read the estimated numbers of coding and non coding drivers

data <- read.table("/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/figures_for_ms/2023_12_25/f5b_n_coding_vs_noncoding_driver_muts.csv", sep = '\t')
colnames(data) <- c("coding", "non_coding", "num_tumours", "bracket") 
data$sum_driver <- data$coding + data$non_coding
data$sum_driver_cohort <- data$sum_driver * data$num_tumours
data$sum_coding_driver_cohort <- data$coding * data$num_tumours
data$sum_noncoding_driver_cohort <- data$non_coding * data$num_tumours

sum_coding_driver_total <- sum(data$sum_coding_driver_cohort)
average_coding <- sum_coding_driver_total/1011

sum_noncoding_driver_total <- sum(data$sum_noncoding_driver_cohort)
average_noncoding <- sum_noncoding_driver_total/1011

# whats the percentage of tumours with these number of drivers

num_tumour_with_coding <- sum(data[data$coding > 0, "num_tumours"])
percentage_tumour_with_coding <- (num_tumour_with_coding/1011)*100
num_tumour_with_noncoding <- sum(data[data$non_coding > 0, "num_tumours"])
percentage_tumour_with_noncoding <- (num_tumour_with_noncoding/1011)*100

# get the other mutation types
driver_table <- read.table(driver_path, head = T, sep = "\t")

CN <- (nrow(driver_table[which(driver_table$CNdriver == "CN"), ]) / 1011) * 100
SV <- (nrow(driver_table[which(driver_table$SVdriver == "SV"), ]) / 1011) * 100
fusion <- (nrow(driver_table[which(driver_table$fusion == "fusion"), ]) / 1011) * 100
germline <- (nrow(driver_table[which(driver_table$germline == "germline"), ]) / 1011) * 100
any <- (nrow(driver_table[!is.na(driver_table$class), ]) / 1011) * 100

driver_counts <- data.frame(class = c("any", "coding", "noncoding", "CN", "SV", "fusion", "germline"),
                            prop = c(any, percentage_tumour_with_coding, percentage_tumour_with_noncoding, CN, SV, fusion, germline))
                            
driver_counts <- driver_counts[order(driver_counts$prop, decreasing = T), ]
driver_counts$class <- factor(driver_counts$class, levels = rev(driver_counts$class))

p_bar <- ggplot(driver_counts, aes(prop, class, label = paste0(round(prop, 1)))) +
  geom_bar(stat = "identity", fill = "black") +
  geom_text(hjust = 1.1, color = "white", size = 2) +
  scale_y_discrete(labels = rev(c("all", "coding", "SCNA", "non-coding", "SV", "fusion", "germline"))) +
  ylab("lung cancer gene mutation type") +
  xlab("% tumours") +
  theme_bw() +
  theme(text = element_text(size = 8))

# make a point range plot for average number of mutations
# for coding and non-coding we need to make the table back into a tumour wise table so we can calculate the SD
SNV_tumour_df <- data[, c("coding", "non_coding", "num_tumours")]
SNV_tumour_df$class <- paste0(SNV_tumour_df$coding, "_", SNV_tumour_df$non_coding)

foo <- rep(SNV_tumour_df$class, SNV_tumour_df$num_tumours)
foo <- data.frame(foo)
foo <- foo %>% separate(foo, into = c("coding", "non_coding"), sep = "_")
foo$coding <- as.numeric(as.character(foo$coding))
foo$non_coding <- as.numeric(as.character(foo$non_coding))

sd_driver_count <- data.frame(apply(foo, 2, function(x){sd(x)}))
sd_driver_count$type <- rownames(sd_driver_count)
colnames(sd_driver_count)[1] <- "sd"

mean_driver_count <- data.frame(apply(foo, 2, function(x){mean(x)}))
mean_driver_count$type <- rownames(mean_driver_count)
colnames(mean_driver_count)[1] <- "mean"

# add this info for the other types

all_files <- list.files(SNV_driver_path, full.names = T)
all_files <- all_files[grep("PerPatient", all_files, invert = T)]
all_files <- all_files[grep("NoncodingDriversOnly", all_files, invert = T)]

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

# how many coding / non-coding driver mutations does each tumour have?
SNV_driver[SNV_driver$gr_id == "CDS", "class"] <- "coding" 
SNV_driver[SNV_driver$gr_id != "CDS", "class"] <- "non_coding" 

SNV_driver[SNV_driver$class == "coding", "coding_SNVdriver"] <- TRUE
SNV_driver[SNV_driver$class == "non_coding", "noncoding_SNVdriver"] <- TRUE

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

# SVs 
SV_driver <- read.table(SV_driver_mut_path, head = T, sep = "\t")
SV_driver <- unique(SV_driver[, c("participant_id", "gene_name")])
SV_driver_count <- data.frame(table(SV_driver$participant_id))
colnames(SV_driver_count) <- c("participant_id", "SV")

# CN
#just use pancan peaks for this
CN_driver_sample <- read.table(CN_driver_mut_path, head = T, sep = "\t")
CN_driver_sample <- CN_driver_sample[grep("PANCAN", CN_driver_sample$new_peak_name), ]
CN_driver_sample <- unique(CN_driver_sample[, c("patient", "new_peak_name")])
CN_driver_sample <- data.frame(table(CN_driver_sample$patient))
colnames(CN_driver_sample) <- c("participant_id", "CN")

# fusion
fusions <- read.table(fusion_path, head = T, sep = "\t")
fusions <- unique(fusions[, c("sample", "symbol1", "symbol2")])
fusions <- fusions[fusions$symbol1 != fusions$symbol2, ]
fusions <- fusions[fusions$symbol2 != "NRG1-IT3", ]
fusions <- unique(fusions[, c("sample", "symbol1")])
fusion_count <- data.frame(table(fusions$sample))
colnames(fusion_count) <- c("participant_id", "fusion")

# germline
germline <- read.table(germline_path, head = T, sep = "\t")
germline <- unique(germline[, c("patient", "germline_driver_present", "germline_driver_LOH_present")])
germline <- germline[germline$germline_driver_LOH_present, ]
germline$germline_driver_present <- NULL
germline$germline_driver_LOH_present <- NULL
germline$germline <- TRUE
colnames(germline) <- c("participant_id", "germline")
germline$participant_id <- as.character(germline$participant_id)
germline_count <- data.frame(table(germline$participant_id))
colnames(germline_count) <- c("participant_id", "germline")

driver_count_mat <- left_join(driver_count_mat, SV_driver_count)
driver_count_mat <- left_join(driver_count_mat, CN_driver_sample)
driver_count_mat <- left_join(driver_count_mat, fusion_count)
driver_count_mat <- left_join(driver_count_mat, germline_count)
driver_count_mat[is.na(driver_count_mat$SV), "SV"] <- 0
driver_count_mat[is.na(driver_count_mat$CN), "CN"] <- 0
driver_count_mat[is.na(driver_count_mat$fusion), "fusion"] <- 0
driver_count_mat[is.na(driver_count_mat$germline), "germline"] <- 0
driver_count_sample <- driver_count_mat
driver_count_mat$participant_id <- NULL

driver_count_mat$any <- driver_count_mat$coding + driver_count_mat$non_coding + driver_count_mat$SV + driver_count_mat$CN + driver_count_mat$fusion + driver_count_mat$germline

driver_count_mat$coding <- NULL
driver_count_mat$non_coding <- NULL


sd_other_driver_count <- data.frame(apply(driver_count_mat, 2, function(x){sd(x)}))
sd_other_driver_count$type <- rownames(sd_other_driver_count)
colnames(sd_other_driver_count)[1] <- "sd"

sd_driver_count <- rbind(sd_driver_count, sd_other_driver_count)

mean_other_driver_count <- data.frame(apply(driver_count_mat, 2, function(x){mean(x)}))
mean_other_driver_count$type <- rownames(mean_other_driver_count)
colnames(mean_other_driver_count)[1] <- "mean"

# calculate quantiles
loquant_other_driver_count <- data.frame(apply(driver_count_mat, 2, function(x){quantile(x, c(0,0.25,0.5,0.75,1))[2]}))
loquant_other_driver_count$type <- rownames(loquant_other_driver_count)
colnames(loquant_other_driver_count)[1] <- "lo_quantile"

upquant_other_driver_count <- data.frame(apply(driver_count_mat, 2, function(x){quantile(x, c(0,0.25,0.5,0.75,1))[4]}))
upquant_other_driver_count$type <- rownames(upquant_other_driver_count)
colnames(upquant_other_driver_count)[1] <- "up_quantile"

x <- left_join(mean_other_driver_count, loquant_other_driver_count)
x <- left_join(x, upquant_other_driver_count)

mean_driver_count <- rbind(mean_driver_count, mean_other_driver_count)

mean_driver_count <- left_join(mean_driver_count, sd_driver_count)
mean_driver_count$min <- mean_driver_count$mean - mean_driver_count$sd
mean_driver_count$max <- mean_driver_count$mean + mean_driver_count$sd

sample_size <- length(unique(sample_table$participant_id))
mean_driver_count$ci <- qnorm(0.975)*mean_driver_count$sd/sqrt(sample_size)

mean_driver_count[mean_driver_count$min <0, "min"] <- 0
mean_driver_count$type <- factor(mean_driver_count$type, levels = c("germline", "fusion", "SV", "non_coding", "CN", "coding", "any"))

p_point <- ggplot(mean_driver_count, aes(x = mean, y = type, xmin = min, xmax = max, label = round(mean, 1))) +
  geom_pointrange(fatten = 2) +
  geom_text(vjust = -0.65, size = 2) +
  ylab("") +
  xlab("# cancer gene mutations") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8))


p <- plot_grid(p_bar, p_point, nrow = 1, rel_widths = c(1.5, 1))

pdf(paste0(output_path, "driver_count_barplot_pointrange.pdf"), width = 4, height = 2)
p
dev.off()


# make a bubble plot for each histology with a bar plot attached
sample_table$participant_id <- as.character(sample_table$participant_id)
driver_count_sample <- left_join(driver_count_sample, sample_table[, c("participant_id", "histology")])
driver_count_sample$any <- driver_count_sample$coding + driver_count_sample$non_coding + driver_count_sample$SV + driver_count_sample$CN + driver_count_sample$fusion

driver_count_sample_means <- driver_count_sample %>% group_by(histology) %>%
  mutate(mean_coding = mean(coding)) %>%
  mutate(median_coding = median(coding)) %>%
  mutate(median_non_coding = median(non_coding)) %>%
  mutate(mean_non_coding = mean(non_coding)) %>%
  mutate(median_SV = median(SV)) %>%
  mutate(mean_SV = mean(SV)) %>%
  mutate(median_CN = median(CN)) %>%
  mutate(mean_CN = mean(CN)) %>%
  mutate(median_fusion = median(fusion)) %>%
  mutate(mean_fusion = mean(fusion)) %>%
  mutate(median_germline = median(germline)) %>%
  mutate(mean_germline = mean(germline)) %>%
  mutate(median_any = median(any)) %>%
  mutate(mean_any = mean(any))

driver_count_sample_means <- driver_count_sample_means %>% pivot_longer(cols = c("mean_coding", "mean_non_coding", "mean_SV", "mean_CN", "mean_fusion", "mean_germline", "mean_any"), names_to = "mean_type", values_to = "mean")
driver_count_sample_means <- driver_count_sample_means %>% pivot_longer(cols = c("median_coding", "median_non_coding", "median_SV", "median_CN", "median_fusion", "median_fusion", "median_any"), names_to = "median_type", values_to = "median")

# read in the driver table

driver_count_histology_means <- unique(driver_count_sample_means[, c("histology", "mean_type", "mean", "median_type", "median")])
driver_count_histology_means <- driver_count_histology_means %>% mutate(Freq_bin = cut(mean, breaks=c(0, 0.1, 0.5, 1, 2, 3, 4, 4.3)))
driver_count_histology_means$mean_type <- factor(driver_count_histology_means$mean_type, levels = c("mean_any", "mean_coding", "mean_non_coding", "mean_CN", "mean_SV", "mean_fusion", "mean_germline"))
driver_count_histology_means <- driver_count_histology_means[order(driver_count_histology_means$mean, decreasing = T), ]
driver_count_histology_means$histology <- factor(driver_count_histology_means$histology, levels = rev(unique(driver_count_histology_means$histology)))

p_bubble <- ggplot() +
  geom_point(data = driver_count_histology_means[driver_count_histology_means$mean > 0,], aes(mean_type, histology, size = mean, color = Freq_bin)) +
  scale_color_manual(name = "mean number of\ncancer gene aberrations", values = c("#41ab5d", "#006837", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")) +
  scale_size_continuous(name = "mean number of\ncancer gene aberrations", range = c(0.4, 8)) +
  scale_y_discrete(labels = rev(c("Large cell",
                                  "Squamous cell",
                                  "Small cell",
                                  "Adenocarcinoma",
                                  "Small cell\nmetastasis",
                                  "Adenosquamous",
                                  "Squamous cell\nmetastasis",
                                  "Adenocarcinoma\nmetastasis",
                                  "Neuroendocrine",
                                  "Other",
                                  "Mesothelioma",
                                  "Carcinoid"))) +
  scale_x_discrete(labels = c("all", "coding", "non\ncoding", "CN", "SV", "fusion", "germline")) +
  xlab("") +
  ylab("") +
  theme_bw()

# and add a bar plot with proportions

driver_count_sample[driver_count_sample$any > 0, "any_driver"] <- TRUE
driver_count_sample[driver_count_sample$any == 0, "any_driver"] <- FALSE

any_driver_histology <- data.frame(table(driver_count_sample[, c("histology", "any_driver")]))
any_driver_histology <- any_driver_histology[any_driver_histology$any_driver == TRUE, ]

hist_count <- data.frame(table(sample_table$histology))
colnames(hist_count) <- c("histology", "hist_count")

any_driver_histology <- left_join(any_driver_histology, hist_count)
any_driver_histology$prop <- (any_driver_histology$Freq / any_driver_histology$hist_count)*100
any_driver_histology$histology <- factor(any_driver_histology$histology, levels = rev(unique(driver_count_histology_means$histology)))

p_bar_any <- ggplot(any_driver_histology, aes(prop, histology, label = paste0(round(prop, 1)))) +
  geom_bar(stat = "identity", fill = "black") +
  geom_text(hjust = 1.1, color = "white", size = 2) +
  xlab("% tumours with\nany cancer gene aberration") +
  theme_bw() +
  theme(text = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

theme_margin <- theme(legend.box.margin = margin(1, 2, 2, 1))
legend  <- cowplot::get_legend(p_bubble + theme_margin)
p_bubble <- p_bubble + theme(legend.position = "none")

p2 <- plot_grid(p_bubble, p_bar_any, legend, nrow = 1, rel_widths = c(2, 1, 1))

pdf(paste0(output_path, "histology_driver_count_bubble_bar.pdf"), width = 7, height = 4)
p2
dev.off()

