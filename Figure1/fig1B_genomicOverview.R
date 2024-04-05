#####################################
###  mutational burden overview   ###
#####################################

#####################################
###      libraries & options      ###
#####################################
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(GenomicRanges)
library(cowplot)
library(gtable)
library(ggbeeswarm)
library(ggpubr)
library(scales)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

#####################################
###           file paths          ###
#####################################
sample_table_path  <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
muttable_path      <- "/re_gecip/cancer_lung/kthol/landscape/input/muttable.RData"
SV_table_path      <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/SV/input_files/input_32matrix_patient.rds"
mutburden_sample_count_table_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
wgd_path            <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/patient_WGD.txt"

output_path        <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/GenomicOverview/"

#####################################
###     load and format data      ###
#####################################
sample_table <- read.table(sample_table_path, head = TRUE, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

# for the purposes of this plot make the "other metastasis" into "other"
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

SV_table     <- readRDS(SV_table_path)

# load(muttable_path)
# 
# # get all the SBSs
# 
# SBS <- muttable[muttable$is_SNV == TRUE, ]
# SBS$histology <- NULL
# SBS <- left_join(SBS, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))
# SBS_patient_hist <- data.frame(table(SBS[, c("sample", "histology")]))
# SBS_patient_hist <- SBS_patient_hist[SBS_patient_hist$Freq != 0, ]
# SBS_patient_hist$histology <- as.character(SBS_patient_hist$histology)
# colnames(SBS_patient_hist)[3] <- "SBScount"
# 
# # get all the DBSs
# DBS <- muttable[muttable$is_Dinuc == TRUE, ]
# DBS$histology <- NULL
# DBS <- left_join(DBS, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))
# DBS_patient_hist <- data.frame(table(DBS[, c("sample", "histology")]))
# DBS_patient_hist <- DBS_patient_hist[DBS_patient_hist$Freq != 0, ]
# DBS_patient_hist$histology <- as.character(DBS_patient_hist$histology)
# colnames(DBS_patient_hist)[3] <- "DBScount"
# 
# # get all the IDs
# ID <- muttable[muttable$is_Indel == TRUE, ]
# ID$histology <- NULL
# ID <- left_join(ID, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))
# ID_patient_hist <- data.frame(table(ID[, c("sample", "histology")]))
# ID_patient_hist <- ID_patient_hist[ID_patient_hist$Freq != 0, ]
# ID_patient_hist$histology <- as.character(ID_patient_hist$histology)
# colnames(ID_patient_hist)[3] <- "IDcount"
# 
# # add SV
# SV_table <- data.frame(SV_table)
# SV_table <- data.frame(rowSums(SV_table))
# SV_table$sample <- rownames(SV_table)
# SV_table$sample <- gsub("X", "", SV_table$sample)
# colnames(SV_table) <- c("SVcount", "sample")
# 
# # make into one dataframe
# patient_hist_df <- left_join(SBS_patient_hist, ID_patient_hist)
# patient_hist_df <- left_join(patient_hist_df, DBS_patient_hist)
# patient_hist_df <- left_join(patient_hist_df, SV_table)
# 
# # delete muttable so R won't crash
# muttable <- NULL
# 
# saveRDS(patient_hist_df, file = mutburden_sample_count_table_path)
patient_hist_df <- readRDS(mutburden_sample_count_table_path)
patient_hist_df[patient_hist_df$histology == "MET_OTHER", "histology"] <- "OTHER"

# order
patient_hist_df <- patient_hist_df[order(patient_hist_df$SBScount,
                                         patient_hist_df$IDcount,
                                         patient_hist_df$DBScount,
                                         patient_hist_df$SVcount,
                                         decreasing = FALSE), ]

# make the mutation burdens into mutation per megabase, assuming that an average whole genome has 2800 Mb with sufficient coverage
patient_hist_df$SBScount  <- patient_hist_df$SBScount/2800
patient_hist_df$IDcount   <- patient_hist_df$IDcount/2800
patient_hist_df$DBScount  <- patient_hist_df$DBScount/2800
patient_hist_df$SVcount   <- patient_hist_df$SVcount/2800

patient_hist_df[is.na(patient_hist_df$SBScount), "SBScount"] <- 0
patient_hist_df[is.na(patient_hist_df$DBScount), "DBScount"] <- 0
patient_hist_df[is.na(patient_hist_df$IDcount), "IDcount"] <- 0
patient_hist_df[is.na(patient_hist_df$SVcount), "SVcount"] <- 0

# add mean
mean_SBS <- sapply(unique(patient_hist_df$histology), function(x){median(patient_hist_df[patient_hist_df$histology == x, "SBScount"], na.rm = TRUE)})
mean_SBS_hist <- data.frame(mean_SBS)
mean_SBS_hist$histology <- rownames(mean_SBS_hist)
patient_hist_df <- left_join(patient_hist_df, mean_SBS_hist)

mean_DBS <- sapply(unique(patient_hist_df$histology), function(x){median(patient_hist_df[patient_hist_df$histology == x, "DBScount"], na.rm = TRUE)})
mean_DBS_hist <- data.frame(mean_DBS)
mean_DBS_hist$histology <- rownames(mean_DBS_hist)
patient_hist_df <- left_join(patient_hist_df, mean_DBS_hist)

mean_ID <- sapply(unique(patient_hist_df$histology), function(x){median(patient_hist_df[patient_hist_df$histology == x, "IDcount"], na.rm = TRUE)})
mean_ID_hist <- data.frame(mean_ID)
mean_ID_hist$histology <- rownames(mean_ID_hist)
patient_hist_df <- left_join(patient_hist_df, mean_ID_hist)

mean_SV <- sapply(unique(patient_hist_df$histology), function(x){median(patient_hist_df[patient_hist_df$histology == x, "SVcount"], na.rm = TRUE)})
mean_SV_hist <- data.frame(mean_SV)
mean_SV_hist$histology <- rownames(mean_SV_hist)
patient_hist_df <- left_join(patient_hist_df, mean_SV_hist)

patient_hist_df$sample  <- factor(patient_hist_df$sample, levels = unique(patient_hist_df$sample))
patient_order           <- patient_hist_df$sample

patient_hist_df$sample  <- factor(patient_hist_df$sample, levels = patient_order)

# order histology by decreasing mean SBS burden
hist_order <- patient_hist_df[, c("histology", "mean_SBS")]
hist_order <- unique(hist_order)
hist_order <- hist_order[order(hist_order$mean_SBS, decreasing = FALSE), ]

hist_order$histology <- factor(hist_order$histology, levels = hist_order$histology)

patient_hist_df$histology <- factor(patient_hist_df$histology, levels = hist_order$histology)

patient_hist_df_long <- melt(patient_hist_df[, c("sample", "histology", "SBScount",
                                                "DBScount", "IDcount", "SVcount")],
                            id.vars = c("sample", "histology"))

means <- patient_hist_df[, c("histology", "mean_SBS", "mean_DBS", "mean_ID", "mean_SV")]
means <- unique(means)
means <- melt(means)
means$variable <- as.character(means$variable)

means[means$variable == "mean_SBS", "variable"] <- "SBScount"
means[means$variable == "mean_DBS", "variable"] <- "DBScount"
means[means$variable == "mean_ID", "variable"] <- "IDcount"
means[means$variable == "mean_SV", "variable"] <- "SVcount"

# also calculate the median across the whole cohort
SBS_cohort_median <- median(patient_hist_df$SBScount)
DBS_cohort_median <- median(patient_hist_df$DBScount)
ID_cohort_median <- median(patient_hist_df$IDcount)
SV_cohort_median <- median(patient_hist_df$SVcount)

cohort_median_df <- data.frame(variable = c("SBScount", "DBScount", "IDcount", "SVcount"),
                                cohort_mean = c(SBS_cohort_median, DBS_cohort_median, ID_cohort_median, SV_cohort_median))

colnames(means) <- c("histology", "variable", "mean")

means$variable <- factor(means$variable, levels = c("SBScount", "DBScount", "IDcount", "SVcount", "Tcell_fraction"))


#####################################
###       plotting function       ###
#####################################
histology_names <- c("Carcinoid",
                    "Mesothelioma",
                    "Adenocarcinoma\nmetastasis",
                    "Neuroendocrine",
                    "Other",
                    "Adenocarcinoma",
                    "Squamous cell\nmetastasis",
                    "Small cell\nmetastasis",
                    "Adenosquamous",
                    "Squamous cell",
                    "Small cell",
                    "Large cell")

names(histology_names) <- unique(hist_order$histology)

mut_names <- c("SBS", "DBS", "ID", "SV")
names(mut_names) <- c("SBScount", "DBScount", "IDcount", "SVcount")

patient_hist_df_long$variable <- factor(patient_hist_df_long$variable, levels = c("SBScount", "DBScount", "IDcount", "SVcount"))
means$variales <- factor(means$variable, levels = c("SBScount", "DBScount", "IDcount", "SVcount"))
cohort_median_df$variable <- factor(cohort_median_df$variable, levels = c("SBScount", "DBScount", "IDcount", "SVcount"))

# change order

patient_hist_df_long$histology <- factor(patient_hist_df_long$histology, levels = hist_order$histology)

plot_burden <- ggplot(data = patient_hist_df_long[patient_hist_df_long$variable %in% c("SBScount", "DBScount", "IDcount", "SVcount"),],
                aes(x = value, fill = histology)) +
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
                geom_density(data = patient_hist_df_long[patient_hist_df_long$variable %in% c("SBScount", "DBScount", "IDcount", "SVcount"),]) +
                geom_vline(data = means[means$variable %in% c("SBScount", "DBScount", "IDcount", "SVcount"),], aes(xintercept = mean), color = "black", linetype = "dashed", size = 1.2) +
                geom_vline(data = cohort_median_df, aes(xintercept = cohort_mean), color = "grey43", size = 1, linetype = "dotted") +
                facet_grid(histology ~ variable, scales = "free_y", labeller = labeller(histology = histology_names,
                                                                                      variable = mut_names)) +
                scale_x_log10(labels = label_comma(), breaks = c(0.01, 1, 5, 10, 100), expand = c(0,0)) +
                scale_y_continuous(position = "right", ) +
                ylab("density") +
                xlab("mutation burden / megabase") +
                theme_classic() +
                theme(axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14),
                      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "none",
                      panel.spacing = unit(0.2, 'lines'),
                      plot.margin = unit(c(0.5,0,0.5,0.5), "cm"),
                      strip.background.y = element_blank(),
                      strip.text.y = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.2))

hist_count <- data.frame(table(sample_table$histology))
hist_count$Var1 <- factor(hist_count$Var1, levels = hist_order$histology)

plot_hist <- ggplot(data = hist_count, aes(x = Freq, y = Var1, fill = Var1)) +
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
                geom_bar(stat = "identity") +
                geom_text(aes(label = Freq), hjust = 0, color = "black", size = 3.5)+
                scale_x_continuous(expand = c(0,0), limits = c(0, 590)) +
                facet_grid(Var1 ~ ., scales = "free", labeller = labeller(Var1 = histology_names), switch = "y") +
                xlab("# samples") +
                ylab("") +
                # xlim(0, 530) +
                theme_classic() +
                theme(axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14),
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      legend.position = "none",
                      panel.spacing = unit(0.2, 'lines'),
                      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                      panel.grid.major.x = element_line(colour = "grey88", size = 0.2),
                      strip.text.y.left = element_text(angle = 0))

# add whole genome doubling
wgd     <- read.table(wgd_path, head = T, sep = "\t")
wgd$patient <- as.character(wgd$patient)
wgd <- left_join(sample_table[, c("participant_id", "histology")], wgd[, c("patient", "first_gd")], by = c("participant_id" = "patient"))
wgd_count <- data.frame(table(wgd[, c("histology", "first_gd")]))

# calculate proportion to all samples in that histology
hist_count <- data.frame(table(sample_table$histology))
colnames(hist_count) <- c("histology", "hist_count")
wgd_count <- left_join(wgd_count, hist_count)
wgd_count$prop <- (wgd_count$Freq/wgd_count$hist_count)*100

wgd_count$histology <- factor(wgd_count$histology, levels = rev(hist_order$histology))
# plot
p_wgd <- ggplot(wgd_count, aes(histology, prop, fill = histology, alpha = first_gd)) +
              geom_bar(stat = "identity") +
              ylab("% samples") +
              xlab(" ") +
              scale_alpha_manual(name = "WGD", 
                                 values = c("TRUE" = 1,
                                            "FALSE" = 0.5),
                                 labels = c("TRUE" = "WGD",
                                            "FALSE" = "no WGD")) +
              scale_fill_manual(name = "histology",
                                values = c("ADENOCARCINOMA" = "#67001f",
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
                                           "OTHER" = "#7a7979",
                                           "MET_OTHER" = "#7a797980"),
                                guide = "none") +
              geom_text(data = wgd_count[wgd_count$first_gd == TRUE, ], aes(label = round(prop)), hjust = 1.1) +
              scale_y_continuous(expand = c(0,0)) +
              theme_bw() +
              theme(axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.position = "right",
                    panel.grid = element_blank()) +
              coord_flip()

p <- plot_grid(plot_hist, plot_burden, p_wgd, align = "h", axis = "bt", rel_widths = c(3, 8, 2.5), nrow = 1)

pdf(paste0(output_path, '/main_fig1A_genomicOverview.pdf'), width = 16, height = 6)
p
dev.off()


means[means$histology == "LARGE_CELL" & means$variable == "SBScount", "mean"]
min(patient_hist_df[patient_hist_df$histology == "LARGE_CELL", "SBScount"])
max(patient_hist_df[patient_hist_df$histology == "LARGE_CELL", "SBScount"])

means[means$histology == "CARCINOID" & means$variable == "SBScount", "mean"]
min(patient_hist_df[patient_hist_df$histology == "CARCINOID", "SBScount"])
max(patient_hist_df[patient_hist_df$histology == "CARCINOID", "SBScount"])

means[means$histology == "MESOTHELIOMA" & means$variable == "SBScount", "mean"]
min(patient_hist_df[patient_hist_df$histology == "MESOTHELIOMA", "SBScount"])
max(patient_hist_df[patient_hist_df$histology == "MESOTHELIOMA", "SBScount"])

SV_count <- data.frame(rowSums(data.frame(SV_table)))
SV_count$sample <- rownames(SV_count)
colnames(SV_count) <- c("SVcount", "sample")
SV_count$sample <- sub("X", "", SV_count$sample)
SV_count <- SV_count[SV_count$sample %in% sample_table$participant_id, ]
length(unique(SV_count$sample))

means[means$histology == "LARGE_CELL" & means$variable == "SVcount", "mean"]
min(patient_hist_df[patient_hist_df$histology == "LARGE_CELL", "SVcount"])
max(patient_hist_df[patient_hist_df$histology == "LARGE_CELL", "SVcount"])

means[means$histology == "MESOTHELIOMA" & means$variable == "SVcount", "mean"]
min(patient_hist_df[patient_hist_df$histology == "MESOTHELIOMA", "SVcount"])
max(patient_hist_df[patient_hist_df$histology == "MESOTHELIOMA", "SVcount"])

means[means$variable == "SVcount", ]
means[means$variable == "SBScount", ]

sum(wgd_count[wgd_count$first_gd == TRUE, "Freq"])
wgd_count[wgd_count$first_gd == TRUE & wgd_count$histology == "LARGE_CELL", "prop"]
wgd_count[wgd_count$first_gd == TRUE & wgd_count$histology == "MESOTHELIOMA", "prop"]
