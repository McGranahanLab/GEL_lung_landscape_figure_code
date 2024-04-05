# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                    clinical characteristics plot                  # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
###################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(dplyr)
library(reshape2)
library(cowplot)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
######################################################################file paths

sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
clin_data_path        <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/sampleTables/v14_lung_patient_clinical_data_2022-07-05_useful.tsv"
CPI_response_path     <- "/re_gecip/cancer_lung/bsimpson/CPI_cohort/Scripts/Analysis_version_lock_v13/Mastersheet_V31.csv"
output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/GenomicOverview/"

################################################################################
######################################################################      MAI

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

# make met other into other
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

# make a plot showing histology, stage, sex, treatment info, age

data <- sample_table[, c("participant_id", "histology", "sex", "patient_age")]

# add stage
clinical <- read.table(clin_data_path, head = T, sep = "\t")
clinical$participant_id <- as.character(clinical$participant_id)
stage <- clinical[, c("participant_id", "stage_stage_best")]
stage$stage <- stage$stage_stage_best
stage[which(stage$stage == ""), "stage"] <- NA
stage[which(stage$stage == "?"), "stage"] <- NA
stage[which(stage$stage == "U"), "stage"] <- NA
stage[which(stage$stage %in% c("1A", "1A_1", "1A3", "1A_3B_", "1A1", "1A2", "1A2", "1B")), "stage"] <- 1
stage[which(stage$stage %in% c("2", "2A", "2B")), "stage"] <- 2
stage[which(stage$stage %in% c("3", "3A", "3B", "3C")), "stage"] <- 3
stage[which(stage$stage %in% c("4", "4A", "4B")), "stage"] <- 4
stage$stage <- as.character(stage$stage)
stage[stage$stage %in% c("1", "2", "3"), "stage_category"] <- "early"
stage[stage$stage %in% c("4"), "stage_category"] <- "late"

data <- left_join(data, stage[, c("participant_id", "stage")])

# which ones had surgery
surg <- clinical[, c("participant_id", "avtreatment_eventdesc")]
surg[grep("surgery", surg$avtreatment_eventdesc, ignore.case = T), "surgery"] <- TRUE
surg[is.na(surg$surgery), "surgery"] <- FALSE

data <- left_join(data, surg[, c("participant_id", "surgery")])

treatment <- clinical[, c("participant_id", "chemo_start_date_of_regimen", "RT_treatmentstartdate")]
treatment[is.na(treatment$chemo_start_date_of_regimen), "hadCHEMO"] <- FALSE
treatment[!is.na(treatment$chemo_start_date_of_regimen), "hadCHEMO"] <- TRUE
treatment[is.na(treatment$RT_treatmentstartdate), "hadRT"] <- FALSE
treatment[!is.na(treatment$RT_treatmentstartdate), "hadRT"] <- TRUE

data <- left_join(data, treatment[, c("participant_id", "hadCHEMO", "hadRT")])

# make into plots
hist <- data.frame(table(data$histology))
colnames(hist) <- c("histology", "count")

hist_order <- c("CARCINOID", "MESOTHELIOMA", "MET_ADENOCARCINOMA", "NEUROENDOCRINE_CARCINOMA",
                "OTHER", "ADENOCARCINOMA", "MET_SQUAMOUS_CELL", "MET_SMALL_CELL", 
                "ADENOSQUAMOUS", "SQUAMOUS_CELL", "SMALL_CELL", "LARGE_CELL")

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

names(histology_names) <- unique(hist_order)

hist$histology <- factor(hist$histology, hist_order)
plot_hist <- ggplot(data = hist, aes(x = count, y = histology, fill = histology)) +
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
                  geom_text(aes(label = count), hjust = 0, color = "black", size = 3.5)+
                  scale_x_continuous(expand = c(0,0), limits = c(0, 590)) +
                  facet_grid(histology ~ ., scales = "free", labeller = labeller(histology = histology_names), switch = "y") +
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
                        plot.margin = unit(c(0,0.1,0,0.1), "cm"),
                        panel.grid.major.x = element_line(colour = "grey88", size = 0.2),
                        strip.text.y.left = element_text(angle = 0))
# plot sex
sex <- data.frame(table(data[, c("histology", "sex")]))
sex <- left_join(sex, hist)
sex$prop <- sex$Freq/ sex$count
sex$prop <- sex$prop*100
sex$histology <- factor(sex$histology, levels = hist_order)

plot_sex <- ggplot(sex, aes(prop, histology, fill = sex)) +
                geom_bar(stat = "identity") +
                scale_fill_manual(values = c("lightblue", "lightpink"), breaks = c("MALE", "FEMALE"), label = c("male", "female"), name = "sex") +
                scale_x_continuous(expand = c(0,0), name = "% tumours") +
                # scale_y_continuous(expand = c(0,0)) +
                facet_grid(histology ~ ., scales = "free") +
                theme_bw() +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      plot.margin = unit(c(0,0.1,0,0.1), "cm"),
                      legend.justification = c(0,1),
                      panel.spacing = unit(0,'lines'),
                      legend.position = "bottom",
                      strip.background.y = element_blank(),
                      strip.text.y = element_blank()) +
                labs(x=NULL)

# plot stage
stage <- data.frame(table(data[, c("histology", "stage")]))
stage <- left_join(stage, hist)
stage$prop <- stage$Freq/ stage$count
stage$prop <- stage$prop*100
stage$histology <- factor(stage$histology, levels = hist_order)
stage$stage <- factor(stage$stage, levels = c("1", "2", "3", "4"))

plot_stage <- ggplot(stage, aes(prop, histology, fill = stage)) +
                    geom_bar(stat = "identity") +
                    scale_fill_manual(values = c("#ffffcc",
                                                  "#a1dab4",
                                                  "#41b6c4",
                                                  "#225ea8")) +
                    scale_x_continuous(expand = c(0,0), name = "% tumours") +
                    # scale_y_continuous(expand = c(0,0)) +
                    facet_grid(histology ~ ., scales = "free") +
                    theme_bw() +
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank(),
                          plot.margin = unit(c(0,0.1,0,0.1), "cm"),
                          legend.justification = c(0,1),
                          panel.spacing = unit(0,'lines'),
                          legend.position = "bottom",
                          strip.background.y = element_blank(),
                          strip.text.y = element_blank()) +
                    labs(x=NULL)

# plot surgery
surgery <- data.frame(table(data[, c("histology", "surgery")]))
surgery <- left_join(surgery, hist)
surgery$prop <- surgery$Freq/ surgery$count
surgery$prop <- surgery$prop*100
surgery$histology <- factor(surgery$histology, levels = hist_order)
surgery$surgery <- factor(surgery$surgery, levels = rev(c("TRUE", "FALSE")))

plot_surgery <- ggplot(surgery, aes(prop, histology, fill = surgery)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b5b5b5" , "#5ab4ac")) +
  scale_x_continuous(expand = c(0,0), name = "% tumours") +
  # scale_y_continuous(expand = c(0,0)) +
  facet_grid(histology ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "bottom",
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  labs(x=NULL)

# plot chemo
chemo <- data.frame(table(data[, c("histology", "hadCHEMO")]))
chemo <- left_join(chemo, hist)
chemo$prop <- chemo$Freq/ chemo$count
chemo$prop <- chemo$prop*100
chemo$histology <- factor(chemo$histology, levels = hist_order)
chemo$hadCHEMO <- factor(chemo$hadCHEMO, levels = rev(c("TRUE", "FALSE")))

plot_chemo <- ggplot(chemo, aes(prop, histology, fill = hadCHEMO)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b5b5b5" , "#ca0020")) +
  scale_x_continuous(expand = c(0,0), name = "% tumours") +
  # scale_y_continuous(expand = c(0,0)) +
  facet_grid(histology ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "bottom",
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  labs(x=NULL)

# plot RT
RT <- data.frame(table(data[, c("histology", "hadRT")]))
RT <- left_join(RT, hist)
RT$prop <- RT$Freq/ RT$count
RT$prop <- RT$prop*100
RT$histology <- factor(RT$histology, levels = hist_order)
RT$hadRT <- factor(RT$hadRT, levels = rev(c("TRUE", "FALSE")))

plot_RT <- ggplot(RT, aes(prop, histology, fill = hadRT)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b5b5b5" , "#fdae61")) +
  scale_x_continuous(expand = c(0,0), name = "% tumours") +
  # scale_y_continuous(expand = c(0,0)) +
  facet_grid(histology ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "bottom",
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  labs(x=NULL)

data$histology <- factor(data$histology, levels = hist_order)
plot_age <- ggplot(data = data, aes(x = patient_age, fill = histology)) +
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
            geom_density(data = data) +
            facet_grid(histology ~ ., scales = "free") +
            scale_x_continuous(expand = c(0,0), name = "age") +
            theme_bw() +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.margin = unit(c(0,0.1,0,0.1), "cm"),
                  legend.justification = c(0,1),
                  panel.spacing = unit(0,'lines'),
                  legend.position = "none",
                  strip.background.y = element_blank(),
                  strip.text.y = element_blank()) +
            labs(x=NULL)

p <- plot_grid(plot_hist, plot_age, plot_sex, plot_stage, plot_surgery, plot_chemo, plot_RT, align = "h", axis = "bt", nrow = 1, rel_widths = c(1.7, 1, 1, 1, 1, 1, 1))

pdf(paste0(output_path, "clincial_characteristics_overview.pdf"), width = 12, height = 6)
p
dev.off()
