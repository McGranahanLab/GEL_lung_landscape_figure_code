# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                   Mutational Signatures Overview figures                # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################
#### libraries & functions ####
###############################
library(ggplot2)
library(dplyr)
library(reshape2)
library(grid)
library(gridExtra)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

source('~/re_gecip/cancer_lung/kthol/IR_signatures/code/functions/fun_alignLegends_ggplot.R')

###############################
####      file paths       ####
###############################
signatures_path        <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"

sample_table_path      <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
output_path            <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"

SBS_matrix_path        <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/SBS/input_files/input_96matrix_patient.rds"
DBS_matrix_path        <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/DBS/input_files/input_78matrix_patient.rds"
ID_matrix_path         <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/ID/input_files/input_83matrix_patient.rds"  
SV_matrix_path         <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/SV/input_files/input_32matrix_patient.rds"

###############################
####      load data        ####
###############################
load(signatures_path)

# make into one data frame
sig_df           <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure)
# sig_df           <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure, MDS_exposure)

sig_df[grep("SBS", sig_df$label), "signature_type"] <- "SBS"
sig_df[grep("DBS", sig_df$label), "signature_type"] <- "DBS"
sig_df[grep("ID", sig_df$label), "signature_type"] <- "ID"
sig_df[grep("CN", sig_df$label), "signature_type"] <- "CN"
sig_df[grep("SV", sig_df$label), "signature_type"] <- "SV"
# sig_df[grep("MDS", sig_df$label), "signature_type"] <- "MDS"

sig_df$signature_label <- sig_df$label

# add histology
sample_table           <- read.table(sample_table_path, head = TRUE, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

sig_df                 <- left_join(sig_df, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))

# order
# order based on histology in same order as in overview plot
hist_order             <- c("ADENOCARCINOMA", "MET_ADENOCARCINOMA", "SQUAMOUS_CELL", "MET_SQUAMOUS_CELL", "MESOTHELIOMA",
                            "SMALL_CELL", "MET_SMALL_CELL", "LARGE_CELL", "OTHER", "CARCINOID", "ADENOSQUAMOUS",
                            "NEUROENDOCRINE_CARCINOMA")
sig_type_order         <- c("SBS", "DBS", "ID", "CN", "SV")
signature_order        <- c("SBS1:clock-like", "SBS2:APOBEC", "SBS3:HRd", "SBS4:smoking", "SBS5:clock-like", "SBS13:APOBEC", "SBS15:MMRd", "SBS17b:5FU-chemotherapy", "SBS18:ROS", "SBS21:MMRd", "SBS25:chemotherapy", "SBS44:MMRd", "SBS92:smoking",
                            "DBS2:smoking", "DBS4:unknown", "DBS11:unknown/APOBEC", "DBS_N1:unknown", "DBS_N3:unknown/platinum-chemo", "DBS_N4:unknown/MMRd",
                            "ID2:DNA-slippage", "ID3:smoking", "ID4:unknown", "ID6:HRd", "ID7:MMRd", "ID_N1:unknown", "ID_N2:unknown/NHEJ", "ID_N3:unknown",
                            "CN1:diploid", "CN2:tetraploid", "CN8:chromothripsis", "CN9:CIN", "CN17:HRd", "CN18:unknown/HRd", "CN20:unknown", "CN_N1:unknown", "CN_N2:unknown", "CN_N6:unknown", "CN_N8:unknown", "CN_N11:unknown", "CN_N12:unknown", "CN_N14:unknown", "CN_N15:unknown",
                            "SV2:unknown", "SV4:unknown", "SV6:unknown", "SV_N1:unknown", "SV_N2:unknown", "SV_N3:unknown", "SV_N4:unknown", "SV_N8:unknown", "SV_N9:unknown", "SV_N10:unknown", "SV_N11:unknown")
# order patients based on decreasing mutational burden
SBS_burden        <- readRDS(SBS_matrix_path)
SBS_burden        <- data.frame(rowSums(SBS_burden))
SBS_burden$sample <- rownames(SBS_burden)
colnames(SBS_burden)[1] <- "SBS_burden"

sig_df <- left_join(sig_df, SBS_burden)

sig_df            <- sig_df[order(sig_df$SBS_burden, decreasing = F), ]

sig_df$histology       <- factor(sig_df$histology, levels = hist_order)
sig_df$patient         <- factor(sig_df$sample, levels = unique(sig_df$sample))

sig_df$signature       <- sig_df$label
sig_df$signature_type  <- factor(sig_df$signature_type, levels = sig_type_order)
sig_df$signature       <- factor(sig_df$signature, levels = signature_order)


###############################
####      bar plot         ####
###############################

colour_pallette_SBS   <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#543005")
colour_pallette_DBS   <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")
colour_pallette_ID    <- c("#d53e4f", "#f46d43", "#fdae61", "#e6f598", "#66bd63", "#3288bd", "#053061", "#0c6624")
colour_pallette_CN    <- c("#67000d","#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffff99", "#80cdc1", "#a8ddb5", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695", "#08306b")
colour_pallette_SV    <- c("#bebada", "#fb8072", "#b3de69", "#8073ac", "#fccde5", "#97b1d1", "#bc80bd", "#ccebc5", "#fdc086", "#4eb3d3", "#238443")
colour_pallette_MDS   <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",  "#ffed6f")

plot_sigs_df  <- sig_df[, c("patient", "signature", "weight", "exposure", "signature_type", "signature_label", "histology")]
plot_sigs_df <- unique(plot_sigs_df)

histology_names <- c("ADENOCARCINOMA" = "Adenocarcinoma",
                     "MET_ADENOCARCINOMA" = "Adenocarcinoma\nmetastasis",
                     "SQUAMOUS_CELL" = "Squamous cell",
                     "MET_SQUAMOUS_CELL" = "Squamous cell\nmetastasis",
                     "MESOTHELIOMA" = "Mesothelioma",
                     "SMALL_CELL" = "Small cell",
                     "MET_SMALL_CELL" = "Small cell\nmetastasis",
                     "LARGE_CELL" = "Large cell",
                     "OTHER" = "Other",
                     "CARCINOID" = "Carcinoid",
                     "ADENOSQUAMOUS" = "Adenosquamous",
                     "NEUROENDOCRINE_CARCINOMA" = "Neuroendocrine")

names(histology_names) <- unique(hist_order)


# plot mutational burden burden
sample_table$participant_id <- as.character(sample_table$participant_id)
SBS_burden <- left_join(SBS_burden, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))
SBS_burden$sample <- factor(SBS_burden$sample, levels = unique(sig_df$sample))
SBS_burden$histology <- factor(SBS_burden$histology, hist_order)

# and add all the other burdens as well
DBS_burden <- readRDS(DBS_matrix_path)
DBS_burden <- data.frame(rowSums(DBS_burden))
DBS_burden$sample <- rownames(DBS_burden)
colnames(DBS_burden)[1] <- "DBS_burden"

ID_burden <- readRDS(ID_matrix_path)
ID_burden <- data.frame(rowSums(ID_burden))
ID_burden$sample <- rownames(ID_burden)
colnames(ID_burden)[1] <- "ID_burden"

SV_burden <- readRDS(SV_matrix_path)
SV_burden <- data.frame(SV_burden)
SV_burden <- data.frame(rowSums(SV_burden))
SV_burden$sample <- sub("X", "", rownames(SV_burden))
colnames(SV_burden)[1] <- "SV_burden"

SBS_burden <- left_join(SBS_burden, DBS_burden)
SBS_burden <- left_join(SBS_burden, ID_burden)
SBS_burden <- left_join(SBS_burden, SV_burden)

SBS_burden <- SBS_burden[order(SBS_burden$SBS_burden,
                               SBS_burden$DBS_burden,
                               SBS_burden$ID_burden,
                               SBS_burden$SV_burden, decreasing = F), ]
pat_order <- unique(SBS_burden$sample)

SBS_burden <- melt(SBS_burden, id.vars = c("sample", "histology"))
SBS_burden$sample <- factor(SBS_burden$sample, levels = pat_order)

p_SBS_burden <- ggplot(SBS_burden, aes(x = sample, y = value, color = variable)) + 
  geom_point(size = 0.3) + 
  theme_bw() + 
  xlab("") +
  ylab("mutation burden") +
  scale_y_log10() +
  scale_color_manual(values = c("SBS_burden" = "#d7191c",
                                "DBS_burden" = "#fdae61",
                                "ID_burden" = "#80cdc1",
                                "SV_burden" = "#5e3c99")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 17),
        strip.text.x.top = element_text(size = 12),
        strip.background = element_rect(fill = 'white'),
        panel.spacing = unit(0.2, 'lines')) +
  facet_grid(.~histology, scales = "free", labeller = labeller(histology = histology_names))

plot_sigs_df$patient <- factor(plot_sigs_df$patient, levels = pat_order)

p_sig_weight <- ggplot(plot_sigs_df, aes(x = patient, y = weight, fill = signature)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + 
  scale_fill_manual(values = c(colour_pallette_SBS, colour_pallette_DBS, colour_pallette_ID, colour_pallette_CN, colour_pallette_SV)) +
  ylab("% mutations") +
  xlab("tumour") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(size = 17),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 17),
        strip.text.y = element_text(size = 15),
        strip.background.y = element_rect(fill = 'white'),
        strip.background.x = element_blank(),
        panel.spacing = unit(0.2, 'lines')) +
  facet_grid(signature_type~histology, scales = "free", labeller = labeller(histology = histology_names)) +
  guides(fill=guide_legend(nrow=4))


#exclude legends from original plots
p_SBS_burden            <- p_SBS_burden + theme(legend.position = 'top')
p_sig_weight            <- p_sig_weight + theme(legend.position = 'bottom')

#combine both
p_SBS_burden          <- ggplotGrob(p_SBS_burden        + theme(plot.margin = unit(c(0.1, -10, 0.01, 0.1), "cm")))
p_sig_weight          <- ggplotGrob(p_sig_weight        + theme(plot.margin = unit(c(0, -10, 0.1, 0.1), "cm")))

p_SBS_burden$widths         <- grid::unit.pmax(p_SBS_burden$widths, p_sig_weight$widths)
p_sig_weight$widths         <- grid::unit.pmax(p_SBS_burden$widths, p_sig_weight$widths)


g <- gtable_rbind(p_SBS_burden, p_sig_weight)

#change colours and labels of strips
colours <- c("ADENOCARCINOMA" = "#67001f",
             "MET_ADENOCARCINOMA" = "#67001f80",
             "SQUAMOUS_CELL" = "#053061",
             "MET_SQUAMOUS_CELL" = "#05306180",
             "MESOTHELIOMA" = "#e3309e",
             "SMALL_CELL" = "#b89704",
             "MET_SMALL_CELL" = "#b8970480",
             "LARGE_CELL" = "#54278f",
             "OTHER" = "#7a7979",
             "CARCINOID" = "#fc8d62",
             "ADENOSQUAMOUS" = "#66c2a5",
             "NEUROENDOCRINE_CARCINOMA" = "#0b8ca3")

strip_name_colours <- c("white", "black", "white", "black", "black",  "black", "black", "white", "black", "black", "black", "black")

stript      <- which(grepl('strip-t', g$layout$name))[1:length(unique(SBS_burden$histology))]
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
  
  t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
  
  k <- k+1
}

id_panels_h <- unique(g$layout[grep("panel", g$layout$name), "t"])
g$heights[id_panels_h] <- grid::unit(c(2, 3, 3, 3, 3, 3), "null")

# also adjust widths
id_panels_w <- unique(g$layout[grep("panel", g$layout$name), "l"])
g$widths[id_panels_w] <- grid::unit(c(9, 2, 9, 2, 2, 2, 2, 2, 2, 2, 2, 2), "null")

# change the legend size
# is_legend <- which(g$layout$name == "guide-box")
# legend <- g$grobs[is_legend][[2]]
# legend <- legend$grobs[legend$layout$name == "guides"][[1]]
# 
# # Set widths in guide gtable
# width <- as.numeric(legend$widths[4]) # save bar width (assumes 'cm' unit) 
# legend$widths[4] <- unit(20, "null") # replace bar width
# 
# # Set width/x of bar/labels/ticks. Assumes everything is 'cm' unit.
# legend$grobs[[2]]$width <- unit(20, "npc")
# legend$grobs[[3]]$children[[1]]$x <- unit(
#   as.numeric(legend$grobs[[3]]$children[[1]]$x) / width, "npc"
# )
# legend$grobs[[5]]$x0 <- unit(as.numeric(legend$grobs[[5]]$x0) / width, "npc")
# legend$grobs[[5]]$x1 <- unit(as.numeric(legend$grobs[[5]]$x1) / width, "npc")

# Replace legend
# g$grobs[[is_legend]] <- legend



pdf(paste0(output_path, "hdp_primary_met_MutationalSignatures_consensus_overview_weights.pdf"), width = 26, height = 14)
grid::grid.draw(g)
dev.off()

png(paste0(output_path, "hdp_primary_met_MutationalSignatures_consensus_overview_weights.png"), width = 40, height = 26, unit = "cm", res = 900)
grid::grid.draw(g)
dev.off()

 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # weights small
# p_SBS_burden <- ggplot(SBS_burden, aes(x = sample, y = value, color = variable)) + 
#   geom_point(size = 0.3) + 
#   theme_bw() + 
#   xlab("") +
#   ylab("log10 mutation burden") +
#   scale_y_log10() +
#   scale_color_manual(values = c("SBS_burden" = "#d7191c",
#                                 "DBS_burden" = "#fdae61",
#                                 "ID_burden" = "#80cdc1",
#                                 "SV_burden" = "#5e3c99")) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size = 15),
#         axis.title.y = element_text(size = 17),
#         strip.text.x = element_text(size = 17),
#         strip.background = element_rect(fill = 'white'),
#         panel.spacing = unit(0.2, 'lines')) +
#   facet_grid(.~histology, scales = "free", labeller = labeller(histology = histology_names))
# 
# p_sig_weight <- ggplot(plot_sigs_df, aes(x = patient, y = weight, fill = signature)) + 
#   geom_bar(stat = 'identity', position = position_stack(reverse = T)) + 
#   scale_fill_manual(values = c(colour_pallette_SBS, colour_pallette_MDS, colour_pallette_DBS, colour_pallette_ID, colour_pallette_CN, colour_pallette_SV)) +
#   ylab("signature weight") +
#   xlab("tumour") +
#   theme_bw() + 
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.x = element_text(size = 17),
#         axis.text.y = element_text(size = 15),
#         axis.title.y = element_text(size = 17),
#         strip.text.y = element_text(size = 15),
#         strip.background.y = element_rect(fill = 'white'),
#         strip.background.x = element_blank(),
#         panel.spacing = unit(0.2, 'lines')) +
#   facet_grid(signature_type~histology, scales = "free", labeller = labeller(histology = histology_names))
# 
# #exclude legends from original plots
# p_SBS_burden            <- p_SBS_burden + theme(legend.position = 'none')
# p_sig_weight            <- p_sig_weight + theme(legend.position = 'bottom')
# 
# #combine both
# p_SBS_burden          <- ggplotGrob(p_SBS_burden        + theme(plot.margin = unit(c(0.1, -10, 0.1, 0.1), "cm")))
# p_sig_weight          <- ggplotGrob(p_sig_weight        + theme(plot.margin = unit(c(0.1, -10, 0.1, 0.1), "cm")))
# 
# p_SBS_burden$widths         <- grid::unit.pmax(p_SBS_burden$widths, p_sig_weight$widths)
# p_sig_weight$widths         <- grid::unit.pmax(p_SBS_burden$widths, p_sig_weight$widths)
# 
# 
# g <- gtable_rbind(p_SBS_burden, p_sig_weight)
# 
# #change colours and labels of strips
# colours <- c("ADENOCARCINOMA" = "#67001f",
#              "MET_ADENOCARCINOMA" = "#67001f80",
#              "SQUAMOUS_CELL" = "#053061",
#              "MET_SQUAMOUS_CELL" = "#05306180",
#              "MESOTHELIOMA" = "#e3309e",
#              "SMALL_CELL" = "#b89704",
#              "MET_SMALL_CELL" = "#b8970480",
#              "LARGE_CELL" = "#54278f",
#              "OTHER" = "#7a7979",
#              "CARCINOID" = "#fc8d62",
#              "ADENOSQUAMOUS" = "#66c2a5",
#              "NEUROENDOCRINE_CARCINOMA" = "#0b8ca3",
#              "MET_OTHER" = "#7a797980")
# 
# strip_name_colours <- c("white", "black", "white", "black", "black",  "black", "black", "white", "black", "black", "black", "black", "black")
# 
# stript      <- which(grepl('strip-t', g$layout$name))[1:length(unique(SBS_burden$histology))]
# k <- 1
# for (i in stript) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
#   
#   t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
#   
#   k <- k+1
# }
# 
# id_panels_h <- unique(g$layout[grep("panel", g$layout$name), "t"])
# g$heights[id_panels_h] <- grid::unit(c(1, 3, 3, 3, 3, 3, 3), "null")
# 
# pdf(paste0(output_path, "MutationalSignatures_consensus_overview_weights_small.pdf"), width = 36, height = 20)
# grid::grid.draw(g)
# dev.off()
# 
# png(paste0(output_path, "MutationalSignatures_consensus_overview_weights_small.png"), width = 50, height = 24, unit = "cm", res = 150)
# grid::grid.draw(g)
# dev.off()
# 
