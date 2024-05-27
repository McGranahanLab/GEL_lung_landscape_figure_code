# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                        signature cluster plot                         # # #    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
###################################################################   file paths

################################################################################
####################################################################   libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggpubr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(survival)
library(survminer)
library(ggbeeswarm)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(plotrix)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
####################################################################  file paths

signatures_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"

output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/clustering/ALL_new/"
cluster_output_path   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/"

num_clusters          <- 12

################################################################################
####################################################################        MAIN
signature_order <- c("SBS1:clock-like", "SBS2:APOBEC", "SBS3:HRd", "SBS4:smoking", "SBS5:clock-like", "SBS13:APOBEC", "SBS15:MMRd", "SBS17b:5FU-chemotherapy", "SBS18:ROS", "SBS21:MMRd", "SBS25:chemotherapy", "SBS44:MMRd", "SBS92:smoking",
                     "DBS2:smoking", "DBS4:unknown", "DBS11:unknown/APOBEC", "DBS_N1:unknown", "DBS_N3:unknown/platinum-chemo", "DBS_N4:unknown/MMRd",
                     "ID2:DNA-slippage", "ID3:smoking", "ID4:unknown", "ID6:HRd", "ID7:MMRd", "ID_N1:unknown", "ID_N2:unknown/NHEJ", "ID_N3:unknown",
                     "CN1:diploid", "CN2:tetraploid", "CN8:chromothripsis", "CN9:CIN", "CN17:HRd", "CN18:unknown/HRd", "CN20:unknown", "CN_N1:unknown", "CN_N2:unknown", "CN_N6:unknown", "CN_N8:unknown", "CN_N11:unknown", "CN_N12:unknown", "CN_N14:unknown", "CN_N15:unknown",
                     "SV2:unknown", "SV4:unknown", "SV6:unknown", "SV_N1:unknown", "SV_N2:unknown", "SV_N3:unknown", "SV_N4:unknown", "SV_N8:unknown", "SV_N9:unknown", "SV_N10:unknown", "SV_N11:unknown", "other_cluster")

colour_pallette_SBS   <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#543005")
colour_pallette_DBS   <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")
colour_pallette_ID    <- c("#d53e4f", "#f46d43", "#fdae61", "#e6f598", "#66bd63", "#3288bd",  "#08519c", "#053061")
colour_pallette_CN    <- c("#67000d","#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffff99", "#80cdc1", "#a8ddb5", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695", "#08306b")
colour_pallette_SV    <- c("#bebada", "#fb8072", "#b3de69", "#8073ac", "#fccde5", "#97b1d1", "#bc80bd", "#ccebc5", "#fdc086", "#4eb3d3", "#238443", "#bababa")
colour_pallette_MDS   <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",  "#ffed6f")

# make a colour palette
names(colour_pallette_SBS) <- c("SBS1:clock-like", "SBS2:APOBEC", "SBS3:HRd", "SBS4:smoking", "SBS5:clock-like", "SBS13:APOBEC", "SBS15:MMRd", "SBS17b:5FU-chemotherapy", "SBS18:ROS", "SBS21:MMRd", "SBS25:chemotherapy", "SBS44:MMRd", "SBS92:smoking")
names(colour_pallette_DBS) <- c("DBS2:smoking", "DBS4:unknown", "DBS11:unknown/APOBEC", "DBS_N1:unknown", "DBS_N3:unknown/platinum-chemo", "DBS_N4:unknown/MMRd")
names(colour_pallette_ID)  <- c("ID2:DNA-slippage", "ID3:smoking", "ID4:unknown", "ID6:HRd", "ID7:MMRd", "ID_N1:unknown", "ID_N2:unknown/NHEJ", "ID_N3:unknown")
names(colour_pallette_CN)  <- c("CN1:diploid", "CN2:tetraploid", "CN8:chromothripsis", "CN9:CIN", "CN17:HRd", "CN18:unknown/HRd", "CN20:unknown", "CN_N1:unknown", "CN_N2:unknown", "CN_N6:unknown", "CN_N8:unknown", "CN_N11:unknown", "CN_N12:unknown", "CN_N14:unknown", "CN_N15:unknown")
names(colour_pallette_SV)  <- c("SV2:unknown", "SV4:unknown", "SV6:unknown", "SV_N1:unknown", "SV_N2:unknown", "SV_N3:unknown", "SV_N4:unknown", "SV_N8:unknown", "SV_N9:unknown", "SV_N10:unknown", "SV_N11:unknown", "other_cluster")

colors <- c(colour_pallette_SBS, colour_pallette_DBS, colour_pallette_ID, colour_pallette_CN, colour_pallette_SV)

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

cluster_assignment <- read.table(paste0(cluster_output_path, "new_signature_cluster_patient_assignment.txt"), sep = "\t", head = T)
cluster_assignment$variable <- as.character(cluster_assignment$variable)

# load signatures
load(signatures_path)

all_sigs <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure)
all_sigs$weight <- as.numeric(as.character(all_sigs$weight))
all_sigs$exposure <- as.numeric(as.character(all_sigs$exposure))

mat_exposure <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure)
mat_exposure$weight <- as.numeric(as.character(mat_exposure$weight))
mat_exposure$exposure <- as.numeric(as.character(mat_exposure$exposure))
mat_exposure <- dcast(mat_exposure, sample ~ label, value.var = "exposure")
mat_exposure <- melt(mat_exposure)
colnames(mat_exposure) <- c("sample", "label", "exposure")
mat_exposure[is.na(mat_exposure$weight), "exposure"] <- 0
mat_exposure$label <- as.character(mat_exposure$label)

# which cluster are these tumours in
mat_exposure <- left_join(mat_exposure, cluster_assignment, by = c("sample" = "variable"))

# for each cluster which are the dominant signatures
sig_list <- list("1" = c("SBS2:APOBEC", "SBS13:APOBEC"),
                 "2" = c("SBS5:clock-like"),
                 "3" = c("SBS18:ROS", "ID_N3:unknown"),
                 "4" = c("SBS3:HRd", "ID6:HRd"),
                 "5" = c("CN_N8:unknown"),
                 "6" = c("CN18:unknown/HRd"),
                 "7" = c("CN8:chromothripsis"),
                 "8" = c("CN20:unknown"),
                 "9" = c("SV_N3:unknown", "CN_N2:unknown"),
                 "10" = c("SBS4:smoking", "ID3:smoking", "DBS2:smoking"),
                 "11" = c("CN_N6:unknown"),
                 "12" = c("SBS4:smoking", "ID3:smoking", "DBS2:smoking"))

# make bar plots with the number of mutations of these signatures in each cluster


for(i in 1:num_clusters){
  print(i)
  
  sigs <- sig_list[[i]]
  
  df <- mat_exposure[mat_exposure$label %in% sigs, ]
  df[df$cluster_assignment == i, "in_cluster"] <- TRUE
  df[df$cluster_assignment != i, "in_cluster"] <- FALSE
  df$in_cluster <- factor(df$in_cluster, levels = c(TRUE, FALSE))
  
  # calculate mean and sd
  df_stat <- df %>% group_by(label, in_cluster) %>%
    mutate(mean = mean(exposure, na.rm = T)) %>%
    mutate(sample_size = n()) %>%
    mutate(sd = sd(exposure, na.rm = T)) %>%
    mutate(ci = qnorm(0.975)*sd/sqrt(sample_size)) %>%
    mutate(ymin = mean - ci) %>%
    mutate(ymax = mean + ci)
  
  df_stat <- unique(df_stat[, c("label", "in_cluster", "mean", "ymin", "ymax")])
  
  df_stat$in_cluster <- as.character(df_stat$in_cluster)
  df_stat[df_stat$in_cluster == TRUE, "colour_label"] <- df_stat[df_stat$in_cluster == TRUE, "label"]
  df_stat[df_stat$in_cluster == FALSE, "colour_label"] <- "other_cluster"
  
  p <- ggplot(df_stat, aes(x = label, y = mean, ymin = ymin, ymax = ymax, color = colour_label)) +
    geom_pointrange(fatten = 2) +
    scale_y_log10() + 
    scale_color_manual(values = colors, name = "") +
    ylab("# mutations") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "top",
          plot.margin = unit(c(2, 3, 2, 3), "cm"))
  
  
  
  n_sigs <- length(unique(df$label))
  
  h <- 4
  if(n_sigs == 3){w <- 4}
  if(n_sigs == 2){w <- 3.5}
  if(n_sigs == 1){w <- 3}
  
  pdf(paste0(output_path, "signature_attributions_in_cluster_", i, ".pdf"), width = w, height = h)
  print(p)
  dev.off()
  
}


