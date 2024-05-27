################################################################################
###               correlate all signatures with each other                   ###
################################################################################

################################################################################
#################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))


library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(rtracklayer)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
####################################################################  file paths

signatures_path        <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
sample_table_path      <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"

output_path            <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"

################################################################################
####################################################################   functions

cor_plot_fun <- function(data = signatures_df,
                         variable_order = signature_order,
                         sig_lables = signature_labels,
                         outfile = "siganture_exposure_correlation"){
  
  cormat            <- data.frame(matrix(nrow = ncol(data), ncol = ncol(data)))
  rownames(cormat)  <- colnames(data)
  colnames(cormat)  <- colnames(data)
  
  for (i in 1:nrow(cormat)) {
    for (j in 1:ncol(cormat)) {
      tryCatch(
        cormat[i,j] <- cor.test(x = data[, rownames(cormat)[i]], y = data[, colnames(cormat)[j]])$estimate
        , error = function(e)
          print(paste0("error in ", i, " and ", j, " moving on")))
    }
  }
  
  cormat <- cormat[variable_order, variable_order]
  
  # compute p values
  pmat          <- data.frame(matrix(nrow = ncol(data), ncol = ncol(data)))
  rownames(pmat) <- colnames(data)
  colnames(pmat) <- colnames(data)
  
  for (i in 1:nrow(pmat)) {
    for (j in 1:ncol(pmat)) {
      tryCatch(
        pmat[i,j] <- cor.test(x = data[, rownames(pmat)[i]], y = data[, colnames(pmat)[j]])$p.value
        , error = function(e)
          print(paste0("error in ", i, " and ", j, " moving on")))
      
    }
  }
  
  pmat <- pmat[variable_order, variable_order]
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  upper_tri_p <- get_upper_tri(pmat)
  
  # adjust p values for multiple testing
  pmat_long <- upper_tri_p
  pmat_long$id <- rownames(pmat_long)
  pmat_long <- melt(pmat_long, na.rm = TRUE)
  pmat_long <- pmat_long[pmat_long$id != pmat_long$variable, ]
  
  pmat_long$variable <- as.character(pmat_long$variable)
  
  for(i in 1:nrow(pmat_long)){
    
    sig_type_one <- gsub("[0-9]", "", strsplit(pmat_long[i, "id"], ":")[[1]][1])
    sig_type_two <- gsub("[0-9]", "", strsplit(pmat_long[i, "variable"], ":")[[1]][1])
    
    sig_type_one <- gsub("[[:lower:]]", "", sig_type_one)
    sig_type_two <- gsub("[[:lower:]]", "", sig_type_two)
    
    sig_type_one <- gsub("_N", "", sig_type_one)
    sig_type_two <- gsub("_N", "", sig_type_two)
    
    both_types <- c(sig_type_one, sig_type_two)
    
    if(length(unique(both_types)) == 1){
      pmat_long[i, "same_sig_sig"] <- TRUE
    }
  }
  
  pmat_long <- pmat_long[is.na(pmat_long$same_sig_sig), ]
  pmat_long$same_sig_sig <- NULL
  
  
  pmat_long$value <- p.adjust(pmat_long$value, method = "fdr")
  pmat <- dcast(pmat_long, variable ~ id, value.var = "value")
  rownames(pmat) <- pmat$variable
  pmat$variable <- NULL
  pmat <- pmat[variable_order[variable_order %in% rownames(pmat)], variable_order[variable_order %in% colnames(pmat)]]
  
  # assemble information for plotting
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  upper_tri <- get_upper_tri(cormat)
  
  cormat_long <- upper_tri
  cormat_long$id <- rownames(cormat_long)
  cormat_long <- melt(cormat_long, na.rm = TRUE)
  cormat_long <- cormat_long[cormat_long$id != cormat_long$variable, ]
  
  cormat_long$variable <- as.character(cormat_long$variable)
  
  for(i in 1:nrow(cormat_long)){
    
    sig_type_one <- gsub("[0-9]", "", strsplit(cormat_long[i, "id"], ":")[[1]][1])
    sig_type_two <- gsub("[0-9]", "", strsplit(cormat_long[i, "variable"], ":")[[1]][1])
    
    sig_type_one <- gsub("[[:lower:]]", "", sig_type_one)
    sig_type_two <- gsub("[[:lower:]]", "", sig_type_two)
    
    sig_type_one <- gsub("_N", "", sig_type_one)
    sig_type_two <- gsub("_N", "", sig_type_two)
    
    both_types <- c(sig_type_one, sig_type_two)
    
    if(length(unique(both_types)) == 1){
      cormat_long[i, "same_sig_sig"] <- TRUE
    }
  }
  
  cormat_long <- cormat_long[is.na(cormat_long$same_sig_sig), ]
  cormat_long$same_sig_sig <- NULL
  
  cormat_long <- dcast(cormat_long, variable ~ id, value.var = "value")
  rownames(cormat_long) <- cormat_long$variable
  cormat_long <- cormat_long[, -1]
  cormat <- cormat_long[variable_order[variable_order %in% rownames(cormat_long)], variable_order[variable_order %in% colnames(cormat_long)]]
  
  # get p values in shape to plot on the cor plot
  pmat$variable <- rownames(pmat)
  
  # Melt the p matrix
  melted_pmat <- melt(pmat, na.rm = TRUE)
  colnames(melted_pmat) <- c("variable_a", "variable_b", "adjust_p_value")
  
  # put starts for signifcance levels
  melted_pmat$significance <- NULL
  melted_pmat[which(melted_pmat$adjust_p_value < 0.0001), "significance"] <- "****"
  melted_pmat[which(melted_pmat$adjust_p_value < 0.001 & melted_pmat$adjust_p_value >= 0.0001), "significance"] <- "***"
  melted_pmat[which(melted_pmat$adjust_p_value < 0.01 & melted_pmat$adjust_p_value >= 0.001), "significance"] <- "**"
  melted_pmat[which(melted_pmat$adjust_p_value < 0.05 & melted_pmat$adjust_p_value >= 0.01), "significance"] <- "*"
  melted_pmat[which(melted_pmat$adjust_p_value > 0.05), "significance"] <- ""
  
  # Melt the correlation matrix
  cormat$variable <- rownames(cormat)
  
  melted_cormat <- melt(cormat, na.rm = TRUE)
  colnames(melted_cormat) <- c("variable_a", "variable_b", "correlation")
  
  # put into one dataframe
  cor_df <- left_join(melted_cormat, melted_pmat)
  
  cor_df$variable_a <- factor(cor_df$variable_a, levels = variable_order[variable_order %in% cor_df$variable_a])
  cor_df$variable_b <- factor(cor_df$variable_b, levels = variable_order[variable_order %in% cor_df$variable_b])
  
  
  # Heatmap
  
  p <- ggplot(data = cor_df, aes(variable_a, variable_b, fill = correlation)) +
    geom_tile(color = "grey") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    scale_x_discrete(labels = sig_lables[grep(paste(c("DBS", "ID", "CN", "SV"), collapse = "|"), unique(cor_df$variable_a), value = T)]) +
    scale_y_discrete(labels = sig_lables[grep(paste(c("SBS", "DBS", "ID", "CN", "SV"), collapse = "|"), unique(cor_df$variable_b), value = T)]) +
    geom_text(aes(label = significance), color = "black", size = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text.y = element_text(size =12),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),) +
    xlab("") +
    ylab("")
  
  pdf(paste0(output_path, outfile, ".pdf"), width = 16, height = 10)
  print(p)
  dev.off()
  
}

################################################################################
####################################################################        MAIN
load(signatures_path)

# make into one data frame
signatures_df           <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure)

# make into matrix
signatures_df <- dcast(signatures_df, sample ~ label, value.var = "weight")
signatures_df[is.na(signatures_df)] <- 0

rownames(signatures_df)           <- signatures_df$sample
signatures_df$sample              <- NULL

###############################
####        MAIN           ####
###############################

# make an order of signatures
signature_order        <- c("SBS1:clock-like", "SBS2:APOBEC", "SBS3:HRd", "SBS4:smoking", "SBS5:clock-like", "SBS13:APOBEC", "SBS15:MMRd", "SBS17b:5FU-chemotherapy", "SBS18:ROS", "SBS21:MMRd", "SBS25:chemotherapy", "SBS44:MMRd", "SBS92:smoking",
                            "DBS2:smoking", "DBS4:unknown", "DBS11:unknown/APOBEC", "DBS_N1:unknown", "DBS_N3:unknown/platinum-chemo", "DBS_N4:unknown/MMRd",
                            "ID2:DNA-slippage", "ID3:smoking", "ID4:unknown", "ID6:HRd", "ID7:MMRd", "ID_N1:unknown", "ID_N2:unknown/NHEJ", "ID_N3:unknown",
                            "CN1:diploid", "CN2:tetraploid", "CN8:chromothripsis", "CN9:CIN", "CN17:HRd", "CN18:unknown/HRd", "CN20:unknown", "CN_N1:unknown", "CN_N2:unknown", "CN_N6:unknown", "CN_N8:unknown", "CN_N11:unknown", "CN_N12:unknown", "CN_N14:unknown", "CN_N15:unknown",
                            "SV2:unknown", "SV4:unknown", "SV6:unknown", "SV_N1:unknown", "SV_N2:unknown", "SV_N3:unknown", "SV_N4:unknown", "SV_N8:unknown", "SV_N9:unknown", "SV_N10:unknown", "SV_N11:unknown")

signature_labels        <- c("SBS1:clock-like", "SBS2:APOBEC", "SBS3:HRd", "SBS4:smoking", "SBS5:clock-like", "SBS13:APOBEC", "SBS15:MMRd", "SBS17b:5FU-chemotherapy", "SBS18:ROS", "SBS21:MMRd", "SBS25:chemotherapy", "SBS44:MMRd", "SBS92:smoking",
                            "DBS2:smoking", "DBS4:unknown", "DBS11:unknown/APOBEC", "DBS_N1:unknown", "DBS_N3:unknown/platinum-chemo", "DBS_N4:unknown/MMRd",
                            "ID2:DNA-slippage", "ID3:smoking", "ID4:unknown", "ID6:HRd", "ID7:MMRd", "ID_N1:unknown", "ID_N2:unknown/NHEJ", "ID_N3:unknown",
                            "CN1:diploid", "CN2:tetraploid", "CN8:chromothripsis", "CN9:CIN", "CN17:HRd", "CN18:unknown/HRd", "CN20:unknown", "CN_N1:unknown", "CN_N2:unknown", "CN_N6:unknown", "CN_N8:unknown", "CN_N11:unknown", "CN_N12:unknown", "CN_N14:unknown", "CN_N15:unknown",
                            "SV2:unknown", "SV4:unknown", "SV6:unknown", "SV_N1:unknown", "SV_N2:unknown", "SV_N3:unknown", "SV_N4:unknown", "SV_N8:unknown", "SV_N9:unknown", "SV_N10:unknown", "SV_N11:unknown")

names(signature_labels) <- signature_order

cor_plot_fun(data = signatures_df,
             variable_order = signature_order,
             sig_lables = signature_labels,
             outfile = "primary_met_timepoint_patient_signature_exposure_correlation_newCNsigs")
