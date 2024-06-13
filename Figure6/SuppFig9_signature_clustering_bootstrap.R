# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                    signature clustering clean                       # # #   
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

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
library(stringr)
library(readr)
library(parallel)
library(dendextend)
library(fpc)
library(data.table)
library(circlize)
library(ComplexHeatmap)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

set.seed(123)

################################################################################
####################################################################  file paths

signatures_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"

output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"

################################################################################
####################################################################        MAIN
signature_order <- c("SBS1:clock-like", "SBS2:APOBEC", "SBS3:HRd", "SBS4:smoking", "SBS5:clock-like", "SBS13:APOBEC", "SBS15:MMRd", "SBS17b:5FU-chemotherapy", "SBS18:ROS", "SBS21:MMRd", "SBS25:chemotherapy", "SBS44:MMRd", "SBS92:smoking",
                     "DBS2:smoking", "DBS4:unknown", "DBS11:unknown/APOBEC", "DBS_N1:unknown", "DBS_N3:unknown/platinum-chemo", "DBS_N4:unknown/MMRd",
                     "ID2:DNA-slippage", "ID3:smoking", "ID4:unknown", "ID6:HRd", "ID7:MMRd", "ID_N1:unknown", "ID_N2:unknown/NHEJ", "ID_N3:unknown",
                     "CN1:diploid", "CN2:tetraploid", "CN8:chromothripsis", "CN9:CIN", "CN17:HRd", "CN18:unknown/HRd", "CN20:unknown", "CN_N1:unknown", "CN_N2:unknown", "CN_N6:unknown", "CN_N8:unknown", "CN_N11:unknown", "CN_N12:unknown", "CN_N14:unknown", "CN_N15:unknown",
                     "SV2:unknown", "SV4:unknown", "SV6:unknown", "SV_N1:unknown", "SV_N2:unknown", "SV_N3:unknown", "SV_N4:unknown", "SV_N8:unknown", "SV_N9:unknown", "SV_N10:unknown", "SV_N11:unknown")

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

# load signatures
load(signatures_path)

all_sigs <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure)
all_sigs$weight <- as.numeric(as.character(all_sigs$weight))
all_sigs$exposure <- as.numeric(as.character(all_sigs$exposure))

# make NAs into 0s
mat0 <- dcast(all_sigs, sample ~ label, value.var = "weight")
mat0 <- melt(mat0)
colnames(mat0) <- c("sample", "label", "weight")

mat0[is.na(mat0$weight), "weight"] <- 0
all_sigs <- mat0
all_sigs$label <- as.character(all_sigs$label)

# for each signature calculate the z scores

all_sigs$sig_type <- sapply(all_sigs$label, function(x){strsplit(x, ":")[[1]][1]})
all_sigs$sig_type <- sub("_N.", "", all_sigs$sig_type )
all_sigs$sig_type <- gsub("[[:digit:]]", "", all_sigs$sig_type)
all_sigs$sig_type <- gsub("b", "", all_sigs$sig_type)

all_sigs_z <- list()
for(sig in unique(all_sigs$label)){
  
  df <- all_sigs[all_sigs$label == sig, ]
  df$weight_z <- scale(df$weight)
  # df$exposure_z <- scale(log10(df$exposure +1))
  
  all_sigs_z[[sig]] <- df
  
}
all_sigs_z <- do.call(rbind, all_sigs_z)

# remove signatures present in less than 3% of tumours
signatures_to_keep <- c()
for(sig in unique(all_sigs_z$label)){

  df <- all_sigs_z[all_sigs_z$label == sig, ]
  num_sample_with_sig <- length(unique(df[df$weight >= 0.05, "sample"]))

  # in how many samples were we able to get signatures for this sig_type
  sig_type <- unique(df$sig_type)

  if(sig_type == "SBS"){
    num_sample_sig_type <- length(unique(SBS_exposure$sample))
  }

  if(sig_type == "DBS"){
    num_sample_sig_type <- length(unique(DBS_exposure$sample))
  }

  if(sig_type == "ID"){
    num_sample_sig_type <- length(unique(ID_exposure$sample))
  }

  if(sig_type == "CN"){
    num_sample_sig_type <- length(unique(CN_exposure$sample))
  }

  if(sig_type == "SV"){
    num_sample_sig_type <- length(unique(SV_exposure$sample))
  }

  # what is the porportion of samples that have this signature out of all of them
  percentage <- (num_sample_with_sig/1011)*100
  if(percentage >= 3){
    signatures_to_keep[[sig]] <- sig
  }

}
signatures_to_keep <- as.character(signatures_to_keep)

all_sigs_a <- all_sigs_z[all_sigs_z$label %in% signatures_to_keep, ]
all_sigs_a <- all_sigs_z

##### 
##### 
##### now we want to check if the clusters stay the same when we remove random samples
# original_assignment <- unique(all_sigs_inf_clust_order[, c("sample", "cluster_assignment")])

cluster.assignment.bootstrap <- function(signatures_df = all_sigs_a,
                                         iterations = 2,
                                         percentage = 5,
                                         value = "weight_z",
                                         num_cluster = 17){
  
  # perform orginial clustering
  og_mat <- dcast(signatures_df, sample ~ label, value.var = "weight_z")
  rownames(og_mat) <- og_mat$sample
  og_mat$sample <- NULL
  
  og_distance <- dist(og_mat, method = "euclidean")
  og_dendro <- hclust(og_distance, method = "ward.D2")
  
  original_assignment <- data.frame(cutree(og_dendro, k = num_cluster))
  original_assignment$sample <- rownames(original_assignment)
  colnames(original_assignment)[1] <- "cluster_a"
  
  dendrogram_comparison <- list()
  sample_cluster_iterations <- list()
  
  out <- lapply(1:iterations, function(i){
    
    print(i)
    samples_to_remove <- sample(original_assignment$sample, round(1011*(percentage/100)), replace = FALSE)
    
    # now cluster samples again but remove these ones
    mat_b <- dcast(signatures_df[-which(signatures_df$sample %in% samples_to_remove), ], sample ~ label, value.var = value)
    rownames(mat_b) <- mat_b$sample
    mat_b$sample <- NULL
    distance <- dist(mat_b, method = "euclidean")
    cluster <- hclust(distance, method = "ward.D2")
    
    new_cluster_assignment <- data.frame(cutree(cluster, k = num_cluster))
    new_cluster_assignment$sample <- rownames(new_cluster_assignment)
    colnames(new_cluster_assignment)[1] <- "cluster_b"
    
    new_dendro <- cluster
    
    # prune leaves off of og dendrgram
    pruned_og_dendro <- prune(og_dendro, samples_to_remove)
    
    # now compare
    mean_relative_difference_merge <- readr::parse_number(all.equal(new_dendro, pruned_og_dendro)[2])
    mean_relative_difference_order <- readr::parse_number(all.equal(new_dendro, pruned_og_dendro)[4])
    cophenetic_distance            <- cor_cophenetic(new_dendro, pruned_og_dendro)
    FM_index                       <- FM_index(cutree(new_dendro, k = num_cluster), cutree(pruned_og_dendro, k = num_cluster))[1]
    
    dendro_compare_df              <- data.frame(number_cluster = num_cluster,
                                                 iteration = i,
                                                 mean_relative_difference_merge = mean_relative_difference_merge,
                                                 mean_relative_difference_order = mean_relative_difference_order,
                                                 cophenetic_distance = cophenetic_distance,
                                                 FM_index = FM_index)
    
    # now for each sample, with which other samples does this sample cluster together and is this the same as in the original clustering
    sample_df <- mclapply(unique(new_cluster_assignment$sample), function(s){
      s_cluster <- new_cluster_assignment[new_cluster_assignment$sample == s, "cluster_b"]
      
      # which are the other samples that are in this cluster
      s_cluster_friends <- new_cluster_assignment[new_cluster_assignment$cluster_b == s_cluster, "sample"]
      
      # did s have the same friends in the original cluster assignment? let's find out
      s_old_cluster <- original_assignment[original_assignment$sample == s, "cluster_a"]
      s_old_cluster_friends <- original_assignment[original_assignment$cluster_a == s_old_cluster, "sample"]
      
      # how high is the overlap
      friend_numbers <- data.frame(table(s_cluster_friends %in% s_old_cluster_friends))
      
      if(nrow(friend_numbers) == 1){
        df <- data.frame(matching = friend_numbers[friend_numbers$Var1 == TRUE, "Freq"],
                         not_matching = 0)
      } else {df <- data.frame(matching = friend_numbers[friend_numbers$Var1 == TRUE, "Freq"],
                               not_matching = friend_numbers[friend_numbers$Var1 == FALSE, "Freq"])}
      
      
      friend_numbers <- df
      colnames(friend_numbers) <- c("num_samples_in_new_cluster_and_old_cluster", "num_samples_NOT_in_new_cluster_and_old_cluster")
      
      # add the total number of sampls in new cluster and in old cluster
      friend_numbers$num_samples_new_cluster <- length(s_cluster_friends)
      friend_numbers$num_samples_old_cluster <- length(s_old_cluster_friends)
      friend_numbers$sample <- s
      friend_numbers$percent_samples_in_new_cluster_and_old_cluster <- (as.numeric(friend_numbers$num_samples_in_new_cluster_and_old_cluster)/as.numeric(friend_numbers$num_samples_new_cluster))*100
      friend_numbers$old_cluster_id <- s_old_cluster
      
      return(friend_numbers)
      
    }, mc.cores = 3)
    
    sample_df <- do.call(rbind, sample_df)
    sample_df$iteration <- i
    sample_df$number_of_clusters <- num_cluster
    
    return(list(sample_df, dendro_compare_df))
    
  })
  
  
  # sort out the output
  
  sample_cluster_iterations_all <- lapply(1:iterations, function(x){return(out[[x]][[1]])})
  sample_cluster_iterations_all <- do.call(rbind, sample_cluster_iterations_all)
  
  dendrogram_comparison_all <- lapply(1:iterations, function(x){return(out[[x]][[2]])})
  dendrogram_comparison_all <- do.call(rbind, dendrogram_comparison_all)
  
  return(list(sample_cluster_iterations_all, dendrogram_comparison_all))
  
}


# iterate this function over a range of clusters
boostrap_dendro_compare <- list()
boostrap_sample_compare <- list()

for(c in c(1:20)){
  print(paste0("k= ", c))
  
  bootstrap_test <- cluster.assignment.bootstrap(signatures_df = all_sigs_a,
                                                 iterations = 100,
                                                 percentage = 5,
                                                 value = "weight_z",
                                                 num_cluster = c)
  
  sample_compare <- bootstrap_test[[1]]
  dendro_compare <- bootstrap_test[[2]]
  
  boostrap_dendro_compare[[c]] <- dendro_compare
  boostrap_sample_compare[[c]] <- sample_compare
  
  
}

boostrap_dendro_compare_all <- do.call(rbind, boostrap_dendro_compare)
boostrap_sample_compare_all <- do.call(rbind, boostrap_sample_compare)

# save(boostrap_sample_compare_all, boostrap_dendro_compare_all , file = paste0(output_path, "signature_cluster_bootstrap_with_zeros_all_signatures.RData"))
# save(boostrap_sample_compare_all, boostrap_dendro_compare_all , file = paste0(output_path, "signature_cluster_bootstrap_with_zeros.RData"))
load(paste0(output_path, "signature_cluster_bootstrap_with_zeros.RData"))
################################################################################
# for each number of clusters plot their dendrogram statistics

boostrap_dendro_compare_all_long <- melt(boostrap_dendro_compare_all, id.vars = c("number_cluster", "iteration"))
boostrap_dendro_compare_all_long$number_cluster <- as.character(boostrap_dendro_compare_all_long$number_cluster)
boostrap_dendro_compare_all_long$number_cluster <- factor(boostrap_dendro_compare_all_long$number_cluster, levels = c(1:20))

plot_dendro_stat <- ggplot(boostrap_dendro_compare_all_long, aes(number_cluster, value))+
  geom_quasirandom(alpha = 0.5) + 
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~variable, ncol = 1, scale = "free_y")

pdf(paste0(output_path, "signature_cluster_robustness_bootstrap_dendrogram_stats.pdf"), width = 8, height = 6)
plot_dendro_stat
dev.off()

# plot the percentages of samples within the same cluster for all iterations
plot_hist <- ggplot(boostrap_sample_compare_all, aes(x = percent_samples_in_new_cluster_and_old_cluster)) + 
  geom_histogram(bins = 100) +
  facet_wrap(~number_of_clusters, scale = "free_y") +
  theme_bw()

# for each cluster plot their the percentage of each sample within that cluster
boostrap_sample_compare_all$old_cluster_id <- as.character(boostrap_sample_compare_all$old_cluster_id)
boostrap_sample_compare_all$old_cluster_id <- factor(boostrap_sample_compare_all$old_cluster_id, levels = c(1:20))

boostrap_sample_compare_all <- boostrap_sample_compare_all %>% group_by(number_of_clusters) %>% mutate(med = median(percent_samples_in_new_cluster_and_old_cluster))
boostrap_sample_compare_all <- boostrap_sample_compare_all %>% group_by(number_of_clusters) %>% mutate(mea = mean(percent_samples_in_new_cluster_and_old_cluster))

plot_cluster <- ggplot(boostrap_sample_compare_all, aes(old_cluster_id, percent_samples_in_new_cluster_and_old_cluster))+
  geom_boxplot() +
  geom_hline(data = boostrap_sample_compare_all, aes(yintercept = med), linetype = "dashed", colour = "#eb34cf") +
  geom_hline(data = boostrap_sample_compare_all, aes(yintercept = mea), linetype = "dashed", colour = "#45d921") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~number_of_clusters, scale = "free_x")

pdf(paste0(output_path, "signature_cluster_robustness_bootstrap_cluster_friends.pdf"), width = 14, height = 6)
plot_cluster
dev.off()

################################################################################
# also use Mark Hill's method of using the clusterboot function and calculating
# jaccard similarity

og_mat <- dcast(all_sigs_a, sample ~ label, value.var = "weight_z")
rownames(og_mat) <- og_mat$sample
og_mat$sample <- NULL

cluster.boot.fun <- function(matrix = og_mat,
                             iters = 500,
                             num_cluster = 1){
  
  sig_boot <- fpc::clusterboot(matrix, 
                               B = iters,
                               bootmethod = "boot",
                               clustermethod = hclustCBI,
                               method = "ward.D2",
                               k = num_cluster)
  
  plot_boot <- sig_boot$bootresult %>% 
    as.data.table() %>% 
    reshape2::melt()
  
  plot_boot$cluster <- rep.int( c(1:num_cluster), (nrow(plot_boot) / num_cluster)) 
  return(plot_boot)
}


all_boot_out <- list()
for(nclust in c(1:20)){
  print(nclust)
  out <- cluster.boot.fun(matrix = og_mat,
                          iters = 500,
                          num_cluster = nclust)
  
  out$num_cluster <- nclust
  all_boot_out[[nclust]] <- out
}

all_boot_out <- do.call(rbind, all_boot_out)

all_boot_out$cluster <- factor(all_boot_out$cluster, levels = str_sort(c(1:20), numeric = T))
all_boot_out$num_cluster <- factor(all_boot_out$num_cluster, levels = str_sort(c(1:20), numeric = T))

all_boot_out_p <- all_boot_out %>% group_by(num_cluster) %>% mutate(med = median(value))
all_boot_out_p <- all_boot_out_p %>% group_by(num_cluster) %>% mutate(mea = mean(value))

p_cluster_m <- ggplot(all_boot_out_p, aes(x = cluster, y = value)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 0.6, width = 0.1) +
  geom_hline(data = all_boot_out_p, aes(yintercept = med), linetype = "dashed", colour = "#eb34cf") +
  geom_hline(data = all_boot_out_p, aes(yintercept = mea), linetype = "dashed", colour = "#45d921") +
  theme_bw() +
  xlab("Cluster") +
  ylab("Jaccard Similarity") +
  theme(legend.position = "none") +
  facet_wrap(~num_cluster, scale = "free_x")

pdf(paste0(output_path, "signature_cluster_robustness_bootstrap_jaccard_index.pdf"), width = 14, height = 6)
p_cluster_m
dev.off()

