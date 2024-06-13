# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # #                      copy number cancer genes                       # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
###################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggplot2)
library(reshape2)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(cowplot)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
##################################################################### file paths

sample_table_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
copy_number_path              <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/combined_ascatTable_segments_joined.rds"
cytoband_path                 <- "/re_gecip/cancer_lung/kthol/general/cytobands_hg38_UCSC_sigminer.txt" #from UCSC via sigminer github
amp_bin_path                  <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/amplification_bin_copynumber.Rdata"
loss_bin_path                 <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/loss_bin_copynumber.Rdata"
blacklist_path                <- "/re_gecip/cancer_lung/code/Annotation/blacklist_hg38.bed"
olfactory_gene_path           <- "/re_gecip/cancer_lung/shared/olfactory_barnes_2020.csv"
tcga_expression_gene_path     <- "/re_gecip/cancer_lung/shared/TCGA_expression_in_tumorsubtypes.csv"
gtex_expression_gene_path     <- "/re_gecip/cancer_lung/shared/GTEx_expression_in_tumorsubtypes.csv"
cn_gene_path                  <- "/re_gecip/cancer_lung/kthol/general/cancer_driver_alex_subset_for_driver_only.txt"
centromer_path                <- "/re_gecip/cancer_lung/kthol/CopyNumberSignatures/functions/Macintyre/data/hg38_centromere_locs_ucsc_single.txt"


output_path                   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/"
plot_output_path              <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/genomic_disruption/"
driver_output_path            <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/"

bin_size <- 50000

###############################################################################
####################################################################  functions
# function to make genome continuous
genome_continous_fun <- function(data,
                                 genome){
  
  data                       <- data.frame(data)
  data$bin                   <- rownames(data)
  data$chr                   <- sub(":.*", "", data$bin)
  data$startpos              <- as.numeric(gsub(".*:(.+):.*", "\\1", data$bin))
  data$endpos                <- as.numeric(sub(".*:", "", data$bin))
  
  chr_sizes                 <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_sizes                 <- data.frame(chr_sizes)
  # chr_sizes$chr             <- gsub("chr", "", rownames(chr_sizes))
  chr_sizes[chr_sizes$chr == "chrX", "chr"] <- 23
  chr_sizes$chr             <- as.numeric(chr_sizes$chr)
  chr_sizes                 <- chr_sizes[order(chr_sizes$chr, decreasing = F), ]
  chr_sizes                 <- chr_sizes$chr_sizes
  
  #adjust pos to continous variable
  cumsum_chr_size               <- cumsum(as.numeric(chr_sizes))
  
  data$adjust_start_pos         <- data$startpos
  data$adjust_end_pos           <- data$endpos
  
  # change chr X and Y in 23 and 24
  data[data$chr == "chrX", "chr"] <- 23
  data[data$chr == "chrY", "chr"] <- 24
  data$chr <- sub("chr", "", data$chr)
  data$chr <- as.numeric(data$chr)
  
  for(i in c(2:24)){
    data$adjust_start_pos[data$chr == i] <- data$adjust_start_pos[data$chr == i] + cumsum_chr_size[i-1]
  }
  
  for(i in c(2:24)){
    data$adjust_end_pos[data$chr == i] <- data$adjust_end_pos[data$chr == i] + cumsum_chr_size[i-1]
  }
  
  data$chr <- paste0("chr", data$chr)
  
  return(data)
}

midpoint_continous_fun <- function(data,
                                   genome){
  
  chr_sizes                 <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_sizes                 <- data.frame(chr_sizes)
  chr_sizes$chr             <- gsub("chr", "", rownames(chr_sizes))
  chr_sizes$chr             <- as.numeric(chr_sizes$chr)
  chr_sizes                 <- chr_sizes[order(chr_sizes$chr, decreasing = F), ]
  chr_sizes                 <- chr_sizes$chr_sizes
  
  #adjust pos to continous variable
  cumsum_chr_size               <- cumsum(as.numeric(chr_sizes))
  
  data$adjust_mid         <- data$mid
  
  # change chr X and Y in 23 and 24
  data$seqnames <- sub("chr", "", data$seqnames)
  data$seqnames <- as.numeric(as.character(data$seqnames))
  
  for(i in c(2:22)){
    data$adjust_mid[data$seqnames == i] <- data$adjust_mid[data$seqnames == i] + cumsum_chr_size[i-1]
  }
  
  data$seqnames <- paste0("chr", data$seqnames)
  
  return(data)
}

anno_midpoint_continous_fun <- function(data,
                                        genome){
  
  chr_sizes                 <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_sizes                 <- data.frame(chr_sizes)
  chr_sizes$chr             <- gsub("chr", "", rownames(chr_sizes))
  chr_sizes$chr             <- as.numeric(chr_sizes$chr)
  chr_sizes                 <- chr_sizes[order(chr_sizes$chr, decreasing = F), ]
  chr_sizes                 <- chr_sizes$chr_sizes
  
  #adjust pos to continous variable
  cumsum_chr_size               <- cumsum(as.numeric(chr_sizes))
  
  data$adjust_mid         <- data$gene_midpoint
  
  # change chr X and Y in 23 and 24
  data$seqnames <- sub("chr", "", data$seqnames)
  data$seqnames <- as.numeric(as.character(data$seqnames))
  
  for(i in c(2:22)){
    data$adjust_mid[data$seqnames == i] <- data$adjust_mid[data$seqnames == i] + cumsum_chr_size[i-1]
  }
  
  data$seqnames <- paste0("chr", data$seqnames)
  
  return(data)
}

# function to bin segments
bin_fun <- function(data,
                    bin_size,
                    genome){
  
  
  #divide genome into bins
  chr_length <- seqlengths(genome)[paste0("chr", c(1:22))]
  bins_df <- lapply(names(chr_length), function(x){
    seq_bin <- c(seq(from = 1, to = chr_length[x], by = bin_size), as.numeric(chr_length[x]))
    data.frame(chr = x, start = seq_bin[1:(length(seq_bin)-1)], end = seq_bin[2:length(seq_bin)]-1)
  })
  
  bins_df     <- Reduce(rbind, bins_df)
  bins_gr     <- makeGRangesFromDataFrame(bins_df)
  
  seg_bin_matrix <- matrix(0, nrow = length(bins_gr) ,ncol = length(unique(data$patient)),
                           dimnames = list(paste(bins_df[,1], bins_df[,2], bins_df[,3], sep = ':'),
                                           unique(data$patient)))
  
  # check for each segment of each sample if they are overlapping a bin
  for(sample in unique(data$patient)){
    sub_seg               <- data[data$patient == sample,]
    sub_seg_gr            <- GRanges(seqnames = sub_seg$chr, IRanges(start = sub_seg$startpos, end = sub_seg$endpos))
    overlap               <- findOverlaps(bins_gr, sub_seg_gr)
    intersect_gr          <- pintersect(bins_gr[queryHits(overlap)], sub_seg_gr[subjectHits(overlap)])
    values(intersect_gr)  <- data.frame('bin' = paste(bins_df[,1], bins_df[,2], bins_df[,3], sep = ':')[queryHits(overlap)])
    
    intersect_df          <- data.frame(intersect_gr)
    
    sample                <- as.character(sample)
    seg_bin_matrix[intersect_df$bin, sample] <- 1
    
  }
  
  return(seg_bin_matrix)
  
}

genome_plot_fun <- function(sample_bins,
                            sample_size,
                            genome){
  
  # segment matrix for one of the timings
  
  # for each bin, how many samples have a segment in this bin
  bin_sum <- data.frame(rowSums(sample_bins))
  colnames(bin_sum) <- "samples_per_bin"
  
  # calculate proportion
  bin_sum$sample_per_bin_prop <- bin_sum$samples_per_bin/sample_size
  bin_sum$bin                 <- rownames(bin_sum)
  
  # make continuous
  bin_sum_cont  <- genome_continous_fun(bin_sum,
                                        genome)
  
  chr_length    <- seqlengths(genome)[paste0('chr', c(1:22, "X", "Y"))]
  chr_lines     <- cumsum(as.numeric(chr_length))
  chr_midpoint  <- chr_lines - chr_length/2
  
  # plot
  
  plot_objects <- list(bin_sum_cont, chr_lines, chr_midpoint)
  
  return(plot_objects)
}


################################################################################
######################################################################      MAIN

# load sample table
sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

# load regions to exclude later
olfactory_genes <- read.table(olfactory_gene_path, head = T)
olfactory_genes <- unique(olfactory_genes$gene_name)

tcga_expression_gene <- read.table(tcga_expression_gene_path, head = T)
gtex_expression_gene <- read.table(gtex_expression_gene_path, head = T)

# load copy number data
copynumber <- readRDS(copy_number_path)

# make an indicator saying whether something is an amplification (double ploidy + 1)
copynumber[copynumber$cnTotal >= ((copynumber$Ploidy*2) + 1), "amplification"] <- TRUE
copynumber[copynumber$cnTotal < ((copynumber$Ploidy*2) + 1), "amplification"] <- FALSE

# and indicate homozygous deletions
copynumber[copynumber$cnTotal == 0, "homdel"] <- TRUE
copynumber[copynumber$cnTotal > 0 , "homdel"] <- FALSE

# only keep the amplifications
amps <- copynumber[copynumber$amplification, ]

# add histology to the data frame
amps <- left_join(amps, sample_table[, c("participant_id", "histology")], by = c("patient" = "participant_id"))

# count number of samples per histology
hist_count <- data.frame(table(sample_table$histology))
colnames(hist_count) <- c("histology", "hist_count")
all_histologies <- c(unique(sample_table$histology), "PANCAN")
all_histologies <- grep("OTHER", all_histologies, value = T, invert = T)

################################################################################
# get the number of amplifications and losses per bin

# amplifications

# amp_plot_data <- list()
# amp_plot_chr_lines <- list()
# amp_plot_chr_midpoint <- list()
# 
# for(hist in all_histologies){
# 
#   if(hist == "PANCAN"){
#     cn_df <- amps
#     sample_size <- length(unique(sample_table$participant_id))
#   } else {
#     cn_df <- amps[amps$histology == hist, ]
#     sample_size <- hist_count[hist_count$histology == hist, "hist_count"]
#   }
# 
#   cn_bins <- bin_fun(cn_df,
#                      50000,
#                      BSgenome.Hsapiens.UCSC.hg38)
# 
#   genome_plot_data <- genome_plot_fun(cn_bins,
#                                       sample_size,
#                                       BSgenome.Hsapiens.UCSC.hg38)
# 
#   amp_plot_data[[hist]] <- genome_plot_data[[1]]
#   amp_plot_chr_lines[[hist]] <- genome_plot_data[[2]]
#   amp_plot_chr_midpoint[[hist]] <- genome_plot_data[[3]]
# 
# }
# 
# save(amp_plot_data,amp_plot_chr_lines,amp_plot_chr_midpoint, file = paste0(output_path, "amplification_bin_copynumber.Rdata"))

load(paste0(output_path, "amplification_bin_copynumber.Rdata"))
amp_bin_sum_cont <- do.call(rbind, amp_plot_data)
amp_bin_sum_cont$histology <- rownames(amp_bin_sum_cont)
amp_bin_sum_cont$histology <- sub("\\..*", "", amp_bin_sum_cont$histology)

##### losses

# only keep the homozygous deletions
loss <- copynumber[copynumber$homdel, ]

# add histology to the data frame
loss <- left_join(loss, sample_table[, c("participant_id", "histology")], by = c("patient" = "participant_id"))

# loss_plot_data <- list()
# loss_plot_chr_lines <- list()
# loss_plot_chr_midpoint <- list()
# 
# for(hist in all_histologies){
# 
#   if(hist == "PANCAN"){
#     cn_df <- loss
#     sample_size <- length(unique(sample_table$participant_id))
#   } else {
#     cn_df <- loss[loss$histology == hist, ]
#     sample_size <- hist_count[hist_count$histology == hist, "hist_count"]
#   }
# 
#   cn_bins <- bin_fun(cn_df,
#                      50000,
#                      BSgenome.Hsapiens.UCSC.hg38)
# 
# 
#   genome_plot_data <- genome_plot_fun(cn_bins,
#                                       sample_size,
#                                       BSgenome.Hsapiens.UCSC.hg38)
# 
#   loss_plot_data[[hist]] <- genome_plot_data[[1]]
#   loss_plot_chr_lines[[hist]] <- genome_plot_data[[2]]
#   loss_plot_chr_midpoint[[hist]] <- genome_plot_data[[3]]
# 
# }
# 
# save(loss_plot_data,loss_plot_chr_lines,loss_plot_chr_midpoint, file = paste0(output_path, "loss_bin_copynumber.Rdata"))

load(paste0(output_path, "loss_bin_copynumber.Rdata"))
loss_bin_sum_cont <- do.call(rbind, loss_plot_data)
loss_bin_sum_cont$histology <- rownames(loss_bin_sum_cont)
loss_bin_sum_cont$histology <- sub("\\..*", "", loss_bin_sum_cont$histology)

################################################################################

# calculate a p value for each bin

chr_length <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0("chr", c(1:22))]
bins_df <- lapply(names(chr_length), function(x){
  seq_bin <- c(seq(from = 1, to = chr_length[x], by = bin_size), as.numeric(chr_length[x]))
  data.frame(chr = x, start = seq_bin[1:(length(seq_bin)-1)], end = seq_bin[2:length(seq_bin)]-1)
})

bins_df     <- Reduce(rbind, bins_df)

threshold <- list()
amp_bin_sum_cont_stat <- list()
loss_bin_sum_cont_stat <- list()

for(hist in all_histologies){
  
  print(hist)
  
  if(hist == "PANCAN"){
    num_amps    <- nrow(amps)
    num_losses  <- nrow(loss)
    
  } else{
    num_amps    <- nrow(amps[amps$histology == hist, ])
    num_losses  <- nrow(loss[loss$histology == hist, ])
  }
  
  num_bins    <- nrow(bins_df)
  
  alpha <- 0.05
  alpha_adjusted <- alpha/num_bins
  
  # GAIN significance threshold
  poisson_amps_lambda <- num_amps/num_bins
  amps_threshold <- qpois(1 - alpha_adjusted, poisson_amps_lambda)
  amps_threshold_p <- ppois(max(amp_bin_sum_cont[amp_bin_sum_cont$histology == hist, "samples_per_bin"]), poisson_amps_lambda,lower.tail = FALSE)
  
  # LOSS significance threshold
  poisson_loss_lambda <- num_losses/num_bins
  loss_threshold <- qpois(1 - alpha_adjusted, poisson_loss_lambda)
  loss_threshold_p <- ppois(max(loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, "samples_per_bin"]), poisson_loss_lambda,lower.tail = FALSE)

  df_threshold <- data.frame(histology = hist, 
                             amplification_threshold = amps_threshold,
                             amplification_threshold_p = amps_threshold_p,
                             losses_threshold = loss_threshold,
                             losses_threshold_p = loss_threshold_p)
  
  threshold[[hist]] <- df_threshold
  
  # now get the p value for every bin
  
  # AMPS
  hist_amp_bin_sum_cont <- amp_bin_sum_cont[amp_bin_sum_cont$histology == hist, ]
  hist_amp_bin_sum_cont <- hist_amp_bin_sum_cont %>%
    mutate(lambda = poisson_amps_lambda) %>%
    mutate(p_val = ppois(samples_per_bin -1 , lambda, lower.tail = FALSE)) 
  
  hist_amp_bin_sum_cont$p_adjust <- p.adjust(hist_amp_bin_sum_cont$p_val, method = "fdr")
  
  # LOSSES
  loss_amp_bin_sum_cont <- loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, ]
  loss_amp_bin_sum_cont <- loss_amp_bin_sum_cont %>%
    mutate(lambda = poisson_loss_lambda) %>%
    mutate(p_val = ppois(samples_per_bin -1 , lambda, lower.tail = FALSE)) 
  
  loss_amp_bin_sum_cont$p_adjust <- p.adjust(loss_amp_bin_sum_cont$p_val, method = "fdr")
  
  amp_bin_sum_cont_stat[[hist]] <- hist_amp_bin_sum_cont
  loss_bin_sum_cont_stat[[hist]] <- loss_amp_bin_sum_cont
}

threshold <- do.call(rbind, threshold)

write.table(threshold, paste0(output_path, "CN_bins_thresholds.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

amp_bin_sum_cont_stat <- do.call(rbind, amp_bin_sum_cont_stat)
loss_bin_sum_cont_stat <- do.call(rbind, loss_bin_sum_cont_stat)

################################################################################

# get the bins / genes which are significant

# get peaks to annotate
cn_genes <- read.table(cn_gene_path, head = T, sep = "\t")
cancer_genes <- unique(cn_genes$gene_name)

# gain_peaks <- list()
# for(hist in unique(amp_bin_sum_cont_stat$histology)){
#   df <- amp_bin_sum_cont_stat[amp_bin_sum_cont_stat$histology == hist, ]
# 
#   peaks <- list()
#   for(chr in unique(df$chr)){
#     chr_df <- df[df$chr == chr, ]
# 
#     # overlap the bins with genes
#     # load gene data
#     txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
#     genes                 <- genes(txdb, single.strand.genes.only=FALSE)
#     genes                 <- data.frame(genes)
# 
#     # extract HUGO symbols
#     gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)
# 
#     # merge with gene locations
#     gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))
# 
#     # gene_names <- gene_names[gene_names$symbol %in% driver$gene, ]
#     gene_names <- gene_names[!is.na(gene_names$symbol), ]
# 
#     # make into genomic ranges
#     gene_gr <- GRanges(seqnames = gene_names$seqnames, IRanges(start = gene_names$start, end = gene_names$end))
#     mcols(gene_gr) <- data.frame(gene_names[, c("symbol")])
# 
#     bin_gr <- GRanges(seqnames = chr_df$chr, IRanges(start = chr_df$startpos, end = chr_df$endpos))
#     mcols(bin_gr) <- data.frame(chr_df[, c("samples_per_bin", "sample_per_bin_prop", "bin", "histology", "lambda", "p_val", "p_adjust")])
# 
#     overlap         <-  findOverlaps(gene_gr, bin_gr)
#     intersect_gr    <-  pintersect(bin_gr[subjectHits(overlap)], gene_gr[queryHits(overlap)])
#     mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(gene_gr[queryHits(overlap)]))
# 
#     intersect_df <- data.frame(intersect_gr)
#     intersect_df$strand <- NULL
#     intersect_df$width <- NULL
#     intersect_df$hit <- NULL
#     colnames(intersect_df)[ncol(intersect_df)] <- "gene"
# 
#     # for each gene take the highest p value
#     intersect_df <- intersect_df %>% group_by(gene) %>% mutate(highest_adjusted_p = max(p_adjust))
#     intersect_df <- data.frame(intersect_df)
# 
#     # only take significant bin/genes forward
#     top <- intersect_df[intersect_df$highest_adjusted_p < 0.05, ]
# 
#     # make this unique per gene
# 
#     if(nrow(top) > 0){
# 
#       # take only one value for each gene
#       one_bin_gene_df <- list()
#       for(gen in unique(top$gene)){
#         gene_df <- top[top$gene == gen, ]
# 
#         gene_df <- gene_df[gene_df$sample_per_bin_prop == max(gene_df$sample_per_bin_prop), ]
# 
#         # what if some still have the same sample per bin prop? let's take the midpoint of those
#         mid_point <- min(gene_df$start) + ((max(gene_df$end) - min(gene_df$start))/2)
#         gene_df$start <- mid_point
#         gene_df$end <- mid_point
# 
#         # we need to keep a bin identifier but now can be for several bins
#         gene_df$binID <- paste(unique(gene_df$bin), collapse = "_")
# 
#         gene_df$bin <- NULL
#         gene_df <- unique(gene_df)
# 
#         one_bin_gene_df[[gen]] <- gene_df
# 
#       }
# 
#       one_bin_gene_df <- do.call(rbind, one_bin_gene_df)
# 
#       # let's cut this off based on our histology threshold
#       top_peak <- one_bin_gene_df[one_bin_gene_df$samples_per_bin > threshold[threshold$histology == hist, "amplification_threshold"], ]
# 
#       if(nrow(top_peak) > 0){
#         # now we can say wheter there are cancer genes in here
#         top_peak <- top_peak[order(top_peak$sample_per_bin_prop, decreasing = T), ]
#         top_peak[top_peak$gene %in% cancer_genes, "is_cancer_gene"] <- TRUE
#         top_peak[is.na(top_peak$is_cancer_gene), "is_cancer_gene"] <- FALSE
# 
#         # let's leave it at this, we can cut things off more for plotting later
#         if(nrow(top_peak) > 0){
#           top_peak$sample_per_bin_prop <- as.numeric(as.character(top_peak$sample_per_bin_prop))
#           top_peak <- top_peak[order(top_peak$sample_per_bin_prop, decreasing = T),]
# 
#           # make this continous
#           top_peak$mid <- top_peak$start
#           top_peak <- midpoint_continous_fun(top_peak, BSgenome.Hsapiens.UCSC.hg38)
# 
#           peaks[[chr]] <- top_peak
#         }
#       }
#     }
#   }
#   peaks <- do.call(rbind, peaks)
#   gain_peaks[[hist]] <- peaks
# }
# 
# gain_peaks <- do.call(rbind, gain_peaks)
# save(gain_peaks, file = paste0(output_path, "gain_peaks.Rdata"))
load(paste0(output_path, "gain_peaks.Rdata"))

# loss_peaks <- list()
# for(hist in unique(loss_bin_sum_cont_stat$histology)){
#   df <- loss_bin_sum_cont_stat[loss_bin_sum_cont_stat$histology == hist, ]
# 
#   peaks <- list()
#   for(chr in unique(df$chr)){
#     chr_df <- df[df$chr == chr, ]
# 
#     # overlap the bins with genes
#     # load gene data
#     txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
#     genes                 <- genes(txdb, single.strand.genes.only=FALSE)
#     genes                 <- data.frame(genes)
# 
#     # extract HUGO symbols
#     gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)
# 
#     # merge with gene locations
#     gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))
# 
#     # gene_names <- gene_names[gene_names$symbol %in% driver$gene, ]
#     gene_names <- gene_names[!is.na(gene_names$symbol), ]
# 
#     # make into genomic ranges
#     gene_gr <- GRanges(seqnames = gene_names$seqnames, IRanges(start = gene_names$start, end = gene_names$end))
#     mcols(gene_gr) <- data.frame(gene_names[, c("symbol")])
# 
#     bin_gr <- GRanges(seqnames = chr_df$chr, IRanges(start = chr_df$startpos, end = chr_df$endpos))
#     mcols(bin_gr) <- data.frame(chr_df[, c("samples_per_bin", "sample_per_bin_prop", "bin", "histology", "lambda", "p_val", "p_adjust")])
# 
#     overlap         <-  findOverlaps(gene_gr, bin_gr)
#     intersect_gr    <-  pintersect(bin_gr[subjectHits(overlap)], gene_gr[queryHits(overlap)])
#     mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(gene_gr[queryHits(overlap)]))
# 
#     intersect_df <- data.frame(intersect_gr)
#     intersect_df$strand <- NULL
#     intersect_df$width <- NULL
#     intersect_df$hit <- NULL
#     colnames(intersect_df)[ncol(intersect_df)] <- "gene"
# 
#     # for each gene take the highest p value
#     intersect_df <- intersect_df %>% group_by(gene) %>% mutate(highest_adjusted_p = max(p_adjust))
#     intersect_df <- data.frame(intersect_df)
# 
#     # only take significant bin/genes forward
#     top <- intersect_df[intersect_df$highest_adjusted_p < 0.05, ]
# 
#     # make this unique per gene
# 
#     if(nrow(top) > 0){
# 
#       # take only one value for each gene
#       one_bin_gene_df <- list()
#       for(gen in unique(top$gene)){
#         gene_df <- top[top$gene == gen, ]
# 
#         gene_df <- gene_df[gene_df$sample_per_bin_prop == max(gene_df$sample_per_bin_prop), ]
# 
#         # what if some still have the same sample per bin prop? let's take the midpoint of those
#         mid_point <- min(gene_df$start) + ((max(gene_df$end) - min(gene_df$start))/2)
#         gene_df$start <- mid_point
#         gene_df$end <- mid_point
# 
#         # we need to keep a bin identifier but now can be for several bins
#         gene_df$binID <- paste(unique(gene_df$bin), collapse = "_")
#         gene_df$bin <- NULL
#         gene_df <- unique(gene_df)
# 
#         one_bin_gene_df[[gen]] <- gene_df
# 
#       }
#       one_bin_gene_df <- do.call(rbind, one_bin_gene_df)
# 
#       # let's cut this off based on our histology threshold
#       top_peak <- one_bin_gene_df[one_bin_gene_df$samples_per_bin > threshold[threshold$histology == hist, "losses_threshold"], ]
# 
#       if(nrow(top_peak) > 0){
#         # now we can say wheter there are cancer genes in here
#         top_peak <- top_peak[order(top_peak$sample_per_bin_prop, decreasing = T), ]
#         top_peak[top_peak$gene %in% cancer_genes, "is_cancer_gene"] <- TRUE
#         top_peak[is.na(top_peak$is_cancer_gene), "is_cancer_gene"] <- FALSE
# 
#         # let's leave it at this, we can cut things off more for plotting later
#         if(nrow(top_peak) > 0){
#           top_peak$sample_per_bin_prop <- as.numeric(as.character(top_peak$sample_per_bin_prop))
#           top_peak <- top_peak[order(top_peak$sample_per_bin_prop, decreasing = T),]
# 
#           # make this continous
#           top_peak$mid <- top_peak$start
#           top_peak <- midpoint_continous_fun(top_peak, BSgenome.Hsapiens.UCSC.hg38)
# 
#           peaks[[chr]] <- top_peak
#         }
#       }
#     }
#   }
#   peaks <- do.call(rbind, peaks)
#   loss_peaks[[hist]] <- peaks
# }
# loss_peaks <- do.call(rbind, loss_peaks)
# save(loss_peaks, file = paste0(output_path, "loss_peaks.Rdata"))
load(paste0(output_path, "loss_peaks.Rdata"))

# get a list of copy number drivers, if there is a known driver gene in a bin then only take that gene
gain_peak_genes <- gain_peaks

cancer_gain_df <- list()
for(hist in unique(gain_peak_genes$histology)){
  hist_df <- gain_peak_genes[gain_peak_genes$histology == hist, ]
  
  all_bin_df <- list()
  for(bin in unique(gain_peak_genes$binID)){
    
    df <- hist_df[grepl(gsub("\\+", "\\\\\\+", bin), hist_df$binID), ]
    
    if(any(df$is_cancer_gene)){
      df_out <- df[df$is_cancer_gene, ]
    } else { df_out <- df}
    
    all_bin_df[[bin]] <- df_out 
    
  }
  all_bin_df_conc <- do.call(rbind, all_bin_df)
  cancer_gain_df[[hist]] <- all_bin_df_conc
}
cancer_gain_df <- do.call(rbind, cancer_gain_df)
cancer_gain_df <- unique(cancer_gain_df)

# and losses 
loss_peak_genes <- loss_peaks

cancer_loss_df <- list()
for(hist in unique(loss_peak_genes$histology)){
  hist_df <- loss_peak_genes[loss_peak_genes$histology == hist, ]
  
  all_bin_df <- list()
  for(bin in unique(loss_peak_genes$binID)){
    
    df <- hist_df[grepl(gsub("\\+", "\\\\\\+", bin), hist_df$binID), ]
    
    if(any(df$is_cancer_gene)){
      df_out <- df[df$is_cancer_gene, ]
    } else { df_out <- df}
    
    all_bin_df[[bin]] <- df_out 
    
  }
  all_bin_df_conc <- do.call(rbind, all_bin_df)
  cancer_loss_df[[hist]] <- all_bin_df_conc
}
cancer_loss_df <- do.call(rbind, cancer_loss_df)
cancer_loss_df <- unique(cancer_loss_df)

# i need to figure out if the segments of the same samples overlap those genes in adjacent bins where one gene is a cancer gene
# get centromere locations
centromers <- read.table(centromer_path, head = TRUE, sep = "\t")
centromers$chrom <- sub("chr", "", centromers$chrom)

# make a data frame with the coordinates of p and q arms
# add chromosome length data
chr_length    <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0('chr', c(1:22, "X", "Y"))]
chr_length    <- data.frame(chr_length)
chr_length$chr <- rownames(chr_length)

arm_locs <- chr_length[, c("chr", "chr_length")]
arm_locs$chr <- sub("chr", "", arm_locs$chr)

# add centromere locations
colnames(centromers)  <- c("chr", "centromere_start", "centromere_end")
arm_locs              <- left_join(arm_locs, centromers)
arm_locs$p_arm_start  <- 1
arm_locs$p_arm_end    <- arm_locs$centromere_start
arm_locs$q_arm_start  <- arm_locs$centromere_end
arm_locs$q_arm_end    <- arm_locs$chr_length

# add general chromosome start and stop
arm_locs$chr_start    <- 1
arm_locs$chr_end      <- arm_locs$chr_length

# make into long format
arm_locs_long <- arm_locs[, c("chr", "p_arm_start", "p_arm_end", "q_arm_start", "q_arm_end")]
arm_locs_long <- melt(arm_locs_long, id.vars = c("chr", "p_arm_end", "q_arm_end"))
colnames(arm_locs_long) <- c("chr", "p_arm_end", "q_arm_end", "arm", "start")
arm_locs_long <- melt(arm_locs_long, id.vars = c("chr", "arm", "start"))
arm_locs_long$arm <- sub("_.*", "", arm_locs_long$arm)
arm_locs_long$variable <- sub("_.*", "", arm_locs_long$variable)
arm_locs_long <- arm_locs_long[arm_locs_long$arm == arm_locs_long$variable, ]
arm_locs_long$variable <- NULL
colnames(arm_locs_long)[ncol(arm_locs_long)] <- "end"
arm_locs_long$chr_arm <- paste0(arm_locs_long$chr, arm_locs_long$arm)
arm_locs_long$chr <- paste0("chr", arm_locs_long$chr)
arm_locs_long <- arm_locs_long[arm_locs_long$chr != "chrY", ]
arm_locs_long <- arm_locs_long[arm_locs_long$chr != "chrX", ]

# amplifications first
# do this separately for each histology
amplification_peak_genes_patient <- list()
for(hist in all_histologies){
  
  print(hist)
  
  peak_df <- cancer_gain_df[cancer_gain_df$histology == hist, ]
  
  peak_gr <- GRanges(seqnames = peak_df$seqnames, IRanges(start = peak_df$start, end = peak_df$end))
  mcols(peak_gr) <- data.frame(peak_df[, c("samples_per_bin", "sample_per_bin_prop", "histology", "lambda", "p_val", "p_adjust", "gene", "highest_adjusted_p", "binID", "is_cancer_gene", "mid", "adjust_mid")])
  
  cn_gr <- GRanges(seqnames = amps$chr, IRanges(start = amps$startpos, end = amps$endpos))
  mcols(cn_gr) <- data.frame(amps)
  
  overlap         <-  findOverlaps(cn_gr, peak_gr)
  intersect_gr    <-  pintersect(peak_gr[subjectHits(overlap)], cn_gr[queryHits(overlap)])
  mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(cn_gr[queryHits(overlap)]))
  
  intersect_df <- data.frame(intersect_gr)
  intersect_df$strand <- NULL
  intersect_df$width <- NULL
  intersect_df$hit <- NULL
  
  colnames(intersect_df) <- c("overlap_chr", "overlap_start", "overlap_end", "samples_per_bin", 
                              "sample_per_bin_prop", "histology", "lambda", "p_val", "p_adjust", 
                              "gene", "highest_adjusted_p", "binID", "is_cancer_gene", "mid", 
                              "adjust_mid", "patient", "cn_chr", "cn_start", "cn_end", "cnTotal", 
                              "nMajor", "nMinor", "Ploidy", "ACF", "amplification", "homdel", "cn_histology")
  
  # only keep the copy number of the patients with the same histology as the peaks come from
  if(hist == "PANCAN"){
    peak_cn <- intersect_df
  } else {
    peak_cn <- intersect_df[intersect_df$histology == intersect_df$cn_histology, ]
  }
  
  # let's be more stringent about this and only take those genes which are on the tip of the peak
  hist_bin_prop_mean <- quantile(amp_bin_sum_cont[amp_bin_sum_cont$histology == hist & amp_bin_sum_cont$samples_per_bin > threshold[threshold$histology == hist, "amplification_threshold"], "sample_per_bin_prop"], probs = c(0.75))
  
  peak_cn <- peak_cn[which(peak_cn$sample_per_bin_prop >= hist_bin_prop_mean), ]
  
  if(nrow(peak_cn) > 0){
    # now for each copy number segments check:
    # does this segment overlap a cancer gene? if so do not count the other genes in this segment
    # but tis will still keep a lot of passenger genes
    # if there is a cancer gene on a segment but no cancer genes on the adjacent segments which are on the same chromosome arm, also remove those
    # so let's do this for each patient separately
    patient_peak_df <- list()
    for(samp in unique(peak_cn$patient)){
      
      patient_df <- peak_cn[peak_cn$patient == samp, ]
      # get all the segments of this patient on a given chromosome arm 
      cancer_arm_segs <- list()
      for(arm in unique(arm_locs_long$chr_arm)){
        
        arm_info <- arm_locs_long[arm_locs_long$chr_arm == arm, ]
        
        arm_df <- patient_df[patient_df$cn_chr == arm_info$chr, ]
        
        if(nrow(arm_df) > 0){
          arm_df <- arm_df[arm_df$cn_end >= arm_info$start & arm_df$cn_start <= arm_info$end, ]
          
          # is there a cancer gene on this chromosome arm?
          if(any(arm_df$is_cancer_gene)){
            cancer_seg <- arm_df[arm_df$is_cancer_gene, ]
          } else { cancer_seg <- arm_df}
          if(nrow(cancer_seg) > 0){
            cancer_seg$chr_arm <- arm
          }
          cancer_arm_segs[[arm]] <- cancer_seg
        }
      }
      cancer_arm_segs <- do.call(rbind, cancer_arm_segs)
      patient_peak_df[[samp]] <-cancer_arm_segs 
      
    }
    patient_peak_df <- do.call(rbind, patient_peak_df)
    
    # how often do we see copy number amplifications in the genes which are left now
    patient_peak_genes <- patient_peak_df[, c("patient", "cn_chr", "cn_start", "cn_end", "cnTotal", "nMajor", "nMinor", "chr_arm", "gene", "is_cancer_gene", "samples_per_bin", "sample_per_bin_prop")]
    patient_peak_genes <- unique(patient_peak_genes)
    
    gene_count <- unique(patient_peak_genes[, c("patient", "gene")])
    gene_count <- data.frame(table(gene_count$gene))
    
    # lot's of these genes actually only have ampliffied copy number in one or two patients, so let's kick them out again based on our threshold
    genes_to_keep <- gene_count[gene_count$Freq > threshold[threshold$histology == hist, "amplification_threshold"], "Var1"]
    
    
    if(length(genes_to_keep) > 0){
      patient_peak_genes <- patient_peak_genes[patient_peak_genes$gene %in% genes_to_keep, ]
      other_genes_to_keep <- unique(patient_peak_genes[, "gene"])
      # for some chromosome arms there is still rubbish in here so let's investigate those a bit further
      chr_arms_false <- unique(patient_peak_genes[patient_peak_genes$is_cancer_gene == FALSE, "chr_arm"])
      
      if(length(chr_arms_false) > 0){
        other_genes_to_keep <- unique(patient_peak_genes[-which(patient_peak_genes$chr_arm %in% chr_arms_false), "gene"])
        new_genes_to_keep <- list()
        for(this_arm in chr_arms_false){
          
          false_chr_arm <- patient_peak_genes[patient_peak_genes$chr_arm == this_arm, ]
          
          if(length(unique(false_chr_arm$gene)) > 1){
            # how many tumours have amplifications in these genes?
            gene_count <- data.frame(table(false_chr_arm$gene))
            
            # are any of them cancer genes?
            gene_count<- left_join(gene_count, unique(patient_peak_genes[, c("gene", "is_cancer_gene")]), by = c("Var1" = "gene"))
            
            if(any(gene_count$is_cancer_gene)){
              new_genes_to_keep[[this_arm]] <- gene_count[gene_count$is_cancer_gene, "Var1"]
            } else { 
              count_median <- median(gene_count$Freq)
              new_genes_to_keep[[this_arm]] <- gene_count[gene_count$Freq > count_median, "Var1"]
            }
          } else {new_genes_to_keep[[this_arm]] <- unique(false_chr_arm$gene)}
          
          
        }
        new_genes_to_keep <- unlist(new_genes_to_keep)
        all_genes_to_keep <- c(new_genes_to_keep, other_genes_to_keep)
      } else {all_genes_to_keep <- other_genes_to_keep}
      
      patient_peak_genes_final <- patient_peak_genes[patient_peak_genes$gene %in% all_genes_to_keep, ]
      
      if(nrow(patient_peak_genes_final) > 0){
        # last thing I need to do is filter out some genes
        patient_peak_genes_final_filter <- patient_peak_genes_final %>% filter(!gene %in% olfactory_genes)
        
        hist_tcga_expression_genes <- tcga_expression_gene[, c("gene_name", hist)]
        colnames(hist_tcga_expression_genes) <- c("gene_name", "expression")
        hist_tcga_expression_genes <- unique(hist_tcga_expression_genes[hist_tcga_expression_genes$expression == 0, ])
        patient_peak_genes_final_filter <- patient_peak_genes_final_filter %>% filter(!gene %in% hist_tcga_expression_genes)
        
        hist_gtex_expression_genes <- gtex_expression_gene[, c("gene_name", hist)]
        colnames(hist_gtex_expression_genes) <- c("gene_name", "expression")
        hist_gtex_expression_genes <- unique(hist_gtex_expression_genes[hist_gtex_expression_genes$expression == 0, "gene_name"])
        patient_peak_genes_final_filter <- patient_peak_genes_final_filter %>% filter(!gene %in% hist_gtex_expression_genes)
        
        patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep("LOC", patient_peak_genes_final_filter$gene, invert = T), ]
        patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep(".*AS1", patient_peak_genes_final_filter$gene, invert = T), ]
        patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep("NEAT1", patient_peak_genes_final_filter$gene, invert = T), ]
        patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep("MALAT", patient_peak_genes_final_filter$gene, invert = T), ]
        
        patient_peak_genes_final_filter$histology <- hist
        amplification_peak_genes_patient[[hist]] <- patient_peak_genes_final_filter
      }
    }
  }
}

amplification_peak_genes_patient <- do.call(rbind, amplification_peak_genes_patient)

###### now do the same for losses

loss_peak_genes_patient <- list()
for(hist in all_histologies){
  
  print(hist)
  
  peak_df <- cancer_loss_df[cancer_loss_df$histology == hist, ]
  
  peak_gr <- GRanges(seqnames = peak_df$seqnames, IRanges(start = peak_df$start, end = peak_df$end))
  mcols(peak_gr) <- data.frame(peak_df[, c("samples_per_bin", "sample_per_bin_prop", "histology", "lambda", "p_val", "p_adjust", "gene", "highest_adjusted_p", "binID", "is_cancer_gene", "mid", "adjust_mid")])
  
  cn_gr <- GRanges(seqnames = loss$chr, IRanges(start = loss$startpos, end = loss$endpos))
  mcols(cn_gr) <- data.frame(loss)
  
  overlap         <-  findOverlaps(cn_gr, peak_gr)
  intersect_gr    <-  pintersect(peak_gr[subjectHits(overlap)], cn_gr[queryHits(overlap)])
  mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(cn_gr[queryHits(overlap)]))
  
  intersect_df <- data.frame(intersect_gr)
  
  if(nrow(intersect_df) > 0){
    intersect_df$strand <- NULL
    intersect_df$width <- NULL
    intersect_df$hit <- NULL
    
    colnames(intersect_df) <- c("overlap_chr", "overlap_start", "overlap_end", "samples_per_bin", 
                                "sample_per_bin_prop", "histology", "lambda", "p_val", "p_adjust", 
                                "gene", "highest_adjusted_p", "binID", "is_cancer_gene", "mid", 
                                "adjust_mid", "patient", "cn_chr", "cn_start", "cn_end", "cnTotal", 
                                "nMajor", "nMinor", "Ploidy", "ACF", "amplification", "homdel", "cn_histology")
    
    # only keep the copy number of the patients with the same histology as the peaks come from
    # only keep the copy number of the patients with the same histology as the peaks come from
    if(hist == "PANCAN"){
      peak_cn <- intersect_df
    } else {
      peak_cn <- intersect_df[intersect_df$histology == intersect_df$cn_histology, ]
    }
    
    # let's be more stringent about this and only take those genes which are on the tip of the peak
    hist_bin_prop_mean <- quantile(loss_bin_sum_cont[loss_bin_sum_cont$histology == hist & loss_bin_sum_cont$samples_per_bin > threshold[threshold$histology == hist, "losses_threshold"], "sample_per_bin_prop"], probs = c(0.75))
    
    peak_cn <- peak_cn[which(peak_cn$sample_per_bin_prop >= hist_bin_prop_mean), ]
    
    if(nrow(peak_cn) > 0){
      # now for each copy number segments check:
      # does this segment overlap a cancer gene? if so do not count the other genes in this segment
      # but tis will still keep a lot of passenger genes
      # if there is a cancer gene on a segment but no cancer genes on the adjacent segments which are on the same chromosome arm, also remove those
      # so let's do this for each patient separately
      patient_peak_df <- list()
      for(samp in unique(peak_cn$patient)){
        
        patient_df <- peak_cn[peak_cn$patient == samp, ]
        # get all the segments of this patient on a given chromosome arm 
        cancer_arm_segs <- list()
        for(arm in unique(arm_locs_long$chr_arm)){
          
          arm_info <- arm_locs_long[arm_locs_long$chr_arm == arm, ]
          
          arm_df <- patient_df[patient_df$cn_chr == arm_info$chr, ]
          
          if(nrow(arm_df) > 0){
            arm_df <- arm_df[arm_df$cn_end >= arm_info$start & arm_df$cn_start <= arm_info$end, ]
            
            # is there a cancer gene on this chromosome arm?
            if(any(arm_df$is_cancer_gene)){
              cancer_seg <- arm_df[arm_df$is_cancer_gene, ]
            } else { cancer_seg <- arm_df}
            if(nrow(cancer_seg) > 0){
              cancer_seg$chr_arm <- arm
            }
            cancer_arm_segs[[arm]] <- cancer_seg
          }
        }
        cancer_arm_segs <- do.call(rbind, cancer_arm_segs)
        patient_peak_df[[samp]] <- cancer_arm_segs 
        
      }
      patient_peak_df <- do.call(rbind, patient_peak_df)
      
      # how often do we see copy number amplifications in the genes which are left now
      patient_peak_genes <- patient_peak_df[, c("patient", "cn_chr", "cn_start", "cn_end", "cnTotal", "nMajor", "nMinor", "chr_arm", "gene", "is_cancer_gene", "samples_per_bin", "sample_per_bin_prop")]
      patient_peak_genes <- unique(patient_peak_genes)
      
      gene_count <- unique(patient_peak_genes[, c("patient", "gene")])
      gene_count <- data.frame(table(gene_count$gene))
      
      # lot's of these genes actually only have ampliffied copy number in one or two patients, so let's kick them out again based on our threshold
      genes_to_keep <- gene_count[gene_count$Freq > threshold[threshold$histology == hist, "amplification_threshold"], "Var1"]
      
      if(length(genes_to_keep) > 0){
        patient_peak_genes <- patient_peak_genes[patient_peak_genes$gene %in% genes_to_keep, ]
        other_genes_to_keep <- unique(patient_peak_genes[, "gene"])
        # for some chromosome arms there is still rubbish in here so let's investigate those a bit further
        chr_arms_false <- unique(patient_peak_genes[patient_peak_genes$is_cancer_gene == FALSE, "chr_arm"])
        
        if(length(chr_arms_false) > 0){
          other_genes_to_keep <- unique(patient_peak_genes[-which(patient_peak_genes$chr_arm %in% chr_arms_false), "gene"])
          new_genes_to_keep <- list()
          for(this_arm in chr_arms_false){
            
            false_chr_arm <- patient_peak_genes[patient_peak_genes$chr_arm == this_arm, ]
            
            if(length(unique(false_chr_arm$gene)) > 1){
              # how many tumours have amplifications in these genes?
              gene_count <- data.frame(table(false_chr_arm$gene))
              
              # are any of them cancer genes?
              gene_count<- left_join(gene_count, unique(patient_peak_genes[, c("gene", "is_cancer_gene")]), by = c("Var1" = "gene"))
              
              if(any(gene_count$is_cancer_gene)){
                new_genes_to_keep[[this_arm]] <- gene_count[gene_count$is_cancer_gene, "Var1"]
              } else { 
                count_median <- median(gene_count$Freq)
                new_genes_to_keep[[this_arm]] <- gene_count[gene_count$Freq > count_median, "Var1"]
              }
            } else {new_genes_to_keep[[this_arm]] <- unique(false_chr_arm$gene)}
            
            
          }
          new_genes_to_keep <- unlist(new_genes_to_keep)
          all_genes_to_keep <- c(new_genes_to_keep, other_genes_to_keep)
        } else {all_genes_to_keep <- other_genes_to_keep}
        
        patient_peak_genes_final <- patient_peak_genes[patient_peak_genes$gene %in% all_genes_to_keep, ]
        
        if(nrow(patient_peak_genes_final) > 0){
          
          # last thing I need to do is filter out some genes
          patient_peak_genes_final_filter <- patient_peak_genes_final %>% filter(!gene %in% olfactory_genes)
          
          hist_tcga_expression_genes <- tcga_expression_gene[, c("gene_name", hist)]
          colnames(hist_tcga_expression_genes) <- c("gene_name", "expression")
          hist_tcga_expression_genes <- unique(hist_tcga_expression_genes[hist_tcga_expression_genes$expression == 0, ])
          patient_peak_genes_final_filter <- patient_peak_genes_final_filter %>% filter(!gene %in% hist_tcga_expression_genes)
          
          hist_gtex_expression_genes <- gtex_expression_gene[, c("gene_name", hist)]
          colnames(hist_gtex_expression_genes) <- c("gene_name", "expression")
          hist_gtex_expression_genes <- unique(hist_gtex_expression_genes[hist_gtex_expression_genes$expression == 0, "gene_name"])
          patient_peak_genes_final_filter <- patient_peak_genes_final_filter %>% filter(!gene %in% hist_gtex_expression_genes)
          
          patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep("LOC", patient_peak_genes_final_filter$gene, invert = T), ]
          patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep(".*AS1", patient_peak_genes_final_filter$gene, invert = T), ]
          patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep("NEAT1", patient_peak_genes_final_filter$gene, invert = T), ]
          patient_peak_genes_final_filter <- patient_peak_genes_final_filter[grep("MALAT", patient_peak_genes_final_filter$gene, invert = T), ]
          
          patient_peak_genes_final_filter$histology <- hist
          loss_peak_genes_patient[[hist]] <- patient_peak_genes_final_filter
        }
      }
    }
  }
}

loss_peak_genes_patient <- do.call(rbind, loss_peak_genes_patient)

# create a file with all the gain and loss drivers

amp_driver_genes <- unique(amplification_peak_genes_patient[, c("gene", "is_cancer_gene", "histology", "samples_per_bin", "sample_per_bin_prop")])
amp_driver_genes$type <- "amplification"

loss_driver_genes <- unique(loss_peak_genes_patient[, c("gene", "is_cancer_gene", "histology", "samples_per_bin", "sample_per_bin_prop")])
loss_driver_genes$type <- "homozygous_deletion"

copy_number_driver <- rbind(amp_driver_genes, loss_driver_genes)
write.table(copy_number_driver, paste0(driver_output_path, "copy_number_driver_genes.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

# are there genes specific to any of the smaller subtypes?
amp_small_subtype_genes <- unique(amp_driver_genes[-which(amp_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")), "gene"])
amp_small_subtype <- unique(amp_driver_genes[-which(amp_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")), ])
amp_small_subtype_genes <- amp_small_subtype_genes[-which(amp_small_subtype_genes %in% unique(amp_driver_genes[amp_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), "gene"]))]
amp_small_subtype <- amp_small_subtype[-which(amp_small_subtype$gene %in% unique(amp_driver_genes[amp_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), "gene"])), ]

# which genes are specific to luad or lusc?
amp_luad_genes <- unique(amp_driver_genes[amp_driver_genes$histology == "ADENOCARCINOMA", "gene"])
amp_notluad_genes <- unique(amp_driver_genes[amp_driver_genes$histology != "ADENOCARCINOMA" & amp_driver_genes$histology != "PANCAN", "gene"])
amp_luad_specific_genes <- amp_luad_genes[-which(amp_luad_genes %in% amp_notluad_genes)]
luad_specific_amp_driver_genes <- amp_driver_genes[amp_driver_genes$gene %in% amp_luad_specific_genes & amp_driver_genes$histology == "ADENOCARCINOMA", ]

amp_lusc_genes <- unique(amp_driver_genes[amp_driver_genes$histology == "SQUAMOUS_CELL", "gene"])
amp_notlusc_genes <- unique(amp_driver_genes[amp_driver_genes$histology != "SQUAMOUS_CELL" & amp_driver_genes$histology != "PANCAN", "gene"])
amp_lusc_specific_genes <- amp_lusc_genes[-which(amp_lusc_genes %in% amp_notlusc_genes)]
lusc_specific_amp_driver_genes <- amp_driver_genes[amp_driver_genes$gene %in% amp_lusc_specific_genes & amp_driver_genes$histology == "SQUAMOUS_CELL", ]

# check losses as well

# are there genes specific to any of the smaller subtypes?
loss_small_subtype_genes <- unique(loss_driver_genes[-which(loss_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")), "gene"])
loss_small_subtype <- unique(loss_driver_genes[-which(loss_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")), ])
loss_small_subtype_genes <- loss_small_subtype_genes[-which(loss_small_subtype_genes %in% unique(loss_driver_genes[loss_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), "gene"]))]
loss_small_subtype <- loss_small_subtype[-which(loss_small_subtype$gene %in% unique(loss_driver_genes[loss_driver_genes$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), "gene"])), ]

# which genes are specific to luad or lusc?
loss_luad_genes <- unique(loss_driver_genes[loss_driver_genes$histology == "ADENOCARCINOMA", "gene"])
loss_notluad_genes <- unique(loss_driver_genes[loss_driver_genes$histology != "ADENOCARCINOMA" & loss_driver_genes$histology != "PANCAN", "gene"])
loss_luad_specific_genes <- loss_luad_genes[-which(loss_luad_genes %in% loss_notluad_genes)]
luad_specific_amp_driver_genes <- loss_driver_genes[loss_driver_genes$gene %in% loss_luad_specific_genes & loss_driver_genes$histology == "ADENOCARCINOMA", ]

loss_lusc_genes <- unique(loss_driver_genes[loss_driver_genes$histology == "SQUAMOUS_CELL", "gene"])
loss_notlusc_genes <- unique(loss_driver_genes[loss_driver_genes$histology != "SQUAMOUS_CELL" & loss_driver_genes$histology != "PANCAN", "gene"])
loss_lusc_specific_genes <- loss_lusc_genes[-which(loss_lusc_genes %in% loss_notlusc_genes)]
lusc_specific_loss_driver_genes <- loss_driver_genes[loss_driver_genes$gene %in% loss_lusc_specific_genes & loss_driver_genes$histology == "SQUAMOUS_CELL", ]

# make a list also including patient
amp_driver_patient <- amplification_peak_genes_patient
amp_driver_patient$type <- "amplification"

loss_driver_patient <- loss_peak_genes_patient
loss_driver_patient$type <- "homozygous_deletion"

driver_patient <- rbind(amp_driver_patient, loss_driver_patient)
write.table(driver_patient, paste0(driver_output_path, "copy_number_driver_genes_per_patient.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

# count how many driver genes there are per patient
driver_count <- unique(driver_patient[, c("patient", "cn_chr", "cn_start", "cn_end", "type")])
driver_count <- data.frame(table(driver_count[,c("patient" , "type")]))
driver_count <- dcast(driver_count, patient ~ type)
driver_count$total_driver_count <- driver_count$amplification + driver_count$homozygous_deletion

# add all missing samples to this as well
missing_samples <- sample_table$participant_id[-which(sample_table$participant_id %in% driver_count$patient)]

missing_df <- data.frame(patient = missing_samples, amplification = 0, homozygous_deletion = 0, total_driver_count = 0)

driver_count <- rbind(driver_count, missing_df)
write.table(driver_count, paste0(driver_output_path, "copy_number_driver_genes_per_patient_count.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

# make a file with drivers and p values
loss_driver_genes_p <- left_join(loss_driver_genes, cancer_loss_df)
amp_driver_genes_p <- left_join(amp_driver_genes, cancer_gain_df)

write.table(amp_driver_genes_p, paste0(output_path, "all_hist_amplification_driver_genes_with_p_values.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(loss_driver_genes_p, paste0(output_path, "all_hist_homdel_driver_genes_with_p_values.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# get the location midpoints of the genes

# load gene data
txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes                 <- genes(txdb, single.strand.genes.only=FALSE)
genes                 <- data.frame(genes)

# extract HUGO symbols
gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)

# merge with gene locations
gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))

# gene_names <- gene_names[gene_names$symbol %in% driver$gene, ]
gene_names <- gene_names[!is.na(gene_names$symbol), ]

amp_driver_genes <- left_join(amp_driver_genes, gene_names, by = c("gene" = "symbol"))
amp_driver_genes$gene_midpoint <- (amp_driver_genes$end - (amp_driver_genes$width/2))
amp_driver_genes <- amp_driver_genes[grep("_alt",amp_driver_genes$seqnames, invert = TRUE), ]

loss_driver_genes <- left_join(loss_driver_genes, gene_names, by = c("gene" = "symbol"))
loss_driver_genes$gene_midpoint <- (loss_driver_genes$end - (loss_driver_genes$width/2))
loss_driver_genes <- loss_driver_genes[grep("alt", loss_driver_genes$seqnames, invert = T), ]

# make continous

amp_driver_genes <- anno_midpoint_continous_fun(amp_driver_genes,
                                                genome = BSgenome.Hsapiens.UCSC.hg38)

loss_driver_genes <- anno_midpoint_continous_fun(loss_driver_genes,
                                                 genome = BSgenome.Hsapiens.UCSC.hg38)

chr_lines <- amp_plot_chr_lines[[1]]
chr_midpoint <- amp_plot_chr_midpoint[[1]]

write.table(amp_driver_genes, paste0(output_path, "all_hist_amplification_driver_genes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(loss_driver_genes, paste0(output_path, "all_hist_homdel_driver_genes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")


# we want LUAD and LUSC to be on the same axis
luad_max <- max(max(loss_bin_sum_cont[loss_bin_sum_cont$histology == "ADENOCARCINOMA", "sample_per_bin_prop"] * 100), max(amp_bin_sum_cont[amp_bin_sum_cont$histology == "ADENOCARCINOMA", "sample_per_bin_prop"] * 100))
lusc_max <- max(max(loss_bin_sum_cont[loss_bin_sum_cont$histology == "SQUAMOUS_CELL", "sample_per_bin_prop"] * 100), max(amp_bin_sum_cont[amp_bin_sum_cont$histology == "SQUAMOUS_CELL", "sample_per_bin_prop"] * 100))
ll_max <- max(luad_max, lusc_max)

# hist_peaks <- list()
for(hist in all_histologies){
  
  print(hist)
  
  gain_peaks <- amp_driver_genes[amp_driver_genes$histology == hist,]
  gain_peaks$is_cancer_gene <- factor(gain_peaks$is_cancer_gene, levels = c(TRUE, FALSE))
  
  if(hist %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")){
    max <- ll_max
  } else {
    max <- max(max(loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, "sample_per_bin_prop"] * 100), max(amp_bin_sum_cont[amp_bin_sum_cont$histology == hist, "sample_per_bin_prop"] * 100))
  }
  
  max <- max + (0.05*max)
  
  p_gain <- ggplot() +
    geom_line(data = amp_bin_sum_cont[amp_bin_sum_cont$histology == hist, ], aes(x = adjust_end_pos, y = sample_per_bin_prop*100, col = histology)) +
    geom_vline(xintercept = chr_lines[1:22], linetype = 2, color = "#bababa") +
    # geom_text(data = og, aes(adjust_mid, y = max_y_gain, label = gene_name), angle = 45) +
    # geom_label(data = gain_peaks, aes(adjust_mid, y = sample_per_bin_prop, label = gene, fill = is_cancer_gene)) +
    ggrepel::geom_text_repel(data = gain_peaks, aes(x = adjust_mid, y = sample_per_bin_prop*100, label = gene, color = is_cancer_gene), direction = "y", hjust = -1, vjust = 0, max.overlaps = Inf, fontface = "italic", size = 2.5) +
    scale_color_manual(values = c("ADENOCARCINOMA" = "#67001f",   
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
                                  "MET_OTHER" = "#7a797980",
                                  "PANCAN" = "black",
                                  "TRUE" = "black",
                                  "FALSE" = "#ef3b2c"),
                       labels = c(hist,
                                  "novel cancer gene",
                                  "COSMIC cancer gene"),
                       name = "") +
    # scale_color_manual(name = "", breaks = c(TRUE, FALSE), labels = c("COSMIC cancer gene", "novel cancer gene"), values = c("white", "#ef3b2c")) +
    scale_x_continuous(expand = c(0,max(loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, "sample_per_bin_prop"]) * 100)) +
    ylim(0, max) +
    # annotate("text", x = chr_midpoint, y = -0.005, label = c(paste0("chr", c(1:22, "X", "Y"))), size = 4.5) +
    ylab("% tumours\nwith amplififcation") +
    theme_classic() +
    theme(text = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.1,0,0,0), "cm"))
  
  loss_peaks <- loss_driver_genes[loss_driver_genes$hist == hist, ]
  loss_peaks$is_cancer_gene <- factor(loss_peaks$is_cancer_gene, levels = c(TRUE, FALSE))
  
  p_loss <- ggplot() +
    geom_line(data = loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, ], aes(adjust_end_pos, sample_per_bin_prop*100, col = histology)) +
    # geom_line(data = loh_bin_sum_cont[loh_bin_sum_cont$histology == hist, ], aes(adjust_end_pos, sample_per_bin_prop*100, col = histology), linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("ADENOCARCINOMA" = "#67001f",   
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
                                  "MET_OTHER" = "#7a797980",
                                  "PANCAN" = "black",
                                  "TRUE" = "black",
                                  "FALSE" = "#ef3b2c"),
                       labels = c(" ",
                                  "novel cancer gene",
                                  "COSMIC cancer gene"),
                       name = "") +
    # ylim(0, 25) +
    geom_vline(xintercept = chr_lines[1:22], linetype = 2, color = "#bababa") +
    ggrepel::geom_text_repel(data = loss_peaks, aes(x = adjust_mid, y = sample_per_bin_prop*100, label = gene, color = is_cancer_gene), direction = "y", hjust = 2, vjust = 0, max.overlaps = Inf, fontface = "italic", size = 2.5) +
    scale_x_continuous(expand = c(0,0), breaks = chr_midpoint[1:22], labels = sub("chr", "", names(chr_midpoint[1:22]))) +
    ylab("% tumours\nwith homozygous deletion\\LOH") +
    scale_y_reverse(limits = c(max, 0)) +
    theme_classic() +
    theme(text = element_text(size = 8),
          axis.title.x = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0,0,0.1,0), "cm"))
  
  p_combined <- plot_grid(p_gain, p_loss, ncol = 1, align = "v", axis = "rlbt", rel_heights = c(5, 5))
  
  pdf(paste0(plot_output_path, hist, "_genome_plot_FIX_CN.pdf"), width = 6, height = 3)
  print(p_combined)
  dev.off()
  
}

################################################################################
# get the peaks

# run peak detection

# amplifications
amp_all_hist <- list()
for(hist in unique(amp_bin_sum_cont$histology)){
  print(hist)
  
  data <- amp_bin_sum_cont[amp_bin_sum_cont$histology == hist, ]
  
  min_samps <- mean(data$sample_per_bin_prop) * 2
  
  # identify the consecutive bins that have more than the min_samps value
  data_high <- data[data$sample_per_bin_prop >= min_samps, ]
  
  if(nrow(data_high) > 1){
    
    # which ones are the consecutive peaks
    
    data_high[1, "peak"] <- paste0(hist, "_amplification_peak_", 1)
    for(i in 2:nrow(data_high)){
      
      previous_end <- data_high[i-1, "endpos"]
      
      # is the current start cust 1+ previous end?
      current_start <- data_high[i, "startpos"]
      
      if(current_start == previous_end+1){
        data_high[i, "peak"] <- data_high[i-1, "peak"]
      } else {data_high[i, "peak"] <- paste0(hist, "_amplification_peak_", i)}
      
    }
    
    data <- left_join(data, data_high)
    amp_all_hist[[hist]] <- data
    
  }
}

amp_all_hist <- do.call(rbind, amp_all_hist)

write.table(amp_all_hist, paste0(output_path, "amplification_peaks.txt"), row.names = F, col.names = T, sep = "\t", quote = F)

# losses
loss_all_hist <- list()
for(hist in unique(loss_bin_sum_cont$histology)){
  print(hist)
  
  data <- loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, ]
  
  min_samps <- mean(data$sample_per_bin_prop) * 3
  
  # identify the consecutive bins that have more than the min_samps value
  data_high <- data[data$sample_per_bin_prop >= min_samps, ]
  
  if(nrow(data_high) > 1){
    
    # which ones are the consecutive peaks
    
    data_high[1, "peak"] <- paste0(hist, "_loss_peak_", 1)
    for(i in 2:nrow(data_high)){
      
      previous_end <- data_high[i-1, "endpos"]
      
      # is the current start cust 1+ previous end?
      current_start <- data_high[i, "startpos"]
      
      if(current_start == previous_end+1){
        data_high[i, "peak"] <- data_high[i-1, "peak"]
      } else {data_high[i, "peak"] <- paste0(hist, "_loss_peak_", i)}
      
    }
    
    data <- left_join(data, data_high)
    loss_all_hist[[hist]] <- data
    
  }
}

loss_all_hist <- do.call(rbind, loss_all_hist)

write.table(loss_all_hist, paste0(output_path, "loss_peaks.txt"), row.names = F, col.names = T, sep = "\t", quote = F)


# now only the the stop and start coordinates of those peaks with genes from our list of copy number genes

# get locations of genes
txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes                 <- genes(txdb, single.strand.genes.only=FALSE)
genes                 <- data.frame(genes)

# extract HUGO symbols
gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)

# merge with gene locations
gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))
gene_names            <- gene_names[!is.na(gene_names$symbol), ]

CN_driver_genes <- read.table(paste0(driver_output_path, "copy_number_driver_genes.txt"), head = T, sep = "\t")

all_hist_copynumber_peaks <- list()
for(hist in c("ADENOCARCINOMA", "SQUAMOUS_CELL", "PANCAN")){
  
  print(hist)
  # amplifications first
  hist_amp_drivers <- CN_driver_genes[CN_driver_genes$histology == hist & CN_driver_genes$type == "amplification", ]
  
  # what are the genomic coordinates of these genes
  hist_amp_drivers <- left_join(hist_amp_drivers, gene_names[, c("seqnames", "start", "end", "symbol")], by = c("gene" = "symbol"))
  
  # overlap these coordinates with the peak coordinates
  hist_amp_drivers_GR <- GRanges(seqnames = hist_amp_drivers$seqnames, IRanges(start = hist_amp_drivers$start, end = hist_amp_drivers$end))
  mcols(hist_amp_drivers_GR) <- data.frame(hist_amp_drivers[, c("gene", "is_cancer_gene", "histology", "samples_per_bin", "sample_per_bin_prop", "type")])
  
  amp_hist <- amp_all_hist[amp_all_hist$histology == hist, ]
  amp_peak_gr <- GRanges(seqnames = amp_hist$chr, IRanges(start = amp_hist$startpos, end = amp_hist$endpos))
  mcols(amp_peak_gr) <- data.frame(amp_hist[, c("samples_per_bin", "sample_per_bin_prop", "bin", "adjust_start_pos", "adjust_end_pos", "histology", "peak")])
  
  overlap         <-  findOverlaps(hist_amp_drivers_GR, amp_peak_gr)
  intersect_df    <-  cbind(hist_amp_drivers[queryHits(overlap), ], amp_hist[subjectHits(overlap), ])
  intersect_df    <- intersect_df[!is.na(intersect_df$peak), ]
  
  peaks_with_genes <- unique(intersect_df$peak)
  
  # where do these peaks start and stop and which genes are on the peaks
  peaks_with_genes_coordinates <- amp_hist[amp_hist$peak %in% peaks_with_genes, ]
  
  # get start /stop coordinates and some important info for each peak
  amp_hist_peak_info <- list()
  for(this_peak in unique(peaks_with_genes_coordinates$peak)){
    
    df <- amp_hist[which(amp_hist$peak == this_peak), ]
    
    peak_start <- min(df$startpos)
    peak_end <- max(df$endpos)
    peak_median_samples_per_bin <- median(df$samples_per_bin)
    peak_median_sample_per_bin_prop <- median(df$sample_per_bin_prop)
    
    # which genes are on this peak
    significant_CN_genes_on_peak <- unique(intersect_df[intersect_df$peak == this_peak, "gene"])
    significant_CN_genes_on_peak <- paste(significant_CN_genes_on_peak, collapse = "_")
    
    peak_info <- data.frame(peak_name = this_peak,
                            peak_type = "amplification",
                            peak_histology = hist,
                            peak_chr = unique(df$chr), 
                            peak_start = peak_start, 
                            peak_end = peak_end, 
                            peak_median_samples_per_bin = peak_median_samples_per_bin,
                            peak_median_sample_per_bin_prop = peak_median_sample_per_bin_prop,
                            significant_CN_genes_on_peak = significant_CN_genes_on_peak)
    
    amp_hist_peak_info[[this_peak]] <- peak_info
    
  }
  amp_hist_peak_info <- do.call(rbind, amp_hist_peak_info)
  
  # # # # # # # # # now do loss peaks
  hist_loss_drivers <- CN_driver_genes[CN_driver_genes$histology == hist & CN_driver_genes$type == "homozygous_deletion", ]
  
  # what are the genomic coordinates of these genes
  hist_loss_drivers <- left_join(hist_loss_drivers, gene_names[, c("seqnames", "start", "end", "symbol")], by = c("gene" = "symbol"))
  
  # overlap these coordinates with the peak coordinates
  hist_loss_drivers_GR <- GRanges(seqnames = hist_loss_drivers$seqnames, IRanges(start = hist_loss_drivers$start, end = hist_loss_drivers$end))
  mcols(hist_loss_drivers_GR) <- data.frame(hist_loss_drivers[, c("gene", "is_cancer_gene", "histology", "samples_per_bin", "sample_per_bin_prop", "type")])
  
  loss_hist <- loss_all_hist[loss_all_hist$histology == hist, ]
  loss_peak_gr <- GRanges(seqnames = loss_hist$chr, IRanges(start = loss_hist$startpos, end = loss_hist$endpos))
  mcols(loss_peak_gr) <- data.frame(loss_hist[, c("samples_per_bin", "sample_per_bin_prop", "bin", "adjust_start_pos", "adjust_end_pos", "histology", "peak")])
  
  overlap         <-  findOverlaps(hist_loss_drivers_GR, loss_peak_gr)
  intersect_df    <-  cbind(hist_loss_drivers[queryHits(overlap), ], loss_hist[subjectHits(overlap), ])
  intersect_df    <- intersect_df[!is.na(intersect_df$peak), ]
  
  peaks_with_genes <- unique(intersect_df$peak)
  
  # where do these peaks start and stop and which genes are on the peaks
  peaks_with_genes_coordinates <- loss_hist[loss_hist$peak %in% peaks_with_genes, ]
  
  # get start /stop coordinates and some important info for each peak
  loss_hist_peak_info <- list()
  for(this_peak in unique(peaks_with_genes_coordinates$peak)){
    
    df <- loss_hist[which(loss_hist$peak == this_peak), ]
    
    peak_start <- min(df$startpos)
    peak_end <- max(df$endpos)
    peak_median_samples_per_bin <- median(df$samples_per_bin)
    peak_median_sample_per_bin_prop <- median(df$sample_per_bin_prop)
    
    # which genes are on this peak
    significant_CN_genes_on_peak <- unique(intersect_df[intersect_df$peak == this_peak, "gene"])
    significant_CN_genes_on_peak <- paste(significant_CN_genes_on_peak, collapse = "_")
    
    peak_info <- data.frame(peak_name = this_peak,
                            peak_type = "homozygous_deletion",
                            peak_histology = hist,
                            peak_chr = unique(df$chr), 
                            peak_start = peak_start, 
                            peak_end = peak_end, 
                            peak_median_samples_per_bin = peak_median_samples_per_bin,
                            peak_median_sample_per_bin_prop = peak_median_sample_per_bin_prop,
                            significant_CN_genes_on_peak = significant_CN_genes_on_peak)
    
    loss_hist_peak_info[[this_peak]] <- peak_info
    
  }
  loss_hist_peak_info <- do.call(rbind, loss_hist_peak_info)
  
  hist_peak_info <- rbind(amp_hist_peak_info, loss_hist_peak_info)
  all_hist_copynumber_peaks[[hist]] <- hist_peak_info
}

all_hist_copynumber_peaks <- do.call(rbind, all_hist_copynumber_peaks)

# # # # # #  which samples have amplified

peak_copynumber_events <- list()
for(hist in unique(all_hist_copynumber_peaks$peak_histology)){
  
  print(hist)
  
  # which samples have amplifications in these peaks
  
  hist_amplification_peaks <- all_hist_copynumber_peaks[all_hist_copynumber_peaks$peak_histology == hist & all_hist_copynumber_peaks$peak_type == "amplification", ]
  
  if(hist %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")){
    hist_amps <- amps[amps$histology == hist, ]
  } else {hist_amps <- amps }
  
  hist_amps_GR <- GRanges(seqnames = hist_amps$chr, IRanges(start = hist_amps$startpos, end = hist_amps$endpos))
  mcols(hist_amps_GR) <- data.frame(hist_amps)
  
  hist_amplification_peaks_GR <- GRanges(seqnames = hist_amplification_peaks$peak_chr, IRanges(start = hist_amplification_peaks$peak_start, end = hist_amplification_peaks$peak_end))
  mcols(hist_amplification_peaks_GR) <- data.frame(hist_amplification_peaks)
  
  overlap         <-  findOverlaps(hist_amps_GR, hist_amplification_peaks_GR)
  intersect_gr    <-  pintersect(hist_amplification_peaks_GR[subjectHits(overlap)], hist_amps_GR[queryHits(overlap)])
  mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(hist_amps_GR[queryHits(overlap)]))
  amp_intersect_df    <-  data.frame(intersect_gr)
  amp_intersect_df$strand <- NULL
  amp_intersect_df$hit <- NULL
  colnames(amp_intersect_df) <- c("overlap_chr", "overlap_start", "overlap_end", "overap_width", "peak_name", "peak_type", "peak_histology", "peak_chr", "peak_start", "peak_end",
                                  "peak_median_samples_per_bin", "peak_median_sample_per_bin_prop", "significant_CN_genes_on_peak", "patient", "seg_chr", "seg_startpos", "seg_endpos",
                                  "cnTotal", "nMajor", "nMinor", "Ploidy", "ACF", "amplification", "homdel", "patient_histology")
  
  # do the same for losses
  hist_loss_peaks <- all_hist_copynumber_peaks[all_hist_copynumber_peaks$peak_histology == hist & all_hist_copynumber_peaks$peak_type == "homozygous_deletion", ]
  
  if(hist %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")){
    hist_loss <- loss[loss$histology == hist, ]
  } else {hist_loss <- loss }
  
  hist_loss_GR <- GRanges(seqnames = hist_loss$chr, IRanges(start = hist_loss$startpos, end = hist_loss$endpos))
  mcols(hist_loss_GR) <- data.frame(hist_loss)
  
  hist_loss_peaks_GR <- GRanges(seqnames = hist_loss_peaks$peak_chr, IRanges(start = hist_loss_peaks$peak_start, end = hist_loss_peaks$peak_end))
  mcols(hist_loss_peaks_GR) <- data.frame(hist_loss_peaks)
  
  overlap         <-  findOverlaps(hist_loss_GR, hist_loss_peaks_GR)
  intersect_gr    <-  pintersect(hist_loss_peaks_GR[subjectHits(overlap)], hist_loss_GR[queryHits(overlap)])
  mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(hist_loss_GR[queryHits(overlap)]))
  loss_intersect_df    <-  data.frame(intersect_gr)
  loss_intersect_df$strand <- NULL
  loss_intersect_df$hit <- NULL
  colnames(loss_intersect_df) <- c("overlap_chr", "overlap_start", "overlap_end", "overap_width", "peak_name", "peak_type", "peak_histology", "peak_chr", "peak_start", "peak_end",
                                   "peak_median_samples_per_bin", "peak_median_sample_per_bin_prop", "significant_CN_genes_on_peak", "patient", "seg_chr", "seg_startpos", "seg_endpos",
                                   "cnTotal", "nMajor", "nMinor", "Ploidy", "ACF", "amplification", "homdel", "patient_histology")
  
  peak_samples <- rbind(amp_intersect_df, loss_intersect_df)
  
  peak_copynumber_events[[hist]] <- peak_samples
  
}
peak_copynumber_events <- do.call(rbind, peak_copynumber_events)

# reformat this for analysis purposes
patient_copynumber_events <- peak_copynumber_events

# how much overlap does the copy number segment have with the peak?
patient_copynumber_events$peak_length <- patient_copynumber_events$peak_end - patient_copynumber_events$peak_start
patient_copynumber_events$peak_overlap_proportion <- patient_copynumber_events$overap_width / patient_copynumber_events$peak_length

# let's say the copy number must overlap the peak by 10%
patient_copynumber_events[patient_copynumber_events$peak_overlap_proportion >= 0.05, "event"] <- TRUE
patient_copynumber_events[patient_copynumber_events$peak_overlap_proportion < 0.05, "event"] <- FALSE

# I need to rename the peaks to reflect where they are in the genome
# get cytoband locations
cytobands <- read.table(cytoband_path, head = T)
cytobands$stain <- NULL

# overlap the peaks with the cytobands
all_hist_copynumber_peaks_GR <- GRanges(seqnames = all_hist_copynumber_peaks$peak_chr, IRanges(start = all_hist_copynumber_peaks$peak_start, end = all_hist_copynumber_peaks$peak_end))
cytobands_GR <- GRanges(seqnames = cytobands$chrom, IRanges(start = cytobands$start, end = cytobands$end))

overlap         <-  findOverlaps(all_hist_copynumber_peaks_GR, cytobands_GR)
intersect_df    <-  cbind(all_hist_copynumber_peaks[queryHits(overlap), ], cytobands[subjectHits(overlap), ])

# for each peak get the start and end cytoband if it covers more than one cytobands
peaks_cytoband_names <- list()
for(peak in unique(intersect_df$peak_name)){
  
  df <- intersect_df[intersect_df$peak_name == peak, ]
  
  if(nrow(df) == 1){
    new_peak_name <- df$band
  }
  
  if(nrow(df) > 1){
    min_peak <- min(df$band)
    max_peak <- max(df$band)
    new_peak_name <- paste0(min_peak, "_", max_peak)
  }
  
  hist <- unique(df$peak_histology)
  type <- unique(df$peak_type)
  chr  <- unique(df$peak_chr)
  
  df$new_peak_name <- paste0(hist, "_", type, "_", chr, "_", new_peak_name)
  peaks_cytoband_names[[peak]] <- df
  
}

peaks_cytoband_names <- do.call(rbind, peaks_cytoband_names)

patient_copynumber_events <- left_join(patient_copynumber_events, unique(peaks_cytoband_names[, c("peak_name", "new_peak_name")]))


# now just make a list with patient and peak name and if they have the event
patient_copynumber_events_final <- patient_copynumber_events[, c("patient", "new_peak_name", "event", "significant_CN_genes_on_peak")]
patient_copynumber_events_final <- patient_copynumber_events_final[patient_copynumber_events_final$event == TRUE, ]
patient_copynumber_events_final <- unique(patient_copynumber_events_final)

write.table(patient_copynumber_events_final, paste0(output_path, "copynumber_peak_driver_per_patient.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# plot the peaks across the genome
#  load peak data
load(amp_bin_path)
amp_bin_sum_cont <- do.call(rbind, amp_plot_data)
amp_bin_sum_cont$histology <- rownames(amp_bin_sum_cont)
amp_bin_sum_cont$histology <- sub("\\..*", "", amp_bin_sum_cont$histology)

load(loss_bin_path)
loss_bin_sum_cont <- do.call(rbind, loss_plot_data)
loss_bin_sum_cont$histology <- rownames(loss_bin_sum_cont)
loss_bin_sum_cont$histology <- sub("\\..*", "", loss_bin_sum_cont$histology)

amp_driver <- amp_driver_genes
loss_driver <- loss_driver_genes

# gotta modify the peak data for plotting purposes
# losses 
loss_peaks_cytoband_names <- peaks_cytoband_names[peaks_cytoband_names$peak_type == "homozygous_deletion", ]
loss_peaks_cytoband_names$chr_hist <- paste0(loss_peaks_cytoband_names$peak_chr, ":", loss_peaks_cytoband_names$peak_histology)
loss_peaks_cytoband_names <- unique(loss_peaks_cytoband_names[, c("chr_hist", "new_peak_name", "peak_start", "peak_end")])
loss_all_hist$chr_hist <- paste0(loss_all_hist$chr, ":", loss_all_hist$histology)

loss_all_hist_GR <- GRanges(seqnames = loss_all_hist$chr_hist, IRanges(start = loss_all_hist$startpos, end = loss_all_hist$endpos))
peaks_cytoband_names_GR <- GRanges(seqnames = loss_peaks_cytoband_names$chr_hist, IRanges(start = loss_peaks_cytoband_names$peak_start, end = loss_peaks_cytoband_names$peak_end))

overlap         <-  findOverlaps(loss_all_hist_GR, peaks_cytoband_names_GR)
intersect_df    <-  cbind(loss_all_hist[queryHits(overlap), ], loss_peaks_cytoband_names[subjectHits(overlap), ])

loss_all_hist <- left_join(loss_all_hist, unique(intersect_df[, c("bin", "chr_hist", "new_peak_name")]))

# amplifications
amp_peaks_cytoband_names <- peaks_cytoband_names[peaks_cytoband_names$peak_type == "amplification", ]
amp_peaks_cytoband_names$chr_hist <- paste0(amp_peaks_cytoband_names$peak_chr, ":", amp_peaks_cytoband_names$peak_histology)
amp_peaks_cytoband_names <- unique(amp_peaks_cytoband_names[, c("chr_hist", "new_peak_name", "peak_start", "peak_end")])
amp_all_hist$chr_hist <- paste0(amp_all_hist$chr, ":", amp_all_hist$histology)

amp_all_hist_GR <- GRanges(seqnames = amp_all_hist$chr_hist, IRanges(start = amp_all_hist$startpos, end = amp_all_hist$endpos))
peaks_cytoband_names_GR <- GRanges(seqnames = amp_peaks_cytoband_names$chr_hist, IRanges(start = amp_peaks_cytoband_names$peak_start, end = amp_peaks_cytoband_names$peak_end))

overlap         <-  findOverlaps(amp_all_hist_GR, peaks_cytoband_names_GR)
intersect_df    <-  cbind(amp_all_hist[queryHits(overlap), ], amp_peaks_cytoband_names[subjectHits(overlap), ])

amp_all_hist <- left_join(amp_all_hist, unique(intersect_df[, c("bin", "chr_hist", "new_peak_name")]))

# plot this to verify

chr_lines <- amp_plot_chr_lines[[1]]
chr_midpoint <- amp_plot_chr_midpoint[[1]]

for(hist in unique(amp_driver$histology)){
  
  print(hist)
  
  gain_peaks <- amp_driver[amp_driver$histology == hist,]
  gain_peaks$is_cancer_gene <- factor(gain_peaks$is_cancer_gene, levels = c(TRUE, FALSE))
  
  p_gain <- ggplot() +
    geom_line(data = amp_bin_sum_cont[amp_bin_sum_cont$histology == hist, ], aes(x = adjust_end_pos, y = sample_per_bin_prop, col = histology)) +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    # geom_text(data = og, aes(adjust_mid, y = max_y_gain, label = gene_name), angle = 45) +
    # geom_label(data = gain_peaks, aes(adjust_mid, y = sample_per_bin_prop, label = gene, fill = is_cancer_gene)) +
    ggrepel::geom_text_repel(data = gain_peaks, aes(x = adjust_mid, y = sample_per_bin_prop, label = gene, color = is_cancer_gene), direction = "y", hjust = -1, vjust = 0, max.overlaps = Inf, fontface = "italic", size = 2) +
    scale_color_manual(values = c("ADENOCARCINOMA" = "#67001f",   
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
                                  "MET_OTHER" = "#7a797980",
                                  "PANCAN" = "black",
                                  "TRUE" = "black",
                                  "FALSE" = "#ef3b2c"),
                       labels = c(hist,
                                  "COSMIC cancer gene", 
                                  "novel cancer gene"),
                       name = "") +
    # scale_color_manual(name = "", breaks = c(TRUE, FALSE), labels = c("COSMIC cancer gene", "novel cancer gene"), values = c("white", "#ef3b2c")) +
    scale_x_continuous(expand = c(0,0)) +
    # ylim(0, 1) +
    # annotate("text", x = chr_midpoint, y = -0.005, label = c(paste0("chr", c(1:22, "X", "Y"))), size = 4.5) +
    ylab("proportion of samples") +
    theme_classic() +
    theme(text = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.1,0,0,0), "cm"))
  
  amp_peaks <- amp_all_hist[amp_all_hist$histology == hist, ] 
  
  p_gain_peaks <- ggplot() +
    geom_tile(data = amp_peaks, aes(adjust_end_pos, 0)) +
    geom_tile(data = amp_peaks[!is.na(amp_peaks$new_peak_name), ], aes(adjust_end_pos, 1, fill = new_peak_name)) +
    # scale_fill_manual(values = c("black", "white")) +
    scale_x_continuous(expand = c(0,0)) +
    # scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    ylab("") +
    xlab("") +
    theme_classic() +
    theme(text = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none")
  
  
  # max_y_loss <- max(loss_bin_sum_cont$sample_per_bin_prop)
  loss_peaks <- loss_driver[loss_driver$hist == hist, ]
  loss_peaks$is_cancer_gene <- factor(loss_peaks$is_cancer_gene, levels = c(TRUE, FALSE))
  
  p_loss <- ggplot() +
    geom_line(data = loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, ], aes(adjust_end_pos, sample_per_bin_prop, col = histology)) +
    scale_color_manual(values = c("ADENOCARCINOMA" = "#67001f",   
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
                                  "MET_OTHER" = "#7a797980",
                                  "PANCAN" = "black",
                                  "TRUE" = "black",
                                  "FALSE" = "#ef3b2c"),
                       labels = c(" ",
                                  "COSMIC cancer gene", 
                                  "novel cancer gene"),
                       name = "") +
    # geom_text(data = tsg, aes(adjust_mid, y = max_y_loss, label = gene_name), angle = 45) +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    # geom_label(data = loss_peaks, aes(adjust_mid, y = sample_per_bin_prop, label = gene, fill = is_cancer_gene)) +
    ggrepel::geom_text_repel(data = loss_peaks, aes(x = adjust_mid, y = sample_per_bin_prop, label = gene, color = is_cancer_gene), direction = "y", hjust = 2, vjust = 0, max.overlaps = Inf, fontface = "italic", size = 2) +
    # scale_fill_manual(name = "", breaks = c(TRUE, FALSE), labels = c("COSMIC cancer gene", "novel cancer gene"), values = c("white", "#ef3b2c")) +
    scale_x_continuous(expand = c(0,0)) +
    # ylim(0, 1) +
    annotate("text", x = chr_midpoint, y = max(loss_bin_sum_cont[loss_bin_sum_cont$histology == hist, "sample_per_bin_prop"]) + 0.05, label = c(paste0("chr", c(1:22, "X", "Y"))), size = 2, angle = 45) +
    ylab("proportion of samples") +
    scale_y_reverse() +
    theme_classic() +
    theme(text = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0,0,0.1,0), "cm"))
  
  
  hist_homdel_peaks <- loss_all_hist[loss_all_hist$histology == hist, ] 
  
  p_loss_peaks <- ggplot() +
    geom_tile(data = hist_homdel_peaks, aes(adjust_end_pos, 0)) +
    geom_tile(data = hist_homdel_peaks[!is.na(hist_homdel_peaks$new_peak_name), ], aes(adjust_end_pos, 1, fill = new_peak_name)) +
    # scale_fill_manual(values = c("black", "white")) +
    scale_x_continuous(expand = c(0,0)) +
    # scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    ylab("") +
    xlab("") +
    theme_classic() +
    theme(text = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none")
  
  
  p_combined <- plot_grid(p_gain, p_gain_peaks, p_loss_peaks, p_loss, ncol = 1, align = "v", axis = "rlbt", rel_heights = c(5, 0.5, 0.5, 5))
  
  pdf(paste0(output_path, hist, "_genome_plot_with_peaks.pdf"), width = 24, height = 12)
  print(p_combined)
  dev.off()
  
}







