# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                 fishhook putput analysis                            # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
####################################################################     library

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(cowplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape2)
library(data.table)
library(tidyverse)
library(ggbeeswarm)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(org.Hs.eg.db)
library(poolr, lib.loc = "/re_gecip/cancer_lung/R_packages_4_1/")
library(pacman, lib.loc = "/re_gecip/cancer_lung/kthol/R_packages/R4/")
library(grid)
library(ggh4x)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
####################################################################   functions

source("/re_gecip/cancer_lung/CBailey1/sv_lung_gecip_project/fishhook/scripts/fishhook_functions.R")
source("/re_gecip/cancer_lung/CBailey1/pancan_ecdna_analysis/scripts/get_overlap_function.R")

################################################################################
####################################################################   file path

fh_output_dir                 <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/model_outputs/"
gene_def_path                 <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/data/regions_for_driver_search.txt"
SV_path                       <- "/re_gecip/cancer_lung/CBailey1/sv_lung_gecip_project/sv_analysis/scripts/prepare_sv_data/sv_data_tables/sv_data.RDS"
sample_table_path             <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/sampleTables/LungLandscape_sampleList_20221223.txt"
olfactory_gene_path           <- "/re_gecip/cancer_lung/shared/olfactory_barnes_2020.csv"
tcga_expression_gene_path     <- "/re_gecip/cancer_lung/shared/TCGA_expression_in_tumorsubtypes.csv"
gtex_expression_gene_path     <- "/re_gecip/cancer_lung/shared/GTEx_expression_in_tumorsubtypes.csv"

output_path                   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/"

################################################################################
####################################################################        MAIN
sample_table <- read.table(sample_table_path, head = T, sep = "\t")

all_subtypes <- c("ALL_SUBTYPES", "ADENOCARCINOMA", "SQUAMOUS_CELL", "LARGE_CELL", "SMALL_CELL", "MESOTHELIOMA", "ADENOSQUAMOUS", "MET_ADENOCARCINOMA")
# all_events <- c("ALL_EVENTS", "DUP", "DEL", "h2hINV", "t2tINV", "TRA")
# all_events <- c("DUP", "DEL", "h2hINV", "t2tINV")
all_events <- c("ALL_EVENTS", "DUP", "DEL", "h2hINV", "t2tINV", "TRA")

# load and format the SVs
sv_database <- readRDS(SV_path)
svdt_chr <- sv_database[ , chr := as.character(chr)][ chr == 23 , chr := "X"][chr==24, chr := "Y"]
svdt_chr <- svdt_chr[, Start_position := as.numeric(POS)][, End_position := as.numeric(POS)][, strand := '+']

svdt_chr <- svdt_chr %>% dplyr::rename(participant_id = sample_number)

# svdt_chr <- svdt_chr[svdt_chr$SV_CLASS != "TRA", ]

events_seqlengths <- dt2gr(svdt_chr)
seqlevelsStyle(events_seqlengths) <- "UCSC"

#load updated histology list
updated_histology_list <- fread(sample_table_path)
updated_histology_list$participant_id <- as.character(updated_histology_list$participant_id)

#update_sv_dt
svdt_chr <- inner_join(svdt_chr, updated_histology_list)
svdt_chr <- unique(svdt_chr)
svdt_chr$sv_key <- paste0(svdt_chr$participant_id, ":", svdt_chr$SV_CLASS, ":", svdt_chr$CHROM, ":", svdt_chr$POS)
svdt_chr$POS <- as.numeric(svdt_chr$POS)

sv_gr <- GRanges(seqnames = svdt_chr$CHROM, IRanges(start = svdt_chr$POS, end = svdt_chr$POS))
mcols(sv_gr) <- DataFrame(svdt_chr[, c("participant_id", "SV_CLASS", "ID", "histology", "sv_key", "SV_ID")])

# also pair the SVs up
sv_pair <- svdt_chr[, c("participant_id", "SV_CLASS", "CHROM", "POS", "SV_ID")]
sv_pair <- as.data.table(sv_pair)[, dcast(.SD, SV_ID ~ rowid(SV_ID), value.var = names(sv_pair)[-5])]
sv_pair <- data.frame(sv_pair)
sv_pair$participant_id_2 <- NULL
colnames(sv_pair)[2] <- "participant_id"
sv_pair$SV_CLASS_2 <- NULL
colnames(sv_pair)[3] <- "SV_CLASS" 
sv_pair <- left_join(sv_pair, updated_histology_list[, c("participant_id", "histology")])

# sv_pair <- sv_pair[sv_pair$SV_CLASS != "TRA", ]
sv_pair$sv_key <- paste0(sv_pair$participant_id, ":", sv_pair$SV_CLASS, ":", sv_pair$CHROM_1, ":", sv_pair$POS_1, ":", sv_pair$POS_2)

# what is the length of our SVs?
# SV_length <- data.frame(abs(sv_pair$POS_2 - sv_pair$POS_1))
# colnames(SV_length) <- "length"
# p <- ggplot(SV_length, aes(x=length)) + 
#   geom_histogram(bins = 1000) +
#   geom_vline(xintercept = 50000, color = "red") +
#   xlim(0, 1000000) +
#   xlab("SV length")

# make SVs into genomic ranges, but need to change translocations
non_TRA <- sv_pair[sv_pair$SV_CLASS != "TRA", ]
TRA <- sv_pair[sv_pair$SV_CLASS == "TRA", ]

TRA1 <- TRA[, c("SV_ID", "participant_id", "SV_CLASS", "CHROM_1", "POS_1", "histology", "sv_key")]
TRA1$CHROM_2 <- TRA1$CHROM_1
TRA1$POS_2 <- TRA1$POS_1
TRA2 <- TRA[, c("SV_ID", "participant_id", "SV_CLASS", "CHROM_2", "POS_2", "histology", "sv_key")]
TRA2$CHROM_1 <- TRA2$CHROM_2
TRA2$POS_1 <- TRA2$POS_2

TRA <- rbind(TRA1, TRA2)

sv_pair <- unique(rbind(non_TRA, TRA))

sv_pair_gr <- GRanges(seqnames = sv_pair$CHROM_1, IRanges(start = sv_pair$POS_1, end = sv_pair$POS_2))
mcols(sv_pair_gr) <- DataFrame(sv_pair[, c("participant_id", "SV_CLASS", "SV_ID", "histology", "sv_key")])

# get the gene definitions
gene_def <- read.table(gene_def_path, head = T, sep = "\t")
gene_def$seqnames <- paste0("chr", gene_def$seqnames)
gene_def$gene_key <- paste0(gene_def$seqnames, ":", gene_def$start, ":", gene_def$end, ":", gene_def$name)

# remove HLA genes
gene_def <- gene_def[grep("HLA", gene_def$name, invert = T), ]

# separate the genes from enhancers
gene_def_no_enhancer <- gene_def[gene_def$gr_id != "enhancer", ]
gene_def_enhancer <- gene_def[gene_def$gr_id == "enhancer", ]

enhancer_separate <- list()
for(this_name in unique(gene_def_enhancer$gene_id)){
  
  df <- gene_def_enhancer[gene_def_enhancer$gene_id == this_name, ]
  
  names <- str_split(unique(df$gene_id), "__")[[1]]
  names <- grep("and", names, invert = T, value = T)
  
  data <- df[, c("seqnames", "start", "end", "width", "strand", "target_genome_version", "gr_id", "name", "gene_key")]
  nrow_og_data <- nrow(data)
  
  # multiply this by how many names there are 
  data <- data[rep(seq_len(nrow(data)), each = length(names)), ]
  data$gene_id <- rep(names, each=nrow_og_data)
  data$gene_name <- rep(names, each=nrow_og_data)
  
  enhancer_separate[[this_name]] <- data
}
enhancer_separate <- do.call(rbind, enhancer_separate)

gene_def <- rbind(gene_def_no_enhancer, enhancer_separate)


# # read fishhook output
# all_subtype_list <- list()
# SVs_in_samples_subtype <- list()
# for(subtype in all_subtypes){
#   print(subtype)
# 
#   SVs_in_samples_event <- list()
#   all_event_list <- list()
#   for(event in all_events){
#     print(event)
# 
#     this_path <- paste0(fh_output_dir, "adjusted_model_bin_summary", event, "_", subtype, ".tsv")
#     fh_data <- read.table(this_path, head = T)
#     fh_data$bin_key <- paste0(fh_data$seqnames, ":", fh_data$start, ":", fh_data$end, ":", fh_data$bin_number)
# 
#     # overlap with gene element definitions
#     fh_gr <- GRanges(seqnames = fh_data$seqnames, IRanges(start = fh_data$start, end = fh_data$end))
#     mcols(fh_gr)  <- DataFrame(fh_data[, c("bin_number", "hid", "p", "fdr", "count", "effectsize", "count.pred", "count.density", "count.pred.density", "eligible", "p.neg", "fdr.neg", "log10_p_observed", "bin_key")])
# 
#     gene_def_gr <- GRanges(seqnames = gene_def$seqnames, IRanges(start = gene_def$start, end = gene_def$end))
#     mcols(gene_def_gr) <- DataFrame(gene_def[, c("name", "gene_key", "gene_name", "gr_id")])
# 
#     overlaps  <- findOverlaps(fh_gr, gene_def_gr)
#     intersect <- pintersect(fh_gr[queryHits(overlaps)], gene_def_gr[subjectHits(overlaps)])
#     mcols(intersect) <- data.frame(mcols(intersect), mcols(gene_def_gr[subjectHits(overlaps)]))
# 
#     intersect_df <- data.frame(intersect)
# 
#     bin_gene_df <- intersect_df[, c("bin_key", "p", "fdr", "count", "effectsize", "count.pred", "count.density", "count.pred.density", "eligible", "p.neg", "fdr.neg", "log10_p_observed", "name")]
#     bin_gene_df <- unique(bin_gene_df)
# 
#     # for each gene element combine the p values using fisher method
#     bin_gene_df <- bin_gene_df %>% group_by(name) %>% mutate(gene_element_fisher_p = poolr::fisher(p, adjust = "none")$p)
# 
#     gene_df <- bin_gene_df[, c("name", "gene_element_fisher_p")]
#     gene_df <- unique(gene_df)
#     gene_df$p_adjust <- p.adjust(gene_df$gene_element_fisher_p, method = "fdr")
# 
#     # take only those gene elements with a fdr p value < 0.05
#     gene_sig <- gene_df[gene_df$p_adjust < 0.05, ]
# 
#     if(nrow(gene_sig) > 0){
#       # now I need to get the number of samples with SVs in these gene elements
#       # if(subtype == "ALL_SUBTYPES"){
#          this_sv_gr <- sv_pair_gr
#       # } else { this_sv_gr <- sv_pair_gr[sv_pair_gr$histology == subtype, ] }
#       #
#       if(event == "ALL_EVENTS"){
#         this_sv_gr <- this_sv_gr
#       } else { this_sv_gr <- this_sv_gr[this_sv_gr$SV_CLASS == event, ] }
# 
#       # also get the unpaired SVs
#       # if(subtype == "ALL_SUBTYPES"){
#         this_sv_sng_gr <- sv_gr
#       # } else { this_sv_sng_gr <- sv_gr[sv_gr$histology == subtype, ] }
# 
#       if(event == "ALL_EVENTS"){
#         this_sv_sng_gr <- this_sv_sng_gr
#       } else { this_sv_sng_gr <- this_sv_sng_gr[this_sv_sng_gr$SV_CLASS == event, ] }
# 
#       # get those SVs which have start and stop on the same bin
#       overlap_SV_bin <- findOverlaps(this_sv_sng_gr, fh_gr)
#       intersect_SV_bin <- pintersect(this_sv_sng_gr[queryHits(overlap_SV_bin)], fh_gr[subjectHits(overlap_SV_bin)])
#       mcols(intersect_SV_bin) <- data.frame(mcols(intersect_SV_bin), mcols(fh_gr[subjectHits(overlap_SV_bin)]))
#       intersect_SV_bin_df <- data.frame(intersect_SV_bin)
# 
#       SV_ids_to_keep <- intersect_SV_bin_df[, c("SV_ID", "bin_key")]
#       SV_ids_to_keep <- unique(SV_ids_to_keep)
#       SV_ids_to_keep <- data.frame(table(SV_ids_to_keep$SV_ID))
#       SV_ids_to_keep <- SV_ids_to_keep[SV_ids_to_keep$Freq == 1, "Var1"]
#       SV_ids_to_keep <- as.character(SV_ids_to_keep)
# 
#       # get the translocations back if there are translocations involved
#       # if(any(grepl("TRA", SV_ids_to_keep))){
#         TRAs_to_keep <- intersect_SV_bin_df[, c("SV_ID", "bin_key")]
#         TRAs_to_keep <- unique(TRAs_to_keep[grep("TRA", TRAs_to_keep$SV_ID), ])
#         TRAs_to_keep <- unique(as.character(TRAs_to_keep$SV_ID))
# 
#         SV_ids_to_keep <- c(SV_ids_to_keep, TRAs_to_keep)
#       # }
# 
#       this_sv_gr <- this_sv_gr[this_sv_gr$SV_ID %in% SV_ids_to_keep,]
# 
#       # overlap the positions of the SVs with the gene elements
#       overlap_SV_gene <- findOverlaps(this_sv_gr, gene_def_gr)
#       intersect_SV_gene <- pintersect(this_sv_gr[queryHits(overlap_SV_gene)], gene_def_gr[subjectHits(overlap_SV_gene)])
#       mcols(intersect_SV_gene) <- data.frame(mcols(intersect_SV_gene), mcols(gene_def_gr[subjectHits(overlap_SV_gene)]))
# 
#       intersect_SV_gene_df <- data.frame(intersect_SV_gene)
# 
#       # save this info for later
#       SVs_in_samples_event[[event]] <- intersect_SV_gene_df
# 
#       if(subtype != "ALL_SUBTYPES"){
#         SV_gene_df <- intersect_SV_gene_df[intersect_SV_gene_df$histology == subtype, ]
#       } else { SV_gene_df <- intersect_SV_gene_df}
# 
#       SV_gene_df <-SV_gene_df[, c("participant_id", "gene_name", "gr_id")]
#       SV_gene_df <- unique(SV_gene_df)
# 
#       # mush all elements that are not enhancer into one
#       SV_gene_df[SV_gene_df$gr_id != "enhancer", "gr_id"] <- "gene"
#       SV_gene_df <- unique(SV_gene_df)
# 
#       # count how many tumours have SVs in these gene elements
#       SV_gene_count <- data.frame(table(SV_gene_df[, c("gene_name", "gr_id")]))
#       colnames(SV_gene_count) <- c("name", "gr_id", "number_tumours_with_SV_in_gene")
#       SV_gene_count <- SV_gene_count[SV_gene_count$number_tumours_with_SV_in_gene != 0, ]
# 
#       # add this to the significant gene list
#       full_gene_sig <- gene_sig
#       full_gene_sig$gene_name <- sapply(full_gene_sig$name, function(x){str_split(x, "--")[[1]][4]})
#       full_gene_sig <- unique(full_gene_sig[, c("gene_name", "gene_element_fisher_p", "p_adjust")])
# 
#       gene_sig <- left_join(full_gene_sig, SV_gene_count, by = c("gene_name" = "name"))
# 
#       # those genes where the value is NA do not have any SVs in this gene element
#       # they were probably just significant because the bin was significant
#       gene_sig[is.na(gene_sig$number_tumours_with_SV_in_gene), "number_tumours_with_SV_in_gene"] <- 0
# 
#       # get the bins back that these come from
#       bin_genes <- unique(intersect_df[, c("bin_key", "gene_name", "count", "effectsize", "count.pred")])
#       gene_sig <- left_join(gene_sig, bin_genes)
#       gene_sig <- data.frame(gene_sig)
# 
#       gene_sig$subtype <- subtype
#       gene_sig$event_type <- event
# 
#       all_event_list[[event]] <- gene_sig
#     }
# 
#   }
#    all_event_list <- do.call(rbind, all_event_list)
#    SVs_in_samples_event <- do.call(rbind, SVs_in_samples_event)
# 
#    all_subtype_list[[subtype]] <- all_event_list
#    SVs_in_samples_subtype[[subtype]] <- SVs_in_samples_event
# 
# }
# 
# all_subtype_list <- do.call(rbind, all_subtype_list)
# SVs_in_samples_subtype_list <- do.call(rbind, SVs_in_samples_subtype)
# 
# write.table(all_subtype_list, paste0(output_path, "significant_genes_50kb_with_TRA.txt"), sep = "\t", row.names = F, quote = F, col.names = T)
# write.table(SVs_in_samples_subtype_list, paste0(output_path, "significant_genes_50kb_with_TRA_SVs_in_samples.txt"), sep = "\t", row.names = F, quote = F, col.names = T)

all_subtype_list <- read.table(paste0(output_path, "significant_genes_50kb_with_TRA.txt"), head = T, sep = "\t")
SVs_in_samples_subtype_list <- read.table(paste0(output_path, "significant_genes_50kb_with_TRA_SVs_in_samples.txt"), head = T, sep = "\t")
######################################################################    PLOT
# plot all genes which have at least 5 SV in them 
eligible_genes <- unique(all_subtype_list[all_subtype_list$number_tumours_with_SV_in_gene >= 5, "gene_name"])

one_sv_genes <- unique(all_subtype_list[all_subtype_list$gene_name %in% eligible_genes, c("gene_name", "p_adjust", "number_tumours_with_SV_in_gene", "subtype", "event_type", "bin_key", "gr_id")])
one_sv_genes <- one_sv_genes[!is.na(one_sv_genes$gr_id), ]

# need to exclude some genes
olfactory_genes <- read.table(olfactory_gene_path, head = T)
olfactory_genes <- unique(olfactory_genes$gene_name)

tcga_expression_gene <- read.table(tcga_expression_gene_path, head = T)
tcga_expression_gene <- tcga_expression_gene[tcga_expression_gene$not_expressed_pantissue == TRUE, "gene_name"]

gtex_expression_gene <- read.table(gtex_expression_gene_path, head = T)
gtex_expression_gene <- gtex_expression_gene[gtex_expression_gene$ADENOCARCINOMA == 0, "gene_name"]

one_sv_genes <- one_sv_genes %>% filter(!gene_name %in% olfactory_genes)
one_sv_genes <- one_sv_genes %>% filter(!gene_name %in% tcga_expression_gene)
one_sv_genes <- one_sv_genes %>% filter(!gene_name %in% gtex_expression_gene)

one_sv_genes <- one_sv_genes[grep("LOC", one_sv_genes$gene, invert = T), ]
one_sv_genes <- one_sv_genes[grep(".*AS1", one_sv_genes$gene, invert = T), ]
one_sv_genes <- one_sv_genes[grep("NEAT1", one_sv_genes$gene, invert = T), ]
one_sv_genes <- one_sv_genes[grep("MALAT", one_sv_genes$gene, invert = T), ]

# check if there are ever several genes in the same bin
max_sv_genes <- one_sv_genes[one_sv_genes$gr_id == "gene", ] %>% group_by(subtype) %>% group_by(event_type) %>% group_by(bin_key) %>% slice_max(number_tumours_with_SV_in_gene)
max_sv_genes <- unique(max_sv_genes$gene_name)

max_sv_enhancer <- one_sv_genes[one_sv_genes$gr_id == "enhancer", ] %>% group_by(subtype) %>% group_by(event_type) %>% group_by(bin_key) %>% slice_max(number_tumours_with_SV_in_gene)
max_sv_enhancer <- unique(max_sv_enhancer$gene_name)

sv_genes <- max_sv_enhancer[max_sv_enhancer %in% max_sv_genes]
sv_genes <- unique(c(sv_genes, max_sv_genes))

sv_genes_info <- one_sv_genes[one_sv_genes$gene_name %in% sv_genes, ]
sv_genes_info$bin_key <- NULL
sv_genes_info <- unique(sv_genes_info)

# p value data
p_data <- unique(sv_genes_info[, c("gene_name", "p_adjust", "subtype", "event_type", "gr_id")])

gene_data <- unique(sv_genes_info[, c("gene_name", "number_tumours_with_SV_in_gene", "event_type", "gr_id")])
gene_data <- na.omit(gene_data)

# if for a gene we have an entry for "ALL_EVENTS" and another category we have to adjust the ALL_EVENTS value by subtractig the other value(s) from it
# gene_data_fix <- list()
# for(this_gene in unique(gene_data$gene_name)){
#   
#   print(this_gene)
#   
#   df <- gene_data[gene_data$gene_name == this_gene, ]
#   
#   # first let's check if we even have 2 entries for this gene
#   if(length(unique(df$event_type)) >= 2){
#     
#     # if we have gene and enhancer for this gene we need to split this up
#     if(length(unique(df$gr_id)) == 2){
#       gene_df <- df[df$gr_id == "gene", ]
#       gene_df[gene_df$event_type == "ALL_EVENTS", "fixed_number_tumours_with_SV_in_gene"] <- gene_df[gene_df$event_type == "ALL_EVENTS", "number_tumours_with_SV_in_gene"] - sum(gene_df[gene_df$event_type != "ALL_EVENTS", "number_tumours_with_SV_in_gene"])
#       gene_df[gene_df$event_type != "ALL_EVENTS", "fixed_number_tumours_with_SV_in_gene"] <- gene_df[gene_df$event_type != "ALL_EVENTS", "number_tumours_with_SV_in_gene"]
#       
#       enhancer_df <- df[df$gr_id == "enhancer", ]
#       enhancer_df[enhancer_df$event_type == "ALL_EVENTS", "fixed_number_tumours_with_SV_in_gene"] <- enhancer_df[enhancer_df$event_type == "ALL_EVENTS", "number_tumours_with_SV_in_gene"] - sum(enhancer_df[enhancer_df$event_type != "ALL_EVENTS", "number_tumours_with_SV_in_gene"])
#       enhancer_df[enhancer_df$event_type != "ALL_EVENTS", "fixed_number_tumours_with_SV_in_gene"] <- enhancer_df[enhancer_df$event_type != "ALL_EVENTS", "number_tumours_with_SV_in_gene"]
#       
#       out_df <- rbind(gene_df, enhancer_df)
#     } else {
#       
#       out_df <- df
#       out_df[out_df$event_type == "ALL_EVENTS", "fixed_number_tumours_with_SV_in_gene"] <- out_df[out_df$event_type == "ALL_EVENTS", "number_tumours_with_SV_in_gene"] - sum(out_df[out_df$event_type != "ALL_EVENTS", "number_tumours_with_SV_in_gene"])
#       out_df[out_df$event_type != "ALL_EVENTS", "fixed_number_tumours_with_SV_in_gene"] <- out_df[out_df$event_type != "ALL_EVENTS", "number_tumours_with_SV_in_gene"]
#     }
#   } else {
#     
#     out_df <- df
#     out_df$fixed_number_tumours_with_SV_in_gene <- out_df$number_tumours_with_SV_in_gene
#   }
#   
#   gene_data_fix[[this_gene]] <- out_df
# }
# 
# gene_data_fix <- do.call(rbind, gene_data_fix)
# gene_data_fix <- gene_data_fix[gene_data_fix$fixed_number_tumours_with_SV_in_gene != 0, ]

# sort stuff for plotting

# gene_data_fix$event_type <- factor(gene_data_fix$event_type, levels = c("DEL", "DUP", "h2hINV", "t2tINV", "TRA", "ALL_EVENTS"))
# gene_data_fix$gr_id <- factor(gene_data_fix$gr_id, levels = c("gene", "enhancer"))
# 
# gene_data_fix <- gene_data_fix[order(gene_data_fix$fixed_number_tumours_with_SV_in_gene, decreasing = T), ]
# 
# 
# gene_order <- gene_data_fix %>% group_by(gene_name) %>% summarise(sum(fixed_number_tumours_with_SV_in_gene))
# gene_order <- gene_order[order(gene_order$`sum(fixed_number_tumours_with_SV_in_gene)`, decreasing = T),]
# gene_order <- gene_order$gene_name
# 
# gene_data_fix <- gene_data_fix[gene_data_fix$fixed_number_tumours_with_SV_in_gene >= 0, ]
# 
# gene_data_fix$gene_name <- factor(gene_data_fix$gene_name, levels = unique(gene_order))
# p1 <- ggplot(gene_data_fix, aes(gr_id, fixed_number_tumours_with_SV_in_gene, fill = event_type)) +
#   geom_bar(stat = "identity", width = 0.98) +
#   scale_fill_manual(name = "SV type", 
#                     breaks = c("DEL", "DUP", "h2hINV", "t2tINV", "TRA", "ALL_EVENTS"), 
#                     labels = c("deletion", "duplication", "head-to-head\ninversion", "tail-to-tail\ninversion", "translocation", "all SV\ntypes"),
#                     values = c("#2c7bb6", "#d7191c", "#80cdc1", "#018571", "orange", "#ab7f20")) +
#   facet_grid(~gene_name, scale = "free_x", space = "free_x") +
#   theme_bw() +
#   theme() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
#         axis.text.y = element_text(size = 8),
#         plot.margin = unit(c(0,0,0.1,0), "cm"),
#         legend.justification = c(0,1),
#         panel.spacing = unit(0.1,'lines'),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   ylab("# tumours with SV") +
#   labs(x=NULL)

# combine p values for p value plot
# gene_data_fix$gene_id <- paste0(gene_data_fix$gene_name, "_", gene_data_fix$gr_id)
# 
# p_data$gene_id <- paste0(p_data$gene_name, "_", p_data$gr_id)
# p_data <- p_data[p_data$gene_id %in% unique(gene_data_fix$gene_id), ]
# p_data <- p_data %>% group_by(gene_id, subtype) %>% mutate(p_combine = poolr::fisher(p_adjust, adjust = "none")$p)
# p_data$gene_name <- factor(p_data$gene_name, levels = unique(gene_order))
# p_data$gr_id <- factor(p_data$gr_id, levels = c("gene", "enhancer"))
# 
# p_data$subtype <- factor(p_data$subtype, levels = c("ALL_SUBTYPES", "ADENOCARCINOMA", "SQUAMOUS_CELL", "LARGE_CELL", "SMALL_CELL", "MESOTHELIOMA", "ADENOSQUAMOUS", "MET_ADENOCARCINOMA"))
# 
# p2 <- ggplot(p_data, aes(gr_id, 1, fill = p_combine)) +
#   scale_fill_gradient(name = "p value", low = "black", high = "#d4d2cd") +
#   geom_tile() +
#   facet_grid(subtype~gene_name, scale = "free_x", space = "free_x", switch = "y") +
#   theme_bw() +
#   theme() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.margin = unit(c(0,0,0.1,0), "cm"),
#         legend.justification = c(0,1),
#         panel.spacing = unit(0.1,'lines'),
#         strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
#         strip.background = element_rect(linewidth = 0.3, fill = "white"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank()) +
#   ylab("") +
#   labs(x=NULL)

# make another plot indicating whetehr enhancer and gene are in the same bin
# bin_sig_genes <- one_sv_genes[, c("gene_name", "bin_key", "gr_id")]
# bin_sig_genes <- bin_sig_genes[bin_sig_genes$gene_name %in% gene_order, ]
# 
# same_bin <- list()
# for(name in unique(bin_sig_genes$gene_name)){
#   
#   df <- unique(bin_sig_genes[bin_sig_genes$gene_name == name, ])
#   
#   if(length(unique(df$gr_id)) > 1){
#     if(any(df[df$gr_id == "enhancer", "bin_key"] %in% df[df$gr_id == "gene", "bin_key"])){
#       df_out <- data.frame(gene_name = name, gr_id = unique(df$gr_id), enhancer_gene_same = TRUE)
#     } else{df_out <- data.frame(gene_name = name, gr_id = unique(df$gr_id))
#           df_out[df_out$gr_id == "enhancer", "enhancer_gene_same"] <- FALSE
#           df_out[df_out$gr_id == "gene", "enhancer_gene_same"] <- TRUE}
#   } else{df_out <- data.frame(gene_name = name, gr_id = unique(df$gr_id), enhancer_gene_same = TRUE)}
#   
#   same_bin[[name]] <- df_out
#   
# }
# 
# same_bin <- do.call(rbind, same_bin)
# 
# same_bin$gene_name <- factor(same_bin$gene_name, levels = unique(gene_order))
# 
# p3 <- ggplot(same_bin, aes(gr_id, 1, fill = enhancer_gene_same)) +
#             scale_fill_manual(name = "same bin", 
#                               breaks = c("TRUE", "FALSE"), 
#                               labels = c("TRUE", "FALSE"),
#                               values = c("beige", "darkblue")) +
#           geom_tile() +
#           facet_grid(~gene_name, scale = "free_x", space = "free_x", switch = "y") +
#           theme_bw() +
#           theme() +
#           scale_y_continuous(expand = c(0,0)) +
#           scale_x_discrete(expand = c(0,0)) +
#           theme_bw() +
#           theme(axis.text.x = element_blank(),
#                 axis.text.y = element_blank(),
#                 axis.ticks.x = element_blank(),
#                 axis.ticks.y = element_blank(),
#                 plot.margin = unit(c(0,0,0.1,0), "cm"),
#                 legend.justification = c(0,1),
#                 panel.spacing = unit(0.1,'lines'),
#                 strip.background = element_blank(),
#                 strip.text.x = element_blank(),
#                 panel.grid.major.x = element_blank(),
#                 panel.grid.minor.x = element_blank()) +
#           ylab("") +
#           labs(x=NULL)


# p_one_sv_genes <- plot_grid(p2, p1, ncol = 1, rel_heights = c(1, 2), align = "v", axis = "rlbt")
# 
# pdf(paste0(output_path, "fishook_genes.pdf"), width = 18, height = 6)
# p_one_sv_genes
# dev.off()


# # # # # # # # # # # # # # # # 
# make summary data frames

# which genes are drivers?
# combine p values for p value plot
# gene_data_fix$gene_id <- paste0(gene_data_fix$gene_name, "_", gene_data_fix$gr_id)

sv_genes_info$gene_id <- paste0(sv_genes_info$gene_name, "_", sv_genes_info$gr_id)

p_data$gene_id <- paste0(p_data$gene_name, "_", p_data$gr_id)
p_data <- p_data[p_data$gene_id %in% unique(sv_genes_info$gene_id), ]
p_data <- p_data %>% group_by(gene_id, subtype) %>% mutate(p_combine = poolr::fisher(p_adjust, adjust = "none")$p)

p_data_fix <- p_data
p_data_fix$p_adjust <- NULL
p_data_fix <- unique(p_data_fix)

write.table(p_data_fix, paste0(output_path, "SV_driver_genes.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

# also include patients
SVs_in_samples_subtype_list[SVs_in_samples_subtype_list$gr_id == "enhancer", "element_type"] <- "enhancer" 
SVs_in_samples_subtype_list[SVs_in_samples_subtype_list$gr_id != "enhancer", "element_type"] <- "gene" 
SVs_in_samples_subtype_list$geneID <- paste0(SVs_in_samples_subtype_list$gene_name, "_", SVs_in_samples_subtype_list$element_type, "_", SVs_in_samples_subtype_list$histology)     

# this way we'll miss the all_subtyps genes, add the 
SVs_in_samples_all <- SVs_in_samples_subtype_list
SVs_in_samples_all$geneID <- paste0(SVs_in_samples_subtype_list$gene_name, "_", SVs_in_samples_subtype_list$element_type, "_ALL_SUBTYPES")     

# SV_driver_genes$geneID <- paste0(SV_driver_genes$gene_id, "_", SV_driver_genes$subtype)
p_data$geneID <- paste0(p_data$gene_name, "_", p_data$gr_id, "_", p_data$subtype)
driver_SVs_in_samples <- SVs_in_samples_subtype_list[SVs_in_samples_subtype_list$geneID %in% unique(p_data$geneID), ]
driver_SVs_in_samples <- unique(driver_SVs_in_samples[, c("participant_id", "SV_CLASS", "histology", "gene_name", "element_type")])

all_driver_SV <- SVs_in_samples_all[SVs_in_samples_all$geneID %in% unique(p_data$geneID), ]
all_driver_SV <- unique(all_driver_SV[, c("participant_id", "SV_CLASS", "histology", "gene_name", "element_type")])

driver_SVs_in_samples <- rbind(driver_SVs_in_samples, all_driver_SV)

# since they are all also pancan drivers we need to add those 

pancan_add <- SVs_in_samples_subtype_list
pancan_add$gene_id <- paste0(pancan_add$gene_name, "_", pancan_add$element_type)
pancan_add <- pancan_add[pancan_add$gene_id %in% unique(p_data$gene_id), ]
pancan_add <- pancan_add[, c("participant_id", "SV_CLASS", "histology", "gene_name", "element_type")]
pancan_add$histology <- "PANCAN"

driver_SVs_in_samples <- rbind(driver_SVs_in_samples, pancan_add)
driver_SVs_in_samples <- unique(driver_SVs_in_samples)
 
write.table(driver_SVs_in_samples, paste0(output_path, "SV_driver_genes_in_samples.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

# finally count how many drivers there are per sample
SV_driver_count <- driver_SVs_in_samples[driver_SVs_in_samples$histology != "PANCAN", ]
SV_driver_count$geneID <- paste0(SV_driver_count$gene_name, "_", SV_driver_count$element_type)
SV_driver_count <- unique(SV_driver_count[, c("participant_id", "gene_name")])
SV_driver_count <- data.frame(table(SV_driver_count$participant_id))

write.table(SV_driver_count, paste0(output_path, "SV_driver_genes_per_sample_count.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

# plot this

data_samples <- unique(driver_SVs_in_samples[, c("participant_id", "histology", "gene_name", "element_type")])
# data_samples <- data_samples[data_samples$histology %in% c("PANCAN", "ADENOCARCINOMA", "SQUAMOUS_CELL"), ]
data_samples <- data.frame(table(data_samples[, c("histology", "gene_name", "element_type")]))

# data_samples[data_samples$histology == "ADENOCARCINOMA", "hist_group"] <- "Adenocarcinoma"
# data_samples[data_samples$histology == "SQUAMOUS_CELL", "hist_group"] <- "Squamouscell"
# data_samples[!data_samples$histology %in% c("SQUAMOUS_CELL", "ADENOCARCINOMA"), "hist_group"] <- "other"


data_samples <- data_samples[data_samples$Freq > 0, ]
data_samples$histology <- factor(data_samples$histology, levels = c("PANCAN", "ADENOCARCINOMA", "MET_ADENOCARCINOMA", "SQUAMOUS_CELL", "MET_SQUAMOUS_CELL", "ADENOSQUAMOUS", "CARCINOID",
                                                                    "LARGE_CELL", "MESOTHELIOMA", "SMALL_CELL", "MET_SMALL_CELL", "NEUROENDOCRINE_CARCINOMA"))

# plot the porportion of tumours to the number of smaples in that histology
hist_count <- data.frame(table(sample_table$histology))
colnames(hist_count) <- c("histology", "hist_count")
hist_count <- rbind(hist_count, data.frame(histology = "PANCAN", hist_count = 1011))

data_samples <- left_join(data_samples, hist_count)
data_samples$prop <- data_samples$Freq / data_samples$hist_count
data_samples$prop <- data_samples$prop*100

# do a prop test to see if there is enrichment of SVs in any of the histologies
# data_samples$geneID <- paste0(data_samples$gene_name, "_", data_samples$element_type)
# 
# SV_enrichment <- list()
# for(gene in unique(data_samples$geneID)){
#   
#   print(gene)
#   
#   hist_df <- list()
#   for(hist in unique(data_samples$histology)){
#     
#     print(hist)
#     
#     n_hist = hist_count[hist_count$histology == hist, "hist_count"]
#     if(length(n_hist) == 0){n_hist <- 0}
#     n_rest = nrow(sample_table) - n_hist
#     if(length(n_rest) == 0){n_rest <- 0}
#     
#     SV_hist = data_samples[data_samples$histology == hist & data_samples$geneID == gene, "Freq"]
#     if(length(SV_hist) == 0){SV_hist <- 0}
#     SV_rest = sum(data_samples[data_samples$histology != hist & data_samples$geneID == gene, "Freq"])
#     if(length(SV_rest) == 0){SV_rest <- 0}
#     
#     test_full_cohort <- prop.test(x = c(SV_hist, SV_rest), n = c(n_hist, n_rest), alternative = "greater")
#     test_within_histology <- prop.test(x = SV_hist, n = n_hist, p = 0.5, alternative = "greater")
#     
#     df_out <- data.frame(histology = hist,
#                          geneID = gene, 
#                          p_value_full_cohort = test_full_cohort$p.value,
#                          p_value_within_histology = test_within_histology$p.value)
#     
#     hist_df[[hist]] <- df_out
#     
#   }
#   
#   hist_df <- do.call(rbind, hist_df)
#   SV_enrichment[[gene]] <- hist_df
# }
# 
# SV_enrichment <- do.call(rbind, SV_enrichment)
# SV_enrichment$p_adjust_full_cohort <- p.adjust(SV_enrichment$p_value_full_cohort, method = "fdr")
# SV_enrichment$p_adjust_within_histology <- p.adjust(SV_enrichment$p_value_within_histology, method = "fdr")
# 
# SV_enrichment$significance <- NULL
# SV_enrichment[which(SV_enrichment$p_adjust_full_cohort < 0.0001), "significance"] <- "****"
# SV_enrichment[which(SV_enrichment$p_adjust_full_cohort < 0.001 & SV_enrichment$p_adjust_full_cohort >= 0.0001), "significance"] <- "***"
# SV_enrichment[which(SV_enrichment$p_adjust_full_cohort < 0.01 & SV_enrichment$p_adjust_full_cohort >= 0.001), "significance"] <- "**"
# SV_enrichment[which(SV_enrichment$p_adjust_full_cohort < 0.05 & SV_enrichment$p_adjust_full_cohort >= 0.01), "significance"] <- "*"
# SV_enrichment[which(SV_enrichment$p_adjust_full_cohort > 0.05), "significance"] <- ""
# 
# data_samples <- left_join(data_samples, SV_enrichment[, c("histology", "geneID", "significance")])

# gene_order <- data_samples %>% group_by(gene_name) %>% summarise(sum(prop))
# gene_order <- gene_order[order(gene_order$`sum(prop)`, decreasing = T), "gene_name"]
# gene_order <- as.character(gene_order$gene_name)

# data_samples$gene_name <- factor(data_samples$gene_name, levels = gene_order)
data_samples$element_type <- factor(data_samples$element_type, levels = c("gene", "enhancer"))
# data_samples$histology <- factor(data_samples$histology, levels = c("ADENOCARCINOMA",
#                                                                     "SQUAMOUS_CELL",
#                                                                     "PANCAN"))

# need to inlcude enhancer of genes where this is missing
# data_samples_fill <- list()
# for(gene in unique(data_samples$gene_name)){
#   
#   df <- data_samples[data_samples$gene_name == gene, ]
#   
#   if(all(!df$element_type %in% "enhancer")){
#     df2 <- df
#     df2$element_type <- "enhancer"
#     df2$Freq <- 0
#     df2$hist_count <- 0
#     df2$prop <- 0
#     df <- rbind(df, df2)
#   }
# 
#   data_samples_fill[[gene]] <- df
#   
# }
# 
# data_samples_fill <- do.call(rbind, data_samples_fill)
# data_samples_fill$prop <- data_samples_fill$prop * 100
# 
# data_samples_fill[data_samples_fill$Freq == 0, "Freq"] <- NA
# data_samples_fill[data_samples_fill$prop == 0, "prop"] <- NA

# histology_names_all <- c("Adenocarcinoma",
#                          "Adenocarcinoma\nmetastasis",
#                          "Squamous cell",
#                          "Squamous cell\nmetastasis",
#                          "Adenosquamous",
#                          "Carcinoid",
#                          "Large cell",
#                          "Mesothelioma",
#                          "Small cell",
#                          "Small cell\nmetastasis",
#                          "Neuroendocrine")
# 
# names(histology_names_all) <- levels(data_samples_fill$histology)

# add a collumn to data indicating if gene is independent of enhancer, ie where enhancer has a different p value to main gene

# main_genes <- list()
# enhancer_genes <- list()
# both_genes <- list()
# both_different_genes <- list()
# 
# for(gene in unique(data_samples_fill$gene_name)){
#   
#   print(gene)
#   
#   df <- data_samples_fill[data_samples_fill$gene_name == gene, ]
#   df <- df[!is.na(df$Freq), ]
#   
#   if(length(unique(df$element_type)) == 2){
#     
#     df_count <- df %>% group_by(histology) %>% summarize(length(unique(prop)))
#     
#     if(length(unique(df_count$`length(unique(prop))`)) == 1 & all(unique(df_count$`length(unique(prop))`) == 1) & nrow(df_count) != length(unique(df_count$histology))){
#       both_genes[[gene]] <- gene
#     } else {
#       both_different_genes[[gene]] <- gene
#     }
#     
#   } else {
#     
#     if(all(unique(df$element_type) %in% c("enhancer"))){
#       enhancer_genes[[gene]] <- gene
#     }
#     
#     if(all(unique(df$element_type) %in% c("gene"))){
#       main_genes[[gene]] <- gene
#     }
#     
#   }
#   
# }
# 
# main_genes <- names(main_genes)
# enhancer_genes <- names(enhancer_genes)
# both_genes <- names(both_genes)
# both_different_genes <- names(both_different_genes)
# 
# data_samples_fill[data_samples_fill$gene_name %in% main_genes, "gene_class"] <- "main_only"
# data_samples_fill[data_samples_fill$gene_name %in% enhancer_genes, "gene_class"] <- "enhancer_only"
# data_samples_fill[data_samples_fill$gene_name %in% both_genes, "gene_class"] <- "both"
# data_samples_fill[data_samples_fill$gene_name %in% both_different_genes, "gene_class"] <- "both_different"
# 

# reorder
# data_samples_fill <- data_samples_fill[order(data_samples_fill$gene_class, 
#                                              data_samples_fill$prop, decreasing = T), ]

data_samples <- data_samples[order(data_samples$element_type, 
                                   data_samples$prop, decreasing = T), ]


# data_samples_fill$gene_name <- factor(data_samples_fill$gene_name, levels = unique(data_samples_fill$gene_name))


# add the p value data
p_data_fill <- list()
for(gene in unique(p_data$gene_name)){
  
  df <- p_data[p_data$gene_name == gene, ]
  
  if(all(!df$gr_id %in% "enhancer")){
    df2 <- df
    df2$gr_id <- "enhancer"
    df2$p_adjust <- NA
    df2$p_combine <- NA
    
    df <- rbind(df, df2)
  }
  
  p_data_fill[[gene]] <- df
  
}

p_data_fill <- do.call(rbind, p_data_fill)

p_data_fill[p_data_fill$subtype == "ALL_SUBTYPES", "subtype"] <- "PANCAN"

# data_samples_fill_p <- left_join(data_samples_fill, unique(p_data_fill[!is.na(p_data_fill$gene_name) & p_data_fill$event_type == "ALL_EVENTS", c("gene_name", "gr_id", "subtype", "p_combine")]), by = c("histology" = "subtype",
#                                                                                                                                                                                                      "element_type" = "gr_id",
#                                                                                                                                                                                                      "gene_name" = "gene_name"))

data_samples_p <- left_join(data_samples, unique(p_data_fill[!is.na(p_data_fill$gene_name) & p_data_fill$event_type == "ALL_EVENTS", c("gene_name", "gr_id", "subtype", "p_combine")]), by = c("histology" = "subtype",
                                                                                                                                                                                                     "element_type" = "gr_id",
                                                                                                                                                                                                     "gene_name" = "gene_name"))


data_samples_p[which(data_samples_p$p_combine <= 0.05), "significance"] <- "significant"
data_samples_p[which(data_samples_p$p_combine > 0.05), "significance"] <- "not_significant"
data_samples_p[is.na(data_samples_p$p_combine), "significance"] <- "not_significant"
data_samples_p$significance <- factor(data_samples_p$significance, levels = c("significant", "not_significant"))
data_samples_p$gene_name <- factor(data_samples_p$gene_name, levels = unique(data_samples_p$gene_name))

data_samples_p$histology <- factor(data_samples_p$histology, levels = c("PANCAN", "ADENOCARCINOMA", "MET_ADENOCARCINOMA", "SQUAMOUS_CELL", "MET_SQUAMOUS_CELL", "ADENOSQUAMOUS", "CARCINOID",
                                                                        "LARGE_CELL", "MESOTHELIOMA", "SMALL_CELL", "MET_SMALL_CELL", "NEUROENDOCRINE_CARCINOMA"))

data_samples_p[is.na(data_samples_p$prop), "prop"] <- 0

# for some histologies enhancer/gene is still missing\
# data_samples_fill_p_en <- list()
# 
# for(gene in unique(data_samples_fill_p$gene_name)){
#   print(gene)
#   
#   df <- data_samples_fill_p[data_samples_fill_p$gene_name == gene, ]
#   
#   comb <- unique(paste0(df$histology, ":", df$element_type))
#   
#   if(length(comb) != 6){
#     
#     # which ones are missing
#     
#     if(!any(comb %in% "ADENOCARCINOMA:gene")){
#       missing1 <- "ADENOCARCINOMA:gene"
#     } else {missing1 <- NA}
#     
#     if(!any(comb %in% "ADENOCARCINOMA:enhancer")){
#       missing2 <- "ADENOCARCINOMA:enhancer"
#     } else {missing2 <- NA}
#     
#     if(!any(comb %in% "SQUAMOUS_CELL:gene")){
#       missing3 <- "SQUAMOUS_CELL:gene"
#     } else {missing3 <- NA}
#     
#     if(!any(comb %in% "SQUAMOUS_CELL:enhancer")){
#       missing4 <- "SQUAMOUS_CELL:enhancer"
#     } else {missing4 <- NA}
#     
#     if(!any(comb %in% "PANCAN:gene")){
#       missing5 <- "PANCAN:gene"
#     } else {missing5 <- NA}
#     
#     if(!any(comb %in% "PANCAN:enhancer")){
#       missing6 <- "PANCAN:enhancer"
#     } else {missing6 <- NA}
#     
#     all_mis <- c(missing1, missing2, missing3, missing4, missing5, missing6)
#     all_mis <- all_mis[!is.na(all_mis)]
#     
#     df_add <- data.frame(histology = c(sapply(all_mis, function(x){strsplit(x, ":")}[[1]][1])),
#                          gene_name = unique(df$gene_name),
#                          element_type = c(sapply(all_mis, function(x){strsplit(x, ":")}[[1]][2])),
#                          Freq = 0, 
#                          hist_count = 0,
#                          prop = 0,
#                          gene_class = NA, 
#                          p_combine = NA,
#                          significance = "not_significant")
#     
#     df_out <- rbind(df, df_add)
#     
#   } else {df_out <- df}
#   
#   data_samples_fill_p_en[[gene]] <- df_out
#   
# }

# data_samples_fill_p_en <- do.call(rbind, data_samples_fill_p_en)

# data_samples_fill_p_en$significance <- factor(data_samples_fill_p_en$significance, levels = c("significant", "not_significant"))
# data_samples_fill_p_en$gene_name <- factor(data_samples_fill_p_en$gene_name, levels = unique(data_samples_fill$gene_name))
# data_samples_fill_p_en$element_type <- factor(data_samples_fill_p_en$element_type, levels = c("gene", "enhancer"))
# data_samples_fill_p_en$histology <- factor(data_samples_fill_p_en$histology, levels = c("ADENOCARCINOMA",
#                                                                                   "SQUAMOUS_CELL",
#                                                                                   "PANCAN"))
 
# p1 <- ggplot(data_samples_fill_p_en, aes(histology, prop, fill = element_type, alpha = significance)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   scale_fill_manual(name = "gene element",
#                     labels = c("main gene", "enhancer"),
#                     values = c("#ef8a62","#67a9cf" )) +
#   scale_alpha_manual(values = c(1, 0.5),
#                      labels = c("significant", "not significant")) +
#   geom_vline(xintercept = c(1.45, 1.51, 2.45, 2.51), linewidth = 0.1) +
#   # geom_text(aes(label = significance), vjust = 1.5, colour = "black") +
#   facet_wrap(.~gene_name, nrow = 1) +
#   theme_bw() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0), labels = c("Adenocarcinoma", "Squamous cell", "Pan-lung")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
#         plot.margin = unit(c(0.1,0,0,0), "cm"),
#         legend.justification = c(0,1),
#         panel.spacing = unit(0.1,'lines'),
#         strip.background.x = element_rect(linewidth = 0.3, fill = "white"),
#         strip.text.x.top = element_text(size = 7),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank()) +
#   ylab("% tumours") +
#   labs(x=NULL)

# let's call gene CDS because in most cases that's what it is
data_samples_p[which(data_samples_p$element_type == "gene"), "element_type"] <- "CDS"
data_samples_p$element_type <- factor(data_samples_p$element_type, levels = c("CDS", "enhancer"))


histology_names <- c("Pan-lung",
                     "Adenocarcinoma",
                     "Adenocarcinoma\nmetastasis",
                     "Squamous cell",
                     "Small cell\nmetastasis",
                     "Adenosquamous",
                     "Carcinoid",
                     "Large cell",
                     "Mesothelioma",
                     "Small cell",
                     "Squamous cell\nmetastasis",
                     "Neuroendocrine")

names(histology_names) <- levels(data_samples_p$histology)

# make an indicator saying whter only ever the CDS is hit or not
CDS_only_genes <- data_samples_p %>% group_by(gene_name) %>% filter(all(element_type != "enhancer"))
CDS_only_genes <- unique(CDS_only_genes[, "gene_name"])
CDS_only_genes <- CDS_only_genes$gene_name
CDS_only_genes <- as.character(CDS_only_genes)
data_samples_p[data_samples_p$gene_name %in% CDS_only_genes, "CDS_only"] <- TRUE
data_samples_p[is.na(data_samples_p$CDS_only), "CDS_only"] <- FALSE

data_samples_p <- data_samples_p[order(data_samples_p$CDS_only,
                                       data_samples_p$Freq, decreasing = T), ]
data_samples_p$gene_name <- factor(data_samples_p$gene_name, levels = unique(data_samples_p$gene_name))

# for plotting purposes only keep LUAD and LUSC and make enhancer and CDS into one
data_samples_p_ll <- data_samples_p[data_samples_p$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), ]
data_samples_p_ll$prop <- NULL

# gotta combine the Freq and significance
data_samples_p_ll$hist_gene <- paste0(data_samples_p_ll$histology, "_", data_samples_p_ll$gene_name)

data_samples_p_ll_comb <- list()
for(this_one in unique(data_samples_p_ll$hist_gene)){
  
  df <- data_samples_p_ll[data_samples_p_ll$hist_gene == this_one, ]
  
  if(nrow(df) >1 ){
    
    df_out <- df
    df_out$Freq_combined <- sum(df$Freq)
    
    if(any(df$significance %in% "significant")){
      df_out$significance_combined <- "significant"
    } else { df_out$significance_combined <- "not_significant" }
    
    
  } else {
    df_out <- df
    df_out$Freq_combined <- df_out$Freq
    df_out$significance_combined <- df_out$significance
  }

  data_samples_p_ll_comb[[this_one]] <- df_out
  
}
data_samples_p_ll_comb <- do.call(rbind, data_samples_p_ll_comb)
data_samples_p_ll_comb$significance <- NULL
data_samples_p_ll_comb$Freq <- NULL
data_samples_p_ll_comb$element_type <- NULL
data_samples_p_ll_comb <- unique(data_samples_p_ll_comb)

data_samples_p_ll_comb$prop <- data_samples_p_ll_comb$Freq_combined / data_samples_p_ll_comb$hist_count
data_samples_p_ll_comb$prop <- data_samples_p_ll_comb$prop*100

data_samples_p_ll_comb$histology <- factor(data_samples_p_ll_comb$histology, levels = c("ADENOCARCINOMA", "SQUAMOUS_CELL"))

data_samples_p_ll_comb <- data_samples_p_ll_comb[order(data_samples_p_ll_comb$histology,
                                                      data_samples_p_ll_comb$Freq, decreasing = T), ]
data_samples_p_ll_comb$gene_name <- factor(data_samples_p_ll_comb$gene_name, levels = unique(data_samples_p_ll_comb$gene_name))
data_samples_p_ll_comb$p_combine <- NULL
data_samples_p_ll_comb <- unique(data_samples_p_ll_comb)

# classify if genes are significant in only LUAD, only LUSC, both or only pancancer
data_samples_p_ll_comb_add <- list()
for(gene in unique(data_samples_p_ll_comb$gene_name)){
  
  df <- data_samples_p_ll_comb[data_samples_p_ll_comb$gene_name == gene, ]
  
  if(length(unique(df[df$significance_combined == "significant", "histology"])) == 1){
    df$hist_specificity <- unique(df[df$significance_combined == "significant", "histology"])
  }
  
  if(length(unique(df[df$significance_combined == "significant", "histology"])) == 2){
    df$hist_specificity <- df$histology
  }
  
  if(length(unique(df[df$significance_combined == "significant", "histology"])) == 0){
    df$hist_specificity <- "panlung"
  }
  
  data_samples_p_ll_comb_add[[gene]] <- df
  
}
data_samples_p_ll_comb_add <- do.call(rbind, data_samples_p_ll_comb_add)

# now if it is not significant in the specific histology remove the label for coloring purposes
data_samples_p_ll_comb_add$histology <- as.character(data_samples_p_ll_comb_add$histology)
data_samples_p_ll_comb_add$hist_specificity <- as.character(data_samples_p_ll_comb_add$hist_specificity)
data_samples_p_ll_comb_add[data_samples_p_ll_comb_add$histology != data_samples_p_ll_comb_add$hist_specificity &
                             data_samples_p_ll_comb_add$significance_combined == "not_significant" &
                             data_samples_p_ll_comb_add$hist_specificity != "panlung", "hist_specificity"] <- ""

data_samples_p_ll_comb_add$hist_specificity <- factor(data_samples_p_ll_comb_add$hist_specificity, levels = c("ADENOCARCINOMA", "SQUAMOUS_CELL", "panlung", ""))
data_samples_p_ll_comb_add$histology <- factor(data_samples_p_ll_comb_add$histology, levels = c("ADENOCARCINOMA", "SQUAMOUS_CELL"))

data_samples_p_ll_comb_add <- data_samples_p_ll_comb_add %>% arrange(histology, desc(prop))
data_samples_p_ll_comb_add$gene_name <- factor(data_samples_p_ll_comb_add$gene_name, levels = unique(data_samples_p_ll_comb_add$gene_name))

p1 <- ggplot(data_samples_p_ll_comb_add, aes(gene_name, y = 1, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(name = "% tumours", low = "#d9d9d9", high = "black") +
  geom_tile(data = data_samples_p_ll_comb_add[data_samples_p_ll_comb_add$hist_specificity != "", ], aes(color = hist_specificity), fill = NA, linewidth = 3) +
  scale_color_manual(name = "histology specificity",
                     breaks = c("ADENOCARCINOMA", "SQUAMOUS_CELL", "panlung"),
                     labels = c("Adenocarcinoma", "Squamous cell", "Panlung"),
                     values = c("#67001f", "#053061", "#fee391")) +
  facet_grid(histology~gene_name, scales = "free_x", space = "free_x", switch = "y", labeller = labeller(histology = histology_names)) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x.top = element_text(vjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0.1,0,0.1,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0.1,'lines'),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.x = element_blank()) +
  ylab("") +
  labs(x=NULL)

# p1 <- ggplot(data_samples_fill, aes(element_type, histology, fill = prop)) +
#   geom_tile() +
#   scale_fill_gradient(low = "#deebf7", high = "#08306b", na.value = "white", name = "% tumours\nper histology\nwith SV") +
#   geom_text(aes(label = significance), vjust = 1.5, colour = "black") +
#   facet_nested(histology~hist_group + gene_class + gene_name , scale = "free", space = "free", switch = "y", labeller = labeller(histology = histology_names_all)) +
#   theme_bw() +
#   theme() +
#   scale_y_discrete(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.margin = unit(c(0.1,0,0,0), "cm"),
#         legend.justification = c(0,1),
#         panel.spacing = unit(0.1,'lines'),
#         strip.background.x = element_rect(linewidth = 0.3, fill = "white"),
#         strip.background.y = element_rect(linewidth = 0, fill = "white"),
#         strip.text.x.top = element_text(size = 7),
#         strip.text.y.left = element_text(size = 7, angle = 0),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank()) +
#   ylab("") +
#   labs(x=NULL)


# add enhancers to p data where this is missing

# p_data_fill <- list()
# for(gene in unique(p_data$gene_name)){
#   
#   df <- p_data[p_data$gene_name == gene, ]
#   
#   if(all(!df$gr_id %in% "enhancer")){
#     df2 <- df
#     df2$gr_id <- "enhancer"
#     df2$p_adjust <- NA
#     df2$p_combine <- NA
# 
#     df <- rbind(df, df2)
#   }
#   
#   p_data_fill[[gene]] <- df
#   
# }
# 
# p_data_fill <- do.call(rbind, p_data_fill)

# p_data_fill$subtype <- factor(p_data_fill$subtype, levels = c("ALL_SUBTYPES", "ADENOCARCINOMA", "SQUAMOUS_CELL", "LARGE_CELL", "SMALL_CELL", "MESOTHELIOMA", "ADENOSQUAMOUS", "MET_ADENOCARCINOMA"))
# p_data_fill$gene_name <- factor(p_data_fill$gene_name, levels = gene_order)
# p_data_fill$gr_id <- factor(p_data_fill$gr_id, levels = c("gene", "enhancer"))
# 
# histology_names <- c(
#   'ALL_SUBTYPES'="Pan-lung",
#   'ADENOCARCINOMA'="Adenocarcinoma",
#   'SQUAMOUS_CELL'="Squamous cell",
#   'SMALL_CELL'="Small cell"
# )
# 
# p2 <- ggplot(p_data_fill, aes(gr_id, 1, fill = p_combine)) +
#   scale_fill_gradient(name = "p value", low = "black", high = "#d4d2cd", na.value = "white", trans = "log10") +
#   geom_tile() +
#   facet_grid(subtype~gene_name, scale = "free_x", space = "free_x", switch = "both", labeller = labeller(subtype = histology_names)) +
#   theme_bw() +
#   theme() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.margin = unit(c(0,0,0.1,0), "cm"),
#         legend.justification = c(0,1),
#         panel.spacing = unit(0.1,'lines'),
#         strip.text.y.left = element_text(angle = 0, size = 7, hjust = 1),
#         strip.background.y = element_blank(),
#         strip.background.x = element_blank(),
#         strip.text.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         strip.placement = "outside") +
#   ylab("") +
#   labs(x=NULL)

# make a plot that shows the proportion of SVTYPES
# if a sample has several types 
class_df <- driver_SVs_in_samples[driver_SVs_in_samples$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), ]
class_df$element_type <- NULL

class_new_class_df <- class_df %>%
  # group_by(participant_id, gene_name, element_type) %>%
  # tally() %>%
  # ungroup() %>%
  # filter(n != 1) %>%
  # left_join(class_df, by = c("participant_id", "gene_name", "element_type")) %>%
  group_by(participant_id, gene_name) %>%
  mutate(new_SV_class = paste0(unique(SV_CLASS), collapse = ":")) %>%
  ungroup() %>%
  select(-SV_CLASS) %>%
  unique() 

class_new_class_df[grep(":", class_new_class_df$new_SV_class), "new_SV_class"] <- "multiple"

class_new_class_df_count <- data.frame(table(class_new_class_df[, c("gene_name", "new_SV_class")]))
class_new_class_df_count <- class_new_class_df_count[class_new_class_df_count$Freq != 0, ]

gene_count <- class_new_class_df_count %>% group_by(gene_name) %>% summarise(total = sum(Freq))

class_new_class_df_count <- left_join(class_new_class_df_count, gene_count)
class_new_class_df_count$proportion <- class_new_class_df_count$Freq / class_new_class_df_count$total

class_new_class_df_count$proportion <- class_new_class_df_count$proportion * 100

# also add enhancers where they are missing

# class_new_class_df_count_fill <- list()
# for(gene in unique(class_new_class_df_count$gene_name)){
#   
#   df <- class_new_class_df_count[class_new_class_df_count$gene_name == gene, ]
#   
#   if(all(!df$element_type %in% "enhancer")){
#     df2 <- df
#     df2$element_type <- "enhancer"
#     df2$Freq <- NA
#     df2$total <- NA
#     df2$proportion <- NA
#     
#     df <- rbind(df, df2)
#   }
#   
#   class_new_class_df_count_fill[[gene]] <- df
#   
# }
# 
# class_new_class_df_count_fill <- do.call(rbind, class_new_class_df_count_fill)

# add missing histology
# class_new_class_df_count_fill_hist <- list()
# 
# for(gene in unique(class_new_class_df_count_fill$gene_name)){
#   
#   print(gene)
#   
#   df <- class_new_class_df_count_fill[class_new_class_df_count_fill$gene_name == gene, ]
#   
#   if(length(unique(df$histology)) != 3){
#     
#     # which histology is missing?
#     if(!any(unique(df$histology) %in% "ADENOCARCINOMA")){
#       missing_hist <- "ADENOCARCINOMA"
#     }
#     
#     if(!any(unique(df$histology) %in% "SQUAMOUS_CELL")){
#       missing_hist <- "SQUAMOUS_CELL"
#     }
#     
#     df_mis <- data.frame(histology = missing_hist, 
#                          gene_name = unique(df$gene_name),
#                          element_type = "gene",
#                          new_SV_class = NA,
#                          Freq = 0,
#                          total = 0,
#                          proportion = 0)
#     
#     df_out <- rbind(df, df_mis)
#     
#   } else {df_out <- df}
#   
#   class_new_class_df_count_fill_hist[[gene]] <- df_out
#   
# }
# 
# class_new_class_df_count_fill_hist <- do.call(rbind, class_new_class_df_count_fill_hist)

# class_new_class_df_count_fill_hist$gene_name <- factor(class_new_class_df_count_fill_hist$gene_name, levels = unique(data_samples_fill$gene_name))
# class_new_class_df_count_fill_hist$element_type <- factor(class_new_class_df_count_fill_hist$element_type, levels = c("gene", "enhancer"))
# class_new_class_df_count_fill_hist$new_SV_class <- factor(class_new_class_df_count_fill_hist$new_SV_class, levels = unique(class_new_class_df_count_fill$new_SV_class))

class_new_class_df_count$gene_name <- factor(class_new_class_df_count$gene_name, levels = unique(data_samples_p_ll_comb_add$gene_name))
class_new_class_df_count$new_SV_class <- factor(class_new_class_df_count$new_SV_class, levels = c("DEL", "DUP", "h2hINV", "t2tINV", "TRA", "multiple"))

p3 <- ggplot(class_new_class_df_count, aes(gene_name, proportion, fill = new_SV_class)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(name = "SV class", values = c("#d53e4f", "#fc8d59", "#fee08b", "#e6f598", "#99d594", "#3288bd")) +
  facet_wrap(~gene_name, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0,0.3,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0.1,'lines'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% SVs") +
  labs(x=NULL)


# g <- ggplotGrob(p1)
# 
# #change colours and labels of strips
# colours <- c("ADENOCARCINOMA" = "#67001f",
#              "MET_ADENOCARCINOMA" = "#67001f80",
#              "SQUAMOUS_CELL" = "#053061",
#              "MET_SQUAMOUS_CELL" = "#05306180",
#              "ADENOSQUAMOUS" = "#66c2a5",
#              "CARCINOID" = "#fc8d62",
#              "LARGE_CELL" = "#54278f",
#              "MESOTHELIOMA" = "#e3309e",
#              "SMALL_CELL" = "#b89704",
#              "MET_SMALL_CELL" = "#b8970480",
#              "NEUROENDOCRINE_CARCINOMA" = "#0b8ca3")
# 
# strip_name_colours <- c("white", "black", "white", "black", "black",  "black", "white", "black", "black", "black", "black")
# 
# stript      <- which(grepl('strip-l', g$layout$name))[1:length(unique(colours))]
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
theme_margin <- theme(legend.box.margin = margin(10, 10, 10, 10))
legend  <- cowplot::get_legend(p1+ theme_margin)
legend1 <- cowplot::get_legend(p3+ theme_margin)
combineLegend <- cowplot::plot_grid(
  legend,
  legend1,
  nrow = 1)

p1 <- p1 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

p_all <- plot_grid(p1, p3, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "rlbt")
p_all_l <- plot_grid(p_all, combineLegend, ncol = 2, rel_widths = c(0.8, 0.2), align = "h", axis = "b")

pdf(paste0(output_path, "fishook_genes_new_LUAD_LUSC.pdf"), width = 10, height = 3)
p_all_l
dev.off()

## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# make a seperate plot showing how many SVs are in gene and in enhancer
data_samples_p_ll$element_type <- factor(data_samples_p_ll$element_type, levels = c("enhancer", "CDS"))
data_samples_p_ll <- data_samples_p_ll %>% arrange(element_type, desc(Freq))

data_samples_p_ll$gene_name <- factor(data_samples_p_ll$gene_name, levels = rev(unique(data_samples_p_ll$gene_name)))

p_type <- ggplot(data_samples_p_ll, aes(Freq, gene_name, fill = element_type)) +
            geom_bar(stat = "identity", position = position_dodge()) +
            scale_fill_manual(values = c("#5ab4ac", "#8e0152")) +
            xlab("# tumours") +
            ylab("") +
            theme_bw() +
            theme(legend.position = "top")
            
pdf(paste0(output_path, "fishook_genes_new_LUAD_LUSC_enhancer_CDS.pdf"), width = 3, height = 5)
p_type
dev.off()





# # # # # # # # # # # # # # # # 
#  investigate LRP1B

subtype = "ALL_SUBTYPES"
event = "ALL_EVENTS"

this_path <- paste0(fh_output_dir, "adjusted_model_bin_summary", event, "_", subtype, ".tsv")
fh_data <- read.table(this_path, head = T)
fh_data$bin_key <- paste0(fh_data$seqnames, ":", fh_data$start, ":", fh_data$end, ":", fh_data$bin_number)

# overlap with gene element definitions
fh_gr <- GRanges(seqnames = fh_data$seqnames, IRanges(start = fh_data$start, end = fh_data$end))
mcols(fh_gr)  <- DataFrame(fh_data[, c("bin_number", "hid", "p", "fdr", "count", "effectsize", "count.pred", "count.density", "count.pred.density", "eligible", "p.neg", "fdr.neg", "log10_p_observed", "bin_key")])

gene_def_gr <- GRanges(seqnames = gene_def$seqnames, IRanges(start = gene_def$start, end = gene_def$end))
mcols(gene_def_gr) <- DataFrame(gene_def[, c("name", "gene_key", "gene_name", "gr_id")])

overlaps  <- findOverlaps(fh_gr, gene_def_gr)
intersect <- pintersect(fh_gr[queryHits(overlaps)], gene_def_gr[subjectHits(overlaps)])
mcols(intersect) <- data.frame(mcols(intersect), mcols(gene_def_gr[subjectHits(overlaps)]))

intersect_df <- data.frame(intersect)

bin_gene_df <- intersect_df[, c("bin_key", "p", "fdr", "count", "effectsize", "count.pred", "count.density", "count.pred.density", "eligible", "p.neg", "fdr.neg", "log10_p_observed", "name")]
bin_gene_df <- unique(bin_gene_df)

# for each gene element combine the p values using fisher method
bin_gene_df <- bin_gene_df %>% group_by(name) %>% mutate(gene_element_fisher_p = poolr::fisher(p, adjust = "none")$p)

gene_df <- bin_gene_df[, c("name", "gene_element_fisher_p")]
gene_df <- unique(gene_df)
gene_df$p_adjust <- p.adjust(gene_df$gene_element_fisher_p, method = "fdr")

# take only those gene elements with a fdr p value < 0.05
gene_sig <- gene_df[gene_df$p_adjust < 0.05, ]

# now I need to get the number of samples with SVs in these gene elements
# if(subtype == "ALL_SUBTYPES"){
this_sv_gr <- sv_pair_gr
# } else { this_sv_gr <- sv_pair_gr[sv_pair_gr$histology == subtype, ] }
#
if(event == "ALL_EVENTS"){
  this_sv_gr <- this_sv_gr
} else { this_sv_gr <- this_sv_gr[this_sv_gr$SV_CLASS == event, ] }

# also get the unpaired SVs
# if(subtype == "ALL_SUBTYPES"){
this_sv_sng_gr <- sv_gr
# } else { this_sv_sng_gr <- sv_gr[sv_gr$histology == subtype, ] }

if(event == "ALL_EVENTS"){
  this_sv_sng_gr <- this_sv_sng_gr
} else { this_sv_sng_gr <- this_sv_sng_gr[this_sv_sng_gr$SV_CLASS == event, ] }

# get those SVs which have start and stop on the same bin
overlap_SV_bin <- findOverlaps(this_sv_sng_gr, fh_gr)
intersect_SV_bin <- pintersect(this_sv_sng_gr[queryHits(overlap_SV_bin)], fh_gr[subjectHits(overlap_SV_bin)])
mcols(intersect_SV_bin) <- data.frame(mcols(intersect_SV_bin), mcols(fh_gr[subjectHits(overlap_SV_bin)]))
intersect_SV_bin_df <- data.frame(intersect_SV_bin)

SV_ids_to_keep <- intersect_SV_bin_df[, c("SV_ID", "bin_key")]
SV_ids_to_keep <- unique(SV_ids_to_keep)
SV_ids_to_keep <- data.frame(table(SV_ids_to_keep$SV_ID))
SV_ids_to_keep <- SV_ids_to_keep[SV_ids_to_keep$Freq == 1, "Var1"]
SV_ids_to_keep <- as.character(SV_ids_to_keep)

# get the translocations back if there are translocations involved
if(any(grepl("TRA", SV_ids_to_keep))){
  TRAs_to_keep <- intersect_SV_bin_df[, c("SV_ID", "bin_key")]
  TRAs_to_keep <- unique(TRAs_to_keep[grep("TRA", TRAs_to_keep$SV_ID), ])
  TRAs_to_keep <- unique(as.character(TRAs_to_keep$SV_ID))
  
  SV_ids_to_keep <- c(SV_ids_to_keep, TRAs_to_keep)
}


this_sv_gr <- this_sv_gr[this_sv_gr$SV_ID %in% SV_ids_to_keep,]

# overlap the positions of the SVs with the gene elements
overlap_SV_gene <- findOverlaps(this_sv_gr, gene_def_gr)
intersect_SV_gene <- pintersect(this_sv_gr[queryHits(overlap_SV_gene)], gene_def_gr[subjectHits(overlap_SV_gene)])
mcols(intersect_SV_gene) <- data.frame(mcols(intersect_SV_gene), mcols(gene_def_gr[subjectHits(overlap_SV_gene)]))

intersect_SV_gene_df <- data.frame(intersect_SV_gene)

stk11 <- unique(intersect_SV_gene_df[intersect_SV_gene_df$gene_name == "LRP1B", c("participant_id", "SV_CLASS", "SV_ID", "histology", "sv_key", "gr_id")])

# check how many genes we find that arent in marias list yet

test <- c("SETDB1", "TERT", "EGFR", "MYC", "CCND1", "MDM2", "FOXA1", "NKX2-1", "U2AF1", "STK11", "CDKN2A", "LRP1B", "SOX2A", "FGFR1", "LZTR1", "MAPK1")

test <- unique(c(test, gene_order))
test <- unique(c(test, "MYH9", "GNAS", "RARA", "CHD4", "MDM2"))

maria_list <- read.delim("/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/driverDiscovery/driver_genetic_elements_2023-03-21.csv")
maria_list <- unique(maria_list$gene_id)

setdiff(test, maria_list)

# double check that we actually don't have any SVs in CDKN2A in mesothelioma

meso_sv_pair_gr <- sv_pair_gr[sv_pair_gr$histology == "MESOTHELIOMA", ]
cdkn2a_gene_def <- gene_def[gene_def$gene_id == "CDKN2A", ]
cdkn2a_gene_def_gr <- GRanges(seqnames = cdkn2a_gene_def$seqnames, IRanges(start = cdkn2a_gene_def$start, end = cdkn2a_gene_def$end))
mcols(cdkn2a_gene_def_gr) <- data.frame(cdkn2a_gene_def[, c("gr_id", "gene_id", "gene_key")])

overlap_meso <- findOverlaps(meso_sv_pair_gr, cdkn2a_gene_def_gr)
intersect_meso_gene <- pintersect(meso_sv_pair_gr[queryHits(overlap_meso)], cdkn2a_gene_def_gr[subjectHits(overlap_meso)])
mcols(intersect_meso_gene) <- data.frame(mcols(intersect_meso_gene), mcols(cdkn2a_gene_def_gr[subjectHits(overlap_meso)]))

# but we are only counting SVs that have start and stop within a bin

meso_sv_gr <- sv_gr[sv_gr$histology == "MESOTHELIOMA", ]

fh_gr <- GRanges(seqnames = fh_data$seqnames, IRanges(start = fh_data$start, end = fh_data$end))
mcols(fh_gr)  <- DataFrame(fh_data[, c("bin_number", "hid", "p", "fdr", "count", "effectsize", "count.pred", "count.density", "count.pred.density", "eligible", "p.neg", "fdr.neg", "log10_p_observed", "bin_key")])

meso_overlap_SV_bin <- findOverlaps(meso_sv_gr, fh_gr)
intersect_meso_SV_bin <- pintersect(meso_sv_gr[queryHits(meso_overlap_SV_bin)], fh_gr[subjectHits(meso_overlap_SV_bin)])
mcols(intersect_meso_SV_bin) <- data.frame(mcols(intersect_meso_SV_bin), mcols(fh_gr[subjectHits(meso_overlap_SV_bin)]))
intersect_meso_SV_bin_df <- data.frame(intersect_meso_SV_bin)

SV_ids_to_keep <- intersect_meso_SV_bin_df[, c("SV_ID", "bin_key")]
SV_ids_to_keep <- unique(SV_ids_to_keep)
SV_ids_to_keep <- data.frame(table(SV_ids_to_keep$SV_ID))
SV_ids_to_keep <- SV_ids_to_keep[SV_ids_to_keep$Freq == 1, "Var1"]
SV_ids_to_keep <- as.character(SV_ids_to_keep)

# get the translocations back if there are translocations involved
# if(any(grepl("TRA", SV_ids_to_keep))){
  TRAs_to_keep <- intersect_meso_SV_bin_df[, c("SV_ID", "bin_key")]
  TRAs_to_keep <- unique(TRAs_to_keep[grep("TRA", TRAs_to_keep$SV_ID), ])
  TRAs_to_keep <- unique(as.character(TRAs_to_keep$SV_ID))
  
  SV_ids_to_keep <- c(SV_ids_to_keep, TRAs_to_keep)
# }


this_meso_sv_gr <- meso_sv_gr[meso_sv_gr$SV_ID %in% SV_ids_to_keep,]

# overlap the positions of the SVs with the gene elements
overlap_SV_gene <- findOverlaps(meso_sv_gr, cdkn2a_gene_def_gr)
intersect_SV_gene <- pintersect(meso_sv_gr[queryHits(overlap_SV_gene)], cdkn2a_gene_def_gr[subjectHits(overlap_SV_gene)])
mcols(intersect_SV_gene) <- data.frame(mcols(intersect_SV_gene), mcols(cdkn2a_gene_def_gr[subjectHits(overlap_SV_gene)]))

intersect_SV_gene_df <- data.frame(intersect_SV_gene)











# # which are the genes which get removed?
# null_sv_genes <- unique(all_subtype_list[all_subtype_list$number_tumours_with_SV_in_gene == 0, c("name", "p_adjust", "number_tumours_with_SV_in_gene", "subtype", "event_type")])
# null_sv_genes$element <- sapply(null_sv_genes$name, function(x){str_split(x, "--")[[1]][2]})
# null_sv_genes$gene <- sapply(null_sv_genes$name, function(x){str_split(x, "--")[[1]][3]})
# null_sv_genes$gene <- gsub("__and__", "_", null_sv_genes$gene)
# null_sv_genes$gene_element <- paste0(null_sv_genes$gene, "_", null_sv_genes$element)
# 
# x <- unique(null_sv_genes[, c("gene", "element")])
# 
# # # # plot removing all genes with SV count less than SV breakpoint count in bin
# # plot all genes which have at least 1 SV in them 
# count_sv_genes <- unique(all_subtype_list[which(all_subtype_list$number_tumours_with_SV_in_gene >= all_subtype_list$count), c("name", "p_adjust", "number_tumours_with_SV_in_gene", "count", "subtype", "event_type")])
# count_sv_genes <- count_sv_genes[count_sv_genes$number_tumours_with_SV_in_gene != 0, ]
# count_sv_genes$element <- sapply(count_sv_genes$name, function(x){str_split(x, "--")[[1]][2]})
# count_sv_genes$gene <- sapply(count_sv_genes$name, function(x){str_split(x, "--")[[1]][3]})
# count_sv_genes$gene <- gsub("__and__", "_", count_sv_genes$gene)
# 
# count_sv_genes$event_type <- factor(count_sv_genes$event_type, levels = c("DEL", "DUP", "h2hINV", "t2tINV"))
# count_sv_genes$element <- factor(count_sv_genes$element, levels = c("CDS", "enhancer", "ss", "3primeUTR", "5primeUTR", "promoter", "shortRNA", "lincRNA", "lincRNA_promoter", "lincRNA_ss", "miRNA"))
# 
# count_sv_genes <- count_sv_genes[order(count_sv_genes$number_tumours_with_SV_in_gene,
#                                        count_sv_genes$element, decreasing = T), ]
# count_gene_order <- count_sv_genes$gene
# 
# count_sv_genes$gene <- factor(count_sv_genes$gene, levels = unique(count_gene_order))
# 
# p1_count <- ggplot(count_sv_genes, aes(gene, number_tumours_with_SV_in_gene, fill = event_type)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(name = "SV type", 
#                     breaks = c("DEL", "DUP", "h2hINV", "t2tINV"), 
#                     labels = c("deletion", "duplication", "head-to-head\ninversion", "tail-to-tail\ninversion"),
#                     values = c("#2c7bb6", "#d7191c", "#80cdc1", "#018571")) +
#   facet_grid(element~., switch = "y") +
#   theme_bw() +
#   theme() +
#   scale_y_continuous(expand = c(0,0), position = "right") +
#   scale_x_discrete(expand = c(0,0)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
#         axis.text.y = element_text(size = 8),
#         plot.margin = unit(c(0,0,0.1,0), "cm"),
#         legend.justification = c(0,1),
#         panel.spacing = unit(0.1,'lines'),
#         strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
#         strip.background = element_rect(linewidth = 0.3, fill = "white")) +
#   ylab("# SVs in gene element") +
#   labs(x=NULL)
# 
# count_sv_genes$subtype <- factor(count_sv_genes$subtype, levels = c("ALL_SUBTYPES", "ADENOCARCINOMA", "SQUAMOUS_CELL", "MESOTHELIOMA"))
# 
# p2_count <- ggplot(count_sv_genes, aes(gene, 1, fill = p_adjust)) +
#   scale_fill_gradient(name = "p value", low = "black", high = "#d4d2cd") +
#   geom_tile() +
#   facet_grid(subtype~., switch = "y") +
#   theme_bw() +
#   theme() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.margin = unit(c(0.1,0,0,0), "cm"),
#         legend.justification = c(0,1),
#         panel.spacing = unit(0.1,'lines'),
#         strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
#         strip.background = element_rect(linewidth = 0.3, fill = "white")) +
#   ylab("") +
#   labs(x=NULL)
# 
# p_count_genes <- plot_grid(p2_count, p1_count, ncol = 1, rel_heights = c(1, 7), align = "v", axis = "rlbt")
# 
# 
