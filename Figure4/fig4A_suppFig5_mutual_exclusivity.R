# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                     mutual exclusivity - plotting                   # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
##################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
library(tidyr)
library(stringr)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) #change to hg19 depending on input data
library(org.Hs.eg.db)


###############################################################################
##################################################################### file paths
results_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/mutual_exclusivity/PANCANpostproc_drivers_cooc_excl_CGConly_ssMerge_tumorSubtypeSpecUPD_with_CN_and_SV.csv"
output_path  <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/mutual_exclusivity/"
cn_drivers   <- '/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copynumber_peak_driver_per_patient.txt'

###############################################################################
#####################################################################       MAIN

# run for only CDS 
results <- read.table(results_path, head = T, sep = "\t")

# only keep the results from testing coding
results <- results %>% filter(comparison == "cod")

# where a pancan pair is not supported by the individual subtype, remove the asterics
results[which(results$support_by_indivTT == FALSE), "plotLab"] <- ""

# now we need to do everthing per subtype
for(subtype in unique(results$tumor_subtype)){
  
  print(subtype)
  subtype_results <- results[results$tumor_subtype == subtype, ]
  
  # need to sort out mode, cause every pair has a result for exclusivity and for co-occurence
  results_long <- unique(subtype_results[, c("gr_id_1", "gene_id_1", "gr_id_2", "gene_id_2", "mode", "discover.p.value", "plotLab")])
  results_long <- results_long %>% tidyr::pivot_wider(names_from = "mode", values_from = c("discover.p.value", "plotLab"))
  results_long <- results_long %>% dplyr::mutate(plotLab = ifelse(plotLab_exclusivity == "", `plotLab_co-occurrence`, plotLab_exclusivity)) %>% dplyr::select(-plotLab_exclusivity, -`plotLab_co-occurrence`)
  # results_long <- reshape2::dcast(results_long, gr_id_1 + gene_id_1 + gr_id_2 + gene_id_2 + plotLab ~ mode,  value.var = "discover.p.value")
  
  # always keep the smallest p value
  colnames(results_long) <- c("gr_id_1", "gene_id_1", "gr_id_2", "gene_id_2", "exclusivity", "co_occurence", "plotLab")
  results_long <- results_long %>% dplyr::rowwise() %>% dplyr::mutate(p_use = min(co_occurence, exclusivity, na.rm = T))
  results_long[results_long$p_use == Inf, "p_use"] <- NA
  
  # calculate the -log10(raw p value)
  results_long$log10p <- -log10(results_long$p_use)
  
  # flip the p values of mutual exclusivity around so that it's visible on the heatmap if its more co-occuring or exclusive
  results_long <- data.frame(results_long)
  results_long[which(results_long$exclusivity == results_long$p_use), "label"] <- "exclusive"
  results_long[which(results_long$co_occurence == results_long$p_use), "label"] <- "co_occurence"
  results_long[which(results_long$label == "exclusive"), "log10p"] <- results_long[which(results_long$label == "exclusive"), "log10p"] *-1
  
  # put a chr in front of the chromosome ones
  results_long <- data.frame(results_long)
  results_long[results_long$gr_id_1 == "CN", "gene_id_1"] <- paste0("chr", results_long[results_long$gr_id_1 == "CN", "gene_id_1"])
  results_long[results_long$gr_id_2 == "CN", "gene_id_2"] <- paste0("chr", results_long[results_long$gr_id_2 == "CN", "gene_id_2"])
  
  results_mat <- results_long
  
  # where a gene is in the same locations as a copynumber segment and there is muttual exclusivity let's exclude that
  cn_peaks <- read.table(cn_drivers, head = T, sep = "\t")
  cn_peaks <- unique(cn_peaks[, c("new_peak_name", "significant_CN_genes_on_peak")])
  
  # reformat the peak names to match the peak names in the other table
  cn_peaks$new_peak_name <- sub("SQUAMOUS_CELL", "SQUAMOUSCELL", cn_peaks$new_peak_name)
  cn_peaks$new_peak_name <- sub("homozygous_deletion", "homozygousdeletion", cn_peaks$new_peak_name)
  cn_peaks <- cn_peaks %>% separate(col = new_peak_name, sep = "_", into = c("histology", "type", "chr", "start", "end"))
  cn_peaks[cn_peaks$type == "amplification", "type"] <- "amp"
  cn_peaks[cn_peaks$type == "homozygousdeletion", "type"] <- "homdel"
  cn_peaks$new_name <- paste0(cn_peaks$chr, "_", cn_peaks$start, "_", cn_peaks$end, ",", cn_peaks$type)
  cn_peaks$new_name <- sub("_NA", "", cn_peaks$new_name)
  
  # get all the significant genes that involve a copy number peak
  cn_sig <- results_mat[results_mat$plotLab != "",]
  cn_sig <- cn_sig[cn_sig$gr_id_1 == "CN" | cn_sig$gr_id_2 == "CN", ]                  
  
  # check if the gene is on the copy number peak
  cn_sig <- left_join(cn_sig, unique(cn_peaks[, c("new_name", "significant_CN_genes_on_peak")]), by = c("gene_id_1" = "new_name"))
  cn_sig <- left_join(cn_sig, unique(cn_peaks[, c("new_name", "significant_CN_genes_on_peak")]), by = c("gene_id_2" = "new_name"))
  cn_sig[is.na(cn_sig$significant_CN_genes_on_peak.x), "significant_CN_genes_on_peak.x"] <- cn_sig[is.na(cn_sig$significant_CN_genes_on_peak.x) & !is.na(cn_sig$significant_CN_genes_on_peak.y), "significant_CN_genes_on_peak.y"]
  cn_sig$significant_CN_genes_on_peak.y <- NULL
  
  # whereever the snv gene appears in the gene on peak column we shall eliminiate the plot label
  cn_sig_removed <- list()
  for(i in 1:nrow(cn_sig)){
    
    df <- cn_sig[i, ]
    df[grepl(df$gene_id_2, df$significant_CN_genes_on_peak.x), "plotLab"] <- ""
    df[grepl(df$gene_id_1, df$significant_CN_genes_on_peak.x), "plotLab"] <- ""
    cn_sig_removed[[i]] <- df
  }
  
  cn_sig_removed <- do.call(rbind, cn_sig_removed)
  
  # go one step further and remove significance for genes which are on the same chromsome as the cn peaks
  cn_sig_removed[grepl("chr", cn_sig_removed$gene_id_1), "chr"] <- sapply(cn_sig_removed$gene_id_1[grepl("chr", cn_sig_removed$gene_id_1)], function(x){str_split(x, "_")[[1]][1]})
  cn_sig_removed[grepl("chr", cn_sig_removed$gene_id_2), "chr"] <- sapply(cn_sig_removed$gene_id_2[grepl("chr", cn_sig_removed$gene_id_2)], function(x){str_split(x, "_")[[1]][1]})
  
  # need to get the chromosomes of the genes
  txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes                 <- genes(txdb, single.strand.genes.only=FALSE)
  genes                 <- data.frame(genes)
  
  # extract HUGO symbols
  gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)
  
  # merge with gene locations
  gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))
  gene_names <- gene_names[!is.na(gene_names$symbol), ]
  
  cn_sig_removed <- left_join(cn_sig_removed, gene_names[, c("symbol", "seqnames")], by = c("gene_id_1" = "symbol"))
  cn_sig_removed <- left_join(cn_sig_removed, gene_names[, c("symbol", "seqnames")], by = c("gene_id_2" = "symbol"))
  
  cn_sig_removed[is.na(cn_sig_removed$seqnames.x) & !is.na(cn_sig_removed$seqnames.y), "seqnames.x"] <- cn_sig_removed[is.na(cn_sig_removed$seqnames.x) & !is.na(cn_sig_removed$seqnames.y), "seqnames.y"]
  cn_sig_removed$seqnames.y <- NULL
  
  cn_sig_removed[which(cn_sig_removed$chr == cn_sig_removed$seqnames.x), "plotLab"] <- ""
  
  # now we need to change the big data frame accordingly
  cn_sig_removed$both_gene_id <- paste0(cn_sig_removed$gene_id_1, ":", cn_sig_removed$gene_id_2)
  plot_lab_remove <- unique(cn_sig_removed[cn_sig_removed$plotLab == "", "both_gene_id"])
  
  results_mat$both_gene_id <- paste0(results_mat$gene_id_1, ":", results_mat$gene_id_2)
  results_mat[results_mat$both_gene_id %in% plot_lab_remove, "plotLab"] <- ""
  results_mat$both_gene_id <- NULL
  
  # gotta flip some things around according to my preferred order gene elements
  results_mat$gene_id_1 <- gsub("chr([0-9])_","chr0\\1_", results_mat$gene_id_1)
  results_mat$gene_id_2 <- gsub("chr([0-9])_","chr0\\1_", results_mat$gene_id_2)
  
  element_order <- c("CDS(SNV)", "CDS(SV)", "CDS(SNV, SV)", "CN")
  
  
  results_mat$gene_gr_1 <- paste0(results_mat$gr_id_1, ":", results_mat$gene_id_1)
  results_mat$gene_gr_2 <- paste0(results_mat$gr_id_2, ":", results_mat$gene_id_2)
  
  
  factorOrder <- tibble(results_mat) %>% mutate(gr_id_1 = factor(gr_id_1, levels = rev(element_order))) %>% arrange(gr_id_1, desc(gene_id_1)) %>% dplyr::select(gr_id_1, gene_id_1, gene_gr_1) %>% unique() %>% pull(gene_gr_1)
  
  tmp <- results_mat %>% 
    dplyr::select(gene_gr_1, gene_gr_2, log10p, plotLab) %>% 
    dplyr::filter(!is.na(log10p)) %>%
    unique() %>%
    dplyr::full_join(results_mat %>% dplyr::select(gene_gr_1 = gene_gr_2, gene_gr_2 = gene_gr_1, log10p, plotLab) %>% dplyr::filter(!is.na(log10p)) %>% unique(), by = c("gene_gr_1", "gene_gr_2", "log10p", "plotLab")) %>%
    unique() %>%
    dplyr::group_by(gene_gr_1, gene_gr_2) %>%
    dplyr::summarize_all(na.omit) %>%
    unique() %>%
    dplyr::mutate(gene_gr_1 = factor(gene_gr_1, levels = factorOrder),
                  gene_gr_2 = factor(gene_gr_2, levels = factorOrder)) %>%
    dplyr::arrange(gene_gr_1, gene_gr_2) %>%
    dplyr::mutate(i = as.numeric(gene_gr_2),
                  j = as.numeric(gene_gr_1)) %>%
    dplyr::filter(j >= i) %>%
    dplyr::select(-i, -j)
  
  tmp <- tmp %>%
    tidyr::separate(gene_gr_1, into = c("gr_id_1", "gene_id_1"), sep = ":") %>% 
    dplyr::mutate(gr_id_1 = factor(gr_id_1, levels = unique(gr_id_1))) %>%
    dplyr::mutate(gene_id_1 = factor(gene_id_1, levels = unique(gene_id_1))) %>%
    tidyr::separate(gene_gr_2, into = c("gr_id_2", "gene_id_2"), sep = ":") %>% 
    dplyr::mutate(gr_id_2 = factor(gr_id_2, levels = unique(gr_id_2))) %>%
    dplyr::mutate(gene_id_2 = factor(gene_id_2, levels = unique(gene_id_2))) %>%
    filter(!is.na(log10p))
  

  coding_sig <- tmp[tmp$plotLab != "", ]
  
  p_cod <- ggplot(tmp, aes(forcats::fct_rev(gene_id_1), gene_id_2, fill = log10p, label = plotLab)) +
    geom_tile() +
    geom_text(size = 3, vjust = 0.8) +
    scale_fill_gradient2(
      low = "darkred",
      mid = "white",
      high = "darkblue",
      midpoint = 0,
      na.value = "grey",
      guide = "colourbar") +
    xlab("") +
    ylab("") +
    facet_grid(forcats::fct_rev(gr_id_2) ~ forcats::fct_rev(gr_id_1), scales = "free", space = "free", switch = "both") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          axis.text.y = element_text(size = ),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text.x.bottom = element_text(size = 7),
          strip.text.y.left = element_text(size = 7),
          strip.background = element_rect(fill="white"))
  
  pdf(paste0(output_path, subtype, "_mutual_exclusivity_coding.pdf"), width = 10, height = 9)
  print(p_cod)
  dev.off()
  
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# run for ALL  gene elements

results <- read.table(results_path, head = T, sep = "\t")

# only keep the results from testing coding
results <- results %>% filter(comparison == "all")

# where a pancan pair is not supported by the individual subtype, remove the asterics
results[which(results$support_by_indivTT == FALSE), "plotLab"] <- ""

# now we need to do everthing per subtype
for(subtype in unique(results$tumor_subtype)){

  print(subtype)
  subtype_results <- results[results$tumor_subtype == subtype, ]
  
  # need to sort out mode, cause every pair has a result for exclusivity and for co-occurence
  results_long <- unique(subtype_results[, c("gr_id_1", "gene_id_1", "gr_id_2", "gene_id_2", "mode", "discover.p.value", "plotLab")])
  results_long <- results_long %>% tidyr::pivot_wider(names_from = "mode", values_from = c("discover.p.value", "plotLab"))
  results_long <- results_long %>% dplyr::mutate(plotLab = ifelse(plotLab_exclusivity == "", `plotLab_co-occurrence`, plotLab_exclusivity)) %>% dplyr::select(-plotLab_exclusivity, -`plotLab_co-occurrence`)
  # results_long <- reshape2::dcast(results_long, gr_id_1 + gene_id_1 + gr_id_2 + gene_id_2 + plotLab ~ mode,  value.var = "discover.p.value")
  
  # always keep the smallest p value
  colnames(results_long) <- c("gr_id_1", "gene_id_1", "gr_id_2", "gene_id_2", "exclusivity", "co_occurence", "plotLab")
  results_long <- results_long %>% dplyr::rowwise() %>% dplyr::mutate(p_use = min(co_occurence, exclusivity, na.rm = T))
  results_long[results_long$p_use == Inf, "p_use"] <- NA
  
  # calculate the -log10(raw p value)
  results_long$log10p <- -log10(results_long$p_use)
  
  # flip the p values of mutual exclusivity around so that it's visible on the heatmap if its more co-occuring or exclusive
  results_long[which(results_long$exclusivity == results_long$p_use), "label"] <- "exclusive"
  results_long[which(results_long$co_occurence == results_long$p_use), "label"] <- "co_occurence"
  results_long[which(results_long$label == "exclusive"), "log10p"] <- results_long[which(results_long$label == "exclusive"), "log10p"] *-1
  
  # put a chr in front of the chromosome ones
  results_long <- data.frame(results_long)
  results_long[results_long$gr_id_1 == "CN", "gene_id_1"] <- paste0("chr", results_long[results_long$gr_id_1 == "CN", "gene_id_1"])
  results_long[results_long$gr_id_2 == "CN", "gene_id_2"] <- paste0("chr", results_long[results_long$gr_id_2 == "CN", "gene_id_2"])
  
  results_mat <- results_long
  
  # where a gene is in the same locations as a copynumber segment and there is muttual exclusivity let's exclude that
  cn_peaks <- read.table(cn_drivers, head = T, sep = "\t")
  cn_peaks <- unique(cn_peaks[, c("new_peak_name", "significant_CN_genes_on_peak")])
  
  # reformat the peak names to match the peak names in the other table
  cn_peaks$new_peak_name <- sub("SQUAMOUS_CELL", "SQUAMOUSCELL", cn_peaks$new_peak_name)
  cn_peaks$new_peak_name <- sub("homozygous_deletion", "homozygousdeletion", cn_peaks$new_peak_name)
  cn_peaks <- cn_peaks %>% separate(col = new_peak_name, sep = "_", into = c("histology", "type", "chr", "start", "end"))
  cn_peaks[cn_peaks$type == "amplification", "type"] <- "amp"
  cn_peaks[cn_peaks$type == "homozygousdeletion", "type"] <- "homdel"
  cn_peaks$new_name <- paste0(cn_peaks$chr, "_", cn_peaks$start, "_", cn_peaks$end, ",", cn_peaks$type)
  cn_peaks$new_name <- sub("_NA", "", cn_peaks$new_name)
  
  # get all the significant genes that involve a copy number peak
  cn_sig <- results_mat[results_mat$plotLab != "",]
  cn_sig <- cn_sig[cn_sig$gr_id_1 == "CN" | cn_sig$gr_id_2 == "CN", ]                  
  
  # check if the gene is on the copy number peak
  cn_sig <- left_join(cn_sig, unique(cn_peaks[, c("new_name", "significant_CN_genes_on_peak")]), by = c("gene_id_1" = "new_name"))
  cn_sig <- left_join(cn_sig, unique(cn_peaks[, c("new_name", "significant_CN_genes_on_peak")]), by = c("gene_id_2" = "new_name"))
  cn_sig[is.na(cn_sig$significant_CN_genes_on_peak.x), "significant_CN_genes_on_peak.x"] <- cn_sig[is.na(cn_sig$significant_CN_genes_on_peak.x) & !is.na(cn_sig$significant_CN_genes_on_peak.y), "significant_CN_genes_on_peak.y"]
  cn_sig$significant_CN_genes_on_peak.y <- NULL
  
  # whereever the snv gene appears in the gene on peak column we shall eliminiate the plot label
  cn_sig_removed <- list()
  for(i in 1:nrow(cn_sig)){
    
    df <- cn_sig[i, ]
    df[grepl(df$gene_id_2, df$significant_CN_genes_on_peak.x), "plotLab"] <- ""
    df[grepl(df$gene_id_1, df$significant_CN_genes_on_peak.x), "plotLab"] <- ""
    cn_sig_removed[[i]] <- df
  }
  
  cn_sig_removed <- do.call(rbind, cn_sig_removed)
  
  # go one step further and remove significance for genes which are on the same chromsome as the cn peaks
  cn_sig_removed[grepl("chr", cn_sig_removed$gene_id_1), "chr"] <- sapply(cn_sig_removed$gene_id_1[grepl("chr", cn_sig_removed$gene_id_1)], function(x){str_split(x, "_")[[1]][1]})
  cn_sig_removed[grepl("chr", cn_sig_removed$gene_id_2), "chr"] <- sapply(cn_sig_removed$gene_id_2[grepl("chr", cn_sig_removed$gene_id_2)], function(x){str_split(x, "_")[[1]][1]})
  
  # need to get the chromosomes of the genes
  txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes                 <- genes(txdb, single.strand.genes.only=FALSE)
  genes                 <- data.frame(genes)
  
  # extract HUGO symbols
  gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)
  
  # merge with gene locations
  gene_names            <- left_join(genes, gene_symbols, by = c("group_name" = "gene_id"))
  gene_names <- gene_names[!is.na(gene_names$symbol), ]
  
  cn_sig_removed <- left_join(cn_sig_removed, gene_names[, c("symbol", "seqnames")], by = c("gene_id_1" = "symbol"))
  cn_sig_removed <- left_join(cn_sig_removed, gene_names[, c("symbol", "seqnames")], by = c("gene_id_2" = "symbol"))
  
  cn_sig_removed[is.na(cn_sig_removed$seqnames.x) & !is.na(cn_sig_removed$seqnames.y), "seqnames.x"] <- cn_sig_removed[is.na(cn_sig_removed$seqnames.x) & !is.na(cn_sig_removed$seqnames.y), "seqnames.y"]
  cn_sig_removed$seqnames.y <- NULL
  
  cn_sig_removed[which(cn_sig_removed$chr == cn_sig_removed$seqnames.x), "plotLab"] <- ""
  
  # now we need to change the big data frame accordingly
  cn_sig_removed$both_gene_id <- paste0(cn_sig_removed$gene_id_1, ":", cn_sig_removed$gene_id_2)
  plot_lab_remove <- unique(cn_sig_removed[cn_sig_removed$plotLab == "", "both_gene_id"])
  
  results_mat$both_gene_id <- paste0(results_mat$gene_id_1, ":", results_mat$gene_id_2)
  results_mat[results_mat$both_gene_id %in% plot_lab_remove, "plotLab"] <- ""
  results_mat$both_gene_id <- NULL
  
  # gotta flip some things around according to my preferred order gene elements
  results_mat$gene_id_1 <- gsub("chr([0-9])_","chr0\\1_", results_mat$gene_id_1)
  results_mat$gene_id_2 <- gsub("chr([0-9])_","chr0\\1_", results_mat$gene_id_2)
  
  element_order <- c("CDS(SNV)", "CDS(SV)", "CDS(SNV, SV)", "CN", "enhancer(SNV)", "ss(SNV)", "enhancer(SV)", "3primeUTR(SNV)", "5primeUTR(SNV)", "promoter(SNV)", "lincRNA_promoter(SNV)", "lincRNA(SNV)")
  
  results_mat$gene_gr_1 <- paste0(results_mat$gr_id_1, ":", results_mat$gene_id_1)
  results_mat$gene_gr_2 <- paste0(results_mat$gr_id_2, ":", results_mat$gene_id_2)
  
  
  factorOrder <- tibble(results_mat) %>% mutate(gr_id_1 = factor(gr_id_1, levels = rev(element_order))) %>% arrange(gr_id_1, desc(gene_id_1)) %>% dplyr::select(gr_id_1, gene_id_1, gene_gr_1) %>% unique() %>% pull(gene_gr_1)
  
  tmp <- results_mat %>% 
    dplyr::select(gene_gr_1, gene_gr_2, log10p, plotLab) %>% 
    dplyr::filter(!is.na(log10p)) %>%
    unique() %>%
    dplyr::full_join(results_mat %>% dplyr::select(gene_gr_1 = gene_gr_2, gene_gr_2 = gene_gr_1, log10p, plotLab) %>% dplyr::filter(!is.na(log10p)) %>% unique(), by = c("gene_gr_1", "gene_gr_2", "log10p", "plotLab")) %>%
    unique() %>%
    dplyr::group_by(gene_gr_1, gene_gr_2) %>%
    dplyr::summarize_all(na.omit) %>%
    unique() %>%
    dplyr::mutate(gene_gr_1 = factor(gene_gr_1, levels = factorOrder),
                  gene_gr_2 = factor(gene_gr_2, levels = factorOrder)) %>%
    dplyr::arrange(gene_gr_1, gene_gr_2) %>%
    dplyr::mutate(i = as.numeric(gene_gr_2),
                  j = as.numeric(gene_gr_1)) %>%
    dplyr::filter(j >= i) %>%
    dplyr::select(-i, -j)
  
  tmp <- tmp %>%
    tidyr::separate(gene_gr_1, into = c("gr_id_1", "gene_id_1"), sep = ":") %>% 
    dplyr::mutate(gr_id_1 = factor(gr_id_1, levels = unique(gr_id_1))) %>%
    dplyr::mutate(gene_id_1 = factor(gene_id_1, levels = unique(gene_id_1))) %>%
    tidyr::separate(gene_gr_2, into = c("gr_id_2", "gene_id_2"), sep = ":") %>% 
    dplyr::mutate(gr_id_2 = factor(gr_id_2, levels = unique(gr_id_2))) %>%
    dplyr::mutate(gene_id_2 = factor(gene_id_2, levels = unique(gene_id_2))) %>%
    filter(!is.na(log10p))
  
  all_sig <- tmp[tmp$plotLab != "", ]
  
  p_all <- ggplot(tmp, aes(forcats::fct_rev(gene_id_1), gene_id_2, fill = log10p, label = plotLab)) +
    geom_tile() +
    geom_text(size = 3, vjust = 0.8) +
    scale_fill_gradient2(
      low = "darkred",
      mid = "white",
      high = "darkblue",
      midpoint = 0,
      na.value = "grey",
      guide = "colourbar") +
    xlab("") +
    ylab("") +
    facet_grid(forcats::fct_rev(gr_id_2) ~ forcats::fct_rev(gr_id_1), scales = "free", space = "free", switch = "both") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          axis.text.y = element_text(size = ),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text.x.bottom = element_text(size = 7),
          strip.text.y.left = element_text(size = 7),
          strip.background = element_rect(fill="white"))
  
  pdf(paste0(output_path, subtype, "_mutual_exclusivity.pdf"), width = 15, height = 13)
  print(p_all)
  dev.off()
  
}  