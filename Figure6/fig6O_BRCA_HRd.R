# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #                 HR signature and mutation analysis                        # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#################################################################################
######################################################################  libraries
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

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

options(stringsAsFactors = F)
options(bitmapType = "cairo")

#################################################################################
######################################################################      paths

sample_table_path  <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
signatures_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
copy_number_path   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/combined_ascatTable_segements_joined.rds"
muttable_base_path <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/2.somaticVariants/C.annotation/"
SV_path            <- "/re_gecip/cancer_lung/kthol/SV/input/SV_list.Rdata"
output_path        <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"
germline_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/germline/germline_drivers.txt"
HRDetect_path      <- "/re_gecip/shared_allGeCIPs/SNZ_for_SwantonGroup/HRDetectScores/HRDetect_Lung_summaryScores.tsv"
mutburden_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
BRCA_path          <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HR_BRCA_mutations.RData"

#################################################################################
######################################################################      MAIN

# # # #
# correlate the HR signatures with each other
load(signatures_path)

sig_df           <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure, MDS_exposure)
hr_sigs          <- grep("HRd", unique(sig_df$label), value = T)

hr_sig_df        <- sig_df[sig_df$label %in% hr_sigs, ]
hr_sig_df$label  <- sub(":.*", "", hr_sig_df$label)
hr_sig_weight    <- reshape2::dcast(hr_sig_df, sample ~ label, value.var = "weight")
colnames(hr_sig_weight) <- sub(":.*", "", colnames(hr_sig_weight))
hr_sig_exposure  <- reshape2::dcast(hr_sig_df, sample ~ label, value.var = "exposure")
colnames(hr_sig_exposure) <- sub(":.*", "", colnames(hr_sig_exposure))

df <- hr_sig_df[hr_sig_df$signature %in% c("SBS3", "ID6"), ]
df <- df[df$weight != 0, ]

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
# 
# # get some gene data
txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes                 <- genes(txdb)
genes                 <- data.frame(genes)

# extract HUGO symbols
gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)

# merge with gene locations
gene_names            <- left_join(genes, gene_symbols)
gene_names_gr         <- GRanges(seqnames = gene_names$seqnames, IRanges(start = gene_names$start, end = gene_names$end))
mcols(gene_names_gr)  <- DataFrame(gene_names[, c("symbol")])

all_samples <- unique(sample_table$participant_id)

# get the mutations in BRCA
HR_genes <- c("BRCA1", "BRCA2")

# # load in and format the structural variants
load(SV_path)

all_SV_concatenated <- all_SV_concatenated[all_SV_concatenated$sample %in% sample_table$participant_id, ]

all_SV_concatenated$chr1 <- sub("chr", "", all_SV_concatenated$chr1)
all_SV_concatenated$chr2 <- sub("chr", "", all_SV_concatenated$chr2)

all_SV_concatenated$pos1 <- as.numeric(as.character(all_SV_concatenated$pos1))
all_SV_concatenated$pos2 <- as.numeric(as.character(all_SV_concatenated$pos2))

# lots of them are duplicated, let's remove

all_SV_concatenated <- all_SV_concatenated %>% rowwise() %>% mutate(ID = paste0(min(chr1, chr2), "_", min(pos1, pos2), "_", max(chr1, chr2), "_", max(pos1, pos2)))
all_SV_concatenated <- data.frame(all_SV_concatenated)

# add sample to the ID column
all_SV_concatenated$ID <- paste0(all_SV_concatenated$ID, "_", all_SV_concatenated$sample)

# kick out duplicates
dups <- duplicated(all_SV_concatenated$ID)

all_SV_concatenated <- all_SV_concatenated[dups == FALSE, ]

# filter them by size a little
all_SV_concatenated$width <- all_SV_concatenated$pos2 - all_SV_concatenated$pos1
all_SV_concatenated$width <- abs(all_SV_concatenated$width)

# remove everything smaller than 50kb
all_SV_concatenated <- all_SV_concatenated[all_SV_concatenated$width >= 50, ]

# get rid of some unnecessary columns
all_SV_concatenated <- all_SV_concatenated[, c("chr1", "pos1", "chr2", "pos2", "sample")]
all_SV_concatenated1 <- all_SV_concatenated[, c("chr1", "pos1", "sample")]
colnames(all_SV_concatenated1) <- c("chr", "pos", "sample")
all_SV_concatenated2 <- all_SV_concatenated[, c("chr2", "pos2", "sample")]
colnames(all_SV_concatenated2) <- c("chr", "pos", "sample")

SV_table_all <- rbind(all_SV_concatenated1, all_SV_concatenated2)
SV_table_all$chr <- paste0("chr", SV_table_all$chr)
# 
# # load the germline variant information
germline_all <- read.table(germline_path, head = T, sep = "\t")

germline_all <- germline_all[germline_all$gene %in% HR_genes, ]

germline_all <- germline_all[, c("patient", "gene", "CHROM", "POS")]
germline_all <- as.data.frame(germline_all)

copy_number_table <- readRDS(copy_number_path)

all_samples <- unique(sample_table$participant_id)

# the sample below keeps failing at the SV gene find overlaps step, I do not know why
all_samples <- as.character(all_samples)

# HR_table <- list()
# for(x in all_samples){
#   
#   print(x)
#   patient  <- x
#   
#   annotated_file  <- paste0(muttable_base_path, patient, '/', patient, '.Genome.SNV.upd3.xls')
#   
#   if(file.exists(annotated_file)){
#     
#     annotation_table <- as.data.frame(data.table::fread(annotated_file, header = T))
#     colnames(annotation_table) <- gsub(patient, "", colnames(annotation_table))
#     colnames(annotation_table) <- gsub("^\\.", "", colnames(annotation_table))
#     
#     annotation_table <- annotation_table[annotation_table$is_Indel == TRUE | annotation_table$is_SNV == TRUE | annotation_table$is_Dinuc == TRUE, ]
#     annotation_table <- annotation_table[annotation_table$Use.For.Plots.Indel | annotation_table$Use.For.Plots == TRUE, ]
#     
#     annotation_table <- annotation_table[annotation_table$Gene.refGene %in% HR_genes, ]
#     
#     mut_table_HR_genes_disr <- annotation_table[which(annotation_table$driverCategory == "1A" |
#                                                         annotation_table$indelDriverCategory == 1 |
#                                                         annotation_table$ExonicFunc.refGene %in% c("frameshift substitution", "stopgain")), ]
#     
#     HR_mutation_genes <- unique(mut_table_HR_genes_disr$Gene.refGene)
#     
#     if(length(HR_mutation_genes) > 0 ){
#       HR_genes_mutated <- paste(HR_mutation_genes, collapse = '_')
#     } else { HR_genes_mutated <- NA}
#     
#     if(length(HR_mutation_genes) > 0 ){
#       HR_gene_mutation <- TRUE
#     } else { HR_gene_mutation <- FALSE}
#     
#     HR_patient_mutation_key <- unique(mut_table_HR_genes_disr$key)
#     
#     # save this information
#     df <- data.frame(patient = patient,
#                      HR_gene_mutation = HR_gene_mutation,
#                      HR_genes_mutated = HR_genes_mutated)
#     
#     # let's check out the LOH status of these mutations
#     if(length(HR_patient_mutation_key) > 0){
#       
#       CN_table <- copy_number_table[copy_number_table$patient == patient, ]
#       
#       # find overlaps with genes
#       CN_gr <- GRanges(seqnames = CN_table$chr, IRanges(start = CN_table$startpos, end = CN_table$endpos))
#       mcols(CN_gr) <- DataFrame(CN_table[, c("patient", "cnTotal", "nMajor", "nMinor", "Ploidy", "ACF")])
#       
#       overlaps <- findOverlaps(CN_gr, gene_names_gr)
#       
#       intersect_CN    <- pintersect(CN_gr[queryHits(overlaps)], gene_names_gr[subjectHits(overlaps)])
#       intersect_genes <- pintersect(gene_names_gr[subjectHits(overlaps)], CN_gr[queryHits(overlaps)], )
#       
#       intersect_CN <- data.frame(intersect_CN)
#       intersect_genes <- data.frame(intersect_genes)
#       
#       CN_genes <- cbind(intersect_CN, intersect_genes)
#       CN_genes <- CN_genes[, c(1, 2, 3, 6, 7, 8, 9, 10, 11, 18)]
#       colnames(CN_genes)[ncol(CN_genes)] <- "gene"
#       
#       CN_genes <- CN_genes[CN_genes$gene %in% HR_mutation_genes, ]
#       
#       #is there LOH here?
#       CN_genes_LOH <- CN_genes[CN_genes$nMinor == 0, ]
#       
#       LOH_mutation_genes <- unique(CN_genes_LOH$gene)
#       
#       if(length(LOH_mutation_genes) > 0 ){
#         HR_genes_LOH_mutated <- paste(LOH_mutation_genes, collapse = '_')
#       } else { HR_genes_LOH_mutated <- NA}
#       
#       if(length(LOH_mutation_genes) > 0 ){
#         HR_gene_LOH_mutation <- TRUE
#       } else { HR_gene_LOH_mutation <- FALSE}
#       
#       df$HR_genes_LOH_mutated <- HR_genes_LOH_mutated
#       df$HR_gene_LOH_mutation <- HR_gene_LOH_mutation
#       
#     } else {
#       df$HR_genes_LOH_mutated <- NA
#       df$HR_gene_LOH_mutation <- FALSE
#     }
#   }
#   
#   # Does this patient have any structural variant breakends in HR genes?
#   
#   SV_table <- SV_table_all[SV_table_all$sample == patient, ]
#   SV_table$pos <- as.character(SV_table$pos)
#   
#   if(nrow(SV_table) > 0){
#     
#     SV_gr <- GRanges(seqnames = SV_table$chr, IRanges(start = as.numeric(SV_table$pos), end = as.numeric(SV_table$pos)))
#     
#     overlaps <- findOverlaps(SV_gr, gene_names_gr)
#     
#     intersect_SV    <- pintersect(SV_gr[queryHits(overlaps)], gene_names_gr[subjectHits(overlaps)])
#     intersect_genes <- pintersect(gene_names_gr[subjectHits(overlaps)], SV_gr[queryHits(overlaps)], )
#     
#     intersect_SV <- data.frame(intersect_SV)
#     intersect_genes <- data.frame(intersect_genes)
#     
#     SV_genes <- cbind(intersect_SV, intersect_genes)
#     SV_genes <- SV_genes[, c(1, 2, 3, 12)]
#     colnames(SV_genes)[ncol(SV_genes)] <- "gene"
#     
#     SV_genes <- SV_genes[SV_genes$gene %in% HR_genes, ]
#     
#     HR_SV_genes <- unique(SV_genes$gene)
#     
#     if(length(HR_SV_genes) > 0 ){
#       HR_genes_SV <- paste(HR_SV_genes, collapse = '_')
#     } else { HR_genes_SV <- NA}
#     
#     if(length(HR_SV_genes) > 0 ){
#       HR_gene_SV_present <- TRUE
#     } else { HR_gene_SV_present <- FALSE}
#     
#     df$HR_genes_SV <- HR_genes_SV
#     df$HR_gene_SV_present <- HR_gene_SV_present
#     
#   } else {
#     df$HR_genes_SV <- NA
#     df$HR_gene_SV_present <- FALSE
#   }
#   
#   # Does this patient have a germline mutation in an HR gene?
#   germline_table <- germline_all[germline_all$patient == patient, ]
#   
#   if(nrow(germline_table) > 0){
#     HR_germline_mut <- TRUE
#   } else { HR_germline_mut <- FALSE }
#   
#   if(nrow(germline_table) > 0){
#     HR_germline_gene <- paste(c(unique(germline_table$gene)), collapse = '_')
#   } else { HR_germline_gene <- NA }
#   
#   df$HR_germline_mut <- HR_germline_mut
#   df$HR_germline_gene <- HR_germline_gene
#   
#   
#   # does this patient have LOH at the germline mutation
#   
#   if(nrow(germline_table) > 0){
#     CN_table <- copy_number_table[copy_number_table$patient == patient, ]
#     
#     # find overlaps with genes
#     CN_gr <- GRanges(seqnames = CN_table$chr, IRanges(start = CN_table$startpos, end = CN_table$endpos))
#     mcols(CN_gr) <- DataFrame(CN_table[, c("patient", "cnTotal", "nMajor", "nMinor", "Ploidy", "ACF")])
#     
#     germline_gr <- GRanges(seqnames = germline_table$CHROM, IRanges(start = germline_table$POS, end = germline_table$POS))
#     mcols(germline_gr) <- DataFrame(germline_table[, c("patient", "gene")])
#     
#     overlaps <- findOverlaps(CN_gr, germline_gr)
#     
#     intersect_CN    <- pintersect(CN_gr[queryHits(overlaps)], germline_gr[subjectHits(overlaps)])
#     intersect_genes <- pintersect(germline_gr[subjectHits(overlaps)], CN_gr[queryHits(overlaps)], )
#     
#     intersect_CN <- data.frame(intersect_CN)
#     intersect_genes <- data.frame(intersect_genes)
#     
#     germline_CN_genes <- cbind(intersect_CN, intersect_genes)
#     germline_CN_genes <- germline_CN_genes[, c(1, 2, 3, 6, 7, 8, 9, 10, 11, 19)]
#     colnames(germline_CN_genes)[ncol(germline_CN_genes)] <- "gene"
#     
#     #is there LOH here?
#     germline_CN_genes <- germline_CN_genes[germline_CN_genes$nMinor == 0, ]
#     
#     LOH_germline_genes <- unique(germline_CN_genes$gene)
#     
#     if(length(LOH_germline_genes) > 0){
#       HR_germline_mut_LOH <- TRUE
#     } else { HR_germline_mut_LOH <- FALSE }
#     
#     if(length(LOH_germline_genes) > 0){
#       LOH_germline_genes <- paste(c(unique(LOH_germline_genes)), collapse = '_')
#     } else { LOH_germline_genes <- NA }
#     
#     df$HR_germline_mut_LOH <- HR_germline_mut_LOH
#     df$LOH_germline_genes <- LOH_germline_genes
#     
#   } else{
#     df$HR_germline_mut_LOH <- FALSE
#     df$LOH_germline_genes <- NA
#   }
#   
#   HR_table[[x]] <- df
#   
# }
# 
# HR_table <- do.call(rbind, HR_table)
# save(HR_table, file = BRCA_path)

HR_table <- get(load(BRCA_path))

# # identify second hit mutations
HR_table[which(HR_table$HR_germline_gene == "BRCA1" & HR_table$LOH_germline_gene == "BRCA1" |
                 HR_table$HR_germline_gene == "BRCA1" & HR_table$HR_genes_SV == "BRCA1" |
                 HR_table$HR_germline_gene == "BRCA1" & HR_table$HR_genes_mutated == "BRCA1" |
                 HR_table$HR_germline_gene == "BRCA1" & HR_table$HR_genes_LOH_mutated== "BRCA1"), "second_hit"] <- TRUE

HR_table[which(HR_table$HR_germline_gene == "BRCA2" & HR_table$LOH_germline_gene == "BRCA2" |
                 HR_table$HR_germline_gene == "BRCA2" & HR_table$HR_genes_SV == "BRCA2" |
                 HR_table$HR_germline_gene == "BRCA2" & HR_table$HR_genes_mutated == "BRCA2" |
                 HR_table$HR_germline_gene == "BRCA2" & HR_table$HR_genes_LOH_mutated== "BRCA2"), "second_hit"] <- TRUE

HR_table[is.na(HR_table$second_hit), "second_hit"] <- FALSE


# add the HR signature weights
ID6_sig <- sig_df[sig_df$signature == "ID6", c("sample", "label", "weight", "exposure")]
HR_table <- left_join(HR_table, ID6_sig, by = c("patient" = "sample"))
HR_table[is.na(HR_table$weight), "weight"] <- 0
HR_table[is.na(HR_table$exposure), "exposure"] <- 0

# split those tumours that have the signature into low and high by upper quantile
mea <- as.numeric(quantile(HR_table[HR_table$weight > 0, "weight"])[4])

HR_table[HR_table$weight < mea & HR_table$weight > 0, "confidence"] <- "low_confidence"
HR_table[HR_table$weight >= mea, "confidence"] <- "high_confidence"
HR_table[HR_table$weight == 0, "confidence"] <- "no_HRd_sig"

HR_table$label <- NULL
colnames(HR_table) <- c("patient", "HR_gene_mutation", "HR_genes_mutated", "HR_genes_LOH_mutated", "HR_gene_LOH_mutation", "HR_genes_SV", "HR_gene_SV_present", "HR_germline_mut", "HR_germline_gene", "HR_germline_mut_LOH", "LOH_germline_genes", "second_hit",  "ID6_weight", "ID6_exposure", "ID6_confidence")

SBS3 <- SBS_exposure[SBS_exposure$signature == "SBS3", ]

HR_table <- left_join(HR_table, SBS3[, c("sample", "weight", "exposure")], by = c("patient" = "sample"))
colnames(HR_table) <- c("patient", "HR_gene_mutation", "HR_genes_mutated", "HR_genes_LOH_mutated", "HR_gene_LOH_mutation", "HR_genes_SV", "HR_gene_SV_present", "HR_germline_mut", "HR_germline_gene", "HR_germline_mut_LOH", "LOH_germline_genes", "second_hit", "ID6_weight", "ID6_exposure", "ID6_confidence", "SBS3_weight", "SBS3_exposure")

# check for SBS signature SBS3

# split those tumours that have the signature into low and high by median
mea <- as.numeric(quantile(HR_table[HR_table$SBS3_weight > 0, "SBS3_weight"])[4])

HR_table[HR_table$SBS3_weight < mea & HR_table$SBS3_weight > 0, "SBS3_confidence"] <- "low_confidence"
HR_table[HR_table$SBS3_weight >= mea, "SBS3_confidence"] <- "high_confidence"
HR_table[HR_table$SBS3_weight == 0, "SBS3_confidence"] <- "no_HRd_sig"

# make a table with mutation burden
murburden <- readRDS(mutburden_path)

# also add the percentage of the genome that has LOH
copy_number_table$seg_length <- copy_number_table$endpos - copy_number_table$startpos
copy_number_table <- copy_number_table %>% group_by(patient) %>%
  mutate(total_length = sum(seg_length))

copy_number_table <- copy_number_table %>% group_by(patient) %>%
  filter(nMinor == 0) %>%
  mutate(loh_length = sum(seg_length))

copy_number_table$loh_percentage <- copy_number_table$loh_length / copy_number_table$total_length

loh_percent <- unique(copy_number_table[, c("patient", "loh_percentage")])

# burden <- left_join(murburden, germline_mutburden)
burden <- murburden
burden$germline_burden <- 0
burden <- left_join(burden, loh_percent, by = c("sample" = "patient"))
burden[is.na(burden$loh_percentage), "loh_percentage"] <- 0
burden$all_mutburden <- burden$SBScount + burden$IDcount + burden$DBScount + burden$SVcount
burden$all_mutburden_germline <- burden$SBScount + burden$IDcount + burden$DBScount + burden$SVcount
burden$MNVburden <- burden$SBScount + burden$IDcount + burden$DBScount

# plot the enrichments of all 3 HRd signatures with germline mutations

enrichment.fun <- function(data = HR_table,
                           mutation_type = c("HR_germline_mut", "second_hit"),
                           signature_confidence = c("ID6_confidence", "SBS3_confidence", "HRdetect_confidence")){
  
  # iterate over each mutation type
  mut_list <- list()
  for(mut in mutation_type){
    print(mut)
    
    sig_list <- list()
    # and within that also iterate over each signature confidence
    for(sig in signature_confidence){
      print(sig)
      
      # make our 4x4 table for high vs low
      dat_high_vs_low <- data.frame(high_samples = c(nrow(data[which(data[, sig] == "high_confidence" & data[, mut] == "TRUE"), ]),
                                                     nrow(data[which(data[, sig] == "high_confidence" & data[, mut] == "FALSE"), ])), 
                                    low_samples = c(nrow(data[which(data[, sig] == "low_confidence" & data[, mut] == "TRUE"), ]),
                                                    nrow(data[which(data[, sig] == "low_confidence" & data[, mut] == "FALSE"), ])))
      rownames(dat_high_vs_low) <- c("has_mutation", "no_mutation")
      
      test_high_vs_low <- fisher.test(dat_high_vs_low)
      
      
      test_high_vs_lowflip <- fisher.test(dat_high_vs_low[, c(2, 1)])
      
      test_high_vs_low_df <- data.frame(test = "high_vs_low",
                                        mutation_type = mut,
                                        signature_confidence = sig,
                                        odds_ratio = test_high_vs_low$estimate,
                                        lowCI = test_high_vs_low$conf.int[1],
                                        hiCI = test_high_vs_low$conf.int[2],
                                        odds_ratio_flip = test_high_vs_lowflip$estimate,
                                        lowCI_flip = test_high_vs_lowflip$conf.int[1],
                                        hiCI_flip = test_high_vs_lowflip$conf.int[2],
                                        p = test_high_vs_low$p.value)
      
      # high vs nothing
      dat_high_vs_nothing <- data.frame(high_samples = c(nrow(data[which(data[, sig] == "high_confidence" & data[, mut] == "TRUE"), ]),
                                                         nrow(data[which(data[, sig] == "high_confidence" & data[, mut] == "FALSE"), ])), 
                                        low_samples = c(nrow(data[which(data[, sig] == "no_HRd_sig" & data[, mut] == "TRUE"), ]),
                                                        nrow(data[which(data[, sig] == "no_HRd_sig" & data[, mut] == "FALSE"), ])))
      rownames(dat_high_vs_nothing) <- c("has_mutation", "no_mutation")
      
      test_high_vs_nothing <- fisher.test(dat_high_vs_nothing)
      test_high_vs_nothingflip <- fisher.test(dat_high_vs_nothing[, c(2, 1)])
      
      test_high_vs_nothing_df <- data.frame(test = "high_vs_noSignature",
                                            mutation_type = mut,
                                            signature_confidence = sig,
                                            odds_ratio = test_high_vs_nothing$estimate,
                                            lowCI = test_high_vs_nothing$conf.int[1],
                                            hiCI = test_high_vs_nothing$conf.int[2],
                                            odds_ratio_flip = test_high_vs_nothingflip$estimate,
                                            lowCI_flip = test_high_vs_nothingflip$conf.int[1],
                                            hiCI_flip = test_high_vs_nothingflip$conf.int[2],
                                            p = test_high_vs_nothing$p.value)
      
      # low vs nothing
      dat_low_vs_nothing <- data.frame(high_samples = c(nrow(data[which(data[, sig] == "low_confidence" & data[, mut] == "TRUE"), ]),
                                                        nrow(data[which(data[, sig] == "low_confidence" & data[, mut] == "FALSE"), ])), 
                                       low_samples = c(nrow(data[which(data[, sig] == "no_HRd_sig" & data[, mut] == "TRUE"), ]),
                                                       nrow(data[which(data[, sig] == "no_HRd_sig" & data[, mut] == "FALSE"), ])))
      rownames(dat_low_vs_nothing) <- c("has_mutation", "no_mutation")
      
      test_low_vs_nothing <- fisher.test(dat_low_vs_nothing)
      test_low_vs_nothingflip <- fisher.test(dat_low_vs_nothing[, c(2, 1)])
      
      test_low_vs_nothing_df <- data.frame(test = "low_vs_noSignature",
                                           mutation_type = mut,
                                           signature_confidence = sig,
                                           odds_ratio = test_low_vs_nothing$estimate,
                                           lowCI = test_low_vs_nothing$conf.int[1],
                                           hiCI = test_low_vs_nothing$conf.int[2],
                                           odds_ratio_flip = test_low_vs_nothingflip$estimate,
                                           lowCI_flip = test_low_vs_nothingflip$conf.int[1],
                                           hiCI_flip = test_low_vs_nothingflip$conf.int[2],
                                           p = test_low_vs_nothing$p.value)
      
      all_test <- rbind(test_high_vs_low_df, test_high_vs_nothing_df, test_low_vs_nothing_df)
      all_test$p_adjust_fisher <- p.adjust(all_test$p, method = "fdr")
      
      all_test[all_test$p_adjust_fisher <= 0.05, "significance"] <- "significant"
      all_test[all_test$p_adjust_fisher > 0.05, "significance"] <- "not_significant"
      
      sig_list[[sig]] <- all_test
      
    }
    sig_list <- do.call(rbind, sig_list)
    mut_list[[mut]] <- sig_list
    
  }
  mut_list <- do.call(rbind, mut_list)
  return(mut_list)
}



enrichment_df <- enrichment.fun(data = HR_table,
                                mutation_type = c("HR_germline_mut", "second_hit"),
                                signature_confidence = c("ID6_confidence", "SBS3_confidence"))

enrichment_df$significance <- factor(enrichment_df$significance, levels = c("not_significant", "significant"))
enrichment_df$mutation_type <- factor(enrichment_df$mutation_type, levels = c("HR_germline_mut", "second_hit"))

enrichment_df[which(enrichment_df$odds_ratio > 1), "odds_ratio_plot"] <- enrichment_df[which(enrichment_df$odds_ratio > 1), "odds_ratio_flip"]
enrichment_df[which(enrichment_df$odds_ratio > 1), "lowCI_plot"] <- enrichment_df[which(enrichment_df$odds_ratio > 1), "lowCI_flip"]
enrichment_df[which(enrichment_df$odds_ratio > 1), "hiCI_plot"] <- enrichment_df[which(enrichment_df$odds_ratio > 1), "hiCI_flip"]

enrichment_df[which(enrichment_df$odds_ratio < 1), "odds_ratio_plot"] <- 0-enrichment_df[which(enrichment_df$odds_ratio < 1), "odds_ratio"]
enrichment_df[which(enrichment_df$odds_ratio < 1), "lowCI_plot"] <- 0-enrichment_df[which(enrichment_df$odds_ratio < 1), "lowCI"]
enrichment_df[which(enrichment_df$odds_ratio < 1), "hiCI_plot"] <- 0-enrichment_df[which(enrichment_df$odds_ratio < 1), "hiCI"]


# plot only high vs no
p_odd_hn <- ggplot(enrichment_df[enrichment_df$test %in% c("high_vs_noSignature"), ], aes(x = mutation_type, y = odds_ratio_plot, colour = signature_confidence, alpha = significance)) + 
  geom_pointrange(aes(ymin = lowCI_plot, ymax = hiCI_plot), position = position_dodge(width = 0.8)) +
  # scale_y_log10() +
  scale_alpha_manual(values = c(1), labels = c("significant"), name = "") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", linewidth = 0.5) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed", color = "grey", linewidth = 0.5) +
  scale_colour_manual(name = "signature", values = c("#1b9e77", "#d95f02", "#7570b3"), labels = c("ID6", "SBS3")) +
  scale_x_discrete(labels = c("pathogenic\ngermline\nmutation", "second\nhit")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0),
        axis.text.y = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.x.top = element_text(size = 6),
        strip.background = element_rect(fill = "white")) +
  xlab("") +
  ylab("odds ratio")

pdf(paste0(output_path, "HRd_signature_mutation_enrichment_high_vs_no.pdf"), width = 4, height = 4)
p_odd_hn
dev.off()

# make a bar plot of the contingency table

contingency.fun <- function(data = HR_table,
                           mutation_type = c("HR_germline_mut", "second_hit"),
                           signature_confidence = c("ID6_confidence", "SBS3_confidence")){
  
  # iterate over each mutation type
  mut_list <- list()
  for(mut in mutation_type){
    print(mut)
    
    sig_list <- list()
    # and within that also iterate over each signature confidence
    for(sig in signature_confidence){
      print(sig)
      
      # make our 4x4 table for high vs low
      dat_high_vs_low <- data.frame(high_samples = c(nrow(data[which(data[, sig] == "high_confidence" & data[, mut] == "TRUE"), ]),
                                                     nrow(data[which(data[, sig] == "high_confidence" & data[, mut] == "FALSE"), ])), 
                                    low_samples = c(nrow(data[which(data[, sig] == "low_confidence" & data[, mut] == "TRUE"), ]),
                                                    nrow(data[which(data[, sig] == "low_confidence" & data[, mut] == "FALSE"), ])),
                                    no_samples = c(nrow(data[which(data[, sig] == "no_HRd_sig" & data[, mut] == "TRUE"), ]),
                                                    nrow(data[which(data[, sig] == "no_HRd_sig" & data[, mut] == "FALSE"), ])))
      rownames(dat_high_vs_low) <- c("has_mutation", "no_mutation")
      dat_high_vs_low$mutation_presence <- rownames(dat_high_vs_low)
      dat_high_vs_low <- melt(dat_high_vs_low, id.vars = "mutation_presence")
      dat_high_vs_low$signature <- sig
      dat_high_vs_low$mutation_type <- mut
      
      
      sig_list[[sig]] <- dat_high_vs_low
      
    }
    sig_list <- do.call(rbind, sig_list)
    mut_list[[mut]] <- sig_list
    
  }
  mut_list <- do.call(rbind, mut_list)
  return(mut_list)
}


contingency_table <- contingency.fun(data = HR_table,
                                     mutation_type = c("HR_germline_mut", "second_hit"),
                                     signature_confidence = c("ID6_confidence", "SBS3_confidence"))

# make this into a bar plot

data_plot <- contingency_table[contingency_table$variable %in% c("high_samples", "no_samples"), ]
data_plot$signature <- sub("_confidence", "", data_plot$signature)

# calculate proportion

data_plot$type <- paste0(data_plot$signature, ":",data_plot$variable, ":", data_plot$mutation_type)
data_plot <- data_plot %>% group_by(type) %>% mutate(group_sum = sum(value))
data_plot$prop <- data_plot$value / data_plot$group_sum

data_plot <- data_plot[data_plot$mutation_type == "second_hit", ]
data_plot$prop <- data_plot$prop * 100
data_plot$mutation_presence <- factor(data_plot$mutation_presence, levels = unique(data_plot$mutation_presence))

p_bar <- ggplot(data_plot, aes(variable, prop, fill = mutation_presence, label = value)) +
          geom_bar(stat = "identity") +
          geom_text(position = position_stack(vjust = 0.5), color = "black") +
          scale_x_discrete(labels = c("high\nsignature", "no\nsignature")) +
          scale_fill_manual(name = "mutation presence", values = c("#225ea8", "#41b6c4"), labels = c("has second hit\nBRCA1/BRCA2 mutation", "no BRCA1/BRCA2\nmutation")) +
          ylab("% tumours") +
          xlab("") +
          facet_wrap(~ signature, scales = "free_y") +
          theme_bw() +
          theme(strip.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))

pdf(paste0(output_path, "HRd_signature_mutation_proportion.pdf"), width = 5, height = 3)
p_bar
dev.off()

