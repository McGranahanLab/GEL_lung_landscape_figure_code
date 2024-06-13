# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #       subset germline variants by cancer germline variants            # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#################################################################################
################################################################### libraries
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(tidyr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ggplot2)

#################################################################################
################################################################### file paths
germline_path     <- "/re_gecip/cancer_lung/shared/germline_clinvar_variants.RDS"
vep_path          <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/germline/vep_output_format.txt"
sample_table_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
rahman_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/data/rahman_germline_gene_list.txt"
output_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/germline/germline_drivers.txt"
output_dir        <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/germline/"
copy_number_path   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/combined_ascatTable_segments_joined.rds"
muttable_base_path <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/2.somaticVariants/C.annotation/"
SV_path            <- "/re_gecip/cancer_lung/kthol/SV/input/SV_list.Rdata"
output_path_plot   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/"

#################################################################################
################################################################### MAIN
germline <- readRDS(germline_path)
sample_table <- read.table(sample_table_path, head = T, sep = "\t")

germline <- germline[germline$patient %in% sample_table$participant_id, ]

# we need to subset these variants to those that are rare in the population
vep <- read.table(vep_path, head = T)
vep <- unique(vep[, c("Uploaded_variation", "Extra")])

#get the allele frequencies across the 1000 genomes
all_split <- list()
for(i in 1:nrow(vep)){
  
  df <- vep[i, ]
  l  <- lapply(df$Extra, function(x){strsplit(x, ";")})[[1]][[1]]
  l2 <- lapply(l, function(x){unlist(strsplit(x, "="))})
  l3 <- data.frame(do.call(rbind, l2))
  l4 <- data.frame(t(l3))
  colnames(l4) <- l4[1,]
  l4 <- l4[-1, ]
  l4$ID <- df$Uploaded_variation
  
  all_split[[i]] <- l4
               
}
all_split <- bind_rows(all_split)

all_split <- all_split[, c(grep("AF", colnames(all_split), value = T), "ID")]
all_split <- all_split[!is.na(all_split$AFR_AF) | !is.na(all_split$AMR_AF) | !is.na(all_split$EAS_AF) | !is.na(all_split$EUR_AF) | !is.na(all_split$SAS_AF), ]
all_split <- unique(all_split)


# also split up the info from the original clinvar table
clinvar <- unique(germline[, c("ID", "INFO")])

clinvar_split <- list()
for(i in 1:nrow(clinvar)){
  
  df <- clinvar[i, ]
  l  <- sapply(df$INFO, function(x){strsplit(x, ";")})[[1]]
  l2 <- sapply(l, function(x){strsplit(x, "=")})
  l3 <- data.frame(do.call(rbind, l2))
  l4 <- data.frame(t(l3))
  colnames(l4) <- l4[1,]
  l4 <- l4[-1, ]
  l4$ID <- df$ID
  
  clinvar_split[[i]] <- l4
  
}
clinvar_split <- bind_rows(clinvar_split)
clinvar_split <- clinvar_split[, c(grep("AF", colnames(clinvar_split), value = T), "ID")]
clinvar_split$GMAF <- NULL
clinvar_split <- clinvar_split[!is.na(clinvar_split$AF1000G), ]
clinvar_split <- unique(clinvar_split)

all_split <- full_join(all_split, clinvar_split)

germline <- left_join(germline, all_split)

# also calculate the allele frequency in GEL
GEL_af <- unique(germline[, c("patient", "ID")])
GEL_af <- data.frame(table(GEL_af$ID))
GEL_af$GEL_AF <- GEL_af$Freq / length(unique(sample_table$participant_id))
colnames(GEL_af) <- c("ID", "Freq", "GEL_AF")

germline <- left_join(germline, GEL_af[, c("ID", "GEL_AF")])

germline <- germline %>% rowwise() %>% mutate(max_Af = max(AFR_AF, AMR_AF, AMR_AF, EUR_AF, SAS_AF, AF1000G, GEL_AF, na.rm = T))

# mark allele fraction smaller than 1%
germline <- data.frame(germline)
germline[which(germline$max_Af <= 0.01), "pop_freq"] <- "rare"
germline[which(germline$max_Af > 0.01), "pop_freq"] <- "ubiquitous"

# now annotate the genes into gain of function or loss of function
rahman <- read.table(rahman_path, head = T, sep = "\t")

germline <- left_join(germline, rahman[, c("Gene..Symbol", "Mechanism.of.action.of.CPG.mutations")], by = c("gene" = "Gene..Symbol"))

# let's apply the TRACERx filters
# https://www.nature.com/articles/s41586-023-05783-5#Sec10
# To identify germline-encoded variants that might act as drivers of cancer development, we analysed a published list of germline predisposition genes12. These were subdivided into those that act through gain-of-function or loss-of-function mutations. Within genes acting through gain-of-function mutations, variants classified as ‘pathogenic’ or ‘likely pathogenic’ by ClinVar (20190305 version) were designated as drivers. Within genes acting through loss-of-function mutations, protein-truncating (stop-gain, frameshift or splice-site) variants (excluding those designated as benign by ClinVar), as well as ClinVar pathogenic or likely pathogenic variants, were designated as drivers. Second hit events were identified in cases with either a somatic mutation affecting the same gene as that containing a germline driver or with a somatic copy number loss affecting the wild-type allele.

germline <- germline %>%
                    mutate(stopgain = grepl('stop_gained',INFO)) %>%
                    mutate(splicesite = grepl('splice',INFO))

germline[germline$Mechanism.of.action.of.CPG.mutations == "gain-of-function" & germline$pathogenic, "isDriver"] <- TRUE
germline[germline$Mechanism.of.action.of.CPG.mutations == "loss-of-function" & germline$pathogenic & germline$frameshift, "isDriver"] <- TRUE
germline[germline$Mechanism.of.action.of.CPG.mutations == "loss-of-function" & germline$pathogenic & germline$splicesite, "isDriver"] <- TRUE
germline[germline$Mechanism.of.action.of.CPG.mutations == "loss-of-function" & germline$pathogenic & germline$stopgain, "isDriver"] <- TRUE
germline[is.na(germline$isDriver), "isDriver"] <- FALSE

# now we need to get rid of some high and low VAF variants
germline <- germline %>% separate(allelic.depths, sep = ",", into = c("REF_count", "VAR_count"))
germline$REF_count <- as.numeric(germline$REF_count)
germline$VAR_count <- as.numeric(germline$VAR_count)
germline$depth <- germline$REF_count + germline$VAR_count
germline$VAF <- germline$VAR_count / germline$depth

germline[germline$VAF < 0.3, "isDriver"] <- FALSE
germline[germline$VAF > 0.7, "isDriver"] <- FALSE

germline$max_Af <- as.numeric(germline$max_Af)

# filter everything out with a max allele frequency greater than 1%
germline_filter <- germline[germline$max_Af < 0.01, ]
germline_filter <- germline_filter[germline_filter$isDriver, ]
germline_filter <- germline_filter[!is.na(germline_filter$ID), ]

write.table(germline_filter, output_path, sep = "\t", col.names = T, row.names = F, quote = F)

# are there second hits in these genes
# # get some gene data
txdb                  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes                 <- genes(txdb, single.strand.genes.only=TRUE)
genes                 <- data.frame(genes)

# extract HUGO symbols
gene_symbols          <- as.data.frame(org.Hs.egSYMBOL)

# merge with gene locations
gene_names            <- left_join(genes, gene_symbols)
gene_names_gr         <- GRanges(seqnames = gene_names$seqnames, IRanges(start = gene_names$start, end = gene_names$end))
mcols(gene_names_gr)  <- DataFrame(gene_names[, c("symbol")])

all_samples <- unique(sample_table$participant_id)

# get the genes of interest

HR_genes <- sample(unique(germline_filter$gene))

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

# get the relevant germline info
germline_all <- germline_filter
germline_all <- germline_all[germline_all$gene %in% HR_genes, ]
germline_all <- germline_all[, c("patient", "gene")]
germline_all <- as.data.frame(germline_all)

copy_number_table <- readRDS(copy_number_path)

all_samples <- unique(sample_table$participant_id)
all_samples <- as.character(all_samples)

second_hit <- list()
for(x in all_samples){
  
  print(x)
  patient  <- x
  
  annotated_file  <- paste0(muttable_base_path, patient, '/', patient, '.Genome.SNV.upd3.xls')
  
  if(file.exists(annotated_file)){
    
    annotation_table <- as.data.frame(data.table::fread(annotated_file, header = T))
    colnames(annotation_table) <- gsub(patient, "", colnames(annotation_table))
    colnames(annotation_table) <- gsub("^\\.", "", colnames(annotation_table))
    
    annotation_table <- annotation_table[annotation_table$is_Indel == TRUE | annotation_table$is_SNV == TRUE | annotation_table$is_Dinuc == TRUE, ]
    annotation_table <- annotation_table[annotation_table$Use.For.Plots.Indel | annotation_table$Use.For.Plots == TRUE, ]
    
    annotation_table <- annotation_table[annotation_table$Gene.refGene %in% HR_genes, ]
    
    mut_table_HR_genes_disr <- annotation_table[which(annotation_table$driverCategory == "1A" |
                                                        annotation_table$indelDriverCategory == 1 |
                                                        annotation_table$ExonicFunc.refGene %in% c("frameshift substitution", "stopgain")), ]
    
    HR_mutation_genes <- unique(mut_table_HR_genes_disr$Gene.refGene)
    
    if(length(HR_mutation_genes) > 0 ){
      HR_genes_mutated <- paste(HR_mutation_genes, collapse = '_')
    } else { HR_genes_mutated <- NA}
    
    if(length(HR_mutation_genes) > 0 ){
      HR_gene_mutation <- TRUE
    } else { HR_gene_mutation <- FALSE}
    
    HR_patient_mutation_key <- unique(mut_table_HR_genes_disr$key)
    
    # save this information
    df <- data.frame(patient = patient,
                     germline_gene_somatic_mutation = HR_gene_mutation,
                     germline_genes_somatically_mutated = HR_genes_mutated)

  }
  
  # Does this patient have any structural variant breakends in germline genes?
  
  SV_table <- SV_table_all[SV_table_all$sample == patient, ]
  SV_table$pos <- as.character(SV_table$pos)
  
  if(nrow(SV_table) > 0){
    
    SV_gr <- GRanges(seqnames = SV_table$chr, IRanges(start = as.numeric(SV_table$pos), end = as.numeric(SV_table$pos)))
    
    overlaps <- findOverlaps(SV_gr, gene_names_gr)
    
    intersect_SV    <- pintersect(SV_gr[queryHits(overlaps)], gene_names_gr[subjectHits(overlaps)])
    intersect_genes <- pintersect(gene_names_gr[subjectHits(overlaps)], SV_gr[queryHits(overlaps)], )
    
    intersect_SV <- data.frame(intersect_SV)
    intersect_genes <- data.frame(intersect_genes)
    
    SV_genes <- cbind(intersect_SV, intersect_genes)
    SV_genes <- SV_genes[, c(1, 2, 3, 12)]
    colnames(SV_genes)[ncol(SV_genes)] <- "gene"
    
    SV_genes <- SV_genes[SV_genes$gene %in% HR_genes, ]
    
    HR_SV_genes <- unique(SV_genes$gene)
    
    if(length(HR_SV_genes) > 0 ){
      HR_genes_SV <- paste(HR_SV_genes, collapse = '_')
    } else { HR_genes_SV <- NA}
    
    if(length(HR_SV_genes) > 0 ){
      HR_gene_SV_present <- TRUE
    } else { HR_gene_SV_present <- FALSE}
    
    df$germline_genes_somatic_SV <- HR_genes_SV
    df$germline_gene_somatic_SV_present <- HR_gene_SV_present
    
  } else {
    df$germline_genes_somatic_SV <- NA
    df$germline_gene_somatic_SV_present <- FALSE
  }
  
  # Does this patient have a germline mutation in an HR gene?
  germline_table <- germline_filter[germline_filter$patient == patient, ]
  
  if(nrow(germline_table) > 0){
    HR_germline_mut <- TRUE
  } else { HR_germline_mut <- FALSE }
  
  if(nrow(germline_table) > 0){
    HR_germline_gene <- paste(c(unique(germline_table$gene)), collapse = '_')
  } else { HR_germline_gene <- NA }
  
  df$germline_driver_present <- HR_germline_mut
  df$germline_driver_gene <- HR_germline_gene
  
  # does this patient have LOH at the germline mutation
  
  if(nrow(germline_table) > 0){
    CN_table <- copy_number_table[copy_number_table$patient == patient, ]
    
    # find overlaps with the germline varinat location
    CN_gr <- GRanges(seqnames = CN_table$chr, IRanges(start = CN_table$startpos, end = CN_table$endpos))
    mcols(CN_gr) <- DataFrame(CN_table[, c("patient", "cnTotal", "nMajor", "nMinor", "Ploidy", "ACF")])
    
    germline_gr <- GRanges(seqnames = germline_table$CHROM, IRanges(start = germline_table$POS, end = germline_table$POS))
    mcols(germline_gr) <- DataFrame(germline_table[, c("patient", "gene")])
    
    overlaps <- findOverlaps(CN_gr, germline_gr)
    
    intersect_CN    <- pintersect(CN_gr[queryHits(overlaps)], germline_gr[subjectHits(overlaps)])
    intersect_genes <- pintersect(germline_gr[subjectHits(overlaps)], CN_gr[queryHits(overlaps)], )
    
    intersect_CN <- data.frame(intersect_CN)
    intersect_genes <- data.frame(intersect_genes)
    
    germline_CN_genes <- cbind(intersect_CN, intersect_genes)
    germline_CN_genes <- germline_CN_genes[, c(1, 2, 3, 6, 7, 8, 9, 10, 11, 19)]
    colnames(germline_CN_genes)[ncol(germline_CN_genes)] <- "gene"
    
    #is there LOH here?
    germline_CN_genes <- germline_CN_genes[germline_CN_genes$nMinor == 0, ]
    
    LOH_germline_genes <- unique(germline_CN_genes$gene)
    
    if(length(LOH_germline_genes) > 0){
      HR_germline_mut_LOH <- TRUE
    } else { HR_germline_mut_LOH <- FALSE }
    
    if(length(LOH_germline_genes) > 0){
      LOH_germline_genes <- paste(c(unique(LOH_germline_genes)), collapse = '_')
    } else { LOH_germline_genes <- NA }
    
    df$germline_driver_LOH_present <- HR_germline_mut_LOH
    df$germline_driver_LOH_genes <- LOH_germline_genes
    
  } else{
    df$germline_driver_LOH_present <- FALSE
    df$germline_driver_LOH_genes <- NA
  }
  
  second_hit[[x]] <- df
  
}

second_hit <- do.call(rbind, second_hit)

write.table(second_hit, paste0(output_dir, "germline_second_hit.txt"), row.names = F, col.names = T, sep = "\t")

# add histology of the tumours
sample_table$participant_id <- as.character(sample_table$participant_id)
second_hit <- left_join(second_hit, sample_table[, c("participant_id", "histology")], by = c("patient" = "participant_id"))

second_hit <- second_hit[second_hit$germline_driver_present, ]

# where one patient hase two genes mutated lets split that up into two rows
germline_gene_split <- list()
for(i in 1:nrow(second_hit)){
  
  df <- second_hit[i, ]
  
  somatic <- strsplit(df$germline_genes_somatically_mutated, "_")[[1]]
  sv <- strsplit(df$germline_genes_somatic_SV, "_")[[1]]
  germline_mut <- strsplit(df$germline_driver_gene, "_")[[1]]
  germline_loh <- strsplit(df$germline_driver_LOH_genes, "_")[[1]]
  
  trues <- c(length(somatic) > 1, length(sv) > 1, length(germline_mut) > 1, length(germline_loh) > 1)
  
  if(any(trues)){
    df_out <- data.frame(patient = df$patient,
                         germline_gene_somatic_mutation = df$germline_gene_somatic_mutation,
                         germline_genes_somatically_mutated = somatic,
                         germline_genes_somatic_SV = sv,
                         germline_gene_somatic_SV_present = df$germline_gene_somatic_SV_present,
                         germline_driver_present = df$germline_driver_present,
                         germline_driver_gene = germline_mut,
                         germline_driver_LOH_present = df$germline_driver_LOH_present,
                         germline_driver_LOH_genes = germline_loh,
                         histology = df$histology)

  } else {df_out <- df}

  germline_gene_split[[i]] <- df_out
}
germline_gene_split <- do.call(rbind, germline_gene_split)

# indicate whether there is actually a second hit
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_driver_LOH_genes), "second_hit"] <- "TRUE"
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_genes_somatically_mutated), "second_hit"] <- "TRUE"
germline_gene_split[which(germline_gene_split$germline_driver_gene == germline_gene_split$germline_genes_somatic_SV), "second_hit"] <- "TRUE"
germline_gene_split[is.na(germline_gene_split$second_hit), "second_hit"] <- "FALSE"

germline_count <- germline_gene_split %>% group_by(histology, germline_driver_gene, second_hit) %>% summarise(n())
colnames(germline_count) <- c("histology", "germline_driver_gene", "second_hit", "count")
germline_count$second_hit <- factor(germline_count$second_hit, levels = c("TRUE", "FALSE"))

germline_count <- germline_count %>% group_by(germline_driver_gene) %>% mutate(gene_n = sum(count))
germline_count <- germline_count[order(germline_count$gene_n, decreasing = T), ]
germline_count$germline_driver_gene <- factor(germline_count$germline_driver_gene, levels = unique(germline_count$germline_driver_gene))

# plot this

p <- ggplot(germline_count, aes(germline_driver_gene, count, fill = histology, alpha = second_hit)) +
          geom_bar(stat = "identity") +
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
                                       "OTHER" = "#7a7979"),
                            labels = c("Adenocarcinoma",
                                       "Carcinoid",
                                       "Large cell",
                                       "Squamous cell\nmetastasis",
                                       "Other",
                                       "Carcinoid",
                                       "Squamous cell")) +
          xlab("germline cancer gene") +
          ylab("# tumours") +
          scale_alpha_manual(values = c(1, 0.3), name = "second hit") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf(paste0(output_path_plot, "germline_driver_genes_second_hit.pdf"), width = 6, height = 4)
p
dev.off()


