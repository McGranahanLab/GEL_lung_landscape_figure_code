# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                              SV N3 Circos plot                      # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
##################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggplot2)
library(dplyr)
library(reshape2)
library(circlize)
library(bedr)

options(stringsAsFactors = F)
options(bitmapType = "cairo")


###############################################################################
#####################################################################file paths

signatures_path        <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
SV_path                <- "~/re_gecip/cancer_lung/kthol/SV/input/SV_list.Rdata"
sample_table_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
output_path            <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/clustering/ALL_new/"

###############################################################################
#####################################################################      MAIN

# sample table
sample_table <- read.table(sample_table_path, head = T, sep = "\t")

# load in SVs
load(SV_path)

all_SV_concatenated <- all_SV_concatenated[all_SV_concatenated$sample %in% sample_table$participant_id, ]

all_SV_concatenated$chr1 <- sub("chr", "", all_SV_concatenated$chr1)
all_SV_concatenated$chr2 <- sub("chr", "", all_SV_concatenated$chr2)

all_SV_concatenated[all_SV_concatenated$chr1 == "X", "chr1"] <- 23
all_SV_concatenated[all_SV_concatenated$chr2 == "X", "chr2"] <- 23
all_SV_concatenated[all_SV_concatenated$chr1 == "Y", "chr1"] <- 24
all_SV_concatenated[all_SV_concatenated$chr2 == "Y", "chr2"] <- 24

all_SV_concatenated$chr1 <- as.numeric(all_SV_concatenated$chr1)
all_SV_concatenated$chr2 <- as.numeric(all_SV_concatenated$chr2)

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
all_SV_concatenated <- all_SV_concatenated[, c("chr1", "pos1", "chr2", "pos2", "classification","sample")]
all_SV_concatenated <- unique(all_SV_concatenated)

# format this for circos plot
# pretend the unpaired ones are paired with itself

all_SV_concatenated[is.na(all_SV_concatenated$chr2), "chr2"] <- all_SV_concatenated[is.na(all_SV_concatenated$chr2 == "NA"), "chr1"] 
all_SV_concatenated[is.na(all_SV_concatenated$pos2 == "NA"), "pos2"] <- all_SV_concatenated[is.na(all_SV_concatenated$pos2 == "NA"), "pos1"] 

SV_paired <- all_SV_concatenated
SV_paired <- SV_paired[which(SV_paired$chr1 != "random"), ]
SV_paired <- SV_paired[which(SV_paired$chr2 != "random"), ]

SV_paired <- SV_paired[which(SV_paired$chr1 != "alt"), ]
SV_paired <- SV_paired[which(SV_paired$chr2 != "alt"), ]

SV_paired <- SV_paired[which(SV_paired$chr1 != "chrUn"), ]
SV_paired <- SV_paired[which(SV_paired$chr2 != "chrUn"), ]

SV_paired$pos1 <- as.numeric(SV_paired$pos1)
SV_paired$pos2 <- as.numeric(SV_paired$pos2)

SV_paired <- SV_paired[which(SV_paired$pos1 != "NA"), ]
SV_paired <- SV_paired[which(SV_paired$pos2 != "NA"), ]

SV_paired <- SV_paired[grep("alt", SV_paired$chr1, invert = T), ]
SV_paired <- SV_paired[grep("random", SV_paired$chr1, invert = T), ]
SV_paired <- SV_paired[grep("Un", SV_paired$chr1, invert = T), ]
SV_paired <- SV_paired[grep("alt", SV_paired$chr2, invert = T), ]
SV_paired <- SV_paired[grep("random", SV_paired$chr2, invert = T), ]
SV_paired <- SV_paired[grep("Un", SV_paired$chr2, invert = T), ]
SV_paired[SV_paired$chr1 == "23", "chr1"] <- "X"
SV_paired[SV_paired$chr2 == "23", "chr2"] <- "X"
SV_paired[SV_paired$chr1 == "24", "chr1"] <- "Y"
SV_paired[SV_paired$chr2 == "24", "chr2"] <- "Y"

# load the SV signatures to find some examples to use
load(signatures_path)

# example samples
# 225000027
# 221004638

for(s in c("225000027", "221004638")){
  
  ww <- which(c("225000027", "221004638") == s)
  
  SVs <- SV_paired[SV_paired$sample == s, ]
  
  # need to separate SVs
  nuc1 <- SVs[, c("chr1", "pos1")]
  nuc1$end <- nuc1$pos1
  colnames(nuc1)<- c("chr", "start", "stop")
  nuc1$chr <- paste0("chr", nuc1$chr)
  
  nuc2 <- SVs[, c("chr2", "pos2")]
  nuc2$end <- nuc2$pos2
  colnames(nuc2)<- c("chr", "start", "stop")
  nuc2$chr <- paste0("chr", nuc2$chr)
  
  # need to sort out the colors
  
  paired_SV_col <- data.frame(region = SVs$classification)
  paired_SV_col$color <- NULL
  paired_SV_col[paired_SV_col$region== "TRA", "color"] <- "#99d594"
    paired_SV_col[paired_SV_col$region == "DEL", "color"] <- "#d53e4f"
      paired_SV_col[paired_SV_col$region == "DUP", "color"] <- "#fc8d59"
        paired_SV_col[paired_SV_col$region == "h2h", "color"] <- "#fee08b"
          paired_SV_col[paired_SV_col$region == "t2t", "color"] <- "#e6f598"
            
          paired_SV_col <- paired_SV_col$color
          
          # need to get the chromososme length
          chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
          
          chr_length <- get.chr.length(species = "human", build = "hg38")
          chr_length <- chr_length[chr_length$chr %in% chrs, ]
          chr_length$start <- 0
          chr_length <- chr_length[, c("start", "length")]
          chr_length <- as.matrix(chr_length)
          
          # make the pdf
          save_path <- paste0(output_path, ww, "_SV_circos_SVTYPE.pdf")
          pdf(save_path, width = 4, height = 5)
          
          # make plot canvas
          circos.clear()
          
          # some plot parameters
          col_text <- "grey40"
            circos.par("track.height"=0.8, gap.degree=3, cell.padding=c(0, 0, 0, 0))
            
            
            # initialize plot
            circos.initialize(factors = chrs, 
                              xlim = chr_length)
            
            
            # make chromosome track
            
            circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
              chr=CELL_META$sector.index
              xlim=CELL_META$xlim
              ylim=CELL_META$ylim
              circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
                          facing="bending.inside")
            }, bg.col="grey90", bg.border=F, track.height=0.06)
            
            # add the SVs
            circos.genomicLink(nuc1, nuc2, col = paired_SV_col, border = NA)
            
            dev.off()
            
            
}





