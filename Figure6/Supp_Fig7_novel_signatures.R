# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                          plot all novel signatures                  # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
##################################################################### libraries
library(ggplot2)
library(cowplot)
library(dplyr)
###############################################################################
##################################################################### functions

# DBS #
plot.dbs <- function(signatures = DBS_signatures,
                      channel_order = channel_order){
  #prepare data to plot
  sigs_df        <- reshape2::melt(signatures)
  colnames(sigs_df) <- c('channel', 'component', 'value')
  sigs_df$group     <- paste0(matrix(unlist(strsplit(as.character(sigs_df$channel), '>')), ncol = 2, byrow = T)[,1], '>NN')
  sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels = channel_order),
           group = factor(group, levels = unique(sigs_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
                       '#e3211d', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'), unique(sigs_df$group))
  strip_name_colours <- c('black','white','black','white','black','white','black','white','black','white')
                                 
 #plot
  
   sigs_df$real_x_label <- substr(sigs_df$channel, 4, 5)
   
   p <- ggplot(sigs_df, aes(x = real_x_label, y = value, fill = group)) +
     geom_bar(stat = 'identity') +
     facet_grid(component~ group, space = 'free_x', scales = 'free') +
     scale_fill_manual(name = '', values = colours, guide = 'none') +
     xlab('') +
     ylab('% DBS') +
     scale_y_continuous(expand = c(0,0,0.05,0)) +
     # scale_x_discrete(labels = xlabels) +
     theme_classic() + 
     theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                             strip.text = element_text(face = 'bold')) 
   
   #change colours and labels of strips
   g       <- ggplot_gtable(ggplot_build(p))
   striprt <- which(grepl('strip-t', g$layout$name))
   k <- 1
   for (i in striprt) {
     j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
     g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
     g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
     
     t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
     g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
     
     k <- k+1
   }
   
   grid::grid.draw(g)
   return(g)
}


# INDELS #
plot.id <- function(signatures = ID_signatures,
                    channel_order = channel_order){
  ##prepare data to plot
  sigs_df        <- reshape2::melt(signatures)
  colnames(sigs_df) <- c('channel', 'component', 'value')
  sigs_df$group     <- paste(matrix(unlist(strsplit(as.character(sigs_df$channel), ':')), ncol = 4, byrow = T)[,1],
                             matrix(unlist(strsplit(as.character(sigs_df$channel), ':')), ncol = 4, byrow = T)[,2],
                             matrix(unlist(strsplit(as.character(sigs_df$channel), ':')), ncol = 4, byrow = T)[,3], sep = ':')
  
  group_order <- c("1:Del:C", "1:Del:T", "1:Ins:C", "1:Ins:T", 
                   "2:Del:R", "3:Del:R", "4:Del:R", "5:Del:R",
                   "2:Ins:R", "3:Ins:R", "4:Ins:R", "5:Ins:R",
                   "2:Del:M", "3:Del:M", "4:Del:M", "5:Del:M")
  
  sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels = channel_order),
           group = factor(group, levels = group_order),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#fdbe6f', '#ff8002', '#b0dd8b', '#36a12e', '#fdcab5', '#fc8a6a', '#f14432','#bc191a', '#d0e1f2', '#94c4df', '#4a98c9', '#1764ab',
                       '#e1e1ef', '#b6b6d8', '#8683bd','#62409b'), group_order)
  
  legend_labels <- setNames(c('1bp C Deletion', '1bp T Deletion','1bp C Insertion', '1bp T Insertion',
                           '2bp Deletion at Repeats', '3bp Deletion at Repeats', '4bp Deletion at Repeats', '5+bp Deletion at Repeats',
                           '2bp Insertion at Repeats', '3bp Insertion at Repeats', '4bp Insertion at Repeats', '5+bp Insertion at Repeats',
                           '2bp Microhomology Deletion', '3bp Microhomology Deletion', '4bp Microhomology Deletion', '5+bp Microhomology Deletion'),
                           group_order)
  
  strip_name_colours <- rep('black', length(colours))
  strip_name_colours[which(names(colours) %in% c("1:Del:T", "1:Ins:T", "5:Del:R", "5:Ins:R", "5:Del:M"))] <- 'white'
  strip_name <- c('C', 'T', 'C', 'T', '2', '3', '4', '5+', '2', '3', '4', '5+', '2', '3', '4', '5+')
                         
  xlabels <- c('1','2','3','4','5','6+', 
               '1','2','3','4','5','6+',
              '0','1','2','3','4','5+',
              '0','1','2','3','4','5+',
              '1','2','3','4','5','6+', 
              '1','2','3','4','5','6+',
              '1','2','3','4','5','6+', 
              '1','2','3','4','5','6+',
              '0','1','2','3','4','5+',
              '0','1','2','3','4','5+',
              '0','1','2','3','4','5+',
              '0','1','2','3','4','5+',
              '1',
              '1','2',
              '1','2','3',
              '1','2','3','4','5+')
  
   p <- ggplot(sigs_df, aes(x = channel, y = value, fill = group)) +
     geom_bar(stat = 'identity') +
     scale_x_discrete(breaks = unique(as.character(levels(sigs_df$channel))), labels = xlabels) +
     facet_grid(component~group, space = 'free_x', scales = 'free') +
     scale_fill_manual(name = '', values = colours, labels = legend_labels) +
     xlab('Homopolymer Length / Number of Repeat Units / Microhomology Length') +
     ylab('% Indels') +
     scale_y_continuous(expand = c(0,0,0.05,0)) +
     theme_classic() + 
     theme(plot.title = element_text(hjust = 0.5, size = 15),
                             legend.position = 'bottom', strip.text = element_text(face = 'bold'))
     
     #change colours and labels of strips
     g       <- ggplot_gtable(ggplot_build(p))
     striprt <- which(grepl('strip-t', g$layout$name))
     k <- 1
     for (i in striprt) {
       j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
       g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
       g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
       
       t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
       g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
       g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$label  <- strip_name[k]
       
       k <- k+1
     }
     
     grid::grid.draw(g)
     return(g)
}



# SVs #
plot.sv <- function(signatures = SV_signatures,
                    channel_order = channel_order){
  ##prepare data to plot
  sigs_df        <- reshape2::melt(signatures)
  colnames(sigs_df) <- c('channel', 'component', 'value')
  sigs_df$group <- sapply(as.character(sigs_df$channel), function(x){ unlist(strsplit(x, ':'))[1] })
  sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels = channel_order),
           group = factor(group, levels = unique(sigs_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#1f78b4', '#33a02c', '#e31a1c','#ff7f00'), unique(sigs_df$group))
  strip_name_colours <- c(rep('black', 4), rep('white', 4))
  strip_name <- matrix(unlist(strsplit(as.character(unique(sigs_df$group)), '_')), ncol = 2, byrow = T)[,2]
  
  legend_labels <- setNames(c('Clustered Deletion', 'Clustered Tandem Duplication','Clustered Inversion', 'Clustered Translocation',
                              'Non-clustered Deletion', 'Non-clustered Tandem Duplication','Non-clustered Inversion', 'Non-clustered Translocation'),
                            unique(sigs_df$group))
  
  xlabels <- as.character(sapply(unique(as.character(sigs_df$channel)), function(x){ unlist(strsplit(x, ':'))[2] }))
  xlabels <- xlabels[!is.na(xlabels)]
  
  #plot
    
    p <- ggplot(sigs_df, aes(x = channel, y = value, fill = group)) +
      geom_bar(stat = 'identity') +
      facet_grid(component~group, space = 'free_x', scales = 'free') +
      scale_fill_manual(name = '', values = colours[c(1,5,2,6,3,7,4,8)], labels = legend_labels[c(1,5,2,6,3,7,4,8)]) +
      xlab('Rearrangment size') +
      ylab('% SVs') +
      scale_y_continuous(expand = c(0,0,0.05,0)) +
      scale_x_discrete(breaks = grep('trans', unique(as.character(sigs_df$channel)), invert = T, value = T), labels = xlabels) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                              strip.text = element_text(face = 'bold'), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), legend.position = 'bottom')
    
    #change colours and labels of strips
    g       <- ggplot_gtable(ggplot_build(p))
    striprt <- which(grepl('strip-t', g$layout$name))
    k <- 1
    for (i in striprt) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
      
      t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
      g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$label  <- strip_name[k]
      
      k <- k+1
    }
    grid::grid.draw(g)
    return(g)
}


# SCNA #
plot.cn <- function(signatures = CN_signatures,
                    channel_order = channel_order){
  ##prepare data to plot
  sigs_df        <- reshape2::melt(signatures)
  colnames(sigs_df) <- c('channel', 'component', 'value')
  sigs_df$CN        <- sapply(as.character(sigs_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[1] })
  sigs_df$type      <- sapply(as.character(sigs_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[2] })
  sigs_df$length    <- sapply(as.character(sigs_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[3] })
  sigs_df$value     <- sigs_df$value * 100
  sigs_df$group     <- paste0(sigs_df$CN, ":", sigs_df$type)
  sigs_df$group     <- factor(sigs_df$group, levels = unique(sigs_df$group))
  
  sigs_df$type      <- factor(sigs_df$type, levels = c("homdel", "LOH", "het"))
  sigs_df$length    <- as.character(sigs_df$length)

  
  #set colours
  colours       <- setNames(c('#fb8072', '#ffd92f','#66c2a5', '#e78ac3',  '#8da0cb', '#fc8d62', '#1b9e77', '#e7298a',  '#7570b3', '#d95f02'), unique(sigs_df$group))
  legend_labels <- setNames(c('homozygous deletion', 'LOH & CN1', 'LOH & CN2', 'LOH & CN3-4', 'LOH & CN5-8', 'LOH & CN9+',
                              'heterozygous & CN2', 'heterozygous & CN3-4', 'heterozygous & CN5-8', 'heterozygous & CN9+'), unique(sigs_df$group))
  
  strip_name_colours <- rep('black', length(colours))
  strip_name_colours[which(names(colours) %in% c("2:het", "3-4:het", "5-8:het", "9+:het"))] <- 'white'
  strip_name <- c('CN 0', 'CN 1', 'CN 2', 'CN 3-4', 'CN 5-8', 'CN 9+', 'CN 2', 'CN 3-4', 'CN 5-8', 'CN 9+')
    
    xlabels <- c("0-100kb", "100kb-1Mb", ">1Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
                 "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb")
    
    channel_order <- c("0:homdel:0-100kb", "0:homdel:100kb-1Mb", "0:homdel:>1Mb",
                       "1:LOH:0-100kb", "1:LOH:100kb-1Mb", "1:LOH:1Mb-10Mb", "1:LOH:10Mb-40Mb", "1:LOH:>40Mb",
                       "2:LOH:0-100kb", "2:LOH:100kb-1Mb", "2:LOH:1Mb-10Mb", "2:LOH:10Mb-40Mb", "2:LOH:>40Mb",
                       "3-4:LOH:0-100kb", "3-4:LOH:100kb-1Mb", "3-4:LOH:1Mb-10Mb", "3-4:LOH:10Mb-40Mb", "3-4:LOH:>40Mb",
                       "5-8:LOH:0-100kb", "5-8:LOH:100kb-1Mb", "5-8:LOH:1Mb-10Mb", "5-8:LOH:10Mb-40Mb", "5-8:LOH:>40Mb",
                       "9+:LOH:0-100kb", "9+:LOH:100kb-1Mb", "9+:LOH:1Mb-10Mb", "9+:LOH:10Mb-40Mb", "9+:LOH:>40Mb",
                       "2:het:0-100kb", "2:het:100kb-1Mb", "2:het:1Mb-10Mb", "2:het:10Mb-40Mb", "2:het:>40Mb",
                       "3-4:het:0-100kb", "3-4:het:100kb-1Mb", "3-4:het:1Mb-10Mb", "3-4:het:10Mb-40Mb", "3-4:het:>40Mb",
                       "5-8:het:0-100kb", "5-8:het:100kb-1Mb", "5-8:het:1Mb-10Mb", "5-8:het:10Mb-40Mb", "5-8:het:>40Mb",
                       "9+:het:0-100kb", "9+:het:100kb-1Mb", "9+:het:1Mb-10Mb", "9+:het:10Mb-40Mb", "9+:het:>40Mb")
    
    sigs_df$channel   <- factor(sigs_df$channel, levels = channel_order)
    
    
    #plot
      p <- ggplot(sigs_df, aes(x = channel, y = value, fill = group)) +
        geom_bar(stat = 'identity') +
        facet_grid(component~group, space = 'free_x', scales = 'free') +
        scale_fill_manual(name = '', values = colours[c(1,2,3,7,4,8,5,9,6,10)], labels = legend_labels[c(1,2,3,7,4,8,5,9,6,10)]) +
        xlab('Segment Length') +
        ylab('% CN segments') +
        scale_y_continuous(expand = c(0,0,0.05,0)) +
        scale_x_discrete(breaks = channel_order, labels = xlabels) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 15),
              axis.text.x = element_text(angle = 90, vjust = 0.5),
              legend.position = 'bottom', strip.text = element_text(face = 'bold'),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
      
      #change colours and labels of strips
      g       <- ggplot_gtable(ggplot_build(p))
      striprt <- which(grepl('strip-t', g$layout$name))
      k <- 1
      for (i in striprt) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
        
        t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
        g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$label  <- strip_name[k]
        
        k <- k+1
      }
      grid::grid.draw(g)
      return(g)

}

###############################################################################
#################################################################### file paths

signatures_path        <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
output_path            <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"

DBS_input_path         <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/DBS/input_files/input_78matrix_patient.rds"
ID_input_path          <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/ID/input_files/input_83matrix_patient.rds"
CN_input_path          <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/CN/input_files/input_48matrix_patient.rds"
SV_input_path          <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/SV/input_files/input_32matrix_patient.rds"

###############################################################################
#####################################################################     MAIN
load(signatures_path)

# only take the novel signatures
DBS_novel <- DBS_signatures[, c(1, grep("_N", colnames(DBS_signatures)))]
ID_novel <- ID_signatures[, c(1, grep("_N", colnames(ID_signatures)))]
CN_novel <- CN_signatures[, c(1, grep("_N", colnames(CN_signatures)))]
SV_novel <- SV_signatures[, c(1, grep("_N", colnames(SV_signatures)))]

# reorder 
CN_novel <- CN_novel[, c("channel", "CN_N1", "CN_N2", "CN_N6", "CN_N8", "CN_N11", "CN_N12", "CN_N14", "CN_N15")]
SV_novel <- SV_novel[, c("channel", "SV_N1", "SV_N2", "SV_N3", "SV_N4", "SV_N8", "SV_N9", "SV_N10", "SV_N11")]

# get the channels in the right order
DBS_input <- readRDS(DBS_input_path)
DBS_order <- colnames(DBS_input)

ID_input <- readRDS(ID_input_path)
ID_order <- colnames(ID_input)

CN_input <- readRDS(CN_input_path)
CN_order <- colnames(CN_input)

SV_input <- readRDS(SV_input_path)
SV_order <- colnames(SV_input)

# make the plots
p_DBS <- plot.dbs(DBS_novel, DBS_order)
p_ID <- plot.id(ID_novel, ID_order)
p_CN <- plot.cn(CN_novel, CN_order)
p_SV <- plot.sv(SV_novel, SV_order)

p_all <- plot_grid(p_DBS, p_ID, p_CN, p_SV, ncol = 2, rel_widths = c(1, 1), rel_heights = c(4, 7))

pdf(paste0(output_path, "novel_signatures_overview.pdf"), width = 24, height = 16)
p_all
dev.off()
