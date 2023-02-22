
#Overlap matrix processing + metagene plotting script
#Authors: Kate LaChance, Danya Martell-Smart, and Hope Merens

#Steps in script:
   #Upload overlap matrices (ie overlap between bedGraphs from Netseq data & bed files containing position information where looking for Pol II pausing ex. Tss)
   #Remove rows with low counts after combining pos/neg strand data
   #add replicates together, normalize data by gene (row), & then calculate mean counts at each genomic position
   # Plot the expanded bed files around loci of interest

#Setup:
library(scales)
library(stringi)
library(data.table)
library(digest)
library(broom)
library(devtools)
library(tidyverse)    
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(dplyr)
library(tidyr)
library(ggpubr)
library("reshape2")
library(zoo)
remove(korg, gene.idtype.list, gene.idtype.bods, cpd.simtypes)

overlapmatrix_processing_2reps_GO <- function(data_positive_strand_1,data_positive_strand_2,data_negative_strand_1,data_negative_strand_2,Day_name){   
   #function to process overlap matrices (input data) for plotting (2 replicates) --> data for one day of erythropoiesis processed at a time 
   #data_positive_strand_1 - positive strand data, replicate 1
   #data_negative_strand_1 - negative strand data, replicate 1
   #data_positive_strand_2 - positive strand data, replicate 2
   #data_negative_strand_2 - negative strand data, replicate 2
   #day_name - string w/ name of column to use of gene TPMs, matching day that's being processed by function
  
   #Combine positive and negative strands together 
   data_strands_combined_1 <- rbind(data_positive_strand_1,data_negative_strand_1)
   data_strands_combined_2 <- rbind(data_positive_strand_2,data_negative_strand_2)
   
   #Remove rows with low counts - only retain expressed genes 
   expressed_genes <- read.csv(file = "2023_02_DM_expressed_genes.csv", row.names = 1) %>% rownames_to_column("ensembl_gene_id")
   
   #filter for respective day
   express_day <- expressed_genes %>% dplyr::select(c('ensembl_gene_id',Day_name))
   express_day_filter <- express_day %>% filter(express_day[,2] > 5) # > 5 TPM for day currently being analyzed 
   
   data_strands_combined_filter_1 <- data_strands_combined_1 %>% filter(ensembl_gene_id %in% express_day_filter$ensembl_gene_id)
   
   data_strands_combined_filter_2 <- data_strands_combined_2 %>% filter(ensembl_gene_id %in% express_day_filter$ensembl_gene_id)

   
   #Add replicates together
   data_reps_combined_pre <- rbind(data_strands_combined_filter_1, data_strands_combined_filter_2)
   data_reps_combined <- data_reps_combined_pre[,-c(1:4)]
   
   #get row mean counts to normalize signal by gene
   data_reps_combined$Avg_score = rowMeans(data_reps_combined)
   data_reps_averaged = data_reps_combined/(data_reps_combined$Avg_score) #Divide by average for normalization
   #remove Avg_score column
   data_reps_averaged_final <- data_reps_averaged %>% dplyr::select(-Avg_score)
   
   ##Average over all genes at each genomic position to get metagene counts 
   #Get # of genes
   Num_genes_ALL = dim(data_reps_averaged_final)[1]/2 #number of ALL genes after filtering taking into account 2 replicates
   print(Num_genes_ALL) 
   
   #for metagene: average over genes first - get average signal across all genes for each position around TSS
   MEAN_data_dist=data.frame(colSums(data_reps_averaged_final,na.rm=TRUE)/Num_genes_ALL) #mean of distribution, na.rm = TRUE
   
   return_list <- list(MEAN_data_dist)
}

######################

plot_metagene <- function(mean_signal_per_position,TSS_total_window_length,Day,window_to_plot_left,window_to_plot_right){   
   #Inputs:
   #mean_signal_per_position is list of outputs from overlap_matrix_processing function
   #TSS total window length is the length (bp) of the window on both sides ie 2kb window on both sides = 4000bp total window length
   #Day given as character vector for labeling plots
   #window_to_plot is subset of total_window to plot (bp)
   
   #set x axis => TSS window
   wnu= TSS_total_window_length/2 #upstream window bp
   wnd= TSS_total_window_length/2 #downstream window bp
   x=-wnu:(wnd-1)
   
   #set up dataframe with position, avg signal, and day information and add dataset from each day to the dataframe for plotting
   DF_x_vs_meansignal_0_pre <- as_tibble(data.frame(x,mean_signal_per_position[[1]],Day[[1]]))
   colnames(DF_x_vs_meansignal_0_pre) <- c("position", "coverage", "Day")
   
   for (i in 1:(length(mean_signal_per_position)-1)){
      nextday <- data.frame(x,mean_signal_per_position[[i+1]],Day[[i+1]])
      colnames(nextday) <- c("position", "coverage", "Day")
      DF_x_vs_meansignal_0_pre <- DF_x_vs_meansignal_0_pre %>% rbind(nextday)
   }
   
   #change the window to plot by filtering for the correct positions 
   DF_x_vs_meansignal_0 <- DF_x_vs_meansignal_0_pre %>% filter(position >= -window_to_plot_left & position <= window_to_plot_right)
   
   #View(DF_x_vs_meansignal_0)
   
   pdf("MetagenePlot1.pdf", 10, 6) #save image to PDF
   
   p1 <- ggplot(DF_x_vs_meansignal_0) +
      geom_line(aes(x = position, y = coverage, color = factor(Day),group=factor(Day)),size=1.5,alpha=0.75) +
      geom_line(data=subset(DF_x_vs_meansignal_0, Day == 4),aes(x = position, y = coverage, color = factor(Day),group=factor(Day)),size=1.5,alpha=0.75) + #move day 4 to front of metagene plot
      xlab("Genomic Position") +
      ylab("Mean Counts") +
      theme_bw() +
      scale_colour_manual(values = c("#a50026","#f46d43","#fee090","#abd9e9","#74add1","#313695",""))  +  #RdYlBu
      theme(axis.text.x = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw() + theme(text = element_text(size=15), plot.background = element_blank() ,
                         panel.grid.major = element_blank() ,
                         panel.grid.minor = element_blank() ,
                         panel.border = element_blank() ,
                         panel.background = element_blank() ) +
      #draws x and y axis line  
      theme(axis.line = element_line(color = 'black') + theme(axis.line = element_line(color = 'black')), axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title = element_text(size = 15))  + ggtitle(NULL) +
      labs(color = "Day") 

   plot(p1)
   
   dev.off()
}  

######################
#File Inputs

setwd("pwd_to_overlap_matrices_folder")

##UPLOAD SENSE TXN FILES (OVERLAP MATRICES)
#Day 0 files
D00a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D00a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D00a_F)[4] <- "ensembl_gene_id"

D00b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D00b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D00b_F)[4] <- "ensembl_gene_id"

D00a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D00a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D00a_R)[4] <- "ensembl_gene_id"

D00b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D00b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D00b_R)[4] <- "ensembl_gene_id"

#Day 2 files
D02a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D02a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D02a_F)[4] <- "ensembl_gene_id"

D02b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D02b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D02b_F)[4] <- "ensembl_gene_id"

D02a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D02a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D02a_R)[4] <- "ensembl_gene_id"

D02b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D02b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D02b_R)[4] <- "ensembl_gene_id"


#Day 3 files
D03a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D03a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D03a_F)[4] <- "ensembl_gene_id"

D03b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D03b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D03b_F)[4] <- "ensembl_gene_id"

D03a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D03a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D03a_R)[4] <- "ensembl_gene_id"

D03b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D03b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D03b_R)[4] <- "ensembl_gene_id"


#Day 4 files
D04a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D04a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D04a_F)[4] <- "ensembl_gene_id"

D04b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D04b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D04b_F)[4] <- "ensembl_gene_id"

D04a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D04a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D04a_R)[4] <- "ensembl_gene_id"

D04b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D04b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D04b_R)[4] <- "ensembl_gene_id"


#Day 5 files
D05a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D05a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D05a_F)[4] <- "ensembl_gene_id"

D05b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D05b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D05b_F)[4] <- "ensembl_gene_id"

D05a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D05a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D05a_R)[4] <- "ensembl_gene_id"

D05b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D05b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D05b_R)[4] <- "ensembl_gene_id"


#Day 6 files
D06a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D06a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D06a_F)[4] <- "ensembl_gene_id"

D06b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D06b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D06b_F)[4] <- "ensembl_gene_id"

D06a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D06a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D06a_R)[4] <- "ensembl_gene_id"

D06b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D06b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D06b_R)[4] <- "ensembl_gene_id"


#Day 7 files
D07a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D07a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D07a_F)[4] <- "ensembl_gene_id"

D07b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D07b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D07b_F)[4] <- "ensembl_gene_id"

D07a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D07a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D07a_R)[4] <- "ensembl_gene_id"

D07b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D07b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D07b_R)[4] <- "ensembl_gene_id"


#Day 8 files
D08a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D08a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D08a_F)[4] <- "ensembl_gene_id"

D08b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D08b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D08b_F)[4] <- "ensembl_gene_id"

D08a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D08a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D08a_R)[4] <- "ensembl_gene_id"

D08b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D08b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D08b_R)[4] <- "ensembl_gene_id"


#Day 9 files
D09a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D09a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D09a_F)[4] <- "ensembl_gene_id"

D09b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D09b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D09b_F)[4] <- "ensembl_gene_id"

D09a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D09a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D09a_R)[4] <- "ensembl_gene_id"

D09b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D09b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D09b_R)[4] <- "ensembl_gene_id"

#Day 10 files
D10a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D10a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D10a_F)[4] <- "ensembl_gene_id"

D10b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D10b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D10b_F)[4] <- "ensembl_gene_id"

D10a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D10a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D10a_R)[4] <- "ensembl_gene_id"

D10b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D10b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D10b_R)[4] <- "ensembl_gene_id"

#Day 11 files
D11a_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D11a_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D11a_F)[4] <- "ensembl_gene_id"

D11b_F <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D11b_F_TSS_2kb_pos_area2000_processed_pos.txt", header = F, skip = 1)
colnames(D11b_F)[4] <- "ensembl_gene_id"

D11a_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D11a_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D11a_R)[4] <- "ensembl_gene_id"

D11b_R <- read.delim(file = "./Overlap_Matrices/OverlapMatrix_D11b_R_TSS_2kb_neg_area2000_processed_neg.txt", header = F, skip = 1)
colnames(D11b_R)[4] <- "ensembl_gene_id"

################
#Compile input information and process overlap matrices

#list of input overlap matrices
DXa_F <- list(D00a_F,D02a_F,D04a_F,D06a_F,D08a_F,D10a_F)
DXb_F <- list(D00b_F,D02b_F,D04b_F,D06b_F,D08b_F,D10b_F)
DXa_R <- list(D00a_R,D02a_R,D04a_R,D06a_R,D08a_R,D10a_R)
DXb_R <- list(D00b_R,D02b_R,D04b_R,D06b_R,D08b_R,D10b_R)

Day_names <- c('D0b','D2b','D4b','D6b','D8b','D10b')

#Run processing script - sense transcripts 
mean_signal_byday_list <- list()
for (i in 1:(length(DXa_F))){
   print(i)
   meansig_byday <- overlapmatrix_processing_2reps_GO(DXa_F[[i]],DXb_F[[i]],DXa_R[[i]],DXb_R[[i]],Day_names[i])
   mean_signal_byday_list[[i]] <- meansig_byday
}


#Make metagene plots

#Plot all days - smaller window 
TSS_total_window_length = 4000
Days <- c(0,2,4,6,8,10)
mean_signal <- c(mean_signal_byday_list[[1]],mean_signal_byday_list[[2]],mean_signal_byday_list[[3]],mean_signal_byday_list[[4]],mean_signal_byday_list[[5]],mean_signal_byday_list[[6]])
plot_metagene(mean_signal,TSS_total_window_length,Days,500,1000)


