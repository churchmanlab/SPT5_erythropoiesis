#Plot output of deeptools matrix counts for HBB gene (+200bp downstream of gene). Processed from unique reads, stt NETseq files
#Author: Hope Merens 

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

##########FUNCTIONS#############

processing_norm <- function(data_a,SF_a,the_day){   
  #Step 1: normalize data_a file with SF_a size factor for one day of experimentation (ie 1 netseq library)
  #Step 2: format data for plotting with additional data for other days of erythropoiesis 
  
  #Normalize data using size factors
  #size factors for normalization
  data_a_norm <- (data_a / SF_a) 
  
  #Prep for plotting with other days...
  data_for_plot <- as_tibble(data_a_norm) %>% mutate(Day=the_day) %>% pivot_longer(cols = contains("V"), names_to = "names", values_to = "Values") %>% mutate(names = str_replace(names, "V", "")) %>% mutate_at(vars(matches("names")), as.numeric)
  
  return_list <- list(data_for_plot)
}


processing_norm_merge <- function(data_a,data_b,SF_a,SF_b,the_day){   
  #Step 1: normalize data_a/b file with SF_a/b size factor for one day of experimentation (ie 1 netseq library)
  #Step 2: Merge 2 donor replicates (note, taking sum not average)
  #Step 3: format data for plotting with additional data for other days of erythropoiesis 
  
  #Normalize data using size factors
  #size factors for normalization
  data_a_norm <- (data_a / SF_a) 
  data_b_norm <- (data_b / SF_b) 
  
  #Merge replicates - get SUM*** (rather than average)
  data_reps_combined <- data_a_norm + data_b_norm 
  
  #Prep for plotting with other days...
  data_for_plot <- as_tibble(data_reps_combined) %>% mutate(Day=the_day) %>% pivot_longer(cols = contains("V"), names_to = "names", values_to = "Values") %>% mutate(names = str_replace(names, "V", "")) %>% mutate_at(vars(matches("names")), as.numeric)
  
  return_list <- list(data_for_plot)
}

########LOADING DATA##############
setwd("pwd")

#Open deeptools computed matrices

D02a_deeptools_computedmatrix1 <-read.delim("./D02a_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D02b_deeptools_computedmatrix1 <-read.delim("./D02b_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D03a_deeptools_computedmatrix1 <-read.delim("./D03a_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D03b_deeptools_computedmatrix1 <-read.delim("./D03b_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D04a_deeptools_computedmatrix1 <-read.delim("./D04a_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D04b_deeptools_computedmatrix1 <-read.delim("./D04b_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D05a_deeptools_computedmatrix1 <-read.delim("./D05a_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D05b_deeptools_computedmatrix1 <-read.delim("./D05b_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D06a_deeptools_computedmatrix1 <-read.delim("./D06a_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D06b_deeptools_computedmatrix1 <-read.delim("./D06b_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D07a_deeptools_computedmatrix1 <-read.delim("./D07a_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D07b_deeptools_computedmatrix1 <-read.delim("./D07b_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D08a_deeptools_computedmatrix1 <-read.delim("./D08a_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D08b_deeptools_computedmatrix1 <-read.delim("./D08b_merged_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D09a_deeptools_computedmatrix1 <-read.delim("./D09a_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D09b_deeptools_computedmatrix1 <-read.delim("./D09b_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D10a_deeptools_computedmatrix1 <-read.delim("./D10a_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D10b_deeptools_computedmatrix1 <-read.delim("./D10b_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

D11a_deeptools_computedmatrix1 <-read.delim("./D11a_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)
D11b_deeptools_computedmatrix1 <-read.delim("./D11b_noRTbias_noPCRdup_noSI_uniq_stt_neg_processed.bedGraph_HBB_stt_200bpDS.tab", header=F, skip=3)

#skip first 3 lines = text
#dimensions = rows are # of genes, columns are positions in genomic region of interest 

#read in size factors
#Note: size factors were derived from TPM calculations (2/2023):

D0a <- 0.008196223      
D2a <- 0.011581701         
D3a <- 0.044621475        
D4a <- 0.010817784         
D5a <- 0.032742367         
D6a <- 0.071769686         
D7a <- 0.054880842  
D8a <- 0.031091986         
D9a <- 0.020591689       
D10a <- 0.023297973        
D11a <- 0.013536743         
D0b <- 0.046996504         
D2b <- 0.024913224         
D3b <- 0.007720441 
D4b <- 0.006235214         
D5b <- 0.026214829          
D6b <- 0.438169822         
D7b <- 0.098162653         
D8b <- 0.042042060         
D9b <- 0.094084031        
D10b <- 0.040291215 
D11b <- 0.090759200

#########PROCESSING DATA###################

#Compile input information 

#size factors
SFa <- c(D2a,D3a,D4a,D5a,D6a,D7a,D8a,D9a,D10a,D11a)
SFb <- c(D2b,D3b,D4b,D5b,D6b,D7b,D8b,D9b,D10b,D11b)

#list of input matrices
DXa <- list(D02a_deeptools_computedmatrix1,D03a_deeptools_computedmatrix1,D04a_deeptools_computedmatrix1,D05a_deeptools_computedmatrix1,D06a_deeptools_computedmatrix1,D07a_deeptools_computedmatrix1,D08a_deeptools_computedmatrix1,D09a_deeptools_computedmatrix1,D10a_deeptools_computedmatrix1,D11a_deeptools_computedmatrix1)
DXb <- list(D02b_deeptools_computedmatrix1,D03b_deeptools_computedmatrix1,D04b_deeptools_computedmatrix1,D05b_deeptools_computedmatrix1,D06b_deeptools_computedmatrix1,D07b_deeptools_computedmatrix1,D08b_deeptools_computedmatrix1,D09b_deeptools_computedmatrix1,D10b_deeptools_computedmatrix1,D11b_deeptools_computedmatrix1)

#List of Days
days <- c('Day 2', 'Day 3','Day 4','Day 5','Day 6','Day 7','Day 8','Day 9', 'Day 10','Day 11')


######################################################################
#Run processing script - normalize by size factors & merge replicates

signal_byday_list <- list()
for (i in 1:(length(DXa))){
  print(i)
  sig_byday <- processing_norm_merge(DXa[[i]],DXb[[i]],SFa[i],SFb[i],days[i]) #normalize and merge biological replicates for each day of erythropoiesis 
  signal_byday_list[[i]] <- sig_byday
}

#Add all days into one dataframe for plotting:

#First day analyzed
DF_x_vs_signal_0 <- as_tibble(data.frame(signal_byday_list[[1]][[1]]))
colnames(DF_x_vs_signal_0) <- c("Day","names", "Values")
#Other days
for (i in 1:(length(signal_byday_list)-1)){
  nextday <- data.frame(signal_byday_list[[i+1]][[1]]) #pick out next day from list of days 
  colnames(nextday) <- c("Day","names", "Values")
  DF_x_vs_signal_0 <- DF_x_vs_signal_0 %>% rbind(nextday)
}

DF_x_vs_signal_0$sorter = factor(DF_x_vs_signal_0$Day, levels=c('Day 2', 'Day 3','Day 4','Day 5','Day 6','Day 7','Day 8','Day 9', 'Day 10','Day 11'))
######################################################################
#Run processing script - normalize & process Donor B only 

signal_byday_list_donorB <- list()
for (i in 1:(length(DXb))){
  print(i)
  sig_byday_donorB <- processing_norm(DXb[[i]],SFb[i],days[i])
  signal_byday_list_donorB[[i]] <- sig_byday_donorB
}

#Add all days into one dataframe for plotting

DF_x_vs_signal_0_donorB<- as_tibble(data.frame(signal_byday_list_donorB[[1]][[1]]))

colnames(DF_x_vs_signal_0_donorB) <- c("Day","names", "Values")

#Other days
for (i in 1:(length(signal_byday_list_donorB)-1)){
  nextday_donorB<- data.frame(signal_byday_list_donorB[[i+1]][[1]])
  colnames(nextday_donorB) <- c("Day","names", "Values")
  DF_x_vs_signal_0_donorB<- DF_x_vs_signal_0_donorB%>% rbind(nextday_donorB)
}

DF_x_vs_signal_0_donorB$sorter = factor(DF_x_vs_signal_0_donorB$Day, levels=c('Day 2', 'Day 3','Day 4','Day 5','Day 6','Day 7','Day 8','Day 9', 'Day 10','Day 11'))

######################################################################
#Run processing script - normalize & process Donor A only 

signal_byday_list_donorA <- list()
for (i in 1:(length(DXa))){
  print(i)
  sig_byday_donorA <- processing_norm(DXa[[i]],SFa[i],days[i])
  signal_byday_list_donorA[[i]] <- sig_byday_donorA
}

#Add all days into one dataframe for plotting

DF_x_vs_signal_0_donorA<- as_tibble(data.frame(signal_byday_list_donorA[[1]][[1]]))

colnames(DF_x_vs_signal_0_donorA) <- c("Day","names", "Values")

#Other days
for (i in 1:(length(signal_byday_list_donorA)-1)){
  nextday_donorA<- data.frame(signal_byday_list_donorA[[i+1]][[1]])
  colnames(nextday_donorA) <- c("Day","names", "Values")
  DF_x_vs_signal_0_donorA<- DF_x_vs_signal_0_donorA%>% rbind(nextday_donorA)
}

DF_x_vs_signal_0_donorA$sorter = factor(DF_x_vs_signal_0_donorA$Day, levels=c('Day 2', 'Day 3','Day 4','Day 5','Day 6','Day 7','Day 8','Day 9', 'Day 10','Day 11'))

########BINNING & PLOTTING##############

bin_plot <- function(input_all,binsize,ymax,title){
  #Input_all is normalized/merged matrices for region of interest, all days
  #binsize = binsize for barplot (~histogram)
  #ymax for barplot
  #title for barplot
  
  all_bin1 <- input_all %>% filter(names <= 1850) %>% mutate(binned1='empty')
  #Note - for binning, adjust cutoff based on downstream region plotted & binning #s (ie 200bp downstream of HBB in matrix - binning by 5 (ie ideal to be divisible by 5) - round down 7 bp to 1850)
  
  #setting up binning
  i_front <- 1
  i_back <- binsize 
  print('the binsize is:')
  print(binsize)
  
  while (i_back <= dim(all_bin1)[1]){
    Mean_value <- mean(all_bin1$Values[i_front:i_back]) #calculate mean value for bin
    all_bin1$binned1[i_back] <- Mean_value
    all_bin1$names[i_back] <- all_bin1$names[i_back] / binsize #adjusting numbering to bin #s - note, bin #s, NOT bp
    #move to next binning position
    i_front <- i_front + binsize
    i_back <- i_back + binsize
  }
    
    all_bin2 <- all_bin1 %>% filter(binned1 != 'empty')
    all_bin2$binned1 <- as.numeric(as.character(all_bin2$binned1))

    #now plot with facet-wrap for binned data 
    #View(all_bin2)
    ggplot(all_bin2,aes(x=names,y=binned1,fill=Day)) +
      geom_bar(stat='identity',width=1) +
      facet_grid(rows = vars(sorter)) +
      xlab("HBB gene binned positions") + #-50 to ~200 bp past TTS
      ylab("Normalized Binned Counts per Position") +
      theme_bw() +
      scale_fill_manual(values = c("#2166ac","#2166ac","#2166ac","#2166ac","#2166ac","#2166ac","#2166ac","#2166ac","#2166ac","#2166ac","#2166ac")) +
      theme(axis.text.x = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw() + theme(text = element_text(size=10), plot.background = element_blank() ,
                         panel.grid.major = element_blank() ,
                         panel.grid.minor = element_blank() ,
                         panel.border = element_blank() ,
                         panel.background = element_blank(),
                         legend.position='none') +
      ggtitle(title) +
      coord_cartesian(ylim = c(0,ymax)) +
      #draws x and y axis line
      theme(axis.line = element_line(color = 'black') + theme(axis.line = element_line(color = 'black')), axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title = element_text(size = 10))  +
      labs(color = "Day")
}

  

######################
bin_plot(DF_x_vs_signal_0,5,1600,'Summed counts, donors a + b')

bin_plot(DF_x_vs_signal_0_donorB,5,400,'Counts, donor b')
bin_plot(DF_x_vs_signal_0_donorA,5,400,'Counts, donor a')
bin_plot(DF_x_vs_signal_0,5,400,'Summed counts, donors a + b')

