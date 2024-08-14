library(data.table)
library(dplyr)
library(tidyverse)
library(future.apply)


nf_fishers <- as.data.table(read_csv("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/nf_fishers_trim20_averaged.csv"))

dir <- "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/Plots//"

################
################
#Aligning reads to age-dependent pause sites

#Age-dependent pause sites that show progressive increase in pausing
nf_fishers_YAO_up <- nf_fishers[padjYO<0.05 & 
                                  odds_ratioYO > 1 &
                                  odds_ratioYO < Inf &
                                  Old_pause>Young_pause &
                                  Old_pause>Adult_pause &
                                  position>=20 &
                                  stopdist<= -21 &
                                  Old > Old_rpc
]

#Function for aligning reads
align_to_ps <- function(x){
  k <- dplyr::select(x, orf, position, stopdist) %>% dplyr::filter(position > 25 & stopdist < -25)
  res <- list()
  res1 <- list()
  p <- list()
  for (i in 1:nrow(k)){
    p[[i]] <- c((k[[i,2]]-40):(k[[i,2]]+40)) #Select positions +/-20 codons from the pause sites
    res[[i]] <- dplyr::filter(nf_fishers, orf %in% k[[i,1]] & position %in% p[[i]]) #select transcript and positions 
    if (nrow(res[[i]])==81 & sum(res[[i]]$Young)/81>=0.5 & sum(res[[i]]$Old)/81>=0.5 & sum(res[[i]]$Adult)/81>=0.5) {
      res1[[i]] <- res[[i]] %>% 
        dplyr::mutate(from_ps=c(-40:40), norm_Young= movingAverage(Young/(sum(Young)/81), n=2, centered = T) , norm_Old= movingAverage(Old/(sum(Old)/81), n=2, centered = T), norm_Adult= movingAverage(Adult/(sum(Adult)/81), n=2, centered = T)  ) %>% #Normalize counts based on local average around each site
        dplyr::select(orf, from_ps, Young_pause, Adult_pause, Old_pause, Young, Adult, Old, norm_Young, norm_Adult, norm_Old)
    }  
  }
  res1
}

YAO_up_aligned <- align_to_ps(nf_fishers_YAO_up)
YAO_up_aligned1 <- bind_rows(YAO_up_aligned)

plot <- ggplot(data = YAO_up_aligned1) + 
  #geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(from_ps, norm_Young, fill = 'Young'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, norm_Adult, fill = 'Adult'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, norm_Old, fill = 'Old'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, norm_Young, color = 'Young'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(aes(from_ps, norm_Adult, color = 'Adult'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(aes(from_ps, norm_Old, color = 'Old'), fun.y = "mean", geom = "line", size = 1.25) +
  coord_cartesian(xlim = c(-25, 25))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from Pause site(codons)", color="condition", fill="condition") +
  theme(legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
plot
ggsave(paste0(dir, "Alignedto_Young_vs_Old_pause_sites_with_Oldpause_greaterthan_Adultpause.pdf"), plot, width = 8, height = 4, dpi = 300, useDingbats = F)




############
############
#Maing LOGO plot around age-dependent pause sites

#Function to get tripeptide of E, P, A site
total_tripeptide_motif <- function(x){
  t <- dplyr::select(x, position, residue)
  t2 <- t %>% dplyr::mutate(position2=position+1) %>% dplyr::select(position2, residue)
  t3 <- t %>% dplyr::mutate(position2=position+2) %>% dplyr::select(position2, residue)
  colnames(t2) <- c("position", "residue2")
  colnames(t3) <- c("position", "residue3")
  t4 <- left_join(x, t2, by="position")
  t5 <- left_join(t4, t3, by="position")
  t8 <- tidyr::unite(t5, "tripeptide", c("residue3","residue2", "residue"), sep="", remove = TRUE)
  dplyr::filter(t8, position!=0 & position!=-1 & position!=max(position) & position!=(max(position)-1) )
}

nf_fishers <- as.data.frame(nf_fishers) %>% as.data.table()
nf_fishers$orf <- as.character(nf_fishers$orf)
nf_fishers_list <- split(nf_fishers, nf_fishers$orf)

nf_fishers_tp_list <- lapply(nf_fishers_list, total_tripeptide_motif)
nf_fishers_tp <- bind_rows(nf_fishers_tp_list)
nf_fishers_tp <- as.data.table(nf_fishers_tp)

nf_fishers_tp_YAO_up <- dplyr::filter(nf_fishers_tp, ID %in% nf_fishers_YAO_up$ID)
nf_fishers_tp_YAO_up_list <- split(nf_fishers_tp_YAO_up, nf_fishers_tp_YAO_up$orf)

nf_fishers_tp_YAO_up_ps6 <- dplyr::filter(nf_fishers_tp_YAO_up, Old_pause>=6)

ce_dt_intern <- ce_dt[position > 21 & stopdist < -21 & residue != "X" ]
proteome_codon <- ce_dt_intern %>% group_by(residue) %>% summarise(count=n()) %>% dplyr::mutate(freq=count/sum(count))
colnames(proteome_codon) <- c("codon", "count", "freq")

AA_bgd <- proteome_codon


###NEED to RUN 'logo_AA_bits' function at the bottom before running this part. 

nf_fishers_YAO_up_logo <- logo_AA_bits(x=nf_fishers_tp_YAO_up$tripeptide, y="wKL_Logo", z = "All")
plot <- ggseqlogo(nf_fishers_YAO_up_logo, method= "custom", seq_type="aa", col_scheme=col_sch)+
  ylab('Bits')+
  scale_x_continuous(breaks=c(1:3), labels = c("E", "P", "A"))+
  labs(title="Young vs. Old with Adult (All 4268 motifs)")
plot
ggsave(paste0(dir, "Young_vs_Old_with_Oldpause_higherthan_Adultpause_All_logo.pdf"), plot, width = 6, height = 8, dpi = 300, useDingbats = F)


nf_fishers_up_YAO_ps6_logo <- logo_AA_bits(x=nf_fishers_tp_YAO_up_ps6$tripeptide, y="wKL_Logo", z = "All")
plot <- ggseqlogo(nf_fishers_up_YAO_ps6_logo, method= "custom", seq_type="aa", col_scheme=col_sch)+
  ylab('Bits')+
  scale_x_continuous(breaks=c(1:3), labels = c("E", "P", "A"))+
  labs(title="Young vs. Old with Adult (ps>6 1129 motifs)")
plot
ggsave(paste0(dir, "Young_vs_Old_with_Oldpause_higherthan_Adultpause_PS6_logo.pdf"), plot, width = 6, height = 8, dpi = 300, useDingbats = F)




#############
#############
#Function for Logo plots
library(ggseqlogo)

col_sch <- make_col_scheme(chars=c("D", "E", 
                                   "K", "R", "H", 
                                   "A", "L", "V", "I", "F", "P", "W", "M", 
                                   "G", "S", "T", "Y", "C", 
                                   "N", "Q"),
                           groups=c("Acidic","Acidic",
                                    "Basic", "Basic","Basic",
                                    "Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic",
                                    "Polar", "Polar", "Polar", "Polar", "Polar", 
                                    "Neutral","Neutral"),
                           cols=c("tomato3", "tomato3",
                                           "dodgerblue4", "dodgerblue4","dodgerblue4",
                                           "gray35", "gray35","gray35","gray35","gray35","gray35","gray35","gray35",
                                           "gray70","gray70","gray70","gray70","gray70",
                                           "gray70","gray70")
)

logo_AA_bits <- function(x, y, z){
  k <- data.frame(sequence=as.character(x))
  k$sequence <- as.character(k$sequence)
  
  #Separate sequences into single characters
  b <- separate(data=k, col = 1, sep =c(1:(str_length(k[1,1])-1)), into = paste0("col", c(1:str_length(k[1,1]))) )
  
  #Get summary of each position
  abc <- list()
  temp <- list()
  temp2 <- list()
  aminoacids <- data.frame(factor=c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'))
  for (i in 1:ncol(b)){
    abc[[i]] <- b %>% dplyr::select(i) %>% dplyr::mutate(count=1) %>% group_by(nt=get(paste0("col",i))) %>% summarise(total=sum(count)) %>% dplyr::filter(nt != "X")
    temp[[i]] <- dplyr::filter(aminoacids, !factor %in% abc[[i]]$nt)
    if (nrow(temp[[i]])>0) {
      temp2[[i]] <- data.frame(nt=temp[[i]]$factor,
                               total=0)
      abc[[i]] <- rbind(abc[[i]], temp2[[i]])
    }
    else {
      abc[[i]] <- abc[[i]]
    }
    
  }
  
  #Calculate frequency at each position
  #bgd <- data.frame(nt=AA_bgd$codon,
  #                  bgd_freq=AA_bgd$freq)
  if (z=="All") {
    bgd <- dplyr::select(AA_bgd, codon, freq)
  }
  if (z=="Cyto"){
    bgd <- dplyr::select(AA_bgd_Cyto, codon, freq)
  }
  if (z=="ER"){
    bgd <- dplyr::select(AA_bgd_ER, codon, freq)
  }
  if (z=="Mito"){
    bgd <- dplyr::select(AA_bgd_Mito, codon, freq)
  }
  
  colnames(bgd) <- c("nt", "bgd_freq")
  
  freq <- function(x){
    k <- dplyr::mutate(x, freq=total/sum(total))
    full_join(k, bgd, by="nt")
  }
  
  abc_freq <- lapply(abc, freq)
  
  #Calculate the height of each residue at each position
  if (y=="Logo") {
    Rvalue <- function(x){
      #Traditional method without background correction
      k <- dplyr::mutate(x, R=log2(20)-(-sum(freq*log2(freq))+((1/logb(2))*((20-1)/(2*nrow(b))))))
      dplyr::mutate(k, height=freq*R)
    }
  }
  
  if (y=="Bgd_Logo") {
    Rvalue <- function(x){
      #Method from Dey KK 2018 (Mini-Review of Sequence Logos)
      k <- dplyr::mutate(x, R=sum(log2((freq/bgd_freq)^freq)))
      dplyr::mutate(k, height=freq*R)
    }
  }
  
  if (y=="wKL_Logo") {
    Rvalue <- function(x){
      k <- dplyr::mutate(x, R=log2((freq/bgd_freq)^freq))
      dplyr::mutate(k, height=R)
    }
  }
  if (y=="wKL_norm_Logo") {
    Rvalue <- function(x){
      k <- dplyr::mutate(x, IC=sum(log2((freq/bgd_freq)^freq)), R=log2((freq/bgd_freq)^freq))
      dplyr::mutate(k, height=R*IC/sum(abs(R)))
    }
  }
  
  if (y=="ED_Logo") {
    Rvalue <- function(x){
      #Method from Saethang T et al. 2019 (PTM-Logo), Douglass J et al. 2012 Am. J Physiol Cell Physiol
      r <- dplyr::mutate(x, R=log2((freq/bgd_freq)))
      k <- dplyr::mutate(r, raw_height=(R-median(R)))
      t <- dplyr::mutate(k, height=raw_height)
    }
  }
  
  
  abc_Rvalue <- lapply(abc_freq, Rvalue)
  
  #Combine the data into single dataframe and reformat for ggseqlogo
  abcd <- dplyr::select(abc_Rvalue[[1]], nt, height)
  colnames(abcd) <- c("nt", "height1") 
  for (i in 1:(ncol(b)-1)){
    abcde <- dplyr::select(abc_Rvalue[[i+1]], nt, height)
    colnames(abcde) <- c("nt", paste0("height", i+1))
    abcd <- full_join(abcd, abcde, by="nt")
  }
  
  abcd <- as.data.frame(abcd)
  rownames(abcd) <- abcd$nt
  abcd <- dplyr::select(abcd, -nt)
  abcd <- as.matrix(abcd)
  
}









