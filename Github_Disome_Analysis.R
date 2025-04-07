library(data.table)
library(dplyr)
library(tidyverse)
library(future.apply)
library(RareVariantVis)


#Import Dataset
ce_fishers <- as.data.table(read_csv("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/nf_fishers_trim20_averaged.csv"))

#Set directory
dir <- "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/Plots//"


#Choose strong pause sites in each condition. In this case, pause score > 6. 
ce_fishers_Young_ps6 <- ce_fishers[position >=20 &
                                      stopdist <=-20 &
                                      Young_pause>=6 &
                                      Young>Young_rpc &
                                      Young>20]

ce_fishers_Old_ps6 <- ce_fishers[position >=20 &
                                    stopdist <=-20 &
                                    Old_pause>=6 &
                                    Old>Old_rpc &
                                    Old>20
]


##Function to identify sites with disome peak. 
##Compare average reads per codon in the entire range (+/-40 from paues site) and position of expected disome peak (-12 ~ -8 from pause site ).
Disome_filter_Young <- function(x){
  a <- dplyr::filter(x, from_ps>=-12 & from_ps<=-8)
  avg <- sum(x$Young)/81
  avg2 <- sum(a$Young)/5
  if (avg2>avg) {
    r <- x
  } else if (avg2<=avg) {
    r <- tibble()
  }
  r
}

##Function to plot normalized ribosome density across the analysis window (+/- 40 codons)
Disome_search_Young <- function(x){
  k <- dplyr::select(x, orf, position )
  res <- list()
  res1 <- list()
  p <- list()
  for (i in 1:nrow(k)){
    p[[i]] <- c((k[[i,2]]-40):(k[[i,2]]+40))
    res[[i]] <- dplyr::filter(ce_fishers, orf %in% k[[i,1]]) %>% dplyr::filter(position %in% p[[i]])
    if (nrow(res[[i]])==81 & sum(res[[i]]$Young)/81>=0.5 & sum(res[[i]]$Adult)/81>=0.5 & sum(res[[i]]$Old)/81>=0.5) {
      res1[[i]] <- res[[i]] %>% 
        dplyr::mutate(from_ps=c(-40:40), Young_norm=res[[i]]$Young/(sum(res[[i]]$Young)/81), Adult_norm=res[[i]]$Adult/(sum(res[[i]]$Adult)/81), Old_norm=res[[i]]$Old/(sum(res[[i]]$Old)/81)) %>% 
        dplyr::select(orf, position, from_ps, Young_pause, Adult_pause, Old_pause, Young, Adult, Old, Young_norm, Adult_norm, Old_norm)
    }
  }
  res2 <- res1[-which(sapply(res1, is.null))]
  res3 <- lapply(res2, Disome_filter_Young)
  res4 <- bind_rows(res3)
}


Disome_filter_Old <- function(x){
  a <- dplyr::filter(x, from_ps>=-12 & from_ps<=-8)
  avg <- sum(x$Old)/81
  avg2 <- sum(a$Old)/5
  if (avg2>avg) {
    r <- x
  } else if (avg2<=avg) {
    r <- tibble()
  }
  r
}

Disome_search_Old <- function(x){
  k <- dplyr::select(x, orf, position )
  res <- list()
  res1 <- list()
  p <- list()
  for (i in 1:nrow(k)){
    p[[i]] <- c((k[[i,2]]-40):(k[[i,2]]+40))
    res[[i]] <- dplyr::filter(ce_fishers, orf %in% k[[i,1]]) %>% dplyr::filter(position %in% p[[i]])
    if (nrow(res[[i]])==81 & sum(res[[i]]$Young)/81>=0.5 & sum(res[[i]]$Adult)/81>=0.5 & sum(res[[i]]$Old)/81>=0.5) {
      res1[[i]] <- res[[i]] %>% 
        dplyr::mutate(from_ps=c(-40:40), Young_norm=res[[i]]$Young/(sum(res[[i]]$Young)/81), Adult_norm=res[[i]]$Adult/(sum(res[[i]]$Adult)/81), Old_norm=res[[i]]$Old/(sum(res[[i]]$Old)/81)) %>% 
        dplyr::select(orf, position, from_ps, Young_pause, Adult_pause, Old_pause, Young, Adult, Old, Young_norm, Adult_norm, Old_norm)
    }
  }
  res2 <- res1[-which(sapply(res1, is.null))]
  res3 <- lapply(res2, Disome_filter_Old)
  res4 <- bind_rows(res3)
}


Young_disome_list <- Disome_search_Young(ce_fishers_Young_ps6)
Old_disome_list <- Disome_search_Old(ce_fishers_Old_ps6)


plot <- ggplot(data = Young_disome_list) + 
  #geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(from_ps, Young_norm, fill = 'Young'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, Old_norm, fill = 'Old'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, Young_norm, color = 'Young'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(aes(from_ps, Old_norm, color = 'Old'), fun.y = "mean", geom = "line", size = 1.25) +
  coord_cartesian(xlim = c(-25, 25))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from Pause site(codons)", color="condition", fill="condition") +
  theme(legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
plot
ggsave(paste0(dir, "Aligned_to_Young_Disome_sites.pdf"), plot, width = 6, height = 4, dpi = 300, useDingbats = F)





###Searching disome sites that are specific for Old condition
Disome_filter_Old_specific <- function(x){
  a <- dplyr::filter(x, from_ps>=-12 & from_ps<=-8)
  Youngavg <- sum(x$Young)/81
  Youngavg2 <- sum(a$Young)/5
  Oldavg <- sum(x$Old)/81
  Oldavg2 <- sum(a$Old)/5
  pause_site <- dplyr::filter(x, from_ps==0)
  disome_site <- dplyr::filter(a, Old_pause==max(Old_pause))
  if (Oldavg2>Oldavg & Youngavg2 <= Youngavg) {
    r <- dplyr::mutate(x, di_mono_ratio=disome_site$Old_pause/Old_pause)
  } else {
    r <- tibble()
  }
  r
}

Disome_search_Old_specific <- function(x){
  k <- dplyr::select(x, orf, position )
  res <- list()
  res1 <- list()
  p <- list()
  for (i in 1:nrow(k)){
    p[[i]] <- c((k[[i,2]]-40):(k[[i,2]]+40))
    res[[i]] <- dplyr::filter(ce_fishers, orf %in% k[[i,1]]) %>% dplyr::filter(position %in% p[[i]])
    if (nrow(res[[i]])==81 & sum(res[[i]]$Young)/81>=0.5 & sum(res[[i]]$Adult)/81>=0.5 & sum(res[[i]]$Old)/81>=0.5) {
      res1[[i]] <- res[[i]] %>% 
        dplyr::mutate(from_ps=c(-40:40), Young_norm=res[[i]]$Young/(sum(res[[i]]$Young)/81), Adult_norm=res[[i]]$Adult/(sum(res[[i]]$Adult)/81), Old_norm=res[[i]]$Old/(sum(res[[i]]$Old)/81)) %>% 
        dplyr::select(orf, position, from_ps, Young_pause, Adult_pause, Old_pause, Young, Adult, Old, Young_norm, Adult_norm, Old_norm)
    }
  }
  res2 <- res1[-which(sapply(res1, is.null))]
  res3 <- lapply(res2, Disome_filter_Old_specific)
  res4 <- bind_rows(res3)
}




Old_specific_disome_list <- Disome_search_Old_specific(ce_fishers_Old_ps6)
Old_specific_disome_list1 <- Old_specific_disome_list %>% dplyr::mutate(ID=str_c(orf, position, sep="_"))
Old_specific_disome_list2 <- Old_specific_disome_list1 %>% dplyr::mutate(ID=str_c(ID, from_ps, sep="_"))

plot <- ggplot(data = Old_specific_disome_list) + 
  #geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(from_ps, Young_norm, fill = 'Young'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, Adult_norm, fill = 'Adult'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, Old_norm, fill = 'Old'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.95)) +
  stat_summary(aes(from_ps, Young_norm, color = 'Young'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(aes(from_ps, Adult_norm, color = 'Adult'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(aes(from_ps, Old_norm, color = 'Old'), fun.y = "mean", geom = "line", size = 1.25) +
  coord_cartesian(xlim = c(-25, 25))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from Pause site(codons)", color="condition", fill="condition") +
  theme(legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
plot
ggsave(paste0(dir, "Aligned_to_Old_specific_Disome.pdf"), plot, width = 8, height = 4, dpi = 300, useDingbats = F)

length(unique(Old_specific_disome_list$orf))





Old_specific_disome1 <- Old_specific_disome_list %>% dplyr::filter(from_ps==0) %>% dplyr::mutate(ID=str_c(orf, position, sep = "_"))
ce_fishers_Old_disome <- ce_fishers %>% dplyr::filter(ID %in% Old_specific_disome1$ID)

ce_fishers_tp <- readRDS(paste0(dir, "ce_fishers_tp.rds"))
ce_fishers_tp_Old_disome <- ce_fishers_tp %>% dplyr::filter(ID %in% Old_specific_disome1$ID)

ce_dt <- readRDS("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/nf_dt.rds")
ce_dt_intern <- ce_dt[position > 21 & stopdist < -21 & residue != "X" ]
proteome_codon <- ce_dt_intern %>% group_by(residue) %>% summarise(count=n()) %>% dplyr::mutate(freq=count/sum(count))
colnames(proteome_codon) <- c("codon", "count", "freq")
AA_bgd <- proteome_codon


Old_disome_logo <- logo_AA_bits(x=ce_fishers_tp_Old_disome$tripeptide, y="wKL_Logo", z = "All")
plot <- ggseqlogo(Old_disome_logo, method= "custom", seq_type="aa", col_scheme=col_sch)+
  ylab('Bits')+
  scale_x_continuous(breaks=c(1:3), labels = c("E", "P", "A"))+
  labs(title="Disome sites(229 motifs)")
plot
ggsave(paste0(dir, "Disome_all_tp.pdf"), plot, width = 6, height = 8, dpi = 300, useDingbats = F)


