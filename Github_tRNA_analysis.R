library(tidyverse)
library(data.table)

#Set directory
dir <- "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon_batch2_all//"

#Import tRNA data
trna_counts <- read_tsv("~/OneDrive - Stanford/GEO_Killifish_new/trnaseq//Anticodon_counts_raw.txt")

##Rename columns for easier manipulation
org_columns <- colnames(trna_counts)

new_columns <- str_remove_all(org_columns, "20240507_fish_mimseq1.3.8_ID0.97_mito/lib2068_ft_tb_")
new_columns <- str_remove_all(new_columns, ".unpaired_uniq.bam")

colnames(trna_counts) <- new_columns

##Calculate RPM for each sample
trna_counts1 <- trna_counts %>% 
  dplyr::mutate(OF_R2_rpm=1000000*OF_R2/sum(trna_counts$OF_R2),
                OF_R1_rpm=1000000*OF_R1/sum(trna_counts$OF_R1),
                OM_R2_rpm=1000000*OM_R2/sum(trna_counts$OM_R2),
                OM_R1_rpm=1000000*OM_R1/sum(trna_counts$OM_R1),
                YF_R2_rpm=1000000*YF_R2/sum(trna_counts$YF_R2),
                YF_R1_rpm=1000000*YF_R1/sum(trna_counts$YF_R1),
                YM_R2_rpm=1000000*YM_R2/sum(trna_counts$YM_R2),
                YM_R1_rpm=1000000*YM_R1/sum(trna_counts$YM_R1)
  ) %>% 
  dplyr::mutate(Old_rpm=(OF_R2_rpm+OF_R1_rpm+OM_R2_rpm+OM_R1_rpm)/4,
                Young_rpm=(YF_R2_rpm+YF_R1_rpm+YM_R2_rpm+YM_R1_rpm)/4
  ) %>% 
  dplyr::mutate(Old_frac=Old_rpm/sum(Old_rpm), 
                Young_frac=Young_rpm/sum(Young_rpm))


##Only filter Cytosolic tRNAs and extract anticodon and residue information
trna_counts_cyto <- trna_counts1 %>% dplyr::filter(str_detect(Anticodon, pattern = "GRZ")==T) %>% 
  dplyr::mutate(anticodon_res = Anticodon) %>% 
  dplyr::filter(str_detect(Anticodon, pattern="iMet")==F)
trna_counts_cyto$Anticodon <- str_remove_all(trna_counts_cyto$Anticodon, "Nothobranchius_furzeri_GRZ_tRNA-")

trna_counts_cyto1 <- separate(trna_counts_cyto, Anticodon, into=c("Residue", "anticodon"), sep = "-")

anticodon_level <- c("TCG", "GCA", "CCT", "ACG", "TCT", "CTC", "CGA", "AGA", "TTC", 
                     "TGA", "GTC", "GCC", "TCC", "CCG", "TAG", "GCT", "CAA", "AAG",
                     "AAT", "CCC", "CTT", "CCA", "TTT", "TAC", "TGC", "AAC", "TCA", 
                     "AGC", "GTG", "CAG", "CAC", "CGC", "TGT", "CAT", "CGG", "TAA", 
                     "TGG", "AGG", "GTA", "CTG", "TAT", "TTG", "GAA", "AGT", "CGT", 
                     "GTT", 
                     "ACA", "AAA")

col_sch <- data.frame(anticodon=c("CTC", "TTC","GTC","CTG", "TTG", "GTT",
                                  "TCG", "CCT", "ACG", "TCT", "CCG", "CTT", "TTT", "GTG", 
                                  "GCA", "GCC", "TCC", "TAG", "CAA", "AAG", "AAT", "CCC", "CCA","ACA",
                                  "TAC", "TGC", "AAC", "TCA", "AGC", "CAG", "CAC", "CGC", "CAT", 
                                  "CGG", "TAA", "TGG", "AGG", "GTA", "TAT", "GAA","AAA",
                                  "CGA", "AGA", "TGA", "GCT", "TGT", "AGT", "CGT" 
),
groups=c("Acidic","Acidic","Acidic","Acidic","Acidic","Acidic",
         "Basic", "Basic","Basic","Basic", "Basic","Basic", "Basic","Basic",
         "Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic",
         "Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic",
         "Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic",
         "Polar", "Polar", "Polar", "Polar", "Polar",  "Polar", "Polar"
),
cols=c("lightsalmon2", "lightsalmon2","lightsalmon2", "lightsalmon2","lightsalmon2", "lightsalmon2",
                     "lightskyblue3", "lightskyblue3","lightskyblue3","lightskyblue3", "lightskyblue3","lightskyblue3","lightskyblue3", "lightskyblue3",
                     "goldenrod", "goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod",
                     "goldenrod", "goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod",
                     "goldenrod", "goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod",
                     "palegreen3","palegreen3","palegreen3","palegreen3","palegreen3",
                     "palegreen3","palegreen3")
)

trna_counts_cyto1 <- left_join(trna_counts_cyto1, col_sch, by="anticodon")

trna_counts_cyto1$anticodon_res <- str_remove_all(trna_counts_cyto1$anticodon_res , "Nothobranchius_furzeri_GRZ_tRNA-")
trna_counts_cyto1 <- as.data.table(trna_counts_cyto1)




##Import tRNA charging data
trna_summary <- read_tsv( "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//trna4_CA1_raw.txt")

trna_young <- trna_summary%>% dplyr::filter(age=="young") %>% 
  group_by(anticodon) %>% 
  summarise(young_avgfrac_charged=mean(end_frac), residue=residue) %>% 
  unique()

trna_old <- trna_summary%>% dplyr::filter(age=="old") %>% 
  group_by(anticodon) %>% 
  summarise(old_avgfrac_charged=mean(end_frac), residue=residue) %>% 
  unique()

trna_all <- full_join(trna_young, trna_old, by=c( "residue", "anticodon"))
trna_all <- dplyr::mutate(trna_all, ratio_charged=old_avgfrac_charged/young_avgfrac_charged) %>% dplyr::filter(residue != "eColiLys")


##Combine tRNA counts and charging information together
trna_counts_cyto_charging <- inner_join(trna_all, trna_counts_cyto1, by="anticodon")

trna_counts_cyto_charging$anticodon <- factor(trna_counts_cyto_charging$anticodon, levels=anticodon_level)



##Import RNAseq and Riboseq data to normalize tRNA counts to demand

RNA_ribo <- read_tsv("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/RNA_Ribo_rpkm.txt")

RNA_ribo_simple <- RNA_ribo %>% dplyr::select(orf, Y_rpkm_rna, O_rpkm_rna)


###Import RIboseq data
ce_dt <- readRDS("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/nf_dt.rds")

ce_dt_cds <- ce_dt[position >= 1 & stopdist <= -2 & residue != "X"]

ce_dt_cds_orf_codon <- ce_dt_cds %>% dplyr::select(orf, codon)

##Combine Riboseq and RNAseq data 
ce_dt_cds_orf_codon_rna_rpkm <- left_join(ce_dt_cds_orf_codon, RNA_ribo_simple, by="orf")

##Summarize each codon abundance based on mRNA abundance (RNAseq data)
cds_codon_summary <- ce_dt_cds_orf_codon_rna_rpkm %>% group_by(codon) %>% summarise(Y_codon_abun=sum(Y_rpkm_rna), O_codon_abun=sum(O_rpkm_rna))


##Anticodon-codon conversion table
anticodon_codon <- c(
  
  "AAC"="GTT",
  "AAG"="CTT",
  "AAT"="ATT",
  "ACG"="CGT",
  "AGA"="TCT",
  "AGC"="GCT",
  "AGG"="CCT",
  "AGT"="ACT",
  "CAA"="TTG",
  "CAC"="GTG",
  "CAG"="CTG",
  "CAT"="ATG",
  "CCA"="TGG",
  "CCC"="GGG",
  "CCG"="CGG",
  "CCT"="AGG",
  "CGT"="ACG",
  "CGC"="GCG",
  "CGA"="TCG",
  "CGG"="CCG",
  "CTC"="GAG",
  "CTG"="CAG",
  "CTT"="AAG",
  "GAA"="TTC",
  "AAT"="ATC",
  "GCA"="TGC",
  "GCT"="AGC",
  "GCC"="GGC",
  "GTT"="AAC",
  "GTC"="GAC",
  "GTA"="TAC",
  "GTG"="CAC",
  "TAT"="ATA",
  "TAC"="GTA",
  "TAA"="TTA",
  "TAG"="CTA",
  "TCT"="AGA",
  "TCC"="GGA",
  "TCA"="TGA",
  "TCG"="CGA",
  "TGT"="ACA",
  "TGC"="GCA",
  "TGA"="TCA",
  "TGG"="CCA",
  "TTC"="GAA",
  "TTG"="CAA",
  "TTT"="AAA",
  
  "TTA"="TAA",
  "CTA"="TAG",
  
  "GTT"="AAT",
  
  "GTC"="GAT",
  
  "GTA"="TAT",
  
  "GTG"="CAT",
  
  "GCT"="AGT",
  
  "GCC"="GGT",
  
  "GCA"="TGT",
  
  "ACG"="CGC",
  
  "GAA"="TTT",
  
  "AAC"="GTC",
  
  "AAG"="CTC",
  
  "AGT"="ACC",
  
  "AGC"="GCC",
  
  "AGA"="TCC",
  
  "AGG"="CCC"
)

anticodon_codon <- as.data.frame(as.matrix(anticodon_codon))
anticodon_codon$anticodon <- rownames(anticodon_codon)
anticodon_codon$codon <- anticodon_codon$V1
anticodon_codon <- anticodon_codon %>% dplyr::select(codon, anticodon)
anticodon_codon$anticodon <- str_remove_all(anticodon_codon$anticodon, pattern = ".1")
anticodon_codon <- as.data.table(anticodon_codon)


##Attach anticodon information to codon abundance 
cds_codon_anticodon <- left_join(cds_codon_summary, anticodon_codon, by="codon")


##Summarize codon abundance (tRNA demand)
cds_codon_anticodon1 <- cds_codon_anticodon %>% group_by(anticodon) %>% summarise(Y_anticodon_abun=sum(Y_codon_abun), 
                                                                                  O_anticodon_abun=sum(O_codon_abun))

cds_codon_anticodon2 <- cds_codon_anticodon1 %>% dplyr::mutate(Y_rel_abun=Y_anticodon_abun/sum(Y_anticodon_abun), 
                                                               O_rel_abun=O_anticodon_abun/sum(O_anticodon_abun)
)

##combine tRNA counts, charging, codon demand into single data table
trna_counts_cyto_charging_anticodon_abun <- left_join(trna_counts_cyto_charging, cds_codon_anticodon2, by="anticodon")



##Plot tRNA abundance vs. tRNA demand 
plot <- ggplot(trna_counts_cyto_charging_anticodon_abun, aes(Young_frac, Y_rel_abun, color=cols))+
  scale_color_identity()+
  geom_smooth( aes(Young_frac, Y_rel_abun, color="grey40"), method="lm", alpha=0.2)+
  geom_point(size=4)+
  labs(x="Rel. tRNA abundance", y="Rel. anticodon Abund.", title="Young tRNA vs. anticodon")+
  theme_classic(base_size = 24)+
  annotate("text", x=0.04, y=0.01,size=8, label=str_c("R=", round(cor(trna_counts_cyto_charging_anticodon_abun$Young_frac, trna_counts_cyto_charging_anticodon_abun$Y_rel_abun, use = "pairwise.complete"), digits = 2)))+
  annotate("text", x=0.04, y=0.00, size=6, label=str_c("p=", cor.test(trna_counts_cyto_charging_anticodon_abun$Young_frac, trna_counts_cyto_charging_anticodon_abun$Y_rel_abun, use = "pairwise.complete")$p.value))
plot
saveRDS(plot, "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNA_vs_ExpNorm_anticodon_Young.rds")
ggsave("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNA_vs_ExpNorm_anticodon_Young.pdf", plot, width = 5, height =5 , dpi = 300, useDingbats = F)

plot <- ggplot(trna_counts_cyto_charging_anticodon_abun, aes(Old_frac, O_rel_abun, color=cols))+
  scale_color_identity()+
  geom_point(size=4)+
  geom_smooth( aes(Old_frac, O_rel_abun, color="grey40"), method="lm", alpha=0.2)+
  labs(x="Rel. tRNA abundance", y="Rel. anticodon Abundance", title="Old tRNA vs. anticodon")+
  theme_classic(base_size = 24)+
  annotate("text", x=0.04, y=0.01,size=8, label=str_c("R=", round(cor(trna_counts_cyto_charging_anticodon_abun$Old_frac, trna_counts_cyto_charging_anticodon_abun$O_rel_abun, use = "pairwise.complete"), digits = 2)))+
  annotate("text", x=0.04, y=0.00, size=6, label=str_c("p=", cor.test(trna_counts_cyto_charging_anticodon_abun$Old_frac, trna_counts_cyto_charging_anticodon_abun$O_rel_abun, use = "pairwise.complete")$p.value))
plot
saveRDS(plot,"~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNA_vs_ExpNorm_anticodon_Old.rds" )
ggsave("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNA_vs_ExpNorm_anticodon_Old.pdf", plot, width = 5, height =5 , dpi = 300, useDingbats = F)


##Plot tRNA abundance change between young and old 
plot <- ggplot(trna_counts_cyto_charging_anticodon_abun, aes(Young_frac, Old_frac, color=cols))+
  scale_color_identity()+
  geom_smooth( aes(Young_frac, Old_frac, color="grey40"), method="lm", alpha=0.2)+
  geom_point(size=4)+
  labs(x="Rel. tRNA Abund.(Young)", y="Rel. tRNA Abund.(Old)", title="tRNA Abund. Young vs. Old")+
  theme_classic(base_size = 24)+
  annotate("text", x=0.04, y=0.01,size=8, label=str_c("R=", round(cor(trna_counts_cyto_charging_anticodon_abun$Young_frac, trna_counts_cyto_charging_anticodon_abun$Old_frac, use = "pairwise.complete"), digits = 4)))+
  annotate("text", x=0.04, y=0.00, size=6, label=str_c("p=", cor.test(trna_counts_cyto_charging_anticodon_abun$Young_frac, trna_counts_cyto_charging_anticodon_abun$Old_frac, use = "pairwise.complete")$p.value))
plot
saveRDS(plot, "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNA_abundance_Young_vs_Old.rds" )
ggsave("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNA_abundance_Young_vs_Old.pdf", plot, width = 5, height =5 , dpi = 300, useDingbats = F)







###Analyze pausing vs. tRNA charging

#import riboseq data
ce_fishers <- as.data.table(read_csv("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/nf_fishers_trim20_averaged.csv"))

#Calculate average pausing on each codon
codon_pause <- ce_fishers %>% 
  dplyr::filter(position >= 1 & stopdist <= -2) %>% 
  group_by(codon) %>% 
  summarise(ps_Old=mean(Old_pause), ps_Young=mean(Young_pause)
  )

##attach anticodon info
codon_summary1 <- left_join(codon_pause, anticodon_codon, by="codon")


##Import tRNA charging data
trna_summary <- read_tsv( "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//trna4_CA1_raw.txt")

trna_summary_OFr1 <- trna_summary %>% dplyr::filter(str_detect(exp, "OF_rep1")==T)
trna_summary_OMr1 <- trna_summary %>% dplyr::filter(str_detect(exp, "OM_rep1")==T) %>% dplyr::mutate(replicate="rep2")
trna_summary_OFr2 <- trna_summary %>% dplyr::filter(str_detect(exp, "OF_rep2")==T) %>% dplyr::mutate(replicate="rep3")
trna_summary_OMr2 <- trna_summary %>% dplyr::filter(str_detect(exp, "OM_rep2")==T) %>% dplyr::mutate(replicate="rep4")
trna_summary_YFr1 <- trna_summary %>% dplyr::filter(str_detect(exp, "YF_rep1")==T)
trna_summary_YMr1 <- trna_summary %>% dplyr::filter(str_detect(exp, "YM_rep1")==T) %>% dplyr::mutate(replicate="rep2")
trna_summary_YFr2 <- trna_summary %>% dplyr::filter(str_detect(exp, "YF_rep2")==T) %>% dplyr::mutate(replicate="rep3")
trna_summary_YMr2 <- trna_summary %>% dplyr::filter(str_detect(exp, "YM_rep2")==T) %>% dplyr::mutate(replicate="rep4")
trna_summary_rest <- trna_summary %>% dplyr::filter(str_detect(condition, "young")==T | str_detect(condition, "old")==T)

trna_summary1 <- rbind(trna_summary_YFr1, trna_summary_OFr1, trna_summary_YMr1, trna_summary_OMr1, 
                       trna_summary_YFr2, trna_summary_OFr2, trna_summary_YMr2, trna_summary_OMr2,
                       trna_summary_rest)

trna_summary1_list <- split(trna_summary1, trna_summary1$replicate)

#Calculate fraction charging for each replicate
get_fraction <- function(x){
  yd <- x %>% dplyr::filter(age=="young") %>% dplyr::mutate(young_frac=end_frac) %>% dplyr::select(replicate,residue, anticodon, res_anti, young_frac)
  od <- x %>% dplyr::filter(age=="old") %>% dplyr::mutate(old_frac=end_frac) %>% dplyr::select(replicate,residue, anticodon, res_anti, old_frac)
  d <- full_join(yd, od, by=c("replicate","residue", "anticodon", "res_anti"))
  d1 <- d %>% dplyr::mutate(ratio_charged=old_frac/young_frac)
  d1
}

trna_summary2_list <- lapply(trna_summary1_list, get_fraction)
trna_summary2 <- bind_rows(trna_summary2_list)

#Summarize ratio of tRNA charging between young and old
trna_all <- trna_summary2 %>% 
  group_by(res_anti) %>% 
  summarise(residue=residue, anticodon=anticodon, avg_ratio=mean(ratio_charged), med_ratio=median(ratio_charged)) %>% 
  unique()

trna_all1 <- left_join(trna_all, col_sch, by="anticodon") %>% as.data.table()


##combined pausing data and tRNA charging data
codon_trna <- full_join(trna_all1, codon_summary1, by="anticodon") %>% 
  dplyr::filter(codon !="TAG" & codon != "TGA") %>% 
  dplyr::mutate(ratio_pause=ps_Old/ps_Young)

plot <- ggplot(data=codon_trna, aes(med_ratio, ratio_pause, color=cols))+
  scale_color_identity()+
  geom_point(size=4, alpha=0.8)+
  theme_classic(base_size=24)+
  theme(legend.position = "none")+
  #coord_cartesian(ylim = c(0.9,1.2))+
  annotate("text", x=0.88, y=0.95,size=8, label=str_c("R=", round(cor(codon_trna$med_ratio, codon_trna$ratio_pause, use = "pairwise.complete"), digits = 2)))+
  annotate("text", x=0.88, y=0.9, size=6, label=str_c("p=", cor.test(codon_trna$med_ratio, codon_trna$ratio_pause, use = "pairwise.complete")$p.value))+
  labs(x="Ratio Frac. tRNA Charged(Old/Young)", y="Ratio Mean Pausing (Old/Young)", title="Average Pausing vs. Frac. tRNA charged")
plot <- plot+ geom_smooth(data=codon_trna, aes(med_ratio, ratio_pause),inherit.aes = F, method="lm", alpha=0.15, color="grey40")
plot
saveRDS(plot, "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNAmed_pausing.rds")
ggsave("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis/tRNA_sum_anticodon//Correlation_tRNAmed_pausing.pdf", plot, width = 5, height =5 , dpi = 300, useDingbats = F)














