library(data.table)
library(dplyr)
library(tidyverse)
library(future.apply)


# Import codon density files
dir <- "~/OneDrive - Stanford/Riboseq_data_deposit/JL014/codon_density/"
files <- list.files(dir, pattern = "*_A.codon")
samples <- str_remove_all(files, "_A.codon")

for (i in 1:length(files)) {
  assign(samples[i], read_tsv(paste0(dir, files[i]), col_names = c("orf", "position", "codon", "counts") ))
  
}



##########Start from here for data table

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

# Create data table with basic ORF properties
nf_dt <- data.table(orf = dfs[[1]][, 1],
                    position = dfs[[1]][, 2],
                    codon = dfs[[1]][, 3]) 
colnames(nf_dt) <- c("orf", "position", "codon")
setkeyv(nf_dt, c("orf"))
orfs <- unique(nf_dt[, 1])
nf_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
nf_dt[, stopdist := position - (length+1)] # stop codon is at stopdist == 0
nf_dt[, ID := as.character(paste(orf, position, sep = "_"))]

# Add residue (codon table built from code at bottom of this script)
codon_table <- as.data.table(read.csv("~/OneDrive - Stanford/Current_analysis/Celegans_Riboseq/Celegans_reference_genome/codon_table_KS.csv"))
i <- cbind(match(nf_dt$codon, codon_table$codon))
nf_dt <- cbind(nf_dt, residue = codon_table[i]$residue)

# Add raw counts
nf_dt <- as.data.frame(nf_dt)
for(i in 1:length(samples)) {
  print(samples[i])
  nf_dt <- full_join(nf_dt, get(samples[i]), by=c("orf", "position", "codon"))
  setnames(nf_dt, "counts", samples[i])
}
nf_dt <- as.data.table(nf_dt)


##Average based on the sequencing depth. 
sum_Young <- sum(Young_1$counts)+sum(Young_2$counts)
sum_Adult <- sum(Adult_1$counts)+sum(Adult_2$counts)
sum_Old <- sum(Old_1$counts)+sum(Old_2$counts)

nf_dt <- nf_dt[, Young:=round((sum_Young/2)*(Young_1/sum(Young_1)+Young_2/sum(Young_2))/2, digits=0)]
nf_dt <- nf_dt[, Adult:=round((sum_Adult/2)*(Adult_1/sum(Adult_1)+Adult_2/sum(Adult_2))/2, digits=0)]
nf_dt <- nf_dt[, Old:=round((sum_Old/2)*(Old_1/sum(Old_1)+Old_2/sum(Old_2))/2, digits=0)]

# Calculate rpc for each orf after excluding first and last 5 sense codons (Rachel Green excludes 100nt in the eIF5A paper)
# Sum reads in order to exclude genes with less than 64 reads in each replicate in downstream analysis
samples <- c(samples, "Young", "Adult", "Old")
temp <- nf_dt[position > 20 & stopdist < -20]
for(i in samples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, coverage1 := length(which(.SD[, ..i] != 0)) / (length - 40), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  gene <- cbind(match(nf_dt$orf, temp$orf))
  nf_dt <- cbind(nf_dt, sum1 = temp[gene]$sum1)
  nf_dt <- cbind(nf_dt, coverage1 = temp[gene]$coverage1)
  nf_dt <- cbind(nf_dt, rpc1 = temp[gene]$rpc1)
  ID1 <- cbind(match(nf_dt$ID, temp$ID))
  setnames(nf_dt, "sum1", paste0(i,"_sum"))
  setnames(nf_dt, "coverage1", paste0(i,"_coverage"))
  setnames(nf_dt, "rpc1", paste0(i,"_rpc"))
}

# Calculate rpm, pause score, and z-score using median absolute deviation
for(i in samples) {
  print(i)
  nf_dt[, paste0(i,"_rpm") := (nf_dt[[i]] / sum(nf_dt[[i]])) * 10^6] # calculate rpm
  nf_dt[, paste0(i,"_pause") := nf_dt[[i]] / nf_dt[[paste0(i,"_rpc")]]] # calculate pause
}

saveRDS(nf_dt, "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis////nf_dt.rds")

nf_dt <- readRDS("~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis///nf_dt.rds")


### Calculate p-value using Fisher's exact test
# Can add filter for genes that have Pearson's correlation > 0.5 based on rpm across transcript
nf_fishers <- nf_dt[Adult_rpc >= 0.5 & 
                      Young_rpc >=0.5 &
                      Old_rpc >= 0.5 &
                      Adult_sum >= 64 & 
                      Old_sum >= 64 & 
                      Young_sum >=64
                      ]

length(unique(nf_fishers$orf))

orf_to_remove <- data.frame(orf=c(unique(nf_fishers[(Young_sum - Young) < 0]$orf),
                                  unique(nf_fishers[(Adult_sum - Adult) < 0]$orf),
                                  unique(nf_fishers[(Old_sum - Old) < 0]$orf)
                                  )
                            )

nf_fishers <- nf_fishers[!orf %in% orf_to_remove$orf]

for (i in 1:nrow(nf_fishers)) {
  print(i)
  counts1 <- matrix(c(nf_fishers[i]$Young, nf_fishers[i]$Adult, 
                      (nf_fishers[i]$Young_sum - nf_fishers[i]$Young), 
                      (nf_fishers[i]$Adult_sum - nf_fishers[i]$Adult)), nrow = 2)
  
  counts2 <- matrix(c(nf_fishers[i]$Young, nf_fishers[i]$Old, 
                      (nf_fishers[i]$Young_sum - nf_fishers[i]$Young), 
                      (nf_fishers[i]$Old_sum - nf_fishers[i]$Old)), nrow = 2)
  
  counts3 <- matrix(c(nf_fishers[i]$Adult, nf_fishers[i]$Old, 
                      (nf_fishers[i]$Adult_sum - nf_fishers[i]$Adult), 
                      (nf_fishers[i]$Old_sum - nf_fishers[i]$Old)), nrow = 2)

  stalling_test1 <- fisher.test(counts1)
  stalling_test2 <- fisher.test(counts2)
  stalling_test3 <- fisher.test(counts3)
  
  nf_fishers[i, pvalueYA := stalling_test1$p.value]
  nf_fishers[i, pvalueYO := stalling_test2$p.value]
  nf_fishers[i, pvalueAO := stalling_test3$p.value]
  
}

nf_fishers[, padjYA := p.adjust(pvalueYA, method = "BH"), by = orf]
nf_fishers[, padjYO := p.adjust(pvalueYO, method = "BH"), by = orf]
nf_fishers[, padjAO := p.adjust(pvalueAO, method = "BH"), by = orf]

nf_fishers[, odds_ratioYA := (Adult/Young)/((Adult_sum-Adult)/(Young_sum-Young)) ]
nf_fishers[, odds_ratioYO := (Old/Young)/((Old_sum-Old)/(Young_sum-Young)) ]
nf_fishers[, odds_ratioAO := (Old/Adult)/((Old_sum-Old)/(Adult_sum-Adult)) ]

write_csv(nf_fishers, "~/OneDrive - Stanford/Current_analysis/Aging_Collaboration/Killifish_Aging_analysis//nf_fishers_trim20_averaged.csv")



