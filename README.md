# riboseq_analysis
This repository contains codes for analyzing changes in translation during N.furzeri (Killifish) aging.<br/> 
Associated publication: Impaired biogenesis of basic proteins impacts multiple hallmarks of the aging brain DOI: [10.1101/2023.07.20.549210](https://doi.org/10.1101/2023.07.20.549210)

Github_Getting_datatables.R<br/>
This file contains codes to combine processed data and perform pausing analysis. 

Github_Pausing_Analysis.R<br/>
This file contains codes to visualize sites that show age-dependent increase in pausing. 

Github_Disome_Analysis.R<br/>
This file contains codes to identify potential disome forming sites specific in old condition. 

Github_tRNA_analysis.R<br/>
This files contains codes for analyzing tRNA abundance and charging data. Please download processed data files (Anticodon_counts_raw.txt and CCAcounts.txt) from GEO. 

For pre-processing of sequencing data, please refer to https://github.com/kcstein/riboSeq_processing_annotation.
The final output from above GEO repository (.codon files) are used as input for the codes. 

