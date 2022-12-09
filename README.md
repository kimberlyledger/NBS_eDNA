# NBS_eDNA
eDNA data from Northern Bering Sea surveys

this repo consists of: 
- "filtering.Rmd" which starts with the trimmed fastq files for eDNA samples taken from the NBS and sequenced using the MiFish primer set and creates an ASV table for taxonomic analyses (includes .csv and .fasta outputs)
- "insect_analysis.Rmd" this script using the insect R package to classify ASVs. outputs include a table with ASVs and taxonomic assingment (asv_taxonomy_insect.csv) and a table with taxonomic assigment for all samples (taxon_table.csv)
- "prelim_analyses.Rmd" this script combines the taxon table from the insect analysis with sample metadata to take a look at some results of NBS samples
- "blastn_taxonomy_and_analysis.Rmd" this script classifies ASV using blastn. outputs a table with ASVs and taxonomic assingment (asv_taxonomy_blastn.csv)
- "blastn_vs_insect_asv_assignment.Rmd" this script compares the outputs of both ASV assignments
