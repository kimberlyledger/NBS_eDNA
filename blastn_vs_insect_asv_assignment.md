blastn_vs_insect_asv_assignment
================
Kimberly Ledger
2022-12-05

This script compares the ASV taxonomic assignment for the 2021 NBS
samples using 12S MiFish primers

## load libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

## read in asv tables

these tables were produced in “filtering.Rmd”,
“blastn_taxonomy_and_analysis.Rmd” and “insect_analysis.Rmd” scripts

``` r
#asv seq and asv id 
asv_table <- read.csv("asv_seq_table.csv") %>%
  select(!X) %>%
  rename(ASV = asv_id) 

#asv id and blastn taxonomic assignment 
taxon_blastn <- read.csv("asv_taxonomy_blastn.csv") %>%
  select(!X) %>%
  rename(taxon_blastn = taxon) %>%
  rename(taxonomic_level_blastn = taxonomic_level)

#asv id and insect taxonomic assingment 
taxon_insect <- read.csv("asv_taxonomy_insect.csv") %>%
  select(!X) %>%
  rename(ASV = representative) %>%
  rename(taxon_insect = taxon) %>%
  rename(taxonomic_level_insect = rank)
```

## join taxonomic assignments

``` r
join <- taxon_insect %>%
  left_join(taxon_blastn, by = "ASV")
```

## join the asv sequence

``` r
join_asv <- join %>%
  left_join(asv_table, by = "ASV")
```

## export for use in excel, etc.

``` r
write.csv(join_asv, "asv_taxonomy_comparison.csv", row.names = FALSE)
```

## summarize asv’s

order of taxonomic rankings: *no rank *class *superorder *order
*suborder *infraorder *superfamily *family *subfamily *genus *species
*subspecies

``` r
insect_asv_summary <- join %>%
  group_by(taxonomic_level_insect) %>%
  summarise(count = n())

insect_asv_summary$list <- c(2, 8, 10, 6, 1, 11, 9, 5, 12, 7, 3)

insect_asv_summary <- insect_asv_summary%>%
  arrange(list) %>%
  select(!list) %>%
  rename(taxonomic_level = taxonomic_level_insect) %>%
  rename(count_insect = count)

#insect_asv_summary
```

``` r
blast_asv_summary <- join %>%
  group_by(taxonomic_level_blastn) %>%
  summarise(count = n())

blast_asv_summary$taxonomic_level_blastn[6] <- "no rank"

blast_asv_summary$list <- c(2, 4, 5, 3, 6 , 1)

blast_asv_summary <- blast_asv_summary %>%
  arrange(list) %>%
  select(!list) %>%
  rename(taxonomic_level = taxonomic_level_blastn) %>%
  rename(count_blast = count)

#blast_asv_summary
```

asv table

``` r
asv_summ <- insect_asv_summary %>%
  left_join(blast_asv_summary, by = "taxonomic_level") 
asv_summ
```

    ## # A tibble: 11 × 3
    ##    taxonomic_level count_insect count_blast
    ##    <chr>                  <int>       <int>
    ##  1 no rank                   39          14
    ##  2 class                      2           5
    ##  3 superorder                 2          NA
    ##  4 suborder                  39          NA
    ##  5 infraorder                 4          NA
    ##  6 superfamily               10          NA
    ##  7 family                    35         117
    ##  8 subfamily                  3          NA
    ##  9 genus                     65          78
    ## 10 species                   29          14
    ## 11 subspecies                 1          NA

taxonomy table

``` r
tax_insect_sum <- taxon_insect %>%
  group_by(taxonomic_level_insect) %>%
  summarise(count_unique = n_distinct(taxon_insect))

tax_insect_sum$list <- c(2, 8, 10, 6, 1, 11, 9, 5, 12, 7, 3)

tax_insect_sum <- tax_insect_sum %>%
  arrange(list) %>%
  select(!list) %>%
  rename(taxonomic_level = taxonomic_level_insect) %>%
  rename(count_insect = count_unique)
```

``` r
tax_blast_sum <- taxon_blastn %>%
  group_by(taxonomic_level_blastn) %>%
  summarise(count_unique = n_distinct(taxon_blastn))

tax_blast_sum$list <-  c(2, 4, 5, 3, 6)

tax_blast_sum <- tax_blast_sum %>%
  arrange(list) %>%
  select(!list) %>%
  rename(taxonomic_level = taxonomic_level_blastn) %>%
  rename(count_blast = count_unique)
```

taxon table

``` r
taxon_summ <- tax_insect_sum %>%
  left_join(tax_blast_sum, by = "taxonomic_level")
taxon_summ
```

    ## # A tibble: 11 × 3
    ##    taxonomic_level count_insect count_blast
    ##    <chr>                  <int>       <int>
    ##  1 no rank                    6          NA
    ##  2 class                      1           1
    ##  3 superorder                 1          NA
    ##  4 suborder                   3          NA
    ##  5 infraorder                 3          NA
    ##  6 superfamily                1          NA
    ##  7 family                     6           6
    ##  8 subfamily                  2          NA
    ##  9 genus                     12          13
    ## 10 species                   26           9
    ## 11 subspecies                 1          NA

## spp of interest for insect

add up the number of ASV’s for species of interest

``` r
asv_per_spp <- taxon_insect %>%
  filter(taxonomic_level_insect == "species") %>%
  group_by(taxon_insect) %>%
  summarize(count = n())
asv_per_spp
```

    ## # A tibble: 26 × 2
    ##    taxon_insect           count
    ##    <chr>                  <int>
    ##  1 Ammodytes hexapterus       1
    ##  2 Anoplopoma fimbria         1
    ##  3 Artedius fenestralis       1
    ##  4 Bos taurus                 1
    ##  5 Clupea pallasii            2
    ##  6 Cymatogaster aggregata     1
    ##  7 Eleginus gracilis          1
    ##  8 Gadus chalcogrammus        3
    ##  9 Gadus macrocephalus        1
    ## 10 Gasterosteus aculeatus     1
    ## # … with 16 more rows

``` r
#filter(asv_per_spp, taxon_insect %in% c("Eleginus gracilis", "Gadus chalcogrammus", "Gadus macrocephalus"))
```

in my quick fishbase search of the species list,some species have
records from Pacific near SE AK, but not from the Bering Sea… and
Ictalurus punctatus (Channel catfish) and Micropterus salmoides
(Largemouth black bass) and the non-fish are contamination(?)

add up the number of ASV’s for genera

``` r
asv_per_genus <- taxon_insect %>%
  filter(taxonomic_level_insect == "genus") %>%
  group_by(taxon_insect) %>%
  summarize(count = n())
asv_per_genus
```

    ## # A tibble: 12 × 2
    ##    taxon_insect  count
    ##    <chr>         <int>
    ##  1 Ammodytes         1
    ##  2 Carassius         1
    ##  3 Gadus            35
    ##  4 Hexagrammos       1
    ##  5 Hippoglossus      1
    ##  6 Ictalurus         2
    ##  7 Limanda           1
    ##  8 Lycodes           1
    ##  9 Myoxocephalus     1
    ## 10 Oncorhynchus     18
    ## 11 Pungitius         2
    ## 12 Tridentiger       1

``` r
#filter(asv_per_genus, taxon_insect %in% c("Gadus"))
```

add up the number of ASV’s for families

``` r
asv_per_family <- taxon_insect %>%
  filter(taxonomic_level_insect == "family") %>%
  group_by(taxon_insect) %>%
  summarize(count = n())
asv_per_family
```

    ## # A tibble: 6 × 2
    ##   taxon_insect   count
    ##   <chr>          <int>
    ## 1 Cottidae           2
    ## 2 Cyprinidae         1
    ## 3 Gadidae           21
    ## 4 Ictaluridae        1
    ## 5 Pleuronectidae     9
    ## 6 Stichaeidae        1

``` r
#filter(asv_per_family, taxon_insect %in% c("Gadidae"))
```

## spp of interest for blastn

add up the number of ASV’s for species of interest

``` r
asv_per_spp <- taxon_blastn %>%
  filter(taxonomic_level_blastn == "species") %>%
  group_by(taxon_blastn) %>%
  summarize(count = n())
asv_per_spp
```

    ## # A tibble: 9 × 2
    ##   taxon_blastn            count
    ##   <chr>                   <int>
    ## 1 Anoplopoma_fimbria          2
    ## 2 Artedius_fenestralis        2
    ## 3 Cymatogaster_aggregata      1
    ## 4 Eleginus_gracilis           1
    ## 5 Gasterosteus_aculeatus      1
    ## 6 Homo_sapiens                2
    ## 7 Mallotus_villosus           3
    ## 8 Nautichthys_pribilovius     1
    ## 9 Osmerus_mordax              1

``` r
#filter(asv_per_spp, taxon_blastn %in% c("Eleginus_gracilis"))
```

add up the number of ASV’s for genera

``` r
asv_per_genus <- taxon_blastn %>%
  filter(taxonomic_level_blastn == "genus") %>%
  group_by(taxon_blastn) %>%
  summarize(count = n())
asv_per_genus
```

    ## # A tibble: 13 × 2
    ##    taxon_blastn  count
    ##    <chr>         <int>
    ##  1 Bos               1
    ##  2 Carassius         2
    ##  3 Hemilepidotus     1
    ##  4 Hexagrammos       1
    ##  5 Ictalurus         4
    ##  6 Lycodes           2
    ##  7 Micropterus      10
    ##  8 Oncorhynchus     46
    ##  9 Pomoxis           2
    ## 10 Pungitius         3
    ## 11 Sus               1
    ## 12 Syngnathus        2
    ## 13 Tridentiger       3

``` r
#filter(asv_per_genus, taxon_insect %in% c("Gadus"))
```

add up the number of ASV’s for families

``` r
asv_per_family <- taxon_blastn %>%
  filter(taxonomic_level_blastn == "family") %>%
  group_by(taxon_blastn) %>%
  summarize(count = n())
asv_per_family
```

    ## # A tibble: 6 × 2
    ##   taxon_blastn   count
    ##   <chr>          <int>
    ## 1 Ammodytidae        8
    ## 2 Clupeidae         35
    ## 3 Cottidae           3
    ## 4 Gadidae           61
    ## 5 Pleuronectidae     7
    ## 6 Stichaeidae        3

``` r
#filter(asv_per_family, taxon_insect %in% c("Gadidae"))
```

## create an output for all ASV’s associated with Gadidae for taking a look at in Geneious

using function from here: <https://github.com/lrjoshi/FastaTabular>

TabularToFasta

``` r
######################### Tabular format to Fasta format###############################

#this is a function to convert tabular fasta into plain fasta file
#first column should be squence names
#second column should be sequence 

TabularToFasta <- function (filename){
  file <- read.csv (file=filename, header=TRUE)
  
  file = as.data.frame(file)
  #delete if any existing file 
  
  unlink(c("dna_fasta.fasta"), force=TRUE)
  
  #give output filename
  
  sink("dna_fasta.fasta")
  
  for (i in 1:nrow(file)){
    name = paste0(">",file[i,1])
    sequence = paste(file[i,2])
    cat(name,sep="\n")
    cat(sequence,sep="\n")
  }
  
  #this is sink all the console output to the file 
  sink()
  
}

#usage 
#TabularToFasta("gene.csv")
```

``` r
gadidae <- join_asv %>%
  filter(taxon_blastn == "Gadidae") %>%
  select(ASV, asv_seq)

write.csv(gadidae, "gadidae_asv.csv")

TabularToFasta("gadidae_asv.csv")
## this saves the fasta as "dna_fasta.fasta"
```
