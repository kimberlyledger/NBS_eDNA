taxonomic analysis of NDS eDNA samples
================
Kimberly Ledger
2022-11-08

this script using the insect R package for taxonomic identification of
amplicon sequence variants and follows this tutorial:
<https://cran.r-project.org/web/packages/insect/vignettes/insect-vignette.html>

load libraries

``` r
library(insect)
```

read in the data (ASV table from DADA2)

``` r
seqs <- read.csv("/home/kimberly.ledger/NBS_eDNA/NBS_eDNA/seqtab.csv")
```

assign taxon IDs to the DADA2 output

``` r
x <- char2dna(colnames(seqs))
## name the sequences sequentially
names(x) <- paste0("ASV", seq_along(x))
```

the MiFish classifier from the insect tutorial…

``` r
download.file("https://www.dropbox.com/s/fv3dpvws6zjvtib/classifier.rds?dl=1", 
              destfile = "/home/kimberly.ledger/NBS_eDNA/classifier.rds", mode = "wb")
```

**doesn’t download :(**

trying another classifier for MiFish…

install the RDP classifier:

- download classifier from here:
  <https://sourceforge.net/projects/rdp-classifier/>
- then copy into VM using the command prompt: scp
  rdp_classifier_2.13.zip
  <kimberly.ledger@161.55.97.134>:/home/kimberly.ledger/NBS_eDNA
- now working in the VM, decompress: unzip rdp_classifier_2.13.zip

load classifier

``` r
classifier <- readRDS("/home/kimberly.ledger/NBS_eDNA/rdp_classifier_2.13/dist/classifier.jar")
```

hmmm… so this classifier is a java file and i’m not sure how to proceed
with using it…
