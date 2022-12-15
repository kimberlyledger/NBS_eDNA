get_BASIS_data
================
Kimberly Ledger
2022-12-15

## this script queries the BASIS database from AKFIN for all 2021 data from stations where eDNA samples were collected

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

## first let’s get a list of the stations from the eDNA metadata

``` r
meta <- read.csv("NBS_2021_metadata_kjl.csv") %>%
  select(!X) %>%
  arrange(Sample_ID) %>%
  filter(Station_ID != "NA")

eDNA_stations <- unique(meta$Station_ID)
eDNA_stations
```

    ##  [1]  2  5  8 11 17 20 23 26 29 35 38 40 43 46 14

## download data

``` r
library(httr) #for accessing web services
library(tidyverse) #for converting data into exportable data frame


#Lookup data from 2021
data2021 <- httr::content(
  httr::GET('https://apex.psmfc.org/akfin/data_marts/akmp/basis_catch_spp_lh_0'), 
  type = "application/json") %>% 
    bind_rows%>%
  filter(SAMPLEYEAR=="2021")
```

## the last three digits of the “STATIONID” corresponds to the station number

``` r
data2021_split <- data2021 %>%
  separate(STATIONID, into= c(NA, "Station"), sep = -3, remove = FALSE, convert = TRUE)
```

## filter by station number

``` r
data2021_stations <- data2021_split %>%
  filter(Station == eDNA_stations)
```

    ## Warning in Station == eDNA_stations: longer object length is not a multiple of
    ## shorter object length

## make a summary table

``` r
catch_summary <- data2021_stations %>%
  group_by(Station, SPECIESNAME) %>%
  summarise(totalcatch = sum(TOTALCATCHNUM),
            totalweight = sum(TOTALCATCHWT)) %>%
  filter(totalcatch != 0)
```

    ## `summarise()` has grouped output by 'Station'. You can override using the
    ## `.groups` argument.

``` r
catch_summary
```

    ## # A tibble: 6 × 4
    ## # Groups:   Station [5]
    ##   Station SPECIESNAME            totalcatch totalweight
    ##     <int> <chr>                       <dbl>       <dbl>
    ## 1       8 Arctic Lamprey                 2           63
    ## 2      14 Chum Salmon                    1         2220
    ## 3      38 Rainbow Smelt                  5          143
    ## 4      38 Threespine stickleback       104.         540
    ## 5      43 Coho Salmon                    2          471
    ## 6      46 Chinook Salmon                64         6710
