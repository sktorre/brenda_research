How to Extract and Clean Data from Brenda
================
November 19, 2020

If you have not yet installed the `brendaDb` package in R, run this code
chuck to get the package installed on your machine.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("brendaDb", dependencies=TRUE)
```

Here are the libraries I use.

``` r
library(brendaDb)
library(data.table)
library(tidyverse)
library(stringr)
library(splitstackshape)
library(janitor)
source("brenda_extract_functions.R")
```

You can use the `brendaDb` package to download the database and read it
into R.

``` r
#downloading the database and reading it into R
brenda.filepath <- DownloadBrenda()
brenda <- ReadBrenda(brenda.filepath)
```

This process can take a while and must be done each time you want to use
the database. Instead of having to constantly download the database each
time you can write the data to a text file and read it in as a data
frame or data table as shown below.

``` r
#writing the db to a text file
fwrite(brenda,file = "brenda.txt")

#reading in the db as a data.table
brenda_data <- fread("brenda.txt")
```

Here is an example query of how to extract meaningful data from the
database. This example pulls out the information for the organism
“Acetobacter pasteurianus” from the EC number 1.1.1.1.

``` r
example <- QueryBrenda(brenda_data, EC = "1.1.1.1", organisms = "Acetobacter pasteurianus")
ace <- example[["1.1.1.1"]][["parameters"]][["km.value"]]

ace %>% filter(fieldInfo == "acetaldehyde")
```

    ## # A tibble: 2 x 5
    ##   proteinID description fieldInfo    commentary                            refID
    ##   <chr>     <chr>       <chr>        <chr>                                 <chr>
    ## 1 61        0.3         acetaldehyde #61# isozyme ADH I, pH 5.5, 25°C <11… 113  
    ## 2 61        24          acetaldehyde #61# isozyme ADH II, pH 6.0, 25°C <1… 113

To be able to extract all the information desired from the brenda
database, we can use the functions detailed [here](https://github.com/sktorre/brenda_research/blob/main/brenda_extract.md). The following is
an example of using those functions to extract the parameters listed
below for EC number 1.1.1.1.

``` r
parameters <- c("km.value", "turnover.number", "ki.value", "temperature.range", 
                "temperature.optimum", "ph.range", "ph.optimum")

brenda <- brenda_EC_data(brenda_data, "1.1.1.1", parameters)

head(brenda)
```

    ##         ec protein_id ref_id                organism uniprot commentary
    ## 1: 1.1.1.1          4     64 Drosophila melanogaster    <NA>       <NA>
    ## 2: 1.1.1.1          4     64 Drosophila melanogaster    <NA>       <NA>
    ## 3: 1.1.1.1          4     64 Drosophila melanogaster    <NA>       <NA>
    ## 4: 1.1.1.1          4     64 Drosophila melanogaster    <NA>       <NA>
    ## 5: 1.1.1.1          4     64 Drosophila melanogaster    <NA>       <NA>
    ## 6: 1.1.1.1          4     64 Drosophila melanogaster    <NA>       <NA>
    ##      substrate km_value km_value_commentary turnover_number
    ## 1:        NAD+     0.26                <NA>            <NA>
    ## 2:  butan-1-ol     3.31                <NA>            <NA>
    ## 3:  butan-2-ol     0.86                <NA>            <NA>
    ## 4:     ethanol     8.92                <NA>            <NA>
    ## 5: propan-2-ol     1.87                <NA>            <NA>
    ## 6: propan-2-ol     3.42                <NA>            <NA>
    ##    turnover_number_commentary ki_value ki_value_commentary temperature_range
    ## 1:                       <NA>     <NA>                <NA>              <NA>
    ## 2:                       <NA>     <NA>                <NA>              <NA>
    ## 3:                       <NA>     <NA>                <NA>              <NA>
    ## 4:                       <NA>     <NA>                <NA>              <NA>
    ## 5:                       <NA>     <NA>                <NA>              <NA>
    ## 6:                       <NA>     <NA>                <NA>              <NA>
    ##    temperature_range_commentary temperature_optimum
    ## 1:                         <NA>                <NA>
    ## 2:                         <NA>                <NA>
    ## 3:                         <NA>                <NA>
    ## 4:                         <NA>                <NA>
    ## 5:                         <NA>                <NA>
    ## 6:                         <NA>                <NA>
    ##    temperature_optimum_commentary ph_range ph_range_commentary ph_optimum
    ## 1:                           <NA>     <NA>                <NA>       <NA>
    ## 2:                           <NA>     <NA>                <NA>       <NA>
    ## 3:                           <NA>     <NA>                <NA>       <NA>
    ## 4:                           <NA>     <NA>                <NA>       <NA>
    ## 5:                           <NA>     <NA>                <NA>       <NA>
    ## 6:                           <NA>     <NA>                <NA>       <NA>
    ##    ph_optimum_commentary
    ## 1:                  <NA>
    ## 2:                  <NA>
    ## 3:                  <NA>
    ## 4:                  <NA>
    ## 5:                  <NA>
    ## 6:                  <NA>
