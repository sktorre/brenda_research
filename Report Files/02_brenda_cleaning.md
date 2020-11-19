Extract and Clean Brenda Data
================
Sarah Torrence
October 21, 2020

In summary, this file shows you how to extract and clean the parameters
and variables we require for one EC number. I created several functions
to do perform different parts of the process and then test them using
the EC of 1.1.1.1.

``` r
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE)
#libraries
library(brendaDb)
library(data.table)
library(tidyverse)
library(stringr)
library(splitstackshape)
library(janitor)
```

I wrote the database into a text file so it was easier to read into R,
but you could also use the code commented out below to read in the data
directly from the website using the brendaDb package.

``` r
#reading in the db as a data.table
brenda_data <- fread("brenda.txt")
```

``` r
#downloading the database and reading it into R
#brenda.filepath <- DownloadBrenda()
#brenda <- ReadBrenda(brenda.filepath)
```

This function, clean\_column, is used to transform the refID and
proteinID columns so that there is only one value in each row for these
variables. For example if a refID = 1,2 there would now be two rows, one
for refID = 1 and one for refID = 2. This is used to clean columns for
all the parameters and the organisms tables.

``` r
clean_column <- function(df, col_name){
  #cleans up a column that is a list of values so each cell has its own row
  
  #splitting up the values of the column
  df <- cSplit(df, splitCols = col_name , sep = ",", direction = "wide", drop = FALSE)
    
    #finding columns to input into the pivot_longer function
    param_columns <- names(df)
    columns <- param_columns[str_detect(param_columns, paste0("^",col_name,"_"))]
      
    #Transforming col_name into one column with one value for each row
    df <- df %>% pivot_longer(cols = columns, 
                                    names_to = "name", 
                                    values_to = paste0("new_",col_name), 
                                    values_drop_na = TRUE) %>% 
                               select(-name, -col_name)
}
```

This function, extract\_parameters, is used to pull all the parameters
for a single EC number. It extracts each parameter, cleans up columns
and columns names and outputs a list with each element being a dataframe
of one parameter.

``` r
extract_parameters <- function(ec, ls, parameters){
  #Extracts specified parameters for one EC, cleans up names and adds them to a new list
  
  #initializing the list
  param <- list()

  #iterate through all parameters needed
  for (i in parameters){
    #extract the parameter and clean up element names
    new_frame <- ls[[ec]][["parameters"]][[i]]
    names(new_frame)[2] <- i
    names(new_frame)[3] <- "substrate"
    names(new_frame)[4] <- paste0(i,"_commentary")
    
    #cleaning refID and proteinID so each row has one ID
    new_frame <- clean_column(new_frame,"refID")
    new_frame <- clean_column(new_frame,"proteinID")
    new_frame <- new_frame %>% select(refID = new_refID, proteinID = new_proteinID, everything())
    
    #append parameter to the list
    param <- append(param, list(new_frame))
  }
  #add names to the list
  names(param) <- parameters
  
  #return the list
  return(param)
}
```

This function, join\_parameters, takes the list of parameter dataframes
from above and joins them creating one larger data table.

``` r
join_parameters <- function(ls){
  #joins a list of parameter tables into one data table
  
  #converts each parameter into a data table
  ls <- lapply(ls, as.data.table)
  for (i in ls){
    setkey(i, proteinID, refID, substrate)
  }
  
  #join parameter data tables
  MergedDT <- Reduce(function(...) merge(..., all = TRUE, allow.cartesian = TRUE), ls)
  MergedDT <- clean_names(MergedDT)
}
```

This function, extract\_organisms, extracts and cleans the organisms
data for one EC number and outputs all the organisms in a tidy data
table.

``` r
extract_organisms <- function(ec, ls){
  #Extracts and cleans the organisms data for one EC
  
  #filtering for organisms
  organisms <- ls[[ec]][["organism"]][["organism"]]
  
  #cleaning variable names
  names(organisms)[2] <- "organism"
  
  #cleaning refID so each row has one ID
  organisms <- clean_column(organisms,"refID")
  organisms <- organisms %>% select(refID = new_refID, everything())
  
  organisms <- clean_names(organisms)
  organisms[, "protein_id"] <- lapply(organisms[, "protein_id"], as.integer)
  organisms <- as.data.table(organisms)
  
  return(organisms)
}
```

This function, join\_param\_organisms, joins the parameters and
organisms for one EC number together outputting on larger data table.

``` r
join_param_organisms <- function(parameters, organisms){
    #joins the parameters and organisms into one data table
  
    EC <- merge(organisms, parameters, by.x = c("protein_id","ref_id") , 
                 by.y = c("protein_id",  "ref_id"), all = TRUE, allow.cartesian = TRUE)
}
```

Here I have tested all of the above functions for the EC 1.1.1.1. The
final output is EC a data table with all the information we require for
EC 1.1.1.1.

``` r
#filtering to just EC 1.1.1.1 for testing
one_ec <- QueryBrenda(brenda_data, EC = "1.1.1.1")

#parameters we want
parameters <- c("km.value", "turnover.number", "ki.value", "temperature.range", "temperature.optimum", "ph.range", "ph.optimum")

#calling all the functions above
param <- extract_parameters("1.1.1.1", one_ec, parameters)

param <- join_parameters(param)

org <- extract_organisms("1.1.1.1", one_ec)

EC <- join_param_organisms(param, org)
```

The next step would be to call these functions on some bigger function
that would perform these operations for all EC numbers.

I am still unsure if there is a way to query all EC numbers at once, but
that would be the next thing I would like to look into. If not, we would
probably have to manually create a list of all the EC numbers in the
database to input. There might be a way of just getting a list of these
values.

The code above works, but I am sure it is not the most efficient way of
doing this process. I tried to add efficiency where I could but I never
did get do.call() to work and I am sure there are plenty of other areas
for improvement as well. If you have any suggestions I am happy to
implement them\! I am sure overtime looking back through this I can make
adjustments as well.
