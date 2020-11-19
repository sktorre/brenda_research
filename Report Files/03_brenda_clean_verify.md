Cleaning, Verify, Testing
================
Sarah Torrence
October 30, 2020

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

``` r
#reading in the db as a data.table
brenda_data <- fread("brenda.txt")
```

I found an error in trying to run the process over several ECs in that
the parameters within some ECs were NA and even further some ECs
themselves were NA with no information included at all.

Here is an examples so you can see for yourself.

``` r
ex_ec <- QueryBrenda(brenda_data, EC = "1.1.1.111")

ex_ec[["1.1.1.111"]][["parameters"]]
```

    ## $km.value
    ## [1] NA
    ## 
    ## $turnover.number
    ## [1] NA
    ## 
    ## $ki.value
    ## [1] NA
    ## 
    ## $pi.value
    ## [1] NA
    ## 
    ## $ph.optimum
    ## [1] NA
    ## 
    ## $ph.range
    ## [1] NA
    ## 
    ## $temperature.optimum
    ## [1] NA
    ## 
    ## $temperature.range
    ## [1] NA
    ## 
    ## $specific.activity
    ## # A tibble: 2 x 5
    ##   proteinID description fieldInfo commentary                               refID
    ##   <chr>     <chr>       <lgl>     <chr>                                    <chr>
    ## 1 1         7.0         NA        #1# in crude extracts, cofactor NADH <1> 1    
    ## 2 1         12.0        NA        #1# in crude extracts, cofactor NADPH <… 1    
    ## 
    ## $ic50
    ## [1] NA
    ## 
    ## attr(,"class")
    ## [1] "brenda.sublist"

For EC 1.1.1.111 all of the parameters of interest are NA. In some cases
only one or two of the parameters of interest are NA so we want to make
sure to include those EC numbers for the parameters that do have values
but not include ones in which there are no parameter values.

I added to the code to not include NA parameters and not include EC
numbers that were inherently NA with no information.

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

extract_parameters <- function(ec, ls, parameters){
  #Extracts specified parameters for one EC, cleans up names and adds them to a new list
  
  #initializing the list
  param <- list()

  #iterate through all parameters needed
  for (i in parameters){
    
    #if the parameter is NA don't include it
    if (length(ls[[ec]][["parameters"]][[i]]) == 0 | is.na(ls[[ec]][["parameters"]][[i]])){
      break
    }
    else {
      #extract the parameter and clean up element names
      new_frame <- ls[[ec]][["parameters"]][[i]]
      names(new_frame)[2:4] <- c(i, "substrate", paste0(i,"_commentary"))
      
      #cleaning refID and proteinID so each row has one ID
      new_frame <- clean_column(new_frame,"refID")
      new_frame <- clean_column(new_frame,"proteinID")
      new_frame <- new_frame %>% select(refID = new_refID, proteinID = new_proteinID, everything())
      
      #append parameter to the list
      param <- append(param, list(new_frame))
    }
  }
  
  #return the list
  return(param)
}

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

join_param_organisms <- function(parameters, organisms){
    #joins the parameters and organisms into one data table
  
    EC <- merge(organisms, parameters, by.x = c("protein_id","ref_id") , 
                 by.y = c("protein_id",  "ref_id"), all = TRUE, allow.cartesian = TRUE)
}
```

I created the function one\_EC that calls all the other functions to
produce the output for one EC. Then I created a function that calls
one\_EC over a list of inputted ECs to produce the desired output for
all the ECs we require. To extract the information in the desired format
we only need to call `brenda_EC_data()`.

``` r
one_EC <- function(ec, ls, parameters){
  #cleans, joins and outputs a data table for one EC number

  param <- extract_parameters(ec, ls, parameters)
  
  #if the ECs has no parameter information assign it to NA
  if (length(param) == 0){
    return(NA)
  }
  else {
    param <- join_parameters(param)
    
    org <- extract_organisms(ec, ls)
    
    EC <- join_param_organisms(param, org)
  }
}

brenda_EC_data <- function(brenda_data, EC_list, parameters){
  #queries data from the brenda db for all desired ECs, processes the data and
  #outputs the data into one data table
  
  s <- Sys.time()
  #querying data from the brenda Db
  ECs <- QueryBrenda(brenda_data, EC = EC_list)
  e <- Sys.time()
  print(e - s)
  
  name <- list()
  out <- list()
  #iterating over all EC numbers and cleaning the data in a data table for each one
  for (i in 1:length(ECs)){
    o <- one_EC(names(ECs[i]),ECs, parameters)
    #if the EC is NA ignore it and move on
    if (is.na(o)){
      next
    }
    #appending each EC data table to a list
    name <- append(name, EC_list[i])
    out <- append(out, list(o))
  }
  names(out) <- name
  #merging all the data tables into one larger data table
  out_table <- rbindlist(out, use.names = TRUE, idcol = "ec", fill = TRUE)
  #filtering out rows where no valuable information is provided
  out_table <- out_table %>% filter((!is.na(km_value) & km_value != "additional information") |
                       (!is.na(turnover_number) & turnover_number != "additional information"))
  return(out_table)
}
```

How many EC’s in Brenda?

``` r
#count of ECs in brenda database
num_ec <- n_distinct(brenda_data$ID)
```

There are 7222 unique EC numbers in the brenda database.

``` r
#creating a vector of all ECs
all_ECs <- unique(brenda_data$ID)
```

By looking at `field` I can see that some of the EC numbers have been
deleted or transferred to a new EC number.

``` r
ECs_to_remove <- brenda_data %>% filter(field == "TRANSFERRED_DELETED")
```

This was causing errors in the code and these EC numbers are not needed
anyways so I removed these in the list of ECs to extract from Brenda.

``` r
brenda_data <- brenda_data %>% filter(field != "TRANSFERRED_DELETED")

all_ECs <- unique(brenda_data$ID)

new_num_ec <- n_distinct(brenda_data$ID)
```

I removed these from the data and now there are 6208 ECs we need to
extract and clean.

I want to test my function over multiple ECs so I will start with just
two and see the results.

``` r
two_ec <- all_ECs[1:2]

parameters <- c("km.value", "turnover.number", "ki.value", "temperature.range", 
                "temperature.optimum", "ph.range", "ph.optimum")

start <- Sys.time()
result <- brenda_EC_data(brenda_data, two_ec, parameters)
```

    ## Time difference of 1.681332 secs

``` r
end <- Sys.time()

runtime <- end - start
runtime
```

    ## Time difference of 2.281268 secs

The results look good and produced the output I expected, but it took a
little over 2 seconds to run which is not very efficient considering the
number of ECs we need to scale out the process to.

I want to see how efficient this function/overall process is and see how
long it would take over 100 of the ECs in the brenda database. I also
added code to calculate the run time for just the `QueryBrenda()`
portion of the code because I believe this function is taking quite a
bit of time.

``` r
hundred_ec <- all_ECs[1:100]

start <- Sys.time()
hundred_result <- brenda_EC_data(brenda_data, hundred_ec, parameters)
```

    ## Time difference of 5.852458 secs

``` r
end <- Sys.time()

runtime <- end - start
runtime
```

    ## Time difference of 12.73707 secs

It took about 13 seconds to run the process over 100 ECs, but almost
half of that time was spent just on the `QueryBrenda()` function.

``` r
#start <- Sys.time()
#final <- brenda_EC_data(brenda_data, all_ECs, parameters)
#end <- Sys.time()
 
#runtime <- end - start
#runtime
```

It took a little less than 13 minutes to run for all EC numbers. It took
5.5 minutes just to run the `QueryBrenda()` in the `brenda_EC_data()`
function. The query function is built into the `brendaDb` package so we
cannot make this portion more efficient but I can work to try to make
the rest of the code more efficient.
