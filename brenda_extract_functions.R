
library(brendaDb)
library(data.table)
library(tidyverse)
library(stringr)
library(splitstackshape)
library(janitor)

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
  
  #querying data from the brenda Db
  ECs <- QueryBrenda(brenda_data, EC = EC_list)
  
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

