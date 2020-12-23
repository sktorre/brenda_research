# brenda_research
Extracting and Cleaning data from the Brenda Database

If you have not yet installed the `brendaDb` package in R, run this code chuck to get the package installed on your machine.

```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("brendaDb", dependencies=TRUE)
```

