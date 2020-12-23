# Extract data from Brenda database
Extracting and Cleaning data from the Brenda Database

The files in this repository contain functions I created to extract and clean data from the [brenda database](https://www.brenda-enzymes.org). This makes the data easier to digest, perform analysis and join with other sources.

You can find an explanation of the functions I created [here](https://github.com/sktorre/brenda_research/blob/main/brenda_extract.md). To find a quick explanation of the `brendadb` package and an example use of my functions take a look at my [tutorial](https://github.com/sktorre/brenda_research/blob/main/brenda_extract_example.md).


If you have not yet installed the `brendaDb` package in R, run this code chuck to get the package installed on your machine.

```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("brendaDb", dependencies=TRUE)
```

