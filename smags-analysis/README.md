# SMAGs analyses for Kjær et al. 2022

Code for the some of the analyses in **A Two-Million-year-old ecosystem in North Greenland uncovered by Environmental DNA** from Kjær et al. 2022
To recreate the analyses follow the instructions below:

  > This has been tested on R 4.1.1

Clone the repo with:

  ```bash
git clone https://github.com/KapCopenhagen/smags-analysis.git
cd smags-analysis
```

Then let's install the packages we used to analyase the data and plot the figures. First start R to get renv installed:

```
R
```

> If you open the project file `samgs-analysis.Rproj` in Rstudio it will perform the same steps.

If everything went well, [renv](https://rstudio.github.io/renv/articles/renv.html) will be installed and you will get a message like:

```
* Installing renv 0.15.4 ... Done!
Successfully installed and loaded renv 0.15.4.
* Project '~/Desktop/repos/smags-analysis' loaded. [renv 0.15.4]
```

And restore the environment:

```r
renv::restore()
q()
```

You will need to download the data from [here](https://doi.org/10.6084/m9.figshare.19404503) and decompress:
```bash
curl -L https://figshare.com/ndownloader/files/XXXXXX -o data.tar.gz
tar -zxvf data.tar.gz
```

Then you can run the code in the folder `scripts`:

- `01-find-ANI-threshold.R` will process the results from filtering the SMAG references after mapping and identify the most suitable read mean ANI using the elbow method
- `02-smags-analysis.R` will produce the tables and phylogenomic tree used for the analyses in the manuscript and that can be imported into anvi'o

The tables and tree can be imported into anvi'o using the following command:

```bash
anvi-interactive --manual-mode -p PROFILE-ebr.db -d smags-cov-dmg-ebr.tsv -t smags-tree.tree
```
