# squirrelBox

Shiny app to browse hibernating ground squirrel brain tissue RNA-seq data, new genome assembly and transcriptome annotation.

Tested on macOS/Win10/Linux, with R3.6/4.0/4.1

Please see manuscripts published in Frontiers in Physiology for [brain](https://www.frontiersin.org/articles/10.3389/fphys.2020.624677) and [liver](https://www.frontiersin.org/articles/10.3389/fphys.2021.662132).

# Getting squirrelBox

## 1. web version

https://raysinensis.shinyapps.io/squirrelBox/?gene=Zfp36l2

## 2. download and run in Rstudio

Alternatively, download all files of this [repo](https://github.com/rnabioco/squirrelbox/archive/master.zip), and run app in RStudio.
`shiny::runApp(folder_location)`

All packages required are listed at the start of `global.R`

## 3. docker

[image on S3](https://squirrelbox.s3-us-west-2.amazonaws.com/squirrelbox_docker.tar)
```
docker load --input squirrelbox_docker.tar
docker run --rm -p 80:80 squirrelbox_docker
```
or find latest repo in [dockerhub](https://hub.docker.com/r/raysinensis/squirrelbox)

# Using squirrelBox

Here is one example of data exploration. 
![ZFPs](www/zfp.gif)

1) Search for zinc finger protein genes. 
2) Move these genes to Genelist. 
3) Make lineplot and heatmap from Genelist.
4) GO term analysis. 
5) Kmer analysis of Genelist 3'UTRs. 
6) Trimming list of genes of interest by intersection to other gene lists. 
7) View display on the genome. 

Note that clicking on a row in any table or any interactive element in the plots readies the gene query process.

# Configuring squirrelBox

squirrelBox is developed and tested on only a handful of machine/OS combinations with common resolution settings. Please check the sidebar "Options" and "Order" options to tailor the presentation. Alternatively, settings can be edited in the `config.R` file.

Libraries are up-to-date as of July 18 2020, however, in case of potential function-breaking updates, please check latest confirmed [sessionInfo](https://github.com/rnabioco/squirrelbox/issues/81).

# BSgenome package (temporary)

A BSgenome package for the new *Ictidomys tridecemlineatus* genome was forged to aid sequence query and manipulation. This is not directly required for squirrelBox browsing. However, for interested parties, a functional draft version is hosted here: [tar.gz on S3](https://squirrelbox.s3-us-west-2.amazonaws.com/BSgenome/BSgenome.SQ1_1.0.tar.gz). Installation as below:

```
install.packages(location to BSgenome.SQ1_1.0.tar.gz file,
                 repo = NULL,
                 type = "source")
```