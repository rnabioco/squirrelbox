# squirrelBox
Shiny app to browse hibernating ground squirrel brain tissue RNA-seq data, new genome assembly and transcriptome annotation.

Tested on macOS/Win10/Linux, with R3.6 and R4.0.

# 1. web version
https://raysinensis.shinyapps.io/squirrelBox/?gene=Zfp36l2

# 2. download and run in Rstudio
Alternatively, download all files of this [repo](https://github.com/rnabioco/squirrelbox/archive/master.zip), and run app in RStudio.
`shiny::runApp(folder_location)`

All packages required are listed at the start of `global.R`

# 3. docker
[image](https://squirrelbox.s3-us-west-2.amazonaws.com/squirrelbox_docker.tar)
```
docker load --input squirrelbox_docker.tar
docker run --rm -p 80:80 squirrelbox_docker
```