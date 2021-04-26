### options
shinyOptions(cache = diskCache("./app-cache", max_size = 100 * 1024^2))
options(readr.num_columns = 0)
options(stringsAsFactors = FALSE)
options(spinner.type = 6)
theme_set(theme_cowplot())
stack_size <- getOption("pandoc.stack.size", default = "1000m")
options(repos = BiocManager::repositories()) # for pushing to shinyapps.io
# options(shiny.reactlog = TRUE) # for checking shiny logic

### folders
rpath <- "R" # additional R code
datapath <- "data" # data, previewed and for download
annotpath <- "annot" # annotations from external sources
listpath <- "data/lists" # csv and txt lists used by Venn and Circos

### general data settings
apptitle_short <- "squirrelBox"
apptitle <- "13-lined ground squirrel hibernating liver and brain sequencing"
url <- "https://github.com/rnabioco/squirrelbox/"
s3 <- "https://squirrelbox.s3-us-west-2.amazonaws.com/zip6/squirrelbox6.tar.gz"
docker <- "https://hub.docker.com/r/raysinensis/squirrelbox"
versionN <- "1.2.0"
geoN <- "GSE106947"
genomeN <- "NCBI Assembly HiC_Itri_2"
genomeL <- "https://www.ncbi.nlm.nih.gov/assembly/GCA_016881025.1"
manuscriptN <- "Dynamic RNA Regulation in the Brain Underlies Physiological Plasticity in a Hibernating Mammal"
manuscriptN2 <- "Liver transcriptome dynamics during hibernation are shaped by a shifting balance between transcription and RNA stability"
manuscriptL <- "https://www.frontiersin.org/articles/10.3389/fphys.2020.624677"
manuscriptL2 <- "https://www.frontiersin.org/articles/10.3389/fphys.2021.662132"
annotN <- "new annotation"
annotL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/881/025/GCA_016881025.1_HiC_Itri_2/GCA_016881025.1_HiC_Itri_2_genomic.gtf.gz"
bsgenomeN <- "BSgenome.Itridecemlineatus.HiC_Itri_2"
bsgenomeL <- "https://squirrelbox.s3-us-west-2.amazonaws.com/BSgenome/BSgenome.SQ1_1.0.tar.gz"
pageN <- 10 # number of lines, for tables
warningN <- 100 # number of genes, for throwing warnings in line and heat plots
plot_width <- 8
plot_height <- 6
set_shinytheme <- "paper"
track_name <- "hub_1512849_KG_HiC"
track_url <- "https://squirrelhub.s3-us-west-1.amazonaws.com/hub/hub.txt"
gmt_file <- "c5.all.v7.1.symbols.gmt"
gmt_short <- "GO_"
sig_cut <- 0.001
ncore <- parallel::detectCores() - 1
start_tutorial <- TRUE # whether to start loaded app with tutorial
start_tabhint <- TRUE # whether to pop up message on first switch to tab
verbose_bench <- FALSE
chrlimit <- 17
qc_report <- "multiqc_report.html"

### choose and order columns
table_cols <- c(
  "gene_id", # comment out to remove from table outputs
  "unique_gene_symbol",
  "clean_gene_symbol",
  "original_gene_name",
  "source"
)

orf_cols_join <- c(
  "gene_id",
  "chrom",
  "source",
  "orf_len",
  "micropeptide_pred",
  "micropeptide_homology",
  "exons",
  "rna_len",
  "transcript_id",
  "novel",
  "min_padj",
  "domains",
  "br_expr",
  "edited",
  "majiq_directed"
)

orf_cols <- c(
  "gene_id",
  "transcript_id",
  "rna_len",
  "orf_len",
  "micropeptide_pred",
  "exons",
  "edited"
)

maj_cols <- c(
  "region",
  "comp_pair",
  "significant_90",
  "significant_99",
  "LSV_ID",
  "A5SS",
  "A3SS",
  "ES"
)

### sample settings, define state colors and order, region order
state_cols <- c(
  SA = rgb(255, 0, 0, maxColorValue = 255),
  IBA = rgb(67, 205, 128, maxColorValue = 255),
  Ent = rgb(155, 48, 255, maxColorValue = 255),
  ET = rgb(120, 10, 130, maxColorValue = 255),
  LT = rgb(25, 25, 112, maxColorValue = 255),
  EAr = rgb(0, 0, 255, maxColorValue = 255),
  Ar = rgb(0, 0, 255, maxColorValue = 255),
  LAr = rgb(0, 102, 102, maxColorValue = 255),
  SpD = rgb(255, 165, 0, maxColorValue = 255)
)
state_order <- c(
  "SA",
  "IBA",
  "Ent",
  "ET",
  "LT",
  "EAr",
  "Ar",
  "LAr",
  "SpD"
)

region_main2 <- c(
  "Forebrain",
  "Hypothalamus",
  "Medulla",
  "Adrenal",
  "Kidney"
)
region_main <- c(
  "Liver",
  "Liver_GRO-seq",
  "Liver_GRO-seq_promoter"
)
region_order <- c(region_main, region_main2)
region_short <- c(
  "liv",
  "gro",
  "gropro",
  "fb",
  "hy",
  "med",
  "adr",
  "kid"
)
region_short_main <- c(
  "liv",
  "gro",
  "gropro"
)
region_one <- c(
  "l",
  "g",
  "p",
  "f",
  "h",
  "m",
  "a",
  "k"
)
hide_region <- c("Adrenal", "Kidney")