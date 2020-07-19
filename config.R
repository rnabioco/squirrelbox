### options
shinyOptions(cache = diskCache("./app-cache", max_size = 100 * 1024^2))
options(readr.num_columns = 0)
options(stringsAsFactors = FALSE)
options(spinner.type = 6)
theme_set(theme_cowplot())
stack_size <- getOption("pandoc.stack.size", default = "1000m")
# options(shiny.reactlog = TRUE) # for checking shiny logic
# options(repos = BiocManager::repositories()) # for pushing to shinyapps.io

### folders
rpath <- "R"
datapath <- "data"
annotpath <- "annot"
listpath <- "data/lists"

### general data settings
apptitle <- "13-lined ground squirrel hibernating brain RNA-seq expression"
url <- "https://github.com/rnabioco/squirrelbox/"
s3 <- "https://s3.console.aws.amazon.com/s3/object/squirrelbox/"
docker <- "https://hub.docker.com/r/raysinensis/squirrelbox"
versionN <- "1.0.0"
geoN <- "GSE106947"
bsgenomeL <- "BSgenome.Itridecemlineatus.v???"
pageN <- 10
warningN <- 100
plot_width <- 8
plot_height <- 6
set_shinytheme <- "paper"
track_name <- "hub_2262371_KG_HiC"
track_url <- "https://squirrelhub.s3-us-west-1.amazonaws.com/hub/hub_brain.txt"
gmt_file <- "c5.all.v7.1.symbols.gmt"
gmt_short <- "GO_"
sig_cut <- 0.001
ncore <- parallel::detectCores() - 1
start_tutorial <- TRUE
verbose_bench <- FALSE
chrlimit <- 20

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
  "LT",
  "EAr",
  "Ar",
  "LAr",
  "SpD"
)
region_order <- c(
  "Forebrain",
  "Hypothalamus",
  "Medulla",
  "Adrenal",
  "Kidney",
  "Liver"
)
region_main <- c(
  "Forebrain",
  "Hypothalamus",
  "Medulla"
)
region_main2 <- c(
  "Adrenal",
  "Kidney",
  "Liver"
)
region_short <- c(
  "fb",
  "hy",
  "med",
  "adr",
  "kid",
  "liv"
)
region_short_main <- c(
  "fb",
  "hy",
  "med"
)
region_one <- c(
  "f",
  "h",
  "m",
  "a",
  "k",
  "l"
)
