# Load libraries

library(igrpah)
library(FELLA)
library(tidyverse)

# Setup working environment

setwd("~/Documents/Projects/Baiqun_metabolomics/")

# Import selected compound KEGG IDs
# Selection was done based on the metabolomics_analysis_pipeline.py script

compounds <- read.table("output/compounds.txt", header = F, sep = "\t")$V1

###############################
####  Build KEGG database  ####
###############################

# Running the database building is run only once, for reusing the database just load it in with loadKEGGdata()

# Create graph from KEGG
graph <- buildGraphFromKEGGREST(organism = 'hsa')

# Initialize directory where the DB will reside
tmpdir <- "~/Documents/Projects/Baiqun_metabolomics/kegg_db/"

# Build the FELLA.DATA object
buildDataFromGraph(
  keggdata.graph = graph,
  databaseDir = tmpdir,
  internalDir = FALSE,
  matrices = "diffusion",
  normality = "diffusion",
  niter = 100)


# Load in the above generated database
fella.data <- loadKEGGdata(
  databaseDir = "~/Documents/Projects/Baiqun_metabolomics/kegg_db/",
  internalDir = FALSE
)


######################
####  ENRICHMENT  ####
######################

# Running the enrichment wrapper
enrich_res <- enrich(
  compounds = compounds,
  data = fella.data,
  method = "diffusion", 
  approx = "normality")

# Which compounds could be mapped?
enrich_res %>%
  getInput %>%
  getName(data = fella.data)

# Which compounds could NOT be mapped?
enrich_res %>%
  getExcluded()


# Visualize the resultant graph
plot(enrich_res,
     data = fella.data,
     nlimit = 100,          #Change this to limit how many nodes are displayed
     method = "diffusion",
     plotLegend = T)


# Extract results
enrich_res_graph <- generateResultsGraph(
  object = enrich_res,
  data = fella.data,
  method = "diffusion")

enrich_res_tbl <- generateResultsTable(
  object = enrich_res,
  data = fella.data,
  method = "diffusion")


# Export results
write.table(enrich_res_tbl, "output/fella_subnetwork_result.tsv", quote = F, row.names = F, col.names = T, sep = "\t")





