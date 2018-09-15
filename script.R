#!/usr/bin/env Rscript

# Loading libraries:
library("raster")
library("XML")
library("snowfall")
library("rJava")
library("rgdal")
library("ggplot2")
library(dplyr)
# source functions in ./fct
## the modeling workflow
#library(devtools)
#install_github("Model-R/modelr_pkg")
library(ModelR)
## several helper functions
source("./fct/read.eval1.R")


# Loading environmental data, study area mask
file <- list.files("./env", pattern = "1K", full.names = T)
#file <- list.files("./env", pattern = "^eig", full.names = T)
predictors <- stack(file[1])
# Cortando pela MataAtlantica:
mascara <- readOGR(dsn = "./data", layer = "Bioma_MA1148")

# FLORA ----
# Loading occurrences:

plants <- read.csv("./data/Flora_endemic_final.csv", row.names = 1)
birds <- read.csv("./data/Birds_endemic_final.csv", row.names = 1)
amphibians <- read.csv("./data/Amphibians_endemic_final.csv", row.names = 1)
occs <- bind_rows(plants, birds) %>% bind_rows(amphibians)
#head(occs)
#count(occs, grupo)
#occs %>% select(sp, grupo) %>% distinct() %>% count(grupo)

#Data on amphibian occurrence was obtained from35, with updates from the authors, and comprised 114 endemic species (3,786 occurrences). Data on bird occurrence was obtained from the Global Biodiversity Information Facility database36, and comprised 223 endemic species (12,085 occurrences). Data on plants occurrence was obtained from NeoTropTree and SpeciesLink37 , and comprised 846 endemic species and 44,024 records . ---
#As a consequence, at the end of this modelling phase 51 amphibian species, 122 bird species, and 612 woody plant species endemic to the Brazilian Atlantic Forest composed the final potential richness maps

#quitar nombres
library(purrr)
occs_names <- occs$sp %>% map(~ flora::remove.authors(.)) %>% simplify2array()
occs2 <- cbind(occs, occs_names)
head(occs2)
# Defining names to be modelled (after taxa and spacial cleaning)
#write.csv(occs2, "./data/occs_final_corrected_names.csv")
lista_locs <- occs2 %>% split(.$sp)
lista_locs <- lista_locs
#using the reproducible example
set.seed(712)
# MODELOS.R----

#iniciar snowfall
{
sfInit(parallel = T, cpus = 24, slaveOutfile = "FLORA_150918_modelagem.log")
# 	#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #
sfLibrary(raster)
sfLibrary(ModelR)
#sfSource("./fct/modelos.R")
#Com buffer ----
tInicial <- Sys.time()
sfClusterApplyLB(lista_locs,
                 fun = function(x) {
                     ModelR::do_enm(species_name = unique(x[,"occs_names"]),
                                    occurrences = x[ ,c("lon", "lat")],
                                    predictors = predictors,
                                    models_dir = "./FLORA_buffermax_15092018",
                                    buffer_type = "max",
                                    seed = 712,
                                    clean_dupl = T,
                                    clean_nas = T,
                                    plot_sdmdata = T,
                                    n_back = 500,
                                    partition_type = "crossvalidation",
                                    cv_n = 1,
                                    cv_partitions = 3,
                                    rf = T,
                                    svm.k = T,
                                    maxent = T,
                                    mask = mascara,
                                    write_png = T)
                     }
)

tFinal <- Sys.time()
tFinal - tInicial
sfStop()
}

# finalModel ----

{
sfInit(parallel = T, cpus = 24, slaveOutfile = "FLORA_150918_final_model.log")
# 	#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #
sfLibrary(raster)
sfLibrary(ModelR)

tInicial <- Sys.time()
sfClusterApplyLB(lista_locs,
                 fun = function(x) {
                     ModelR::final_model(species_name = unique(x[,"occs_names"]),
                                         select_partitions = T,
                                         threshold = "spec_sens",
                                         select_par = "TSS",
                                         select_par_val = 0.7,
                                         which_models = "bin_mean",
                                         write_png = T,
                                         models_dir = "./FLORA_buffermax_15092018")
                     })
tFinal <- Sys.time()
tFinal - tInicial
sfStop()
}

# ENSEMBLE ----
{
sfInit(parallel = T, cpus = 24, slaveOutfile = "FLORA_150918_ensemble_model.log")
# 	#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #
sfLibrary(raster)
sfLibrary(ModelR)

tInicial <- Sys.time()
sfClusterApplyLB(lista_locs,
                 fun = function(x) {
                     ModelR::ensemble_model(
                         species_name = unique(x[, "occs_names"]),
                         occurrences = x[, c("lon", "lat")],
                         models_dir = "./FLORA_buffermax_15092018",
                         which_models = "bin_mean",
                         consensus = T,
                         consensus_level = 0.5,
                         write_png = T)
                     })
tFinal <- Sys.time()
tFinal - tInicial
sfStop()
}
