# Title: Jump around: selecting Markov Chain Monte Carlo parameters and
# diagnostics for improved food web model quality and ecosystem representation

# Authors: Gerber & Scharler (2024)
# Journal: Ecological Informatics
# GitHub repo:

###############################################################################
# Reproducible Code for Section 2 Materials and Methods
###############################################################################

# All code used to simulate and analyse results
# All experimental data used in the manuscript is too large to upload
# To re-generate the outputs, please use the below code
# Remember to set a working directory/use an R project
# There are many steps - please follow in order

###############################################################################
# Required R packages

library(tidyverse)
library(ggh4x)
library(ggpubr)
library(enaR) # v3.0.0 (Lau et al., 2017)
library(autoLIMR) # v3.0.1 (Gerber et al., 2023)
library(effsize) # v0.8.1 Torchiano (2020)
library(rstatix) # v.0.7.1, Kassambara (2022)

###############################################################################
# 2.1 Model Construction
###############################################################################

# Construct network model from .csv input files
# Used together with two .csv input files in folder 'autoGen_input'
# Network input data = 'autoGen_input/4node_input_data.csv'
# Adjacency matrix data = 'autoGen_input/4node_admat_data.csv'

# Construct LIM declaration file
autoGen(
  net_data_input = 'autoGen_input/4node_input_data.csv',
  adj_mat_input = 'autoGen_input/4node_admat_data.csv',
  primary_producer = c('Phytoplankton', 'Microphytobenthos'),
  NLNode = c('susPOC', 'sedPOC'),
  respiration = TRUE,
  force = T
) # Save to working directory, select "1"

###############################################################################
# 2.2 Flow & ecological indicator uncertainty analysis with LIM-MCMC algorithm
#     scenarios
###############################################################################

###############################################################################
# 2.2.1: Algorithm scenario development

# As described in manuscript
# Already included below
# Nothing to do here, simply keeping this section for continuity.

###############################################################################
# 2.2.2 Calculating multiple plausible networks

# One ensemble of networks generated for each of 30 LIM-MCMC algorithm scenarios
# Scenarios differ according to: 1) starting point, 2) jump size, 3) iterations
# This step automatically calculates MCMC diagnostics on the solved flow values
# Outputs for each scenario:
#     1) .rds file with all information (solved flows, packed objects)
#     2) .csv file with all MCMC diagnostics

# First, set up a function to loop over experiments
# The function does a few things
#     - Solves multiple plausible networks
#     - Runs MCMC diagnostics
#     - Packs solved network values into network objects (for analysis with ENA)
#     - Automatically saves to working directory:
#       1) .rds file with all information (solved flows, packed objects)
#       2) .csv file with all MCMC diagnostics

#' @title Function 'experiment' to calculate multiple plausible networks and
#' evaluate with MCMC convergence diagnostics.
#' @description Writes flow networks to file as .rds files.
#' Writes MCMC stats & visual diagnostics table to file.
#' @param file The LIM declaration file
#' @param jmp The jump size (mgC m^-2 d^-1)
#' @param iter The number of iterations to return
#' @param x0 The starting solution, either default LSEI
#' (Haskell and Hansen) or 'central' solution calculated with
#' LIM::Xranges() (van Oevelen et al., 2010)
#' @param pack Pack flows into network objects? Default = TRUE
#' @param replicate set.seed() replicate. Defaults to NULL.

experiment <-
  function (file, jmp, iter, x0, pack, flow, replicate) {
    ## Mini function to create folders if they do not exist
    if (!dir.exists("mcmc_diags")) {
      dir.create("mcmc_diags")
      message("Folder 'mcmc_diags' created. Writing files into this folder.")
    } else {
      message("Writing files into existing folder(s).")
    }
    
    if (!dir.exists("networks_rds")) {
      dir.create("networks_rds")
      message("Folder 'networks_rds' created. Writing files into this folder.")
    } else {
      message("Please wait...")
    }
    
    # Run networks
    set.seed(replicate)
    x <-
      multi_net(
        file = file,
        jmp = jmp,
        iter = iter,
        x0 = x0,
        pack = pack
      )
    
    # Set up metadata
    if (is.null(x0)) {
      x0 <- "NULL"
    } else {
      x0 <- x0
    }
    
    if (is.null(jmp)) {
      jmp <- "NULL"
    } else {
      jmp <- paste0(jmp)
    }
    
    # Save network rds into a file
    objects <-
      paste0("jmp", jmp, "iter", iter, "x0", x0, "rep", replicate)
    saveRDS(object = x,
            file = paste0("networks_rds/", objects, ".rds"))
    
    # MCMC Diags
    diags <- mcmc_diags(x = x)
    diags_table <- do.call(cbind, diags)
    diags_table <- diags_table %>%
      mutate(
        jmp = paste0(jmp),
        iter = paste0(iter),
        x0 = paste0(x0),
        replicate = paste0(replicate)
      )
    write.csv(diags_table, file = paste0("mcmc_diags/", "diags", objects, ".csv"))
    
  }

# Next, using the 'experiment' function, solve multiple plausible networks of
# the weighted 4node network (defined by the LIM declaration file) for 30
# LIM-MCMC algorithm scenarios of differing jump sizes, starting solutions,
# and number of iterations.

lim <-
  (paste0(getwd(), "/weighted_limfiles/Weighted_Network_LIMfile.R"))

# Scenario 1 (LSEI, jmp = 0.01, iter = 5000)
experiment(
  file = lim,
  jmp = 0.01,
  iter = 5000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 2 (LSEI, jmp = 0.01, iter = 10000)
experiment(
  file = lim,
  jmp = 0.01,
  iter = 10000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 3 (LSEI, jmp = 0.01, iter = 20000)
experiment(
  file = lim,
  jmp = 0.01,
  iter = 20000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 4 (LSEI, jmp = 0.1, iter = 5000)
experiment(
  file = lim,
  jmp = 0.1,
  iter = 5000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 5 (LSEI, jmp = 0.1, iter = 10000)
experiment(
  file = lim,
  jmp = 0.1,
  iter = 10000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 6 (LSEI, jmp = 0.1, iter = 20000)
experiment(
  file = lim,
  jmp = 0.1,
  iter = 20000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 7 (LSEI, jmp = 1, iter = 5000)
experiment(
  file = lim,
  jmp = 1,
  iter = 5000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 8 (LSEI, jmp = 1, iter = 10000)
experiment(
  file = lim,
  jmp = 1,
  iter = 10000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 9 (LSEI, jmp = 1, iter = 20000)
experiment(
  file = lim,
  jmp = 1,
  iter = 20000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 10 (LSEI, jmp = 10, iter = 5000)
experiment(
  file = lim,
  jmp = 10,
  iter = 5000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 11 (LSEI, jmp = 10, iter = 10000)
experiment(
  file = lim,
  jmp = 10,
  iter = 10000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 12 (LSEI, jmp = 10, iter = 20000)
experiment(
  file = lim,
  jmp = 10,
  iter = 20000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 13 (LSEI, jmp = 100, iter = 5000)
experiment(
  file = lim,
  jmp = 100,
  iter = 5000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 14 (LSEI, jmp = 100, iter = 10000)
experiment(
  file = lim,
  jmp = 100,
  iter = 10000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 15 (LSEI, jmp = 100, iter = 20000)
experiment(
  file = lim,
  jmp = 100,
  iter = 20000,
  x0 = NULL,
  pack = TRUE,
  replicate = 1
)

# Scenario 16 (Central, jmp = 0.01, iter = 5000)
experiment(
  file = lim,
  jmp = 0.01,
  iter = 5000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 17 (Central, jmp = 0.01, iter = 10000)
experiment(
  file = lim,
  jmp = 0.01,
  iter = 10000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 18 (Central, jmp = 0.01, iter = 20000)
experiment(
  file = lim,
  jmp = 0.01,
  iter = 20000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 19 (Central, jmp = 0.1, iter = 5000)
experiment(
  file = lim,
  jmp = 0.1,
  iter = 5000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 20 (Central, jmp = 0.1, iter = 10000)
experiment(
  file = lim,
  jmp = 0.1,
  iter = 10000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 21 (Central, jmp = 0.1, iter = 20000)
experiment(
  file = lim,
  jmp = 0.1,
  iter = 20000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 22 (Central, jmp = 1, iter = 5000)
experiment(
  file = lim,
  jmp = 1,
  iter = 5000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 23 (Central, jmp = 1, iter = 10000)
experiment(
  file = lim,
  jmp = 1,
  iter = 10000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 24 (Central, jmp = 1, iter = 20000)
experiment(
  file = lim,
  jmp = 1,
  iter = 20000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 25 (Central, jmp = 10, iter = 5000)
experiment(
  file = lim,
  jmp = 10,
  iter = 5000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 26 (Central, jmp = 10, iter = 10000)
experiment(
  file = lim,
  jmp = 10,
  iter = 10000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 27 (Central, jmp = 10, iter = 20000)
experiment(
  file = lim,
  jmp = 10,
  iter = 20000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 28 (Central, jmp = 100, iter = 5000)
experiment(
  file = lim,
  jmp = 100,
  iter = 5000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 29 (Central, jmp = 100, iter = 10000)
experiment(
  file = lim,
  jmp = 100,
  iter = 10000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

# Scenario 30 (Central, jmp = 100, iter = 20000)
experiment(
  file = lim,
  jmp = 100,
  iter = 20000,
  x0 = 'central',
  pack = TRUE,
  replicate = 1
)

###############################################################################
# 2.2.3 Ecological Network Analysis

#' @title Function 'extract' to balance networks, extract flow values, and
#' calculate ENA metrics for each set of multiple plausible networks.
#' @description Another multi-purpose function that inputs the multiple plausible
#' network object (.rds), and then (1) balances each network (for ENA calculations),
#' (2) extracts all flow values and saves into the working directory (.csv),
#' (3) calculates ENA metrics, (4) saves ENA metrics to working directory (.csv).
#' A few outputs are written to the working directory:
# - Folder 'balanced_rds' stores balanced network objects for each scenario (.rds)
# - Folder 'flows' stores flow values for each scenario (.csv)
# - Folder 'ena' stores ENA values for each scenario (.csv)
# - Folder 'ena_rds' stores ENA objects for each scenario (.rds)
#'
#' @param rdsfile The rds object of multiple plausible networks generated with
#' function multi_net()
#' @param jmp The jump size (mgC m^-2 d^-1) used to solve the MPN
#' @param iter The number of iterations used to solve the MPN
#' @param x0 The starting solution, either default LSEI
#' (Haskell and Hansen) or 'central' solution calculated with
#' LIM::Xranges() (van Oevelen et al., 2010) used to solve the MPN
#' @param balance Balance the networks, or not. Default = TRUE.

extract <- function(rdsfile, jmp, iter, x0, replicate, balance) {
  if (is.null(x0)) {
    x0 <- "NULL"
  } else {
    x0 <- x0
  }
  
  if (is.null(jmp)) {
    jmp <- "NULL"
  } else {
    jmp <- paste0(jmp)
  }
  
  objects <-
    paste0("jmp",
           jmp,
           "iter",
           iter,
           "x0",
           x0,
           "bal",
           balance,
           "rep",
           replicate)
  
  ### Create folders if required
  
  # Check if the 'flows' folder exists
  if (!dir.exists("flows")) {
    # Create the folder if it does not exist
    dir.create("flows")
    message("Folder 'flows' created. Writing files into this folder.")
  } else {
    message("Folder 'flows' already exists. Writing files into existing folder.")
  }
  
  # Check if the 'ena' folder exists
  if (!dir.exists("ena")) {
    # Create the folder if it does not exist
    dir.create("ena")
    message("Folder 'ena' created. Writing files into this folder.")
  } else {
    message("Folder 'ena' already exists. Writing files into existing folder.")
  }
  
  # Check if the 'ena_rds' folder exists
  if (!dir.exists("ena_rds")) {
    # Create the folder if it does not exist
    dir.create("ena_rds")
    message("Folder 'ena_rds' created. Writing files into this folder.")
  } else {
    message("Folder 'ena_rds' already exists. Writing files into existing folder.")
  }
  
  # Check if the 'balanced_rds' folder exists
  if (!dir.exists("balanced_rds")) {
    # Create the folder if it does not exist
    dir.create("balanced_rds")
    message("Folder 'balanced_rds' created. Writing files into this folder.")
  } else {
    message("Folder 'balanced_rds' already exists. Writing files into existing folder.")
  }
  
  flow_path <- paste0(getwd(), "/flows/")
  enards_path <- paste0(getwd(), "/ena_rds/")
  bal_path <- paste0(getwd(), "/balanced_rds/")
  ena_path <- paste0(getwd(), "/ena/")
  
  # 1. Read in RDS file
  x <- readRDS(paste0(rdsfile))
  
  # 2. Extract flow table
  flow <- as.data.frame(x[["solved.flow.values"]])
  flow <- flow %>% mutate(
    Iteration = as.factor(1:nrow(flow)),
    jmp = as.factor(paste0(jmp)),
    iter = as.factor(paste0(iter)),
    x0 = as.factor(paste0(x0)),
    replicate = as.factor(paste0(replicate))
  )
  
  # 3. Calculate total herbivory and total detritivory
  # Then calculate D:H ratio
  flow <- flow |>
    mutate(
      Microphtyobenthos_Q_Phytoplankton = NULL,
      Phytoplankton_Q_Microphytobenthos = NULL
    )
  flow$detritivory <-
    flow %>% dplyr::select(starts_with(c("susPOCNLNode_Q_", "sedPOCNLNode_Q_"))) %>% rowSums
  flow$herbivory <-
    flow %>% dplyr::select(starts_with(c(
      "Phytoplankton_Q_", "Q_Microphytobenthos_Q_"
    ))) %>% rowSums
  flow$DH <- flow$detritivory / flow$herbivory
  DHtable <-
    flow %>% dplyr::select(Iteration, jmp, iter, x0, replicate, DH) # Will add to ena list later
  
  # 4. Write flow table to .csv
  write.csv(flow, paste0(flow_path, "flows_", objects, ".csv"))
  rm(flow)
  
  # 5. Perform ENA on balanced networks
  invisible(nets <-
              lapply(x[["packed.nets"]], balance, method = "AVG2")) # balance networks
  rm(x)
  saveRDS(nets, file = paste0(bal_path, "balanced_", objects, ".rds"))
  Flow <- do.call(rbind, sapply(lapply(nets, enaFlow), "[", "ns"))
  Asc <- do.call(rbind, lapply(nets, enaAscendency))
  rm(nets)
  ena <- as.data.frame(cbind(Flow, Asc))
  ena <- ena %>% mutate(
    Iteration = as.factor(1:nrow(ena)),
    jmp = as.factor(paste0(jmp)),
    iter = as.factor(paste0(iter)),
    x0 = as.factor(paste0(x0)),
    replicate = as.factor(paste0(replicate)),
    TSTc = FCI * TST,
    FCIwithTSTp = TSTc / TSTp
  )
  rm(Flow)
  rm(Asc)
  ena_all <- left_join(ena, DHtable)
  rm(ena)
  rm(DHtable)
  
  # 6. write getnslist to rds and .csv
  saveRDS(ena_all, file = paste0(enards_path, "ena_rds_", objects, ".rds"))
  write.csv(ena_all, file = paste0(ena_path, "ena_", objects, ".csv"))
  rm(ena_all)
  
}

## Apply custom function 'extract' to each .rds file (so, each scenario)
file_path <- paste0(getwd(), "/networks_rds/")

# Scenario 1 (LSEI, jmp = 0.01, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp0.01iter5000x0NULLrep1.rds'),
  jmp = 0.01,
  iter = 5000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 2 (LSEI, jmp = 0.01, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp0.01iter10000x0NULLrep1.rds'),
  jmp = 0.01,
  iter = 10000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 3 (LSEI, jmp = 0.01, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp0.01iter20000x0NULLrep1.rds'),
  jmp = 0.01,
  iter = 20000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 4 (LSEI, jmp = 0.1, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp0.1iter5000x0NULLrep1.rds'),
  jmp = 0.1,
  iter = 5000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 5 (LSEI, jmp = 0.1, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp0.1iter10000x0NULLrep1.rds'),
  jmp = 0.1,
  iter = 10000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 6 (LSEI, jmp = 0.1, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp0.1iter20000x0NULLrep1.rds'),
  jmp = 0.1,
  iter = 20000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 7 (LSEI, jmp = 1, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp1iter5000x0NULLrep1.rds'),
  jmp = 1,
  iter = 5000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 8 (LSEI, jmp = 1, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp1iter10000x0NULLrep1.rds'),
  jmp = 1,
  iter = 10000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 9 (LSEI, jmp = 1, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp1iter20000x0NULLrep1.rds'),
  jmp = 1,
  iter = 20000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 10 (LSEI, jmp = 10, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp10iter5000x0NULLrep1.rds'),
  jmp = 10,
  iter = 5000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 11 (LSEI, jmp = 10, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp10iter10000x0NULLrep1.rds'),
  jmp = 10,
  iter = 10000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 12 (LSEI, jmp = 10, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp10iter20000x0NULLrep1.rds'),
  jmp = 10,
  iter = 20000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 13 (LSEI, jmp = 100, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp100iter5000x0NULLrep1.rds'),
  jmp = 100,
  iter = 5000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 14 (LSEI, jmp = 100, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp100iter10000x0NULLrep1.rds'),
  jmp = 100,
  iter = 10000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 15 (LSEI, jmp = 100, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp100iter20000x0NULLrep1.rds'),
  jmp = 100,
  iter = 20000,
  x0 = "LSEI",
  balance = TRUE,
  replicate = 1
)

# Scenario 16 (central, jmp = 0.01, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp0.01iter5000x0NULLrep1.rds'),
  jmp = 0.01,
  iter = 5000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 17 (central, jmp = 0.01, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp0.01iter10000x0NULLrep1.rds'),
  jmp = 0.01,
  iter = 10000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 18 (central, jmp = 0.01, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp0.01iter20000x0NULLrep1.rds'),
  jmp = 0.01,
  iter = 20000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 19 (central, jmp = 0.1, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp0.1iter5000x0NULLrep1.rds'),
  jmp = 0.1,
  iter = 5000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 20 (central, jmp = 0.1, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp0.1iter10000x0NULLrep1.rds'),
  jmp = 0.1,
  iter = 10000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 21 (central, jmp = 0.1, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp0.1iter20000x0NULLrep1.rds'),
  jmp = 0.1,
  iter = 20000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 22 (central, jmp = 1, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp1iter5000x0NULLrep1.rds'),
  jmp = 1,
  iter = 5000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 23 (central, jmp = 1, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp1iter10000x0NULLrep1.rds'),
  jmp = 1,
  iter = 10000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 24 (central, jmp = 1, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp1iter20000x0NULLrep1.rds'),
  jmp = 1,
  iter = 20000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 25 (central, jmp = 10, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp10iter5000x0NULLrep1.rds'),
  jmp = 10,
  iter = 5000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 26 (central, jmp = 10, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp10iter10000x0NULLrep1.rds'),
  jmp = 10,
  iter = 10000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 27 (central, jmp = 10, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp10iter20000x0NULLrep1.rds'),
  jmp = 10,
  iter = 20000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 28 (central, jmp = 100, iter = 5000)
extract(
  rdsfile = paste0(file_path, 'jmp100iter5000x0NULLrep1.rds'),
  jmp = 100,
  iter = 5000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 29 (central, jmp = 100, iter = 10000)
extract(
  rdsfile = paste0(file_path, 'jmp100iter10000x0NULLrep1.rds'),
  jmp = 100,
  iter = 10000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

# Scenario 30 (central, jmp = 100, iter = 20000)
extract(
  rdsfile = paste0(file_path, 'jmp100iter20000x0NULLrep1.rds'),
  jmp = 100,
  iter = 20000,
  x0 = "central",
  balance = TRUE,
  replicate = 1
)

###############################################################################
# 2.2.4 Statistical analysis

# Compare flows and ENA metrics between scenarios with Cliffs Delta
# Read in all flows from folder 'flows', save as one table, tidy, and extract
# the required flows (2 small, 2 medium, 2 large)

myfiles <-
  list.files(
    path = paste0(getwd(), "/flows"),
    pattern = "*.csv",
    full.names = TRUE
  )
flows_tidy <-
  sapply(myfiles, read.csv, simplify = FALSE) %>% bind_rows(.id = "id") |>
  select(
    jmp,
    iter,
    x0,
    Iteration,
    Bacteria_R,
    susPOCNLNode_Q_Oligochaeta,
    Gastropoda_EX,
    susPOCNLNode_Q_Nematoda,
    Platyhelminthes_R,
    Arachnida_R
  ) |>
  pivot_longer(Bacteria_R:Arachnida_R,
               names_to = 'flow',
               values_to = 'value') |>
  mutate(
    jmp = as.factor(jmp),
    iter = as.factor (iter),
    x0 = as.factor(x0),
    Iteration = as.factor(Iteration),
    flow = as.factor(flow),
    flow = str_replace(flow, 'susPOCNLNode', 'susPOC'),
    flow_magnitude = as.factor(
      case_when(
        flow %in% c('Bacteria_R', 'susPOC_Q_Oligochaeta') ~ 'Large flows',
        flow %in% c('Gastropoda_EX', 'susPOC_Q_Nematoda') ~ 'Medium flows',
        flow %in% c('Platyhelminthes_R', 'Arachnida_R') ~ 'Small flows'
      )
    ),
    x0 = case_when(x0 == 'central' ~ 'Central',
                   x0 == 'LSEI' ~ 'LSEI')
  )

# Write clean flow data to file (if it exists, if not - create new folder)
# We will use this file for Cliffs Delta and for figure
# Check if the 'clean_results' folder exists
if (!dir.exists("final_results")) {
  dir.create("final_results")
  message("Folder 'final_results' created. Writing files into this folder.")
} else {
  message("Folder 'final_results' already exists. Writing files into existing folder.")
}
write.csv(flows_tidy, 'final_results/flows_tidy.csv')

# Create custom cliffs delta function
cliff.flow <- function(df, var1, flow) {
  df <- as.data.frame(df)
  labelnames <- levels(df[[var1]])
  cliff3 <- list()
  
  for (i in 1:length(labelnames)) {
    group1 <- paste0(labelnames[i])
    cliff <- function(df, flow, labelnames, var1, group1) {
      df <- as.data.frame(df)
      groups.to.compare <- labelnames[!labelnames %in% group1]
      cliff2 <- list()
      
      for (i in 1:length(groups.to.compare)) {
        # Filter out only group 1 and for each group to compare
        data1 <-
          df[df[[var1]] %in% c(paste0(group1), paste0(groups.to.compare[i])),]
        data1[[var1]] <- as.numeric(data1[[var1]])
        data <- data1[order(data1[[var1]]),]
        data[[var1]] <- as.factor(data[[var1]])
        
        # Run cliffs delta
        cliff <- cliff.delta(
          d = data[[flow]],
          f = data[[var1]],
          conf.level = .95,
          use.unbiased = TRUE,
          use.normal = FALSE,
          return.dm = FALSE
        )
        
        clifftable <- as.data.frame(do.call(cbind, cliff))
        clifftable$flow <- paste0(flow)
        clifftable$Group1 <- paste0(group1)
        clifftable$Group2 <- paste0(groups.to.compare[i])
        
        cliff2[[i]] <- clifftable
        
      }
      
      cliff2df <- do.call(rbind, cliff2)
      return(cliff2df)
      
    }
    
    cliff3[[i]] <-
      cliff(
        df = df,
        var1 = var1,
        flow = flow,
        group1 = group1,
        labelnames = labelnames
      )
    
  }
  
  cliff3df <- do.call(rbind, cliff3)
  return(cliff3df)
  
}

# Read in tidy flow data, if not already loaded
flows_tidy <- read.csv('final_results/flows_tidy.csv')

# Arrange flows_tidy for Cliffs Delta analysis
flows2 <- flows_tidy |>
  mutate(X = NULL,
         flow_magnitude = NULL,
         group = as.factor(paste0(x0, "_", jmp, "_", iter))) |>
  pivot_wider(names_from = "flow", values_from = "value") |>
  arrange(x0, as.numeric(jmp),
          as.numeric(iter)) |>
  group_by(x0, jmp, iter) |>
  mutate(group_id = cur_group_id(),
         group_id = as.factor(group_id))

# Cliffs Delta for each flow (comparing all groups)
c1 <-
  cliff.flow(df = flows2, var1 = "group_id", flow = "Bacteria_R")
c2 <-
  cliff.flow(df = flows2, var1 = "group_id", flow = "susPOC_Q_Oligochaeta")
c3 <-
  cliff.flow(df = flows2, var1 = "group_id", flow = "Gastropoda_EX")
c4 <-
  cliff.flow(df = flows2, var1 = "group_id", flow = "susPOC_Q_Nematoda")
c5 <-
  cliff.flow(df = flows2, var1 = "group_id", flow = "Platyhelminthes_R")
c6 <-
  cliff.flow(df = flows2, var1 = "group_id", flow = "Arachnida_R")

all <- rbind(c1, c2, c3, c4, c5, c6) |>
  rownames_to_column() |>
  filter(!grepl("upper", rowname)) |>
  filter(Group1 < Group2) |>
  mutate(
    estimate = as.numeric(estimate),
    Group1 = as.numeric(Group1),
    Group2 = as.numeric(Group2),
    flow = as.factor(flow)
  )

# Add group names back to numerics
group_names <- flows2 |>
  ungroup() |>
  select(group, group_id) |>
  unique()
lookUp1 <-
  setNames(as.character(group_names$group), group_names$group_id)
all$Group1_name <- lookUp1[as.character(all$Group1)]
all$Group2_name <- lookUp1[as.character(all$Group2)]
write.csv(all, 'final_results/flow_cliff.csv')

# Now the same comparisons, but grouped according to jump size and starting
# solution only i.e., iteration scenarios combined

cliff.flow.big <- function(df, var1, flow) {
  df <- as.data.frame(df)
  labelnames <- levels(df[[var1]])
  cliff3 <- list()
  
  for (i in 1:length(labelnames)) {
    group1 <- paste0(labelnames[i])
    
    cliff <- function(df, flow, labelnames, var1, group1) {
      df <- as.data.frame(df)
      groups.to.compare <- labelnames[!labelnames %in% group1]
      
      cliff2 <- list()
      
      for (i in 1:length(groups.to.compare)) {
        # Filter out only group 1 and for each group to compare
        data1 <-
          df[df[[var1]] %in% c(paste0(group1), paste0(groups.to.compare[i])), ]
        data1[[var1]] <- as.numeric(data1[[var1]])
        data <- data1[order(data1[[var1]]), ]
        data[[var1]] <- as.factor(data[[var1]])
        
        # Run cliffs delta
        cliff <- cliff.delta(
          d = data[[flow]],
          f = data[[var1]],
          conf.level = .95,
          use.unbiased = TRUE,
          use.normal = FALSE,
          return.dm = FALSE
        )
        
        clifftable <- as.data.frame(do.call(cbind, cliff))
        clifftable$flow <- paste0(flow)
        clifftable$Group1 <- paste0(group1)
        clifftable$Group2 <- paste0(groups.to.compare[i])
        cliff2[[i]] <- clifftable
        
      }
      
      cliff2df <- do.call(rbind, cliff2)
      return(cliff2df)
      
    }
    
    cliff3[[i]] <-
      cliff(
        df = df,
        var1 = var1,
        flow = flow,
        group1 = group1,
        labelnames = labelnames
      )
    
  }
  
  cliff3df <- do.call(rbind, cliff3)
  return(cliff3df)
  
  
}

flows_tidy <- read.csv('final_results/flows_tidy.csv')
flows_big <- flows_tidy |>
  mutate(X = NULL,
         flow_magnitude = NULL,
         group = as.factor(paste0(x0, "_", jmp, "_", iter))) |>
  pivot_wider(names_from = "flow", values_from = "value") |>
  arrange(x0, as.numeric(jmp),
          as.numeric(iter)) |>
  group_by(x0, jmp) |>
  mutate(group_id = as.factor(cur_group_id()),
         group_name = paste0(x0, "_jmp", jmp))

c1 <-
  cliff.flow.big(df = flows_big, var1 = "group_id", flow = "Bacteria_R")
c2 <-
  cliff.flow.big(df = flows_big, var1 = "group_id", flow = "susPOC_Q_Oligochaeta")
c3 <-
  cliff.flow.big(df = flows_big, var1 = "group_id", flow = "Gastropoda_EX")
c4 <-
  cliff.flow.big(df = flows_big, var1 = "group_id", flow = "susPOC_Q_Nematoda")
c5 <-
  cliff.flow.big(df = flows_big, var1 = "group_id", flow = "Platyhelminthes_R")
c6 <-
  cliff.flow.big(df = flows_big, var1 = "group_id", flow = "Arachnida_R")

all <- rbind(c1, c2, c3, c4, c5, c6) |>
  rownames_to_column() |>
  filter(!grepl("upper", rowname)) |>
  filter(Group1 < Group2) |>
  mutate(
    estimate = as.numeric(estimate),
    Group1 = as.numeric(Group1),
    Group2 = as.numeric(Group2),
    flow = as.factor(flow)
  )

# Add group names back to numerics
group_names <- flows_big |>
  ungroup() |>
  select(group_id, group_name) |>
  unique()
lookUp1 <-
  setNames(as.character(group_names$group_name), group_names$group_id)
all$Group1_name <- lookUp1[as.character(all$Group1)]
all$Group2_name <- lookUp1[as.character(all$Group2)]
write.csv(all, 'final_results/flow_cliff_iters_combined.csv')

# Now the same, but for ENA metrics
# Can use the same custom cliff.flow function
# Read in all ENA metrics from folder 'ena', save as one table
myfiles <-
  list.files(
    path = paste0(getwd(), "/ena"),
    pattern = "*.csv",
    full.names = TRUE
  )
ena <-
  sapply(myfiles, read.csv, simplify = FALSE) %>% bind_rows(.id = "id")

# Tidy data & extract necessary ENA metrics
ena_tidy <- ena |>
  select(jmp,
         iter,
         x0,
         Iteration,
         TST,
         FCI,
         DH,
         AMI,
         H,
         robustness) |>
  pivot_longer(TST:robustness,
               names_to = 'ena',
               values_to = 'value') |>
  mutate(
    jmp = as.factor(jmp),
    iter = as.factor (iter),
    x0 = as.factor(x0),
    Iteration = as.factor(Iteration),
    ena = as.factor(ena),
    x0 = case_when(x0 == 'central' ~ 'Central',
                   x0 == 'LSEI' ~ 'LSEI')
  )

write.csv(ena_tidy, 'final_results/ena_tidy.csv')

cliff.flow <- function(df, var1, flow) {
  df <- as.data.frame(df)
  labelnames <- levels(df[[var1]])
  cliff3 <- list()
  
  for (i in 1:length(labelnames)) {
    group1 <- paste0(labelnames[i])
    cliff <- function(df, flow, labelnames, var1, group1) {
      df <- as.data.frame(df)
      groups.to.compare <- labelnames[!labelnames %in% group1]
      cliff2 <- list()
      
      for (i in 1:length(groups.to.compare)) {
        # Filter out only group 1 and for each group to compare
        data1 <-
          df[df[[var1]] %in% c(paste0(group1), paste0(groups.to.compare[i])),]
        data1[[var1]] <- as.numeric(data1[[var1]])
        data <- data1[order(data1[[var1]]),]
        data[[var1]] <- as.factor(data[[var1]])
        
        # Run cliffs delta
        cliff <- cliff.delta(
          d = data[[flow]],
          f = data[[var1]],
          conf.level = .95,
          use.unbiased = TRUE,
          use.normal = FALSE,
          return.dm = FALSE
        )
        
        clifftable <- as.data.frame(do.call(cbind, cliff))
        clifftable$flow <- paste0(flow)
        clifftable$Group1 <- paste0(group1)
        clifftable$Group2 <- paste0(groups.to.compare[i])
        
        cliff2[[i]] <- clifftable
        
      }
      
      cliff2df <- do.call(rbind, cliff2)
      return(cliff2df)
      
    }
    
    cliff3[[i]] <-
      cliff(
        df = df,
        var1 = var1,
        flow = flow,
        group1 = group1,
        labelnames = labelnames
      )
    
  }
  
  cliff3df <- do.call(rbind, cliff3)
  return(cliff3df)
  
}

# Read in tidy ena data
ena_tidy <- read.csv('final_results/ena_tidy.csv')

# Arrange ena_tidy for Cliffs Delta analysis
ena2 <- ena_tidy |>
  mutate(group = as.factor(paste0(x0, "_", jmp, "_", iter))) |>
  pivot_wider(names_from = "ena", values_from = "value") |>
  arrange(x0, as.numeric(jmp),
          as.numeric(iter)) |>
  group_by(x0, jmp, iter) |>
  mutate(group_id = cur_group_id(),
         group_id = as.factor(group_id))

# Cliffs Delta for each flow (comparing all groups)
c1 <- cliff.flow(df = ena2, var1 = "group_id", flow = "TST")
c2 <- cliff.flow(df = ena2, var1 = "group_id", flow = "FCI")
c3 <- cliff.flow(df = ena2, var1 = "group_id", flow = "DH")
c4 <- cliff.flow(df = ena2, var1 = "group_id", flow = "AMI")
c5 <- cliff.flow(df = ena2, var1 = "group_id", flow = "H")
c6 <- cliff.flow(df = ena2, var1 = "group_id", flow = "robustness")

all <- rbind(c1, c2, c3, c4, c5, c6) |>
  rownames_to_column() |>
  filter(!grepl("upper", rowname)) |>
  mutate(
    estimate = as.numeric(estimate),
    Group1 = as.numeric(Group1),
    Group2 = as.numeric(Group2),
    flow = as.factor(flow)
  )

all <- all |>
  filter(Group1 < Group2)

# Add group names back to numerics
group_names <- ena2 |>
  ungroup() |>
  select(group, group_id) |>
  unique()
lookUp1 <-
  setNames(as.character(group_names$group), group_names$group_id)
all$Group1_name <- lookUp1[as.character(all$Group1)]
all$Group2_name <- lookUp1[as.character(all$Group2)]

write.csv(all, 'final_results/ena_cliff.csv')

# ENA grouped according to jump size and starting
# solution only i.e., iteration scenarios combined

cliff.flow.big <- function(df, var1, flow) {
  df <- as.data.frame(df)
  labelnames <- levels(df[[var1]])
  cliff3 <- list()
  
  for (i in 1:length(labelnames)) {
    group1 <- paste0(labelnames[i])
    
    cliff <- function(df, flow, labelnames, var1, group1) {
      df <- as.data.frame(df)
      groups.to.compare <- labelnames[!labelnames %in% group1]
      
      cliff2 <- list()
      
      for (i in 1:length(groups.to.compare)) {
        # Filter out only group 1 and for each group to compare
        data1 <-
          df[df[[var1]] %in% c(paste0(group1), paste0(groups.to.compare[i])), ]
        data1[[var1]] <- as.numeric(data1[[var1]])
        data <- data1[order(data1[[var1]]), ]
        data[[var1]] <- as.factor(data[[var1]])
        
        # Run cliffs delta
        cliff <- cliff.delta(
          d = data[[flow]],
          f = data[[var1]],
          conf.level = .95,
          use.unbiased = TRUE,
          use.normal = FALSE,
          return.dm = FALSE
        )
        
        clifftable <- as.data.frame(do.call(cbind, cliff))
        clifftable$flow <- paste0(flow)
        clifftable$Group1 <- paste0(group1)
        clifftable$Group2 <- paste0(groups.to.compare[i])
        cliff2[[i]] <- clifftable
        
      }
      
      cliff2df <- do.call(rbind, cliff2)
      return(cliff2df)
      
    }
    
    cliff3[[i]] <-
      cliff(
        df = df,
        var1 = var1,
        flow = flow,
        group1 = group1,
        labelnames = labelnames
      )
    
  }
  
  cliff3df <- do.call(rbind, cliff3)
  return(cliff3df)
  
  
}

ena_tidy <- read.csv('final_results/ena_tidy.csv')
ena_big <- ena_tidy |>
  mutate(group = as.factor(paste0(x0, "_", jmp, "_", iter))) |>
  pivot_wider(names_from = "ena", values_from = "value") |>
  arrange(x0, as.numeric(jmp), as.numeric(iter)) |>
  group_by(x0, jmp) |>
  mutate(group_id = as.factor(cur_group_id()),
         group_name = paste0(x0, "_jmp", jmp))

c1 <- cliff.flow.big(df = ena_big, var1 = "group_id", flow = "TST")
c2 <- cliff.flow.big(df = ena_big, var1 = "group_id", flow = "FCI")
c3 <- cliff.flow.big(df = ena_big, var1 = "group_id", flow = "DH")
c4 <- cliff.flow.big(df = ena_big, var1 = "group_id", flow = "AMI")
c5 <- cliff.flow.big(df = ena_big, var1 = "group_id", flow = "H")
c6 <-
  cliff.flow.big(df = ena_big, var1 = "group_id", flow = "robustness")

all <- rbind(c1, c2, c3, c4, c5, c6) |>
  rownames_to_column() |>
  filter(!grepl("upper", rowname)) |>
  filter(Group1 < Group2) |>
  mutate(
    estimate = as.numeric(estimate),
    Group1 = as.numeric(Group1),
    Group2 = as.numeric(Group2),
    flow = as.factor(flow)
  )

# Add group names back to numerics
group_names <- flows_big |>
  ungroup() |>
  select(group_id, group_name) |>
  unique()
lookUp1 <-
  setNames(as.character(group_names$group_name), group_names$group_id)
all$Group1_name <- lookUp1[as.character(all$Group1)]
all$Group2_name <- lookUp1[as.character(all$Group2)]
write.csv(all, 'final_results/ena_cliff_iters_combined.csv')

###############################################################################
# 2.3 Algorithm quality assessments
###############################################################################

# MCMC diagnostics already calculated in step 2.2.2
# MCMC diagnostics only presented on BOUNDED flows

# Read in all mcmc_diags for each flow
myfiles <-
  list.files(
    path = paste0(getwd(), "/mcmc_diags"),
    pattern = "*.csv",
    full.names = TRUE
  )
diags <-
  sapply(myfiles, read.csv, simplify = FALSE) %>% bind_rows(.id = "id")

# Remove unbounded flows (1e^30) per group
# Unbounded is where max flow values = 1e^30
bounded <- diags |>
  filter(Summary.Xranges_max != 1e+30) # 200 bounded flows

# Visual assessment with traceplots
# Can do for any flow

flows <- read.csv('final_results/flows_tidy.csv')
flow_plot <- flows |> filter(flow == 'Bacteria_R') |> droplevels()

traceplot <- function(data, flow) {
  data %>%
    mutate(
      iter = as.factor(iter),
      jmp = as.factor(jmp),
      x0 = as.factor(x0),
      flow = as.factor(flow)
    ) |>
    ggplot(aes(x = Iteration, y = value, colour = x0)) +
    geom_line(linewidth = 0.5, position = position_jitter(w = 0.02, h =
                                                            0)) +
    facet_grid2(cols = vars(iter),
                rows = vars(jmp),
                scales = "free") +
    theme_classic2() +
    guides(color = guide_legend(title = "Starting Point")) +
    xlab("Iteration") +
    ylab(bquote(.(paste0(flow)) ~  ~ (mgC ~ m ^ -2 ~ d ^ -1))) +
    theme(legend.position = "top")
  
}

traceplot(data = flow_plot, flow = 'Bacteria_R')



###############################################################################
# 2.3.1 Statistical analysis

# Clean up diags_data for analysis
diags_clean <- bounded |>
  rename(flow = X) |>
  mutate(
    cov = Summary.SD / Summary.Mean,
    cov_test =
      case_when(cov < 1 ~ 1, cov > 1 ~ 0),
    geweke_test = case_when(
      Geweke.geweke.diag > 1.96 | Geweke.geweke.diag < -1.96 ~ 0,
      Geweke.geweke.diag > -1.96 &
        Geweke.geweke.diag < 1.96 ~ 1
    ),
    HW_test = case_when(
      Heidelberger.Welch.HW.Stationarity.Test == "Passed" &
        Heidelberger.Welch.HW.Halfwidth.Test == "Passed" ~ 1,
      Heidelberger.Welch.HW.Stationarity.Test == "Failed" |
        Heidelberger.Welch.HW.Halfwidth.Test == "Failed" ~ 0
    ),
    I_test = case_when(
      Raftery.Lewis.RL.Dependence.Factor < 5 ~ 1,
      Raftery.Lewis.RL.Dependence.Factor > 5 ~ 0,
      Raftery.Lewis.RL.Dependence.Factor == 5 ~ 0
    ),
    AC_test = case_when(
      Autocorrelation.Lag.5 / Autocorrelation.Lag.1 <
        Autocorrelation.Lag.50 / Autocorrelation.Lag.1 |
        Autocorrelation.Lag.1 < 0 |
        Autocorrelation.Lag.5 < 0 |
        Autocorrelation.Lag.10 < 0 |
        Autocorrelation.Lag.50 < 0 ~ 1,
      Autocorrelation.Lag.5 / Autocorrelation.Lag.1 >
        Autocorrelation.Lag.50 / Autocorrelation.Lag.1 ~ 0
    ),
    N_test = case_when(
      Raftery.Lewis.RL.Required.Sample.Size > iter ~ 0,
      Raftery.Lewis.RL.Required.Sample.Size < iter ~ 1
    ),
    ess_test = case_when(
      effective.sample.size / iter > 0.1 ~ 1,
      effective.sample.size / iter < 0.1 ~ 0
    ),
    group = paste0(x0, "_", jmp, "_", iter)
  ) |>
  pivot_longer(cov_test:ess_test, names_to = 'diag_test', values_to = 'test_pass')





# Create table for analysis
summary <- diags_clean |>
  group_by(group, diag_test) |>
  summarise(passed_n = sum(test_pass),
            failed_n = 200 - passed_n) |>
  ungroup()

# Fisher and post hoc tests per diagnostic
# COV
fish <-
  summary %>% filter(diag_test == "cov_test") %>% droplevels() %>% mutate(diag_test = NULL)
f1_cov <-
  as.data.frame(do.call(cbind, fisher.test(table(
    fish$passed_n, fish$failed_n
  ))))
fish2 <- fish |> column_to_rownames(var = "group")
f2_cov <-
  pairwise_fisher_test(as.matrix(fish2), p.adjust.method = "fdr")
f1_cov$diag <- "cov"
f2_cov$diag <- "cov"

# HW
fish <-
  summary %>% filter(diag_test == "HW_test") %>% droplevels() %>% mutate(diag_test = NULL)
f1_hw <-
  as.data.frame(do.call(cbind, fisher.test(table(
    fish$passed_n, fish$failed_n
  ))))
fish2 <- fish |> column_to_rownames(var = "group")
f2_hw <-
  pairwise_fisher_test(as.matrix(fish2), p.adjust.method = "fdr")
f1_hw$diag <- "HW"
f2_hw$diag <- "HW"

# Geweke
fish <-
  summary %>% filter(diag_test == "geweke_test") %>% droplevels() %>% mutate(diag_test = NULL)
f1_g <-
  as.data.frame(do.call(cbind, fisher.test(table(
    fish$passed_n, fish$failed_n
  ))))
fish2 <- fish |> column_to_rownames(var = "group")
f2_g <-
  pairwise_fisher_test(as.matrix(fish2), p.adjust.method = "fdr")
f1_g$diag <- "geweke"
f2_g$diag <- "geweke"

# I
fish <-
  summary %>% filter(diag_test == "I_test") %>% droplevels() %>% mutate(diag_test = NULL)
f1_i <-
  as.data.frame(do.call(cbind, fisher.test(table(
    fish$passed_n, fish$failed_n
  ))))
fish2 <- fish |> column_to_rownames(var = "group")
f2_i <-
  pairwise_fisher_test(as.matrix(fish2), p.adjust.method = "fdr")
f1_i$diag <- "I"
f2_i$diag <- "I"

# Autocorrelation
fish <-
  summary %>% filter(diag_test == "AC_test") %>% droplevels() %>% mutate(diag_test = NULL)
f1_ac <-
  as.data.frame(do.call(cbind, fisher.test(table(
    fish$passed_n, fish$failed_n
  ))))
fish2 <- fish |> column_to_rownames(var = "group")
f2_ac <-
  pairwise_fisher_test(as.matrix(fish2), p.adjust.method = "fdr")
f1_ac$diag <- "autocorrelation"
f2_ac$diag <- "autocorrelation"

# ESS
fish <-
  summary %>% filter(diag_test == "ess_test") %>% droplevels() %>% mutate(diag_test = NULL)
f1_ess <-
  as.data.frame(do.call(cbind, fisher.test(table(
    fish$passed_n, fish$failed_n
  ))))
fish2 <- fish |> column_to_rownames(var = "group")
f2_ess <-
  pairwise_fisher_test(as.matrix(fish2), p.adjust.method = "fdr")
f1_ess$diag <- "ESS"
f2_ess$diag <- "ESS"

# N
fish <-
  summary %>% filter(diag_test == "N_test") %>% droplevels() %>% mutate(diag_test = NULL)
f1_n <-
  as.data.frame(do.call(cbind, fisher.test(table(
    fish$passed_n, fish$failed_n
  ))))
fish2 <- fish |> column_to_rownames(var = "group")
f2_n <-
  pairwise_fisher_test(as.matrix(fish2), p.adjust.method = "fdr")
f1_n$diag <- "N"
f2_n$diag <- "N"


# Write all tables
# Fishers Exact results
x <- rbind(f1_cov, f1_hw, f1_g, f1_i, f1_ac, f1_ess, f1_n)
write.csv(x, 'final_results/fishers_exact_mcmcdiags.csv')

# Posthoc
y <- rbind(f2_cov, f2_hw, f2_g, f2_i, f2_ac, f2_ess, f2_n)
write.csv(y, 'final_results/fishers_pairwise_mcmcdiags.csv')
