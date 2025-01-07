# -------------------------------------------------------------------
# R Script for Calculating tNICE and fNICE Metrics
# -------------------------------------------------------------------

# Load Required Packages
library(cluster)
library(ape)
library(geometry)
library(StatMatch)
library(reshape2)

# Source External Functions
source('./sExtinct/R/OLE.R')
source('./sExtinct/R/OLE.fun.R')

# -------------------------------------------------------------------
# Function: complementary_volume
# -------------------------------------------------------------------

#' Compute Complementary Volumes Between Consecutive Time Points
#'
#' This function calculates the unique volumes occupied by species at each pair of consecutive time points in trait space.
#'
#' @param trait_mat A matrix of species traits with species as rows.
#' @param comm_mat A community matrix with time points as rows and species as columns.
#' @param nb.axes Integer. Number of principal coordinates axes to use. Default is 3.
#' @return A list of data frames containing the unique volumes lost and gained between each pair of years.
complementary_volume <- function(trait_mat, comm_mat, nb.axes = 3) {
  
  require(cluster)
  require(ape)
  require(geometry)
  
  # Compute Gower distance matrix
  dist_mat <- StatMatch::gower.dist(data.x = trait_mat)
  
  # Compute PCoA using the Gower distance matrix
  pcoa <- ape::pcoa(dist_mat, correction = 'cailliez')
  
  # Compute the complementary volumes for each year
  years <- nrow(comm_mat)
  vol <- list()
  for (i in 1:(years - 1)) {
    ch1 <- convhulln(pcoa$vectors[which(comm_mat[i, ] > 0), 1:nb.axes], option = 'FA')
    ch2 <- convhulln(pcoa$vectors[which(comm_mat[i + 1, ] > 0), 1:nb.axes], option = 'FA')
    
    # Intersection
    common_vol <- intersectn(ch1$hull, ch2$hull)
    
    # Unique part of each volume
    if (all(which(comm_mat[i, ] > 0) %in% which(comm_mat[i + 1, ] > 0))) {
      vol1 <- 0
    } else {
      vol1 <- ch1$vol - common_vol$ch$vol
    }
    if (all(which(comm_mat[i + 1, ] > 0) %in% which(comm_mat[i, ] > 0))) {
      vol2 <- 0
    } else {
      vol2 <- ch2$vol - common_vol$ch$vol
    }
    
    # Store the results
    vol[[paste0(i, "-", i - 1)]] <- data.frame(
      year1 = as.numeric(rownames(comm_mat)[i]),
      year2 = as.numeric(rownames(comm_mat)[i + 1]),
      lost_vol = abs(vol1),
      gained_vol = abs(vol2)
    )
  }
  return(vol)
}

# -------------------------------------------------------------------
# Function: compute_fNICE
# -------------------------------------------------------------------

#' Compute Functional Net Immigration-Extinction (fNICE) Metric
#'
#' This function calculates the functional Net Imbalance between Colonisation-Extinction (fNICE) metric 
#' for a given community dataset over time, based on species traits. 
#' It computes the functional space occupied by species at each time point, 
#' calculates the unique volumes lost and gained, and compares them against null models.
#'
#' @param samp A list containing the following elements:
#'   \itemize{
#'     \item \code{com_matrix}: A community matrix with rows as time points (e.g., years) and columns as species names. 
#'                              Values are species abundances or presence (counts) at each time point.
#'     \item \code{traits}: A data frame of species traits with species as rows and trait variables as columns.
#'     \item \code{infos}: Metadata about the study (e.g., study ID, location, year range).
#'   }
#' @param nb.null Integer. The number of null model iterations to perform. Default is 49.
#' @param traits Optional. A character vector specifying the columns of traits to include. If \code{NULL} (default), all traits except the first column are used.
#' @param all_species Logical. Whether to use all species in the traits data. If \code{FALSE} (default), only species present in \code{com_matrix} are used.
#' @param gower_weights Optional. A numeric vector of weights for the traits when computing Gower distance. If \code{NULL} (default), all traits are weighted equally.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{eig}: A vector of eigenvalues from the PCoA, indicating the percentage of variance explained by each axis.
#'     \item \code{fNICE}: A data frame with observed and null model volumes of functional space lost and gained between consecutive time points, and standardized effect sizes (SES).
#'     \item \code{fNICE_cum}: A data frame with cumulative observed and SES values, and the fNICE metric over time.
#'   }
compute_fNICE <- function(samp, nb.null = 49, traits = NULL, all_species = FALSE, gower_weights = NULL) {
  
  print(samp$infos[, 1])
  
  samp_com <- samp$com_matrix
  # Replace spaces with periods in the species names
  colnames(samp_com) <- gsub(' ', '.', colnames(samp_com))
  
  # Prepare trait data
  if (!is.null(traits)) {
    samp_traits <- samp$traits[, traits]
  } else {
    samp_traits <- samp$traits[, -1]
  }
  
  # Formatting species names in the trait table
  if ('Species' %in% colnames(samp$traits)) {
    # Replace spaces with periods in the species names
    rownames(samp_traits) <- gsub(' ', '.', samp$traits$Species, fixed = TRUE)
  } else {
    # Use 'sci_name' column to format row names for samp_traits
    rownames(samp_traits) <- gsub(' ', '.', samp$traits$sci_name, fixed = TRUE)
  }
  
  # Remove species with only NA
  species_to_remove <- rownames(samp_traits)[apply(samp_traits, 1, function(x) all(is.na(x)))]
  
  # Remove these species from both samp_com and samp_traits
  samp_com <- samp_com[, !(colnames(samp_com) %in% species_to_remove)]
  samp_traits <- samp_traits[!rownames(samp_traits) %in% species_to_remove, ]
  
  # Print the number of species removed for debugging purposes
  print(paste(length(species_to_remove), "species removed due to missing trait values"))
  
  # Remove species never sampled
  if (any(apply(samp_com, 2, sum) == 0)) {
    samp_com <- samp_com[, -which(apply(samp_com, 2, sum) == 0)]
  }
  
  # Remove samples with less than 3 species
  if (length(which(apply(ifelse(samp_com > 0, 1, 0), 1, sum) < 4)) > 0) {
    samp_com <- samp_com[-which(apply(ifelse(samp_com > 0, 1, 0), 1, sum) < 4), ]
  }
  
  # Ensure there are enough species and samples for analysis
  if (ncol(samp_com) < 5 | nrow(samp_com) < 5) {
    return(NULL)
  }
  
  # Assign column names if missing (e.g., in simulations)
  if (!all_species) {
    if (is.null(colnames(samp_com))) {
      colnames(samp_com) <- paste('sp', 1:ncol(samp_com), sep = '.')
    }
  }
  
  samp_traits <- samp_traits[rownames(samp_traits) %in% colnames(samp_com), ]
  samp_traits <- samp_traits[match(colnames(samp_com), rownames(samp_traits)), ]
  
  # Compute functional distance matrix using Gower distance
  gower_dist <- StatMatch::gower.dist(data.x = samp_traits, var.weights = gower_weights)
  
  # Compute PCoA on the Gower distance matrix
  samp_pcoa <- ape::pcoa(D = gower_dist, correction = 'cailliez')
  
  # Calculate eigenvalues and extract relevant PCoA vectors
  eig <- round(100 * (samp_pcoa$values[1:min(c(3, ncol(samp_pcoa$vectors))), 1] / sum(samp_pcoa$values$Eigenvalues)), 1)
  samp_pcoa_vectors <- samp_pcoa$vectors[, 1:min(c(3, ncol(samp_pcoa$vectors)))]
  
  # Update the trait matrix with PCoA vectors
  samp_traits_pcoa <- samp_pcoa_vectors
  
  # Calculate complementary volumes
  comp_vols <- try(complementary_volume(trait_mat = samp_traits_pcoa, comm_mat = samp_com), silent = FALSE)
  
  if (inherits(comp_vols, "try-error")) {
    return(NULL)
  } else {
    comp_vols <- do.call(rbind, as.list(comp_vols))
    
    # Run Null Models
    null_out <- list()
    k <- 1
    while (length(null_out) != nb.null) {
      random_traits <- samp_traits_pcoa[sample(1:nrow(samp_traits_pcoa)), ]
      rownames(random_traits) <- rownames(samp_traits_pcoa)
      
      null_k <- try(complementary_volume(trait_mat = random_traits, comm_mat = samp_com), silent = TRUE)
      if (inherits(null_k, "try-error")) {
        next
      } else {
        null_out[[k]] <- do.call(rbind, as.list(null_k))
        k <- k + 1
      }
    }
    
    null_out <- do.call(rbind, null_out)
    null_mean <- do.call(rbind, as.list(by(data = null_out, INDICES = null_out[, c('year1', 'year2')], FUN = function(x) apply(x, 2, mean, na.rm = TRUE))))
    null_sd <- do.call(rbind, as.list(by(data = null_out, INDICES = null_out[, c('year1', 'year2')], FUN = function(x) apply(x, 2, sd, na.rm = TRUE))))
    
    # Formatting the output
    fNICE <- data.frame(
      year1 = comp_vols$year1, 
      year2 = comp_vols$year2, 
      obs.loss = comp_vols$lost_vol, 
      obs.gain = comp_vols$gained_vol, 
      mean.loss = null_mean[, 3], 
      mean.gain = null_mean[, 4], 
      sd.loss = null_sd[, 3], 
      sd.gain = null_sd[, 4], 
      ses.loss = ((comp_vols$lost_vol - null_mean[, 3]) / null_sd[, 3]), 
      ses.gain = ((comp_vols$gained_vol - null_mean[, 4]) / null_sd[, 4])
    )
    
    # Handling NA and infinite values for cumulative sum
    fNICE$ses.loss <- ifelse(is.na(fNICE$ses.loss) | is.infinite(fNICE$ses.loss), 0, fNICE$ses.loss)
    fNICE$ses.gain <- ifelse(is.na(fNICE$ses.gain) | is.infinite(fNICE$ses.gain), 0, fNICE$ses.gain)
    
    # Computing the fNICE metric
    fNICE_cum <- data.frame(
      year = fNICE$year2, 
      cum.obs.loss = cumsum(fNICE$obs.loss), 
      cum.obs.gain = cumsum(fNICE$obs.gain), 
      cum.ses.loss = cumsum(fNICE$ses.loss), 
      cum.ses.gain = cumsum(fNICE$ses.gain), 
      fNICE.obs = (cumsum(fNICE$obs.gain) - cumsum(fNICE$obs.loss)) / (cumsum(fNICE$obs.gain) + cumsum(fNICE$obs.loss)),
      fNICE.ses = (cumsum(fNICE$ses.gain) - cumsum(fNICE$ses.loss)) / (cumsum(fNICE$ses.gain) + cumsum(fNICE$ses.loss))
    )
    
    out_FD <- list(
      eig = eig, 
      fNICE = fNICE, 
      fNICE_cum = fNICE_cum
    )
    
    return(out_FD)
  }
}

# -------------------------------------------------------------------
# Function: compute_tNICE
# -------------------------------------------------------------------

#' Compute Taxonomic Net Immigration-Extinction (tNICE) Metric
#'
#' This function calculates the taxonomic Net Immigration-Extinction (tNICE) metric 
#' for a given community dataset over time. 
#' It estimates the timing of species colonizations and extinctions, 
#' counts events per year, and computes the cumulative tNICE metric over time.
#'
#' @param samp A list containing the following elements:
#'   \itemize{
#'     \item \code{com_matrix}: A community matrix with rows as time points (e.g., years) and columns as species names. Values are species abundances or presence (counts) at each time point.
#'     \item \code{infos}: Metadata about the study (e.g., study ID, location, year range).
#'   }
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item \code{YEAR}: Time points (e.g., years).
#'     \item \code{COL}: Cumulative number of species colonizations up to each year.
#'     \item \code{EXT}: Cumulative number of species extinctions up to each year.
#'     \item \code{NICE}: The tNICE metric calculated as \eqn{(Cumulative COL - Cumulative EXT) / (Cumulative COL + Cumulative EXT)}.
#'     \item \code{SR}: Species richness at each time point.
#'   }
compute_tNICE <- function(samp) {
  
  print(samp$infos[, 1])
  
  samp.com <- samp$com_matrix
  
  # Remove species never sampled
  if (any(apply(samp.com, 2, sum) == 0)) {
    samp.com <- samp.com[, -which(apply(samp.com, 2, sum) == 0)]
  }
  
  # Ensure there are enough species and samples for analysis
  if (ncol(samp.com) < 5 | nrow(samp.com) < 5) {
    return(NULL)
  }
  
  # Assign column names if missing
  if (is.null(colnames(samp.com))) {
    colnames(samp.com) <- paste('sp', 1:ncol(samp.com), sep = '.')
  }
  
  # Time points
  years <- as.numeric(rownames(samp.com))
  
  # Extinction timing estimation
  est.EXT <- apply(samp.com, 2, function(x) OLE(sightingdata = cbind(years, x), alpha = .05))
  est.EXT <- data.frame(Sp = colnames(samp.com), do.call(rbind, lapply(est.EXT, function(x) cbind(x[1], x[2], x[3]))))
  colnames(est.EXT) <- c('Sp', 'Est.E', 'low.E', 'up.E')
  
  # Colonization timing estimation
  samp.COL <- cbind(abs(years - max(years)), samp.com)
  samp.COL <- samp.COL[order(samp.COL[, 1]), ]
  
  est.COL <- apply(samp.COL[, -1], 2, function(x) lapply(OLE(
    sightingdata = cbind(years, x), alpha = .05),
    function(z) abs(z - max(samp.com[, 1]))))
  est.COL <- data.frame(Sp = colnames(samp.COL[, -1]), do.call(rbind, lapply(est.COL, function(x) cbind(x[1], x[2], x[3]))))
  colnames(est.COL) <- c('Sp', 'Est.C', 'low.C', 'up.C')
  
  # Combine extinction and colonization estimates
  sp.infos <- merge(est.COL, est.EXT, by = 'Sp', all = TRUE)
  sp.infos[, -1] <- apply(sp.infos[, -1], 2, function(x) as.numeric(as.character(substr(x, 1, 4))))
  
  # Counting events for each year
  COL <- apply(table(sp.infos$Sp, sp.infos$Est.C), 2, sum)
  COL <- data.frame(YEAR = as.numeric(names(COL)), COL)
  EXT <- apply(table(sp.infos$Sp, sp.infos$Est.E), 2, sum)
  EXT <- data.frame(YEAR = as.numeric(names(EXT)), EXT)
  events <- merge(COL, EXT, all = TRUE)
  
  # Remove years outside the sampling window
  events <- merge(data.frame(YEAR = years), events, all.x = TRUE, all.y = FALSE)
  
  # Replace NA with zero
  events[is.na(events)] <- 0
  
  # Cumulative events over time
  events.cum <- cbind(events$YEAR, as.matrix(apply(events[, -1], 2, cumsum)))
  colnames(events.cum)[1] <- 'YEAR'
  
  # Compute tNICE
  NICE <- (events.cum[, 'COL'] - events.cum[, 'EXT']) / (events.cum[, 'COL'] + events.cum[, 'EXT'])
  
  NICE <- data.frame(YEAR = events.cum[, 1], NICE = NICE)
  NICE[NICE == Inf] <- NA
  NICE[NICE == -Inf] <- NA
  NICE[is.nan(NICE[, 2]), 2] <- NA
  
  out.TD <- merge(events.cum, NICE)
  SR <- data.frame(YEAR = years, SR = apply(ifelse(samp.com > 0, 1, 0), 1, sum))
  out.TD <- merge(out.TD, SR)
  
  return(out.TD)
}

# -------------------------------------------------------------------
# Example Usage
# -------------------------------------------------------------------

# 1. Load Your Data
# Replace 'your_file.RData' with the path to your data file
# load('your_file.RData')

# For this example, we'll assume 'data_list' is following this structure: 

# Example data (Replace with your actual data)
# data_list <- list(
#   site1 = list(
#     com_matrix = your_community_matrix,
#     traits = your_traits_data,
#     infos = your_metadata),
#   site2 = list(
#     com_matrix = your_community_matrix,
#     traits = your_traits_data,
#     infos = your_metadata)
# )
#

# 2. Run the tNICE Function
# tNICE_result <- lapply(data_list, function(x) compute_tNICE)
# save(tNICE_result, file = 'tNICE_result.RData')

# 3. Run the fNICE Function
# fNICE_results <- lapply(input, function(x) {
#   compute_fNICE(
#     samp = x,
#     nb.null = 49,             # Adjust as needed
#     traits = NULL,            # Use all traits
#     all_species = FALSE,      # Use only species in com_matrix
#     gower_weights = NULL      # Equal weighting for traits
#   )
# })
# save(fNICE_result, file = 'fNICE_result.RData')

# -------------------------------------------------------------------
# Notes:
# - Ensure that the species names in 'com_matrix' and 'traits' match.
# - Replace placeholder data and functions with your actual data and implementations.
# - Install any required packages that are not already installed.
# -------------------------------------------------------------------
