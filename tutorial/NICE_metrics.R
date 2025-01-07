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
#' This function calculates the unique volumes occupied by species at 
#' each pair of consecutive time points in trait space.
#'
#' @param trait_mat A matrix (or data frame) of species traits, 
#'   with species as rows and traits as columns.
#' @param comm_mat A community matrix with time points as rows (named by year) 
#'   and species as columns. Values indicate species abundances (or presence) 
#'   at each time point.
#' @param nb.axes Integer. Number of principal coordinate axes to use from the 
#'   PCoA output. Default is 3.
#'
#' @return A list of data frames containing the unique volumes lost (extinctions) 
#'   and gained (colonizations) between each pair of consecutive time points.
#'
#' @details
#' \itemize{
#'   \item \code{geometry::convhulln} is used to compute the convex hull volume 
#'         for the trait space occupied by species in each time point.
#'   \item The intersection (\code{common_vol}) is subtracted to find 
#'         the unique volumes.
#'   \item You may need to handle degenerate cases (e.g., <3 points in trait space).
#' }
#'
#' @examples
#' # Assume you have a trait matrix 'trait_mat' and a community matrix 'comm_mat'
#' # vol_list <- complementary_volume(trait_mat, comm_mat, nb.axes = 3)
#' # str(vol_list)
#'
complementary_volume <- function(trait_mat, comm_mat, nb.axes = 3) {
  
  require(cluster)
  require(ape)
  require(geometry)
  
  # 1. Compute Gower distance on the trait matrix
  dist_mat <- StatMatch::gower.dist(data.x = trait_mat)
  
  # 2. Perform PCoA using the Gower distance matrix
  pcoa <- ape::pcoa(dist_mat, correction = 'cailliez')
  
  # 3. For each consecutive pair of rows in comm_mat (i.e., consecutive years),
  #    compute the convex hull volume in trait space and the intersection.
  years <- nrow(comm_mat)
  vol_list <- list()
  
  for (i in 1:(years - 1)) {
    
    # Species present in time i
    sp_i  <- which(comm_mat[i, ]     > 0)
    # Species present in time i+1
    sp_ip <- which(comm_mat[i + 1, ] > 0)
    
    # If either time point has fewer than 3 species, convhulln can fail
    # We skip or handle gracefully
    if (length(sp_i) < 3 || length(sp_ip) < 3) {
      vol_list[[paste0(i, "-", i+1)]] <- data.frame(
        year1     = as.numeric(rownames(comm_mat)[i]),
        year2     = as.numeric(rownames(comm_mat)[i + 1]),
        lost_vol  = NA,
        gained_vol= NA
      )
      next
    }
    
    # Compute convex hulls in trait space
    ch1 <- convhulln(pcoa$vectors[sp_i, 1:nb.axes],  option = 'FA')
    ch2 <- convhulln(pcoa$vectors[sp_ip, 1:nb.axes], option = 'FA')
    
    # Intersection volume
    # geometry::intersectn() expects hull objects, not just volumes.
    # ch1$hull and ch2$hull are the facets. 
    common_vol <- intersectn(ch1$hull, ch2$hull)
    
    # Unique volumes
    vol1 <- ch1$vol - common_vol$ch$vol  # volume unique to time i
    vol2 <- ch2$vol - common_vol$ch$vol  # volume unique to time i+1
    
    # If all species at time i are also present at time i+1, lost_vol=0, etc.
    if (all(sp_i %in% sp_ip)) vol1 <- 0
    if (all(sp_ip %in% sp_i)) vol2 <- 0
    
    # Store results
    vol_list[[paste0(i, "-", i+1)]] <- data.frame(
      year1     = as.numeric(rownames(comm_mat)[i]),
      year2     = as.numeric(rownames(comm_mat)[i + 1]),
      lost_vol  = abs(vol1),
      gained_vol= abs(vol2)
    )
  }
  
  return(vol_list)
}

# -------------------------------------------------------------------
# Function: compute_fNICE
# -------------------------------------------------------------------

#' Compute Functional Net Immigration-Extinction (fNICE) Metric
#'
#' This function calculates the functional Net Imbalance between Colonisation-Extinction (fNICE) metric
#' for a given community dataset over time, based on species traits. 
#' It computes the functional space occupied by species at each time point 
#' and calculates the unique volumes lost and gained between consecutive time points.
#'
#' @param com A community matrix with rows as time points (e.g., years) and columns 
#'   as species names. Values are species abundances or presence (counts) at each time point.
#' @param infos A data frame or list containing metadata about the study 
#'   (e.g., study ID, location, year range). Not strictly required for calculations,
#'   but useful for logging or QA.
#' @param traits A data frame of species traits, with one row per species 
#'   and traits in columns. Must include a column of species names 
#'   (by default \code{traits[,1]} is used) which will be matched 
#'   to \code{colnames(com)}.
#' @param subset_traits Optional. A character vector specifying a subset of trait columns 
#'   (excluding the first column of species names) to use in Gower distance. 
#'   If \code{NULL} (default), all available traits (except the first column) are used.
#' @param all_species Logical. Whether to use all species in the trait data. 
#'   If \code{FALSE} (default), only species present in \code{com} are retained.
#' @param gower_weights Optional. A numeric vector of weights for the traits 
#'   when computing Gower distance. If \code{NULL} (default), all traits 
#'   are weighted equally.
#'
#' @return A list with:
#'   \item{eig}{A vector of \% variance explained by the first 1-3 PCoA axes.}
#'   \item{fNICE}{A data frame of volumes lost and gained between consecutive time points.}
#'   \item{fNICE_cum}{A data frame tracking cumulative lost and gained volumes, 
#'         plus a simple imbalance metric over time.}
#'
#' @examples
#' # Suppose com, infos, and traits are loaded:
#' # out_fNICE <- compute_fNICE(com = your_com, infos = your_infos, traits = your_traits)
#' # str(out_fNICE)
#'
compute_fNICE <- function(com, 
                          infos, 
                          traits,
                          subset_traits = NULL, 
                          all_species   = FALSE, 
                          gower_weights = NULL) {
  
  # Optional: log the site or study ID
  print(infos[, 1])
  
  # 1. Standardize species names in com
  colnames(com) <- gsub(' ', '.', colnames(com))
  
  # 2. Subset traits if requested
  if (!is.null(subset_traits)) {
    traits <- traits[, subset_traits]
  }
  
  # 3. Convert 'sci_name' in traits to row names (using a safe pattern)
  rownames(traits) <- gsub(' ', '.', traits[, 1], fixed = TRUE)
  # Remove the original species column now that it's in rownames
  traits <- traits[, -1]
  
  # 4. Remove species with only NA traits
  species_to_remove <- rownames(traits)[apply(traits, 1, function(x) all(is.na(x)))]
  if (length(species_to_remove) > 0) {
    com    <- com[, !(colnames(com) %in% species_to_remove), drop = FALSE]
    traits <- traits[!(rownames(traits) %in% species_to_remove), , drop = FALSE]
  }
  print(paste(length(species_to_remove), "species removed due to missing trait values"))
  
  # 5. Remove species never sampled in com
  #    i.e., total abundance = 0 across all time points
  if (any(apply(com, 2, sum) == 0)) {
    com <- com[, -which(apply(com, 2, sum) == 0), drop = FALSE]
  }
  
  # 6. Remove time points with < 4 total species present
  sp_rich <- apply(ifelse(com > 0, 1, 0), 1, sum)
  drop_rows <- which(sp_rich < 4)
  if (length(drop_rows) > 0) {
    com <- com[-drop_rows, , drop = FALSE]
  }
  
  # 7. If not using all_species, ensure that com and traits match up
  if (!all_species) {
    if (is.null(colnames(com))) {
      colnames(com) <- paste('sp', 1:ncol(com), sep = '.')
    }
  }
  
  # 8. Match species order between com and traits
  #    Keep only species that appear in both
  traits <- traits[rownames(traits) %in% colnames(com), , drop = FALSE]
  # Match the order
  traits <- traits[match(colnames(com), rownames(traits)), , drop = FALSE]
  
  # 9. Compute Gower distance for functional space
  gower_dist <- StatMatch::gower.dist(data.x = traits, var.weights = gower_weights)
  
  # 10. Perform PCoA on the Gower distance matrix
  samp_pcoa <- ape::pcoa(D = gower_dist, correction = 'cailliez')
  
  # 11. Extract eigenvalues (percentage of variance explained)
  #     We take up to 3 axes, or fewer if the dataset is smaller
  max_axes <- min(c(3, ncol(samp_pcoa$vectors)))
  eig_raw  <- samp_pcoa$values$Eigenvalues
  eig      <- round(100 * (samp_pcoa$values[1:max_axes, 1] / sum(eig_raw)), 1)
  
  # 12. Extract the PCoA vectors (coordinates in functional space)
  samp_pcoa_vectors <- samp_pcoa$vectors[, 1:max_axes]
  
  # 13. Use complementary_volume() to get volumes lost/gained between consecutive time points
  comp_vols <- try(
    complementary_volume(trait_mat = samp_pcoa_vectors, comm_mat = com),
    silent = FALSE
  )
  
  # If volume calculation fails (e.g. degenerate hull), return NULL
  if (inherits(comp_vols, "try-error")) {
    return(NULL)
  } else {
    # Flatten the list of data frames into one data frame
    comp_vols <- do.call(rbind, comp_vols)
    
    # Prepare the output
    # Observed volumes lost/gained
    fNICE <- data.frame(
      year1    = comp_vols$year1,
      year2    = comp_vols$year2,
      obs.loss = comp_vols$lost_vol,
      obs.gain = comp_vols$gained_vol
    )
    
    # Compute a simple imbalance metric:
    #   fNICE.obs = (cumulative gain - cumulative loss) / (cumulative gain + cumulative loss)
    fNICE_cum <- data.frame(
      year         = fNICE$year2,
      cum.obs.loss = cumsum(fNICE$obs.loss),
      cum.obs.gain = cumsum(fNICE$obs.gain),
      fNICE.obs    = (cumsum(fNICE$obs.gain) - cumsum(fNICE$obs.loss)) /
        (cumsum(fNICE$obs.gain) + cumsum(fNICE$obs.loss))
    )
    
    # Combine into a single output object
    out_FD <- list(
      eig        = eig,
      fNICE      = fNICE,
      fNICE_cum  = fNICE_cum
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
#' for a given community dataset over time. It estimates the timing of species 
#' colonizations and extinctions, counts these events per year, and computes 
#' the cumulative tNICE metric over time.
#'
#' @param com A community matrix with rows as time points (e.g., years) 
#'   and columns as species names. Values are species abundances or presence 
#'   (counts) at each time point.
#' @param infos A data frame or list containing metadata about the study 
#'   (e.g., study ID, location, year range). Not strictly required for calculations,
#'   but useful for logging or QA.
#'
#' @return A data frame with:
#'   \item{YEAR}{The time points (e.g., years) from rownames of \code{com}.}
#'   \item{COL}{Cumulative number of species colonizations up to each time.}
#'   \item{EXT}{Cumulative number of species extinctions up to each time.}
#'   \item{NICE}{\eqn{(Cumulative COL - Cumulative EXT) / (Cumulative COL + Cumulative EXT)} 
#'               at each time.}
#'   \item{SR}{Species richness (number of species present) at each time.}
#'
#' @details
#' \itemize{
#'   \item \strong{Colonization time}: estimated by reversing the time axis 
#'         and applying the OLE method as if searching for "last sightings."
#'   \item \strong{Extinction time}: estimated using the OLE method on the 
#'         forward time axis (years).
#' }
#'
#' @examples
#' # Suppose com and infos are loaded:
#' # out_tNICE <- compute_tNICE(com = your_com, infos = your_infos)
#' # head(out_tNICE)
#'
compute_tNICE <- function(com, infos) {
  
  # Optional: print site/study ID for QA/logging
  print(infos[, 1])
  
  # 1. Remove species never sampled
  #    (i.e., total abundance = 0 across all time points)
  if (any(apply(com, 2, sum) == 0)) {
    com <- com[, -which(apply(com, 2, sum) == 0), drop = FALSE]
  }
  
  # 2. Ensure at least 5 species and 5 sampling points
  if (ncol(com) < 5 | nrow(com) < 5) {
    return(NULL)
  }
  
  # 3. Assign column names if missing
  if (is.null(colnames(com))) {
    colnames(com) <- paste('sp', 1:ncol(com), sep = '.')
  }
  
  # 4. Extract numeric time points (from rownames)
  years <- as.numeric(rownames(com))
  
  # 5. Extinction timing estimation via OLE on forward time
  #    OLE returns (est, lower, upper). We store them as Est.E, low.E, up.E.
  est.EXT <- apply(com, 2, function(x) {
    OLE(sightingdata = cbind(years, x), alpha = .05)
  })
  est.EXT <- data.frame(
    Sp = colnames(com),
    do.call(rbind, lapply(est.EXT, function(x) cbind(x[1], x[2], x[3])))
  )
  colnames(est.EXT) <- c('Sp', 'Est.E', 'low.E', 'up.E')
  
  # 6. Colonization timing estimation via reversed time
  #    invert years with respect to the max(years)
  #    then apply OLE as if it's extinction in reversed time.
  samp.COL <- cbind(abs(years - max(years)), com)
  samp.COL <- samp.COL[order(samp.COL[, 1]), ]
  
  est.COL <- apply(samp.COL[, -1], 2, function(x) {
    # for each species, run OLE on forward direction
    # but the forward direction is "reversed years"
    lapply(
      OLE(sightingdata = cbind(years, x), alpha = .05),
      function(z) abs(z - max(com[, 1]))
    )
  })
  
  est.COL <- data.frame(
    Sp = colnames(samp.COL[, -1]),
    do.call(rbind, lapply(est.COL, function(x) cbind(x[1], x[2], x[3])))
  )
  colnames(est.COL) <- c('Sp', 'Est.C', 'low.C', 'up.C')
  
  # 7. Combine both extinction and colonization info
  sp.infos <- merge(est.COL, est.EXT, by = 'Sp', all = TRUE)
  # Convert to numeric (in case they're factors)
  sp.infos[, -1] <- apply(sp.infos[, -1], 2, function(x) as.numeric(as.character(substr(x, 1, 4))))
  
  # 8. Count colonizations (COL) and extinctions (EXT) by year
  #    We create frequency tables and sum across species for each year.
  COL_counts <- apply(table(sp.infos$Sp, sp.infos$Est.C), 2, sum)
  COL_df <- data.frame(YEAR = as.numeric(names(COL_counts)), COL = COL_counts)
  
  EXT_counts <- apply(table(sp.infos$Sp, sp.infos$Est.E), 2, sum)
  EXT_df <- data.frame(YEAR = as.numeric(names(EXT_counts)), EXT = EXT_counts)
  
  # 9. Merge into one data frame of events
  events <- merge(COL_df, EXT_df, by = 'YEAR', all = TRUE)
  
  # 10. Restrict to sampled years and replace NA with 0
  events <- merge(data.frame(YEAR = years), events, all.x = TRUE, all.y = FALSE)
  events[is.na(events)] <- 0
  
  # 11. Cumulative sum of COL & EXT across the time points
  events.cum <- cbind(
    YEAR = events$YEAR,
    as.matrix(apply(events[, c('COL','EXT')], 2, cumsum))
  )
  
  # 12. Compute tNICE = (cumulative COL - cumulative EXT) / (cumulative COL + cumulative EXT)
  NICE <- (events.cum[, 'COL'] - events.cum[, 'EXT']) /
    (events.cum[, 'COL'] + events.cum[, 'EXT'])
  # Replace Inf, -Inf, NaN with NA
  NICE[is.infinite(NICE)] <- NA
  NICE[is.nan(NICE)]      <- NA
  
  # 13. Combine into final output
  out.TD <- data.frame(
    YEAR = events.cum[, 'YEAR'],
    COL  = events.cum[, 'COL'],
    EXT  = events.cum[, 'EXT'],
    NICE = NICE
  )
  
  # 14. Compute species richness by presence (>0) for each row/year
  SR <- data.frame(
    YEAR = years,
    SR   = apply(ifelse(com > 0, 1, 0), 1, sum)
  )
  
  # Merge SR into output
  out.TD <- merge(out.TD, SR, by = 'YEAR')
  
  return(out.TD)
}

# -------------------------------------------------------------------
# Example Usage
# -------------------------------------------------------------------

# 1. Load Your Data
# Replace 'your_file.RData' with the path to your data file
# load('your_file.RData')

# For this example, assume you have:
#   com    = your_community_matrix
#   infos  = your_metadata
#   traits = your_traits_data

# 2. Run the tNICE Function
# tNICE_result <- compute_tNICE(com = com, infos = infos)
# head(tNICE_result)

# 3. Run the fNICE Function
# fNICE_result <- compute_fNICE(
#   com          = com,
#   infos        = infos,
#   traits       = traits,
#   subset_traits= NULL,        # Use all traits
#   all_species  = FALSE,       # Only species in `com`
#   gower_weights= NULL         # Equal weighting for traits
# )
# str(fNICE_result)

# -------------------------------------------------------------------
# Notes:
# - Ensure that species names in 'com' match row names in 'traits' 
#   (minus the first column if you store species there).
# - Replace placeholder data and paths with your actual data.
# - Install any required packages that are not already installed.
# - If you have fewer than 3 species in a time step, the convex hull 
#   computation may fail; handle carefully or skip those time steps.
# -------------------------------------------------------------------
