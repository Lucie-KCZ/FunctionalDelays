# -------------------------------------------------------------------
# R Script for Estimating Extinction and Colonization Events
# -------------------------------------------------------------------

# Source External Functions
source('./sExtinct/R/OLE.R')
source('./sExtinct/R/OLE.fun.R')

# -------------------------------------------------------------------
# Function: estimated_rates
# -------------------------------------------------------------------

#' Estimate Extinction and Colonization Rates
#'
#' This function calculates species' extinction and colonization timings from 
#' abundance time-series data using the Optimal Linear Estimator (OLE) method, 
#' then derives colonization and extinction rates per unit time.
#'
#' The underlying methodology relies on OLE to estimate when a species 
#' might have gone extinct (from the last time it was observed) or colonized 
#' (from the first time it was observed, inverted by re-indexing time). 
#' The function returns both per-year rates (i.e., number of events divided by 
#' number of sampled years) and a detailed table of estimated event years 
#' (with lower and upper bounds).
#'
#' @param mat A matrix or data frame of species abundances (or presence-absence) 
#'   where rows are sampling years (named by the year) and columns are species.
#'   \strong{Note}: Row names must be numeric and correspond to the sampling year.
#' @return A list containing:
#'   \item{rates}{A data frame with `colonisation` and `extinction` rates.}
#'   \item{events}{A data frame containing estimated colonization (Est.C) 
#'   and extinction (Est.E) years, with lower and upper 95\% confidence bounds 
#'   (low.C, up.C, low.E, up.E).}
#'
#' @details
#' \enumerate{
#'   \item Species never observed (i.e., total abundance == 0 across all years) 
#'         are removed before analysis.
#'   \item The function guards against very small matrices (fewer than 5 rows or columns).
#'   \item Extinction times are estimated by applying \code{OLE()} separately 
#'         to each species' time series in the forward sense (actual year).
#'   \item Colonization times are estimated by applying \code{OLE()} in a reversed 
#'         timescale (i.e., measuring from the most recent year backward).
#'   \item Colonization and extinction rates are the number of species 
#'         that colonized or went extinct \emph{within} the observed time window, 
#'         divided by the number of sampled years.
#' }
#'
#' @examples
#' # Assume `mat` is your community matrix with row names as years
#' # result <- estimated_rates(mat)
#' # result$rates      # to see the colonisation and extinction rates
#' # head(result$events)  # to see the estimated events and their confidence intervals
#'
estimated_rates <- function(mat){
  
  # -----------------------------------------------------------------
  # 1. Preliminary Data Filtering
  # -----------------------------------------------------------------
  # Remove species never observed (i.e., total abundance = 0)
  if(any(apply(mat, 2, sum) == 0)) {
    mat <- mat[, -which(apply(mat, 2, sum) == 0)]
  }
  
  # "Garde fou" / Guard condition: 
  # if fewer than 5 columns or 5 rows remain, return NULL
  if(ncol(mat) < 5 | nrow(mat) < 5) return(NULL)
  
  # -----------------------------------------------------------------
  # 2. Ensure Species Names
  # -----------------------------------------------------------------
  # If species have no column names, assign generic "sp1", "sp2", etc.
  if(is.null(colnames(mat))) colnames(mat) <- paste('sp', 1:ncol(mat), sep = '.')
  
  # Extract the sampling years (row names should be numeric)
  years <- as.numeric(rownames(mat))
  
  # -----------------------------------------------------------------
  # 3. Estimating Extinction Times
  # -----------------------------------------------------------------
  # The OLE function typically requires a 2-column matrix: 
  # [year, presence/abundance]. For each species, we pass a vector x 
  # with length = number of years. cbind(years, x) forms the input for OLE.
  #
  # `OLE(sightingdata, alpha=.05)` returns a list (or numeric vector) 
  # containing the estimated time of extinction plus lower/upper confidence bounds.
  # We apply it to each column (i.e., each species).
  est.EXT <- apply(mat, 2, function(x) {
    OLE(sightingdata = cbind(years, x), alpha = .05)
  })
  
  # We then reshape the list into a data frame: each species gets 3 columns: 
  # (Est.E, low.E, up.E)
  est.EXT <- data.frame(
    Sp = colnames(mat),
    do.call(rbind, lapply(est.EXT, function(x) cbind(x[1], x[2], x[3])))
  )
  colnames(est.EXT) <- c('Sp', 'Est.E', 'low.E', 'up.E')
  
  # -----------------------------------------------------------------
  # 4. Estimating Colonization Times
  # -----------------------------------------------------------------
  # Colonization is conceptually the "first" appearance in chronological order. 
  # Here, we invert time by measuring from the maximum year 
  # (i.e., 'the present') backwards:
  #
  #   samp.COL <- cbind( abs(years - max(years)), mat )
  #
  # So if our years are [1980, 1981, ... 2015], and the last sampling year is 2015,
  # we transform it into [ |1980-2015|, |1981-2015|, ... , |2015-2015| ] = [35, 34, ..., 0].
  #
  # We then apply the OLE in exactly the same way, but on this reversed time scale, 
  # effectively treating what used to be "first sightings" as "last sightings" 
  # in reversed time.
  samp.COL <- cbind(abs(years - max(years)), mat)
  
  # Sort by this new "time" so that the earliest in the reversed timescale 
  # is at the top:
  samp.COL <- samp.COL[order(samp.COL[, 1]), ]
  
  # Then we apply OLE column-wise again
  est.COL <- apply(samp.COL[, -1], 2, function(x) {
    lapply(
      OLE(sightingdata = cbind(years, x), alpha = .05),
      function(z) abs(z - max(mat[, 1]))
    )
  })
  
  # Reshape as a data frame with columns (Est.C, low.C, up.C)
  est.COL <- data.frame(
    Sp = colnames(samp.COL[, -1]),
    do.call(rbind, lapply(est.COL, function(x) cbind(x[1], x[2], x[3])))
  )
  colnames(est.COL) <- c('Sp', 'Est.C', 'low.C', 'up.C')
  
  # -----------------------------------------------------------------
  # 5. Combine Colonization and Extinction Estimates
  # -----------------------------------------------------------------
  # Merge both data frames by species name. 
  # This gives us a table with columns: Sp, Est.C, low.C, up.C, Est.E, low.E, up.E
  events <- merge(est.COL, est.EXT, by = 'Sp', all = TRUE)
  rm(est.COL, est.EXT)
  
  # Convert the merged columns from factors/characters to numeric 
  # (sometimes OLE outputs can be factors)
  events[, -1] <- apply(events[, -1], 2, function(x) as.numeric(as.character(substr(x, 1, 4))))
  
  # -----------------------------------------------------------------
  # 6. Compute Colonization and Extinction Rates
  # -----------------------------------------------------------------
  # We define colonization (or extinction) "within the sampled years" 
  # as: event year is between min(years) and max(years).
  #
  # The rate is: (Number of events) / (Number of years in the sampling window).
  # We add 1 to the denominator to account for inclusive range 
  # if that is the desired interpretation.
  colonisation_rate <- length(which(events$Est.C >= min(years) & events$Est.C <= max(years))) /
    (1 + max(years) - min(years))
  
  extinction_rate <- length(which(events$Est.E >= min(years) & events$Est.E <= max(years))) /
    (1 + max(years) - min(years))
  
  # Return both the rates and the detailed events
  return(list(
    rates = data.frame(colonisation = colonisation_rate, 
                       extinction    = extinction_rate),
    events = events
  ))
}

# -------------------------------------------------------------------
# Technical Notes on OLE for Colonization and Extinction
# -------------------------------------------------------------------
# The 'OLE' function (Optimal Linear Estimator) estimates the most likely 
# time a species went extinct, assuming its last detection is known. 
# Formally, it uses maximum likelihood under the assumption that the sightings 
# (for a species that is truly extinct) follow an exponential distribution 
# tail for the gap after the last sighting.
#
# To adapt the same logic to colonization, we reverse time so that 
# a species' earliest appearance in real time becomes "the last appearance" 
# in the reversed timescale. This allows us to apply OLE's extinction logic 
# to a "colonization" perspective.
#
# Both estimated times (colonization and extinction) come with lower and upper 
# 95% confidence bounds, specified by alpha=0.05. 
# These help capture the uncertainty in exactly which year colonization 
# or extinction occurred, based on the sightings data.

# -------------------------------------------------------------------
# Example Usage
# -------------------------------------------------------------------
# 1. Prepare Your Data
# mat <- read.csv("your_data.csv") 
# row.names = sampling years
#
# 2. Run the Function
# result <- estimated_rates(mat)
#
# 3. Retrieve the Rates and Events
# result$rates
# head(result$events)
#
# -------------------------------------------------------------------
# Notes:
# - Ensure that row names of `mat` are numeric and correspond to sampling years.
# - `OLE.R` and `OLE.fun.R` must be correctly sourced; these scripts contain 
#   the logic for computing the Optimal Linear Estimator.
# - If species or row (year) naming conventions differ, adjust accordingly.
# -------------------------------------------------------------------
