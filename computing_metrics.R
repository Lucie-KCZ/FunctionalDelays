#########################################
###    Functional diversity loss      ###
###       and taxonomic delays        ###
###   of European freshwater fish     ###
### and North American Breeding birds ###
#########################################

# Authors list
# Lucie Kuczynski, Ana Maria Bastidas Urrutia, Helmut Hillebrand

# Script infos
# L. Kuczynski
# lucie.kuczynski@hotmail.com
# Mar, 2023
# last edit: Mar, 2022

dataset.name <- 'RivFish'
# dataset.name <- 'bbs'
# dataset.name <- 'mzb'

dir.create(paste0('./out/', dataset.name))

# Import the community data formatted as a list with the traits
if(dataset.name == 'RivFish') load('~/Dropbox/fNICE/out/RivFish/RivFishxtraits.RData', v = T)
if(dataset.name == 'BBS') load('~/Dropbox/fNICE/out/BBS/BBSxtraits.RData', v = T)
Obs.list <- out ; rm(out)

# Getting the SR distributions for simulations
Obs.infos <- do.call(rbind, lapply(Obs.list, function(x) x$infos))
names(Obs.list) <- Obs.infos$STUDY_ID
SR.distribution <- unlist(lapply(Obs.list, function(x) rowSums(ifelse(x$com.matrix > 0, 1, 0))))

# Final clean up for fish data
if(dataset.name == 'RivFish') {
  # only European sites
  Obs.Europe <- Obs.infos[Obs.infos$BIOME_MAP == 'Palearctic', ]
  
  # at least 5 sampled years
  Obs.Europe10 <- Obs.Europe[Obs.Europe$NB_SAMP > 4, ]
  
  # final data
  Obs.list <- Obs.list[Obs.Europe10$STUDY_ID]  
}

# Simulated data
library(tidyverse)

# Import the functions from island pkg
files <- list.files('~/Dropbox/fNICE/analysis/functions island package/R/', full.names = T)
for(i in 1:17) source(files[i]) ; rm(files, i)

# Import the traits
load(paste0('../out/', dataset.name, '/traits.RData'), v = T)

# Simulation - with sp having trait information
simulating.ts <- function(i, constraints = Obs.infos[, 1:4], autocor.factor = 1, rate.c, rate.e, SR.dist = SR.distribution, traits.df = traits) {
    
  # Initial parameters of the simulated time series
  length.sim <- sample(constraints$TS_LENGTH, size = 1)
  nb.samp.sim <- sample(constraints$NB_SAMP[constraints$NB_SAMP <= length.sim], size = 1)
  nb.sp.sim <- sample(SR.dist, size = 1)
  
  # Sp pool
  sp.pool <- 0
  while(sp.pool < nb.sp.sim){
    sp.pool <- sample(x = constraints$NB.TOT.SP[constraints$NB.TOT.SP >= nb.sp.sim], 1)
  }
  
  initial <- matrix(c(rep(0, sp.pool - nb.sp.sim), rep(1, nb.sp.sim)), ncol = 1)
  
  # Set colonization and extinction rates
  col <- rate.c / autocor.factor
  ext <- rate.e / autocor.factor
  1 - sum(cetotrans(col, ext)) # This is the temporal autocorrelation. See Ontiveros et al 2021 appendix.
  
  # Running the simulation to get the com dynamics
  sim.ts <- t(PA_simulation(initial, 1, cetotrans(col, ext), length.sim))
  
  # Having the 'years'
  rownames(sim.ts) <- 1:nrow(sim.ts)
  
  # Selecting the sampled years based on the nb of samples
  sim.ts <- sim.ts[sort(sample(1:length.sim, nb.samp.sim)), ]
  
  # Getting random species
  # to then match with traits
  sp.list <- traits.df$sci_name
  colnames(sim.ts) <- sample(sp.list, ncol(sim.ts))
  
  # Format
  list.sim <- list(infos = data.frame(STUDY_ID = paste0('sim_', i),
                                      NB_SAMP = nb.samp.sim,
                                      TS_LENGTH = length.sim,
                                      NB.TOT.SP = ncol(sim.ts),
                                      LATITUDE = NA, LONGITUDE = NA, ID_LUCIE = NA,
                                      PROTECTED_AREA = NA, BIOME.MAP = NA, GRAIN_SQ_KM = NA,
                                      AREA_SQ_KM = NA, SUMMARY_METHODS = 0,
                                      ABUNDANCE_TYPE = 0, BIOMASS = NA),
                   com.matrix = sim.ts, 
                   traits = traits.df[traits.df$sci_name %in% colnames(sim.ts), ])
  
  return(list.sim)
  
}

# # For a first time, run the NICE function first on the obs data
# # as the rates are from the estimated COL and EXT events (i.e. by OLEs)
# load(paste0('./out/', dataset.name, '/NICE.RData'))
# rates <- do.call(rbind, lapply(NICE$observed, function(x) x[ncol(x), 3:4] / ncol(x)))
# rates <- rates[-which(apply(rates, 1, sum) == 0), ]
# rates <- apply(rates, 2, mean)
# rm(NICE)

system.time(sim.list <- mapply(FUN = simulating.ts, i = 1:9, 
                               MoreArgs = list(autocor.factor = 1, rate.c = rates[1], rate.e = rates[1], traits.df = traits),
                               SIMPLIFY = F)) 

names(sim.list) <- paste('autocor.', 1)
rm(rates)

# Merging the data
out <- list(observed = Obs.list, 
            simulated = sim.list)
save(out, file = paste0('../out/', dataset.name, '/complete.datasets.RData'))

# cleaning up
rm(list = ls()[-which(ls() %in% c('out', 'dataset.name'))])

# Running tNICE
source('./fNICE.R')
system.time(NICE <- lapply(out, function(x) lapply(x, compute_tNICE)))
save(NICE, file = paste0('../out/', dataset.name, '/tNICE.RData'))

# Running fNICE
if(dataset.name == 'BBS'){
    system.time(fNICE.all.out <- lapply(out, function(x) lapply(x, compute_fNICE, nb.null = 59, traits = c(2:20, 35:44))))
    names(fNICE.all.out$observed) <- names(out$observed)
    save(fNICE.all.out, file = paste0('../out/', dataset.name, '/fNICE.all.RData'))
  }
 
if(dataset.name == 'RivFish'){
    system.time(fNICE.all.out <- lapply(out, function(x) lapply(x, compute_fNICE, nb.null = 59, traits = c(3:16, 18:34))))
    names(fNICE.all.out$observed) <- names(out$observed)
    save(fNICE.all.out, file = paste0('../out/', dataset.name, '/fNICE.all.RData'))
  }

# Formating fNICE output
format.NICE <- function(fNICE){
  # Eigen values
  eig.out <- lapply(fNICE, function(x) do.call(rbind, lapply(x, function(y) y$eig)))
  print(lapply(eig.out, function(x) apply(x, 2, function(x) round(mean(x), 2))))
  
  # # to get the average and sd explained variance by the axes
  # lapply(lapply(eig.out, function(x) apply(x, 1, sum, na.rm = T)), mean)
  # lapply(lapply(eig.out, function(x) apply(x, 1, sum, na.rm = T)), sd)
  
  # DF format
  fNICE.df <- lapply(fNICE, function(x) lapply(x, function(y) y$fNICE.cum)) # VOLUME
  
  # DF merging
  NICEs.obs <- list()
  for(i in 1:length(NICE$observed)){
    if(!is.null(fNICE.df$observed[[i]])){
      NICEs.obs[[i]] <- data.frame(ID = names(NICE$observed)[i], fNICE.df$observed[[i]], NICE$observed[[i]][-1, -1]) # VOLUME
    }
  } ; rm(i)
  NICEs.obs <- do.call(rbind, NICEs.obs)
  NICEs.obs$fNICE.raw <- (NICEs.obs$cum.obs.gain - NICEs.obs$cum.obs.loss) / (NICEs.obs$cum.obs.gain + NICEs.obs$cum.obs.loss)
  
  NICEs.sim <- list()
  for(i in 1:length(NICE$simulated)){
    if(!is.null(fNICE.df$simulated[[i]])) NICEs.sim[[i]] <- data.frame(ID = paste('sim', i, sep = '_'), fNICE.df$simulated[[i]], NICE$simulated[[i]][-1, -1]) # VOLUME
  } ; rm(i)
  NICEs.sim <- do.call(rbind, NICEs.sim)
  NICEs.sim$fNICE.raw <- (NICEs.sim$cum.obs.gain - NICEs.sim$cum.obs.loss) / (NICEs.sim$cum.obs.gain + NICEs.sim$cum.obs.loss)
  
  print(table(sign(NICEs.obs$NICE), sign(NICEs.obs$fNICE))) #tNICE: raw / fNICE: col
  
  NICEs <- list(obs = NICEs.obs, sim = NICEs.sim)
  
  return(NICEs)
  
}

NICE.all <- format.NICE(fNICE = fNICE.all.out)

# Plotting the trends
plot.trends <- function(samp, col = 'deeppink', k = 8, obs = T){
  
  # selecting the right variable of the list
  if(obs) samp <- samp$obs
  if(!obs) samp <- samp$sim
  
  if(max(samp$year) < 1000) samp$year <- samp$year + 1950
  
  # opening pdf
  pdf(file = paste0('../plots/', dataset.name, '/Trends.', 
                    ifelse(obs, 'obs.', 'sim.'), 
                    readline('Name of the plot: '), '.pdf'),   # The directory you want to save the file in
      width = k, # The width of the plot in inches
      height = k) # The height of the plot in inches
  
  par(mar = c(7, 7, 2, 2), pty = 's')
  
  # empty plot
  plot(x = 1, y = 1, xlim = c(1950, 2020), ylim = c(-1, 1), type = 'n', axes = F, xlab = '', ylab = '')
  
  # axes
  axis(side = 1, tcl = -.5, lwd = 2, las = 2, at = seq(1950, 2020, by = 10), cex.axis = 2.5, padj = .5)
  axis(side = 2, tcl = -.5, lwd = 2, las = 1, at = seq(-1, 1, by = .5), cex.axis = 2.5)
  
  # adding obs data
  by(samp, samp$ID, function(x) lines(x$fNICE.raw ~ x$year, col = 'gray75', lwd = 1, lty = 1))
  
  # getting the model
  library(nlme)
  model <- lme(fixed = fNICE.raw ~ year,
               random = ~ 1 | ID, method = "REML",
               data = samp)
  
  # getting estimates
  print(summary(model))
  
  # getting r2
  library(MuMIn)
  print(r.squaredGLMM(model))
  
  # getting IC
  library(ggeffects)
  pre <- as.data.frame(ggpredict(model, "year", nsim = 999, type = 'sim'))
  
  # plot IC + model
  polygon(y = c(pre$conf.low, rev(pre$conf.high), pre$conf.low[1]),
          x = c(pre$x, rev(pre$x), pre$x[1]), col = paste0(col, 1), border = paste0(col, 1))
  lines(pre$x, pre$predicted, lwd = 3, col = paste0(col, 4))
  
  # closing pdf
  dev.off()
  
}

plot.trends(NICE.all)
plot.trends(NICE.all, obs = F)
  

# comparing tNICE and fNICE
plot.NICE <- function(fNICE, k = 8, obs = T, raw = T){
  
  if(obs) final.NICE <- do.call(rbind, as.list(by(data = fNICE$obs, INDICES = fNICE$obs$ID, FUN = function(x) x[nrow(x), ])))
  if(!obs) final.NICE <- do.call(rbind, as.list(by(data = fNICE$sim, INDICES = fNICE$sim$ID, FUN = function(x) x[nrow(x), ])))

  if(raw) {
    var <- final.NICE$fNICE.raw
  } else {
    var <- final.NICE$fNICE
  }

  
  col.NICE <- rep(x = 'snow2', times = nrow(final.NICE))
  # upper right panel
  col.NICE[which(var > final.NICE$NICE & sign(var) >= 0 & sign(final.NICE$NICE) > 0)] <- 'steelblue1'
  col.NICE[which(var < final.NICE$NICE & sign(var) >= 0 & sign(final.NICE$NICE) > 0)] <- 'steelblue4'
  # lower left panel
  col.NICE[which(var > final.NICE$NICE & sign(var) <= 0 & sign(final.NICE$NICE) < 0)] <- 'darkgoldenrod1'
  col.NICE[which(var < final.NICE$NICE & sign(var) <= 0 & sign(final.NICE$NICE) < 0)] <- 'darkgoldenrod4'
  # upper left panel
  col.NICE[which(var > final.NICE$NICE & sign(var) >= 0 & sign(final.NICE$NICE) < 0)] <- 'darkorange3'
  # lower right panel
  col.NICE[which(var < final.NICE$NICE & sign(var) <= 0 & sign(final.NICE$NICE) > 0)] <- 'darkolivegreen3'
  
  print(table(col.NICE))
  
  print(cor.test(var, final.NICE$NICE, method = 'spearman'))
  
  pdf(file = paste0('../plots/', dataset.name, '/', ifelse(obs, 'Obs.', 'Sim.'), readline('Name of the plot: '), '.pdf'),   # The directory you want to save the file in
      width = k, # The width of the plot in inches
      height = k) # The height of the plot in inches
   
  par(mfrow = c(1, 1), pty = 's')
   
  plot(var ~ final.NICE$NICE, 
       xlab = '', ylab = '', type = 'n', las = 1,
       xlim = c(-1, 1), ylim = c(-1, 1),
       frame.plot = F, axes = F)
  axis(side = 1, tcl = -.5, lwd = 2, las = 1, at = seq(-1, 1, by = .5), cex.axis = 3, padj = .5)
  axis(side = 2, tcl = -.5, lwd = 2, las = 1, at = seq(-1, 1, by = .5), cex.axis = 2.75)
  abline(h = 0, v = 0, lwd = 2)
  abline(0, 1, lty = 2, lwd = 2)
  points(y = var, x = final.NICE$NICE, pch = 21, bg = col.NICE, cex = 2)
  
  barplot(sort(table(col.NICE)), las = 2, col = names(sort(table(col.NICE))), lwd = 2)
  dev.off()
  
}

plot.NICE(NICE.all, obs = T, raw = T)

### END OF THE SCRIPT
