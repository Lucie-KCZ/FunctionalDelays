######################
### fNICE function ###
######################

# Authors list
# Lucie Kuczynski, Helmut Hillebrand

# L. Kuczynski
# lucie.kuczynski@hotmail.com
# Apr, 2022
# last edit: Apr, 2022


source('./sExtinct/R/OLE.R')
source('./sExtinct/R/OLE.fun.R')

# # samp <- out$simulated[[80]]
# # samp <- out$observed[[1]]
# nb.null <- 9
# samp <- out$simulated[[5]]
# traits = c(25:33)

# to compute the unique volume occupied buy each pair of yr
complementary_volume <- function(trait_mat, comm_mat, nb.axes = 2){
  
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
  for (i in 1:(years-1)) {
    ch1 <- convhulln(pcoa$vectors[which(comm_mat[i,] > 0), 1:nb.axes], option = 'FA')
    ch2 <- convhulln(pcoa$vectors[which(comm_mat[i + 1,] > 0), 1:nb.axes], option = 'FA')
    
    # intersection
    common_vol <- intersectn(ch1$hull, ch2$hull) 
    
    # unique part of each volume
    # adding the if for subset and making sure we get zeros.
    if(all(which(com.matrix[i,] > 0) %in% which(com.matrix[i + 1,] > 0))) vol1 <- 0 else vol1 <- ch1$vol - common_vol$ch$vol
    if(all(which(com.matrix[i + 1,] > 0) %in% which(com.matrix[i,] > 0))) vol2 <- 0 else vol2 <- ch2$vol - common_vol$ch$vol
    
    
    # what's lost is what only have during year_1
    # what's gained is what only have during year_2
    vol[[paste0(i,"-",i - 1)]] <- data.frame(year1 = as.numeric(rownames(comm_mat)[i]), 
                                             year2 = as.numeric(rownames(comm_mat)[i + 1]), 
                                             lost_vol = abs(vol1), gained_vol = abs(vol2))
    
  }
  return(vol)
}

# compute functional version - volume
compute_fNICE <- function(samp, nb.null = 49, traits = NULL){
  
  print(samp$infos[, 1])
  
  samp.com <- samp$com.matrix
  
  # working a little bit on the traits
  # selecting traits
  if(!is.null(traits)) samp.traits <- samp$traits[, traits] else samp.traits <- samp$traits[, -1]
  # format to have the sp names
  rownames(samp.traits) <- samp$traits$sci_name
  # make sure not undesirable sp
  samp.traits <- samp.traits[rownames(samp.traits) %in% colnames(samp.com), ]
  # gettting the smae order
  samp.traits <- samp.traits[match(colnames(samp.com), rownames(samp.traits)), ]
  # getting everything as numeric
  samp.traits <- apply(samp.traits, 2, function(x) as.numeric(as.character(x)))
  
  # removing sp never sampled
  if(any(apply(samp.com, 2, sum) == 0)) samp.com <- samp.com[, -which(apply(samp.com, 2, sum) == 0)]
  # garde fou
  if(ncol(samp.com) < 5 | nrow(samp.com) < 5) return(NULL)
  
  # format
  if(is.null(colnames(samp.com))) colnames(samp.com) <- paste('sp', 1:ncol(samp.com), sep = '.')
  
  # functional distance matrix
  gower.dist <- cluster::daisy(x = samp.traits, metric = 'gower')
  
  # making sure enough trait data
  if(any(is.na(gower.dist))) {
    out.FD <- NULL
    return(out.FD)
    
  } else { 
    samp.pcoa <- ape::pcoa(D = gower.dist, correction = 'cailliez')
    
    eig <- round(100 * (samp.pcoa$values[1:min(c(4, ncol(samp.pcoa$vectors))), 1] / sum(samp.pcoa$values$Eigenvalues)), 1)
    samp.pcoa <- samp.pcoa$vectors[, 1:min(c(4, ncol(samp.pcoa$vectors)))]
    
    # Calculate the complementary volumes
    comp_vols <- try(exp = complementary_volume(trait_mat = samp.traits, comm_mat = samp.com), silent = T)
    
    if(inherits(comp_vols, "try-error")) {
      out.FD <- NULL
    } else {
      comp_vols <- do.call(rbind, as.list(comp_vols))
      
      # Running NULL MODEL on that
      null.out <- list() ; k <- 1
      while(length(null.out) != nb.null){
        random.traits <- samp.traits[sample(1:nrow(samp.traits)), ]
        rownames(random.traits) <- rownames(samp.traits)
        
        null.k <- try(complementary_volume(trait_mat = random.traits, comm_mat = samp.com), silent = T)
        if(inherits(null.k, "try-error")) {
          next
        } else {
          null.out[[k]] <- do.call(rbind, as.list(null.k)) 
          rm(random.traits) ; k <- k + 1
        }
      } ; rm(k)
      
      null.out <- do.call(rbind, null.out)
      null.mean <- do.call(rbind, as.list(by(data = null.out, INDICES = null.out[, c('year1', 'year2')], FUN = function(x) apply(x, 2, mean, na.rm = T))))
      null.sd <- do.call(rbind, as.list(by(data = null.out, INDICES = null.out[, c('year1', 'year2')], FUN = function(x) apply(x, 2, sd, na.rm = T))))
      
      fNICE <- data.frame(year1 = comp_vols$year1, year2 = comp_vols$year2, 
                          obs.loss = comp_vols$lost_vol, obs.gain = comp_vols$gained_vol, 
                          mean.loss = null.mean[, 3], mean.gain = null.mean[, 4], 
                          sd.loss = null.sd[, 3], sd.gain = null.sd[, 4], 
                          ses.loss = ((comp_vols$lost_vol - null.mean[, 3]) / null.sd[, 3]), 
                          ses.gain = ((comp_vols$gained_vol - null.mean[, 4]) / null.sd[, 4]))
      
      fNICE.cum <- data.frame(year = fNICE$year2, 
                              cum.obs.loss = cumsum(fNICE$obs.loss), 
                              cum.obs.gain = cumsum(fNICE$obs.gain), 
                              cum.ses.loss = cumsum(fNICE$ses.loss), 
                              cum.ses.gain = cumsum(fNICE$ses.gain), 
                              fNICE = (cumsum(fNICE$ses.gain) - cumsum(fNICE$ses.loss)) / (cumsum(fNICE$ses.gain) + cumsum(fNICE$ses.loss)))
      
      out.FD <- list(eig = eig, 
                     fNICE = fNICE, 
                     fNICE.cum = fNICE.cum)
      
    }
    
    return(out.FD)
  }
}

# compute taxonomic version
compute_tNICE <- function(samp){
  
  # to get an idea where we are at...
  print(samp$infos[, 1])
  
  samp.com <- samp$com.matrix
  
  # removing sp never sampled
  if(any(apply(samp.com, 2, sum) == 0)) samp.com <- samp.com[, -which(apply(samp.com, 2, sum) == 0)]
  # garde fou
  if(ncol(samp.com) < 5 | nrow(samp.com) < 5) return(NULL)
  
  # format
  if(is.null(colnames(samp.com))) colnames(samp.com) <- paste('sp', 1:ncol(samp.com), sep = '.')
  
  # str of the ts
  years <- as.numeric(rownames(samp.com))
  
  # extinction timing estimation
  est.EXT <- apply(samp.com, 2, function(x) OLE(sightingdata = cbind(years, x), alpha = .05))
  est.EXT <- data.frame(Sp = colnames(samp.com), do.call(rbind, lapply(est.EXT, function(x) cbind(x[1], x[2], x[3]))))
  colnames(est.EXT) <- c('Sp', 'Est.E', 'low.E', 'up.E')
  
  # colonisation timing estimation
  samp.COL <- cbind(abs(years - max(years)), samp.com)
  samp.COL <- samp.COL[order(samp.COL[, 1]), ]
  
  est.COL <- apply(samp.COL[, -1], 2, function(x) lapply(OLE(
    sightingdata = cbind(years, x), alpha = .05),
    function(z) abs(z - max(samp.com[, 1]))))
  est.COL <- data.frame(Sp = colnames(samp.COL[, -1]), do.call(rbind, lapply(est.COL, function(x) cbind(x[1], x[2], x[3]))))
  colnames(est.COL) <- c('Sp', 'Est.C', 'low.C', 'up.C')
  
  # combining EXT and COL
  sp.infos <- merge(est.COL, est.EXT, by = 'Sp', all = T) ; rm(est.COL, est.EXT)
  sp.infos[, -1] <- apply(sp.infos[, -1], 2, function(x) as.numeric(as.character(substr(x, 1, 4))))
  
  # counting events for each year
  COL <- apply(table(sp.infos$Sp, sp.infos$Est.C), 2, sum)
  COL <- data.frame(YEAR = as.numeric(names(COL)), COL)
  EXT <- apply(table(sp.infos$Sp, sp.infos$Est.E), 2, sum)
  EXT <- data.frame(YEAR = as.numeric(names(EXT)), EXT)
  events <- merge(COL, EXT, all = T) ; rm(COL, EXT)
  
  # removing years from outside the sampling window
  events <- merge(matrix(years, dimnames = list(NULL, 'YEAR')), events, all.x = T, all.y = F)
  
  # adding the zero when no event
  events[is.na(events)] <- 0
  
  # cumulative events over time
  events.cum <- cbind(events$YEAR, as.matrix(apply(events[, -1], 2, cumsum))) ; rm(events)
  colnames(events.cum)[1] <- 'YEAR'
  
  # when the magic happens
  NICE <- (events.cum[, 'COL'] - events.cum[, 'EXT']) / (events.cum[, 'COL'] + events.cum[, 'EXT'])
  
  NICE <- data.frame(YEAR = events.cum[, 1], NICE = NICE)
  NICE[NICE == Inf] <- NA
  NICE[NICE == -Inf] <- NA
  NICE[is.nan(NICE[, 2]), 2] <- NA
  
  out.TD <- merge(events.cum, NICE) ; rm(events.cum, NICE)
  SR <- data.frame(YEAR = years, SR = apply(ifelse(samp.com > 0, 1, 0), 1, sum))
  out.TD <- merge(out.TD, SR) ; rm(SR)
  
  return(out.TD)
}
