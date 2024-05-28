
##' @title Box information
##' @param bgm.file BGM file,  Atlantis input
##' @param depths Cummulative depths (de max depth of each layer)
##' @return a dataframe with information (box-id, area (m), )
##' @author Demiurgo - https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis
##' @edited Holly Perryman (June 2022)
##' 
##' Comments - Holly Perryman - June 2022:
##' While working with the GoM model, the values produced for Vol and Vol2 proved to be problematic:
##'    1) there are NAs in the vol metrics - why not zeros?
##'    2) I could not figure out if the layer ordering matched that of inputs ncs or output ncs,
##'    3) the cum depth with being used to calc vol when it should be the box depth (diff in cum dep layers) - this data is provided by nc out file nominal_dz
##'    4) the sediment layer seemed to be zeroed out for one of the vol metrics 
##' Ultimately, I removed all the code collecting for vol and vol2. I opted
##' to compute vol directly in the functions computing biomass from nc out data as
##' vol data across polygon, layer and time is provided by nc out data. 
##' Also, I removed cum.depth from the function as it can be determined by the bgm
##' file. Less information asked from the user, the better.  
##' 
boxes.prop <- function(bgm.file){
  bgm       <- readLines(bgm.file, warn = FALSE)
  vert      <- text2num(bgm, 'bnd_vert ', Vector = TRUE, lineal = TRUE)
  centroids <- text2num(bgm, '.inside', Vector = TRUE, lineal = TRUE)
  dynamic   <- as.numeric(apply(centroids, 1, function(x) .point_in_polygon(vert, x)))
  boxes     <- text2num(bgm, 'nbox', FG = 'look')
  out       <- NULL
  # figure out depth profile from bgm file: 
  depths <- NULL
  for(b in 1 : boxes$Value){
    depths[b] <- text2num(bgm, paste0('box', b - 1,'.botz'), FG = 'look')$Value
  }
  depths <- sort(unique(abs(depths)))
  # with this code, the following correction is no longer needed... 
  # if(all(depths[1 : 2] == 0)){
  #   depths <- depths[-1]
  # } else {
  #   depths <- depths[ - which(depths == 0)]
  # }
  max.nlyrs <- length(depths)              ## maximum number of water layers
  #vol       <- array(NA, dim = c(boxes$Value, length(depths)))
  for(b in 1 : boxes$Value){
    area        <- text2num(bgm, paste0('box', b - 1,'.area'), FG = 'look')
    area$Value  <- area$Value * dynamic[b] ## removing the non - dynamic boxes
    z           <- text2num(bgm, paste0('box', b - 1,'.botz'), FG = 'look')
    box.lyrs    <- sum(depths <  - z$Value)
    box.lyrs    <- pmin(box.lyrs, max.nlyrs) # bound by maximum depth
    out         <- rbind(out, data.frame(Boxid   = b - 1, 
                                         Area    = area$Value, 
                                         Volumen = area$Value * -z$Value,
                                         Depth   = -z$Value, 
                                         Layers  = box.lyrs))
    #vol[b, 1 : box.lyrs] <- area$Value * depths[1 : box.lyrs]
  }
  #vol                <- cbind(out$Area, vol) # to include the sediment layer for later calculations
  #vol                <- t(vol[, ncol(vol) : 1])
  #vol[1, ]           <- 0 # removing area from the boundary box
  #vol2               <- vol ## arrange not for infauna
  #vol2[nrow(vol2), ] <- 0
  ## if(any(out$Area < 1)) warning('\nOne (or more) of the boxes areas is less than 1m2,  Check if the right BGM file in xy coordinates')
  ## force are in boundary boxes to be zero
  out[c(1, which(out$Depth <= 0)), 2 : ncol(out)] <- 0
  out$dynamic      <-  dynamic
  out$Depth_Bound  <- out$Depth * out$dynamic ## removing boundary boxes
  # out <- list(info   = out,
  #             Vol    = vol2,
  #             VolInf = vol)
  out <- list(info = out) # I dropped the data for Vol and VolInf as I do not think it is accurate wrt GOMatlantis
  return(out)
}

##' @title Parameter file reader
##' @param text Biological parametar file for Atlatnis
##' @param pattern Text that you are looking
##' @param FG Name of the functional groups
##' @param Vector Logic argument, if the data is on vectors or not
##' @param pprey Logic argument, if the data is a pprey matrix or not
##' @param lineal Logic argument,  if the data is in a vector a lineal
##' @return A matrix with the values from the .prm file
##' @author Demiurgo - https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis
##' 
text2num <- function(text, pattern, FG = NULL, Vector = FALSE, pprey = FALSE, lineal = FALSE){
  if(!isTRUE(Vector)){
    text <- text[grep(pattern = pattern, text)]
    if(length(text) == 0) warning(paste0('\n\nThere is no ', pattern, ' parameter in your file.'))
    txt  <- gsub(pattern = '[[:space:]]+' ,  '|',  text)
    col1 <- col2 <- vector()
    for( i in 1 : length(txt)){
      tmp     <- unlist(strsplit(txt[i], split = '|', fixed = TRUE))
      if(grepl('#', tmp[1])) next
      tmp2    <- unlist(strsplit(tmp[1], split = '_'))
      if(FG[1] == 'look') {
        col1[i] <- tmp2[1]
      } else {
        id.co   <- which(tmp2 %in% FG)
        if(sum(id.co) == 0) next
        col1[i] <- tmp2[id.co]
      }
      col2[i] <- as.numeric(tmp[2])
    }
    if(is.null(FG)) col1 <- rep('FG', length(col2))
    out.t <- data.frame(FG = col1, Value = col2)
    if(any(is.na(out.t[, 1]))){
      out.t <- out.t[-which(is.na(out.t[, 1])), ]
    }
    return(out.t)
  } else {
    l.pat <- grep(pattern = pattern, text)
    nam   <- gsub(pattern = '[[:space:]]+' ,  '|',  text[l.pat])
    fg    <- vector()
    pos   <- 1
    for( i in 1 : length(nam)){
      tmp     <- unlist(strsplit(nam[i], split = '|', fixed = TRUE))
      if(grepl('#', tmp[1]) || (!grepl('^pPREY', tmp[1]) && pprey  == TRUE)) next
      fg[pos] <- tmp[1]
      if(isTRUE(lineal)){
        t.text <- gsub('+[[:space:]]+', ' ',  text[l.pat[i]])
      } else {
        t.text <- gsub('+[[:space:]]+', ' ',  text[l.pat[i] + 1])
      }
      oldw <- getOption("warn")
      options(warn = -1)
      pp.tmp <- matrix(as.numeric(unlist(strsplit(t.text, split = ' +', fixed = FALSE))), nrow = 1)
      options(warn = oldw)
      if(pos == 1) {
        pp.mat <- pp.tmp
      } else {
        if(ncol(pp.mat) != ncol(pp.tmp)) stop('\nError: The pPrey vector for ', tmp[1], ' has ', ncol(pp.tmp), ' columns and should have ', ncol(pp.mat))
        pp.mat <- rbind(pp.mat, pp.tmp)
      }
      pos    <- pos + 1
    }
    if(all(is.na(pp.mat[, 1]))) pp.mat <- pp.mat[, - 1]
    row.names(pp.mat)                  <- fg
    return(pp.mat)
  }
}

##' Raycasting Algorithm to find out whether a point is in a given polygon.
##' Performs the even-odd-rule Algorithm to find out whether a point is in a given polygon.
##' This runs in O(n) where n is the number of edges of the polygon.
##' @param polygon an array representation of the polygon where polygon[i,1] is the x Value of the i-th point and polygon[i,2] is the y Value.
##' @param point   an array representation of the point where point[1] is its x Value and point[2] is its y Value
##' @return whether the point is in the polygon (not on the edge, just turn < into <= and > into >= for that)
##' @author Javier Porobic - https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis
.point_in_polygon <- function(polygon, point){
  ## A point is in a polygon if a line from the point to infinity crosses the polygon an odd number of times
  odd = FALSE
  ## For each edge (In this case for each point of the polygon and the previous one)
  i = 0
  j = nrow(polygon) - 1
  while(i < nrow(polygon) - 1){
    i = i + 1
    ## If a line from the point into infinity crosses this edge
    ## One point needs to be above, one below our y coordinate
    ## ...and the edge doesn't cross our Y corrdinate before our x coordinate (but between our x coordinate and infinity)
    if (((polygon[i,2] > point[2]) != (polygon[j,2] > point[2]))
        && (point[1] < ((polygon[j,1] - polygon[i,1]) * (point[2] - polygon[i,2]) / (polygon[j,2] - polygon[i,2])) + polygon[i,1])){
      ## Invert odd
      odd = !odd
    }
    j = i
  }
  ## If the number of crossings was odd, the point is in the polygon
  return (odd)
}

##' @title Biomass for age groups
##' @param age.grp Age groups
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t Miligram to tons scalar
##' @param x.cn Carbon to Nitrogen transformation (Value/scalar)
##' @param inf.box Polygons Information
##' @param Time Time (Date) vector
##' @return a dataframe with the biomass for all the functional groups with age classes
##' @author Demiurgo - https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis
##' @editor Holly Perryman - 05.2022
##' 
bio.age <- function(age.grp, 
                    nc.out, 
                    ctg, 
                    mg2t, 
                    x.cn, 
                    inf.box, 
                    Time,
                    collapse.by.space = TRUE,
                    polygons.sub = c(),
                    path.bio.prm = c(),
                    path.harvest.prm = c()){
  
  if (!length(polygons.sub)==0) {
    bound   <- ifelse(inf.box$info[which(inf.box$info$Boxid %in% polygons.sub),"Depth_Bound"] > 0, 1, 0)
  } else {
    bound   <- ifelse(inf.box$info$Depth_Bound > 0, 1, 0)
  }
  
  grp.bio <- NULL
  
  for(age in 1 : nrow(age.grp)){ # age = 1 ; looping through fids 
    # message to user 
    print(paste("Getting data for group #",age," - ",
                round((age / nrow(age.grp)) * 100,2),
                "% Complete.",
                sep = ""))
    cohort <- NULL
    start.age <- 1
    # If the user wants SSB or harvestable biomass,
    # then re-set start.age based on the age_mat 
    # if bio prm file is provided we are computing SSB - set start.age appropriately 
    if (!is.null(path.bio.prm)) {
      # if harvest prm is also provided then we must error out as the function can only do SSB or Harvestable biomass, not both at the same time 
      if (!is.null(path.harvest.prm)) {stop('Must compute sSB or Harvestable Biomass, not both at the same time. Please select one, and re-run function.')}
      prm.text.age.mat <- grep(paste(age.grp$Name[age],"_age_mat",sep = ""), readLines(con = path.bio.prm, warn = FALSE), value = TRUE) # get age_mat value
      print(prm.text.age.mat)
      prm.text.age.mat <- gsub("^\\S+\\s+|\\s+\\S+$", "", prm.text.age.mat) # get age_mat value
      prm.text.age.mat <- strsplit(prm.text.age.mat," ") # get age_mat value
      start.age <- as.numeric(prm.text.age.mat[[1]][1]) # get age_mat value
    }
    
    # if harvest prm file provided, we are computing harvestable biomass - set start.age appropriately 
    if (!is.null(path.harvest.prm)) {
      index.line.num <- sum(grep(paste(age.grp$Name[age],"_mFC_startage",sep = ""), readLines(con = path.harvest.prm, warn = FALSE)) )
      index.line.num.data <- readLines(con = path.harvest.prm, warn = FALSE)[index.line.num + 1]
      index.line.num.data.numeric <- as.numeric(strsplit(index.line.num.data, "\\s+")[[1]])
      start.age <- (index.line.num.data.numeric[1]) 
      print(paste0(age.grp$Name[age]," - first harvestable age class: ",start.age))
    }
    
    # compute biomass based on start.age and polygons 
    for(coh in start.age : age.grp[age, 'NumCohorts']){ # coh = 1
      name.fg <- paste0(age.grp$Name[age], coh)
      b.coh   <- (ncdf4::ncvar_get(nc.out, paste0(name.fg, '_ResN'))  +
                    ncdf4::ncvar_get(nc.out, paste0(name.fg, '_StructN')))  *
        ncdf4::ncvar_get(nc.out, paste0(name.fg, '_Nums')) * mg2t * x.cn
      # insert break here if you want to subset the nc data by polygon
      # to compute biomass for specific grids
      if (!length(polygons.sub)==0) {
        b.coh <- b.coh[,polygons.sub+1,]# storage locations are boxid plus 1!! EG,box id 0 is stored in place 1 
      } 
      b.coh1 <- apply(b.coh, 2, colSums, na.rm = TRUE) # collapse by depth leaving time:polygon
      if (collapse.by.space) {
        # clean up data
        b.coh2 <- colSums(apply(b.coh1, 1, function(x) x * bound), na.rm = TRUE) # collapse by polygon leaving time
        # store data 
        cohort  <- cbind(cohort, b.coh2)
      } else {
        # clean up data
        b.coh1<-as.data.frame(b.coh1)
        names(b.coh1)<-as.character(c(0:65))
        b.coh1$Time <- Time
        b.coh2 <- melt(b.coh1, id.vars=c("Time"))
        b.coh2$age <- coh-1
        names(b.coh2) <- c("Time", "box_id", "Biomass", "age")
        # store data 
        cohort  <- rbind(cohort, b.coh2)
      }
      remove(b.coh,b.coh1,b.coh2)
    }
    # generate function output
    if (collapse.by.space) {
      grp.bio <- rbind(grp.bio, data.frame(Time = Time,
                                           FG = paste0(as.character(age.grp$Code[age])),
                                           FG.long = paste0(as.character(age.grp$Long.Name[age])),
                                           Biomass  = rowSums(cohort, na.rm  = TRUE), 
                                           Simulation = ctg))
    } else {
      cohort$FG         <- paste0(as.character(age.grp$Code[age]))
      cohort$FG.long    <- paste0(as.character(age.grp$Long.Name[age]))
      cohort$Simulation <- ctg
      #
      grp.bio <- rbind(grp.bio, cohort)
    }
    
  }
  if (!is.null(path.bio.prm)) {names(grp.bio) <- c("Time", "FG", "FG.long", "SSB", "Simulation")}
  #if (!is.null(path.harvest.prm)) {names(grp.bio) <- c("Time", "FG", "FG.long", "Harvestable_B", "Simulation")}
  
  return(grp.bio)
}

##' @title Biomass for age groups
##' @param pol.grp Biomass pool groups
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t Miligrams to tons scalar
##' @param x.cn Ration from nitrogen to carbon
##' @param box.info Information by box and layer
##' @param Time Date Vector
##' @return a dataframe with the biomass for all the functional groups with age classes
##' @author Demiurgo - https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis
##' @editor Holly Perryman - 06.2022
##' 
bio.pool <- function(pol.grp, 
                     nc.out, 
                     ctg, 
                     mg2t, 
                     x.cn, 
                     box.info, # i still need this for area
                     Time, 
                     collapse.by.space = TRUE,
                     inf.include.water.column = FALSE, # biomass data from the log and bioindx.txt file seem to only report what is in the sediment layer
                     polygons.sub = c()){
  
  # pol.grp <- pol.grp
  # nc.out <- nc.out
  # ctg <- "TEST"
  # mg2t <- mg2t
  # x.cn <- x.cn
  # box.info <- inf.box
  # Time <- Time
  # polygons.sub <- polygons.USA
  
  
  #'Get volume data from nc file for converting from per m3
  volume    <- ncdf4::ncvar_get(nc.out, "volume") # vol of each cell at each time step, dim - layer:box:time
  
  #' prep function output df
  grp.bio <- NULL
  
  #' loop through groups and build output data
  for(pool in 1 : nrow(pol.grp)){   # pool<-20
    name.fg <- paste0(pol.grp$Name[pool], '_N')
    N.m3    <- ncdf4::ncvar_get(nc.out, name.fg)
    
    if(length(dim(N.m3))== 3){ #' if the dim of N data is three then the dim is layer:box:time - for debugging : pool<-9
      
      #' Holly note - june 2022
      #' 
      #' The OG script used Vol data gathered by the boxes.prop() function 
      #' BUT!
      #' I found a few errors with this:
      #'     1) it had NAs in it, which did not make much sense... 
      #'     2) the layer order did not make sense, it did not match input nc or output nc layer order
      #'     3) vol = VolInf, except one of the layers was zeroed out and I am not sure why (it was not the sediment layer)
      #'     4) VolInf zeroed out the sediment layer which ended up computing almost 0 biomass for benthic groups
      #'     5) cumdepth was used to compute volume when it should be depth of the individual cells (nominal_dz in the nc out file)
      #'
      #' Rather than fixing that boxes.prop() function, I am making changes here as this is already calling the nc out file 
      #' and this file contains volumne data for converting N per m3 to N. 
      #' see 5.6. Outputs of physics tracers in the output.nc file in manual 
      
      if(pol.grp$InvertType[pool] == 'LG_INF'){
        # data reported in log and bioindx.txt for INF grp reports back just the biomass data in the sediment
        # thus there is a flag to also drop this data, thus matching these outputs, or to include it
        if (inf.include.water.column) {
          N <- N.m3 * volume                      # dim layer:box:time
          N <- apply(N, 2, colSums, na.rm = TRUE) # dim time:box
          N <- t(N)                               # dim box:time
        } else {
          N <- N.m3 * volume # dim layer:box:time
          N <- N[7,,]        # select data from the sediment layer - dim box:time
        }
      }else{
        #N.m3 <- apply(N.m3, 3, '*', box.info$Vol) # Holly: the vol data from boxes.prop() are not correct... 
        N <- N.m3 * volume                      # dim layer:box:time
        N <- apply(N, 2, colSums, na.rm = TRUE) # dim time:box
        N <- t(N)                               # dim box:time
      }
      # Remove polygons before summing!
      if (!length(polygons.sub)==0) {N <- N[polygons.sub+1,]}# storage locations are boxid plus 1!! EG, box id 0 is stored in place 1 
      # convert to biomass
      biom    <- N * mg2t * x.cn 
      # prep data for storage
      if (collapse.by.space) {
        biom <- apply(biom, 2, sum, na.rm = TRUE) # sum across box 
      } else {
        biom        <- t(biom)
        biom        <- as.data.frame(biom)
        names(biom) <- as.character(c(0:65))
        biom$Time   <- Time
        biom        <- melt(biom, id.vars=c("Time"))
        names(biom) <- c("Time", "box_id", "Biomass")
      }
      
    } else {  #' if the dim of N data is box:time
      
      N <- apply(N.m3, 2, function(x) x * box.info$info$Area) 
      # Remove polygons before summing!
      if (!length(polygons.sub)==0) {N <- N[polygons.sub+1,]}# storage locations are boxid plus 1!! EG,box id 0 is stored in place 1 
      # convert to biomass 
      biom    <- N * mg2t * x.cn 
      # prep data for saving to output df 
      if (collapse.by.space) {
        biom <- apply(biom, 2, sum, na.rm = TRUE)
      } else {
        biom        <- t(biom) 
        biom        <- as.data.frame(biom)
        names(biom) <- as.character(c(0:65))
        biom$Time   <- Time
        biom        <- melt(biom, id.vars=c("Time"))
        names(biom) <- c("Time", "box_id", "Biomass")
      }
    }
    #biom[1] <- biom[2] # I have no flipping clue why this is here... 
    #
    # save data to output df:
    if (collapse.by.space) {
      grp.bio <- rbind(grp.bio, data.frame(Time = Time, 
                                           FG = paste0(as.character(pol.grp$Code[pool])),
                                           FG.long = paste0(as.character(pol.grp$Long.Name[pool])),
                                           Biomass  = biom, 
                                           Simulation = ctg))
    } else {
      biom$FG         <- paste0(as.character(pol.grp$Code[pool]))
      biom$FG.long    <- paste0(as.character(pol.grp$Long.Name[pool]))
      biom$Simulation <- ctg
      grp.bio <- rbind(grp.bio, biom)
    }
  }
  return(grp.bio)
}

##' @title Biomass for age groups
##' @param pwn.grp Biomass pools with age classes
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t mg C converted to wet weight in tonnes == 20 / 1000000000
##' @param x.cn Redfield ratio of C:N 5.7
##' @param box.info Information for each box and layer
##' @param Time Dates vector
##' @return a dataframe with the biomass for all the functional groups with Biomass pool age classes
##' @author Demiurgo - https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis
##' @editor Holly Perryman - 05.2022
##' 
bio.pwn <- function(pwn.grp, 
                    nc.out, 
                    ctg, 
                    mg2t, 
                    x.cn, 
                    box.info,
                    Time, 
                    collapse.by.space = TRUE,# collapse.by.space = FALSE
                    polygons.sub = c()){
  
  #'Get volume data from nc file for converting from per m3
  volume    <- ncdf4::ncvar_get(nc.out, "volume") # vol of each cell at each time step, dim - layer:box:time
  
  #' initialize data returned by function:
  grp.bio <- NULL
  
  # loop through groups 
  for(pwn in 1 : nrow(pwn.grp)){ # pwn = 1
    #' prep placeholder for age data 
    cohort <- NULL
    
    #' loop through cohorts
    for(coh in 1 : pwn.grp[pwn, 'NumCohorts']){ # debugging --- coh <- 1
      
      name.fg <- paste0(pwn.grp$Name[pwn], '_N', coh)
      N.m3 <- ncdf4::ncvar_get(nc.out, name.fg) # dim - box:time
      
      if(length(dim(N.m3)) > 2){
        ## units -- N per m3 (volume) 
        # OG code : apply(b.coh, 3, '*', box.info$Vol) --- I had issues with vol data coming back from boxes.prop()
        N <- N.m3 * volume                      # dim layer:box:time
        N <- apply(N, 2, colSums, na.rm = TRUE) # dim time:box
        N <- t(N)                               # dim box:time
      } else {
        ## units -- N per m2 (area) 
        N <- apply(N.m3, 2, '*', box.info$info$Area) # dim --- box:time
      }
      # Remove polygons before summing!
      if (!length(polygons.sub)==0) {N <- N[polygons.sub+1,]}# storage locations are boxid plus 1!! EG,box id 0 is stored in place 1 
      # convert to biomass
      b.coh <- N * mg2t * x.cn 
      # prep data for output 
      if (collapse.by.space) {
        b.coh   <- apply(b.coh, 2, sum, na.rm = TRUE) # dim --- time 
        cohort  <- cbind(cohort, b.coh)
      } else {
        b.coh        <- t(b.coh) 
        b.coh        <- as.data.frame(b.coh)
        names(b.coh) <- as.character(c(0:65))
        b.coh$Time   <- Time
        b.coh        <- melt(b.coh, id.vars=c("Time"))
        names(b.coh) <- c("Time", "box_id", "Biomass")
        b.coh$age    <- coh
        # store data 
        cohort  <- rbind(cohort, b.coh)
      }
    }
    
    if (collapse.by.space) {
      grp.bio <- rbind(grp.bio, data.frame(Time = Time,
                                           FG = paste0(as.character(pwn.grp$Code[pwn])),
                                           FG.long = paste0(as.character(pwn.grp$Long.Name[pwn])),
                                           Biomass  = rowSums(cohort, na.rm  = TRUE), 
                                           Simulation = ctg))
    } else {
      cohort$FG         <- paste0(as.character(pwn.grp$Code[pwn]))
      cohort$FG.long    <- paste0(as.character(pwn.grp$Long.Name[pwn]))
      cohort$Simulation <- ctg
      #
      grp.bio <- rbind(grp.bio, cohort)
    }
    
  }
  
  return(grp.bio)
}
