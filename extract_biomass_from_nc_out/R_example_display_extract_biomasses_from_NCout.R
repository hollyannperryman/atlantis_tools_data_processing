



####################' needed from user running the script 

# path to model run inputs
dir.inputs <- "C:/Users/haperryman/Documents/GOM/Inputs_GOMAtlantis/GOMAtlantis_CIEreview_prep/"
# path to model run outputs
dir.outputs <- "C:/Users/haperryman/Documents/GOM/outputs_from_ubuntu/HCRtesting/HCR_on_totssb/"
# path save location
dir.save <- paste(dir.outputs,"R_output/",sep = "")
# make this directory if it does not already exist:
ifelse(!dir.exists(file.path(dir.save)), dir.create(file.path(dir.save)), FALSE)
# name of groups.csv input file
input.grps <- "GulfGroups_6665.csv"
# name of groups.csv input file
input.bio <- "gom_prm_2021_sg4_pprey3.prm"
# name of groups.csv input file
input.harvest <- "at_harvest_gom2020_V2.prm"
# name of bgm output file
output.bgm <- "GOM_BGM.bgm"

####################' libraries

library(reshape)
library(dplyr)

####################' additional functions

source("R_tools_from_ReactiveAtlantis_new.R") # misc functions

####################' get data

# groups data
dir.grps.csv <- list.files(path = dir.inputs, 
                           full.names = TRUE, 
                           recursive = F, 
                           pattern = input.grps)[1]
grp <- utils::read.csv(dir.grps.csv)

# biology data
dir.bio.prm <- list.files(path = dir.inputs, 
                           full.names = TRUE, 
                           recursive = F, 
                           pattern = input.bio)[1]

# harvest data
dir.harv.prm <- list.files(path = dir.inputs, 
                          full.names = TRUE, 
                          recursive = F, 
                          pattern = input.harvest)[1]
# bgm data
dir.bgm.file <- list.files(path = dir.inputs,
                           full.names = TRUE,
                           recursive = F,
                           pattern = output.bgm)[1]
inf.box <- boxes.prop(dir.bgm.file) # function from ReactiveAtlantis, I have edited it to isolate what is need (see source)

# output data from nc file
dir.out.nc <- list.files(path = dir.outputs,
                         full.names = TRUE,
                         recursive = F,
                         pattern = ".nc")[1]
nc.out <- ncdf4::nc_open(dir.out.nc)
Time   <- nc.out$dim$t$vals / (60*60*24) # convert from seconds to days

####################' 
####################' get biomass for age structured groups
####################' 

#' This is done using the function bio.age()
#' 
#' (1) This function allows for the user to choose which biomass metric is being returned:
#'   (a) total biomass
#'   (b) ssb
#'   (c) mature biomass
#'   (d) harvestable biomass
#'   
#'   The metric being returned depends on the specifications of the flags.
#'   (i) 
#'   flag.return.whole.system.biomass # T - return whole system biomass metric, F - return in system biomass metric
#'   (ii)
#'   path.bio.prm - if provided, the function will compute either mature biomass or SSB
#'   flag.return.SSB # 1 - use FSPB to compute SSB, 0 - use age_mat to return mature biomass
#'   (iii)
#'   path.harvest.prm - if provided, the function will use _mFC_startage to return harvestable biomass
#' 
#' (2) This function allows the user to choose to return
#'   (a) whole system biomass metric
#'   (b) in system biomass metric
#'   (c) sub system biomass metric


#' First, this function depends on the groups.csv data isolating that appropriate groups: 
# age structured groups:
age.grp <- grp[grp$NumCohorts  > 1 & !grp$InvertType %in% c('PWN', 'PRAWNS', 'PRAWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE', 'MANGROVES', 'SPONGE'), ]

#' for the exercise, I am going to isolate one functional group, but the function
#' will process all age-structured groups provided. 

age.grp.sub <- age.grp[which(age.grp$Code %in% c("BSH")),] # isolating group BSH

#####'
#####' First, lets go over how to get different biomass metrics from the function
#####'

#' 
#' (A) In-system total biomass over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time)
head(nc.age.bio)

#' 
#' (B) whole-system total biomass over time (if you group is not migrating, then I think this will be the same as the above)
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       flag.return.whole.system.biomass = TRUE)
head(nc.age.bio)

#' 
#' (C) sub-system total biomass over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       polygons.sub = c(5,6,7,8,9,10)) # !!! VALUES ARE BOX_ID
head(nc.age.bio)

#' 
#' (D) In-system SSB over time (based on FSPB)
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       path.bio.prm = dir.bio.prm)
head(nc.age.bio)

#' 
#' (E) whole-system SSB over time (based on FSPB)
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       flag.return.whole.system.biomass = TRUE,
                       path.bio.prm = dir.bio.prm)
head(nc.age.bio)

#' 
#' (F) sub-system SSB over time (based on FSPB)
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       polygons.sub = c(5,6,7,8,9,10), # !!! VALUES ARE BOX_ID
                       path.bio.prm = dir.bio.prm) 
head(nc.age.bio)

#' 
#' (G) In-system mature biomass over time (based on age_mat)
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       path.bio.prm = dir.bio.prm,
                       flag.use.FSPB = F)
head(nc.age.bio)

#' 
#' (H) whole-system mature biomass over time (based on age_mat)
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       flag.return.whole.system.biomass = TRUE,
                       path.bio.prm = dir.bio.prm,
                       flag.use.FSPB = F)
head(nc.age.bio)

#' 
#' (I) sub-system mature biomass over time (based on age_mat)
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       polygons.sub = c(5,6,7,8,9,10), # !!! VALUES ARE BOX_ID
                       path.bio.prm = dir.bio.prm,
                       flag.use.FSPB = F) 
head(nc.age.bio)

#' 
#' (J) In-system harvestable biomass over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       path.harvest.prm = dir.harv.prm)
head(nc.age.bio)

#' 
#' (K) whole-system harvestable biomass over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       flag.return.whole.system.biomass = TRUE,
                       path.harvest.prm = dir.harv.prm)
head(nc.age.bio)

#' 
#' (L) sub-system harvestable biomass over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       polygons.sub = c(5,6,7,8,9,10), # !!! VALUES ARE BOX_ID
                       path.harvest.prm = dir.harv.prm) 
head(nc.age.bio)


#####'
#####' Now, lets go over how to get age-specific and/or polygon-specific biomass
#####'

#' 
#' (A) biomass over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time)
head(nc.age.bio)

#' 
#' (B) biomass at age over time
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       return.collapse.by.cohort = F)
head(nc.age.bio)

#' 
#' (C) biomass per polygon over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       return.collapse.by.space = F)
head(nc.age.bio)

#' 
#' (D) biomass at age per polygon over time 
#' 
nc.age.bio  <- bio.age(age.grp = age.grp.sub,
                       nc.out  = nc.out,
                       inf.box = inf.box,
                       Time    = Time,
                       return.collapse.by.space = F,
                       return.collapse.by.cohort = F)
head(nc.age.bio)



####################' 
####################' get biomass for biomass pool groups
####################' 

#' First, this function depends on the groups.csv data isolating that appropriate groups: 
# biomass pools
pol.grp <- grp[grp$NumCohorts == 1, ]

#####'
#####' First, lets go over how to get different biomass metrics from the function
#####'

#' 
#' (A) In-system total biomass over time 
#' 
nc.pool.bio  <- bio.pool(pol.grp  = pol.grp,
                         nc.out   = nc.out,
                         box.info = inf.box,
                         Time     = Time)
head(nc.pool.bio)

#' 
#' (B) whole-system total biomass over time (if you group is not migrating, then I think this will be the same as the above)
#' 
nc.pool.bio  <- bio.pool(pol.grp  = pol.grp,
                         nc.out   = nc.out,
                         box.info = inf.box,
                         Time     = Time,
                         flag.return.whole.system.biomass = TRUE)
head(nc.pool.bio)

#' 
#' (C) sub-system total biomass over time 
#' 

nc.pool.bio  <- bio.pool(pol.grp  = pol.grp,
                         nc.out   = nc.out,
                         box.info = inf.box,
                         Time     = Time,
                         polygons.sub = c(5,6,7,8,9,10)) # !!! VALUES ARE BOX_ID
head(nc.pool.bio)

#####'
#####' Now, lets go over how to get age-specific and/or polygon-specific biomass
#####'

#' 
#' (A) biomass over time 
#' 
nc.pool.bio  <- bio.pool(pol.grp  = pol.grp,
                         nc.out   = nc.out,
                         box.info = inf.box,
                         Time    = Time)
head(nc.pool.bio)

#' 
#' (B) biomass per polygon over time
#' 
nc.pool.bio  <- bio.pool(pol.grp  = pol.grp,
                         nc.out   = nc.out,
                         box.info = inf.box,
                         Time     = Time,
                         return.collapse.by.space = F)
head(nc.pool.bio)



####################' 
####################' get biomass for age structured pool groups
####################' 

#' First, this function depends on the groups.csv data isolating that appropriate groups: 
# age structured biomass pools
pwn.grp <- grp[grp$NumCohorts  > 1 & grp$InvertType %in% c('PWN', 'PRAWNS', 'PRAWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE', 'MANGROVES', 'SPONGE'), ]


#####'
#####' First, lets go over how to get different biomass metrics from the function
#####'

#' 
#' (A) In-system total biomass over time 
#' 
nc.pwn.bio  <- bio.pwn(pwn.grp  = pwn.grp,
                       nc.out   = nc.out,
                       box.info = inf.box,
                       Time     = Time)
head(nc.pwn.bio)

#' 
#' (B) whole-system total biomass over time (if you group is not migrating, then I think this will be the same as the above)
#' 
nc.pwn.bio  <- bio.pwn(pwn.grp  = pwn.grp,
                       nc.out   = nc.out,
                       box.info = inf.box,
                       Time     = Time,
                       flag.return.whole.system.biomass = TRUE)
head(nc.pwn.bio)

#' 
#' (C) sub-system total biomass over time 
#' 

nc.pwn.bio  <- bio.pwn(pwn.grp  = pwn.grp,
                       nc.out   = nc.out,
                       box.info = inf.box,
                       Time     = Time,
                       polygons.sub = c(5,6,7,8,9,10)) # !!! VALUES ARE BOX_ID
head(nc.pwn.bio)

#####'
#####' Now, lets go over how to get age-specific and/or polygon-specific biomass
#####'

#' 
#' (A) biomass over time 
#' 
nc.pwn.bio  <- bio.pwn(pwn.grp  = pwn.grp,
                       nc.out   = nc.out,
                       box.info = inf.box,
                       Time     = Time)
head(nc.pwn.bio)

#' 
#' (B) biomass at age over time
#' 
nc.pwn.bio  <- bio.pwn(pwn.grp  = pwn.grp,
                       nc.out   = nc.out,
                       box.info = inf.box,
                       Time     = Time,
                       return.collapse.by.cohort = F)
head(nc.pwn.bio)

#' 
#' (C) biomass per polygon over time 
#' 
nc.pwn.bio  <- bio.pwn(pwn.grp  = pwn.grp,
                       nc.out   = nc.out,
                       box.info = inf.box,
                       Time     = Time,
                       return.collapse.by.space = F)
head(nc.pwn.bio)

#' 
#' (D) biomass at age per polygon over time 
#' 
nc.pwn.bio  <- bio.pwn(pwn.grp  = pwn.grp,
                       nc.out   = nc.out,
                       box.info = inf.box,
                       Time     = Time,
                       return.collapse.by.space = F,
                       return.collapse.by.cohort = F)
head(nc.pwn.bio)




