

#' @title Spatially explicit prey compositions againsy predator distribution   
#' @author Holly Perryman - May 2024 
#' 
#' The following script produces a map of piecharts that indicate the 
#' the biomass composition of prey against the spatial distribution of the 
#' selected predator



#* ------------------------------
#* libraries

library(ncdf4)      # needed for self-made functions for open nc files
library(rbgm)       # to read bgm file
library(dplyr)      # need for left_join()
library(broom)      # need for tidy()
library(ggplot2)    # need for ggplot
library(stringr)    # needed for str_detect()
library(scatterpie) # spatial pie charts 
library(ggnewscale) # new scale in ggplot - very important 
library(reshape2)   # need to shape between wide and long 
library(tidyr )     # gather() 
library(RColorBrewer)


#* ------------------------------
#* source code

#* The following is a source file I developed heavily based on functions from ReactiveAtlantis
#* https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis

source("R_tools_from_ReactiveAtlantis.R")  # this has the function text2num() 


#* ------------------------------ (NEED TO UPDATE FOR YOUR FILES)
#* Paths to input data and output data
#* 
path.input.files <- "C:/Users/haperryman/Documents/GOM/Inputs_GOMAtlantis/GOMAtlantis_CIEreview/"
path.out.files   <- "C:/Users/haperryman/Documents/GOM/outputs_from_ubuntu/piemaps/"


#* ------------------------------ (NEED TO UPDATE FILE NAMES IF NECESSARY)
#* paths to key input and output files 
#* 

#' inputs 
path.grps <- list.files(path = path.input.files,
                        full.names = TRUE,
                        recursive = F,
                        pattern = "GulfGroups_6665.csv")[1]

path.bio.prm <- list.files(path = path.input.files, 
                           full.names = TRUE, 
                           recursive = F, 
                           pattern = "gom_prm_2021_sg5.prm")[1]
#' outputs 
path.bgm <- list.files(path = path.input.files,
                       full.names = TRUE,
                       recursive = F,
                       pattern = ".bgm")[1]

path.nc.file <- list.files(path = path.out.files,
                           full.names = TRUE,
                           recursive = F,
                           pattern = ".nc")[1]
#'
#' the goal was to plot spatial diet interactions, I just did not get there yet
#' If you do, please let me know
#' 
# path.dietcheck <- list.files(path = path.out.files, 
#                              full.names = TRUE, 
#                              recursive = F, 
#                              pattern = "DietCheck.txt")[1]



###* ------------------------------ 
###* Build spatial biomass data frame  
###*  

#' ----- get data

#' bgm data
#' 
#' 1) use boxes.prop() - I have edited the function so that you no longer have to read in cum.depths  # cum.depths <- c(0, 10, 20, 50, 200, 2000, 4000) # define cumulative depths - needed for some atlantis tools 
inf.box <- boxes.prop(path.bgm)


#' get group data - need to subset by InvertType in order to use reactiveatlantis functions to compute biomass
#' 
# read data file
grp     <- utils::read.csv(path.grps)
# age structured groups 
age.grp <- grp[grp$NumCohorts  > 1 & !grp$InvertType %in% c('PWN', 'PRAWNS', 'PRAWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE', 'MANGROVES', 'SPONGE'), ]
# biomass pools 
pol.grp <- grp[grp$NumCohorts  == 1, ]
# age structured biomass pools 
pwn.grp <- grp[grp$NumCohorts  > 1 & grp$InvertType %in% c('PWN', 'PRAWNS', 'PRAWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE', 'MANGROVES', 'SPONGE'), ]


#' nc out data
#' 
nc.cur <- ncdf4::nc_open(path.nc.file) # names(nc.cur$var)
Time   <- nc.cur$dim$t$vals / (60*60*24)  # convert from seconds to days

#' ----- calc biomass

# Holly - 06.2022 
# I used functions from ReactiveAtlantis to get this data.
#
# I have edited these functions in order to extract specific biomass metrics,
# but the script below is currently extracting biomass by cohort and polygon
#
# I doubled check the results to data from biomindx.txt and the log out plots
# and everything lines up - so these functions seem to be working correctly

age.cur.bio.all  <- bio.age(age.grp, 
                            collapse.by.space = FALSE, # set to TRUE, and will return total biomass collapsed by space and cohort
                            nc.cur, 
                            'Whole System', 
                            0.00000002, # mgC converted to wet weight in tonnes = 20/1000000000, 
                            5.7,   	    # Redfield ratio of C:N 5.7 
                            inf.box, 
                            Time)
head(age.cur.bio.all) 


pool.cur.bio.all  <- bio.pool(pol.grp, 
                              collapse.by.space = FALSE, # set to TRUE, and will return total biomass collapsed by space and cohort
                              nc.cur, 
                              'Whole System',  
                              0.00000002, # mgC converted to wet weight in tonnes = 20/1000000000, 
                              5.7,   	    # Redfield ratio of C:N 5.7 
                              #inf.include.water.column = TRUE, # un-comment this to report back INF data that is also in the water column
                              inf.box, 
                              Time)
#' Notes - holly 06.2022
#' 
#' I checked and the data coming out of this matches that from bioindx.txt:
#'
# data.biomass <- read.csv(list.files(path = path.out.files, 
#                                     full.names = TRUE, 
#                                     recursive = F, 
#                                     pattern = "BiomIndx.txt")[1], 
#                          sep="", stringsAsFactors=FALSE)
# index <- 9
# data.biomass[,c("Time",pol.grp[index,"Name"])]
# pool.cur.bio.all[which(pool.cur.bio.all$FG == pol.grp[index,"Name"]),]
#'
#' an important thing to remember is that for INF groups, the data from log and 
#' bioindx.txt are just what is in the sediment layer. This function assumes 
#' does this as well, but there is an option to include data in the water column 
#' 
#' other notes from checking groups:
#' 
# 17 - I think LPP is right, just very dynamic within 1 time step
head(pool.cur.bio.all) 


pwn.cur.bio.all <- bio.pwn(pwn.grp, 
                           collapse.by.space = FALSE,
                           nc.cur, 
                           'Whole System',  
                           0.00000002,  # mgC converted to wet weight in tonnes = 20/1000000000, 
                           5.7,   	    # Redfield ratio of C:N 5.7 
                           inf.box, 
                           Time)
# data.biomass[,c("Time",pwn.grp[1,"Name"])]
head(pwn.cur.bio.all) 



#' ----- merge

pool.cur.bio.all$age <- NA

biomass.space <- rbind(age.cur.bio.all,
                       pool.cur.bio.all,
                       pwn.cur.bio.all)

remove(age.cur.bio.all,
       pool.cur.bio.all,
       pwn.cur.bio.all)


###* ------------------------------ 
###* For the specified predator and time step,
###* make a data frame indicating predator biomass by polygon
###* 
###* this DF has been revise so it now indicates:
###*    total predator biomass, total juv. predator biomass, and total ad. predator biomss
###* by polygon
###* 
###*  NOTE!!!!!
###*  the distinction for juv vs ad was made using age_mat
###*  
###*  if you want SSB, consider using FSPB instead
###*  !!!!!
###*  

#' ----- specify fid for predator and timestep of focus (NEED TO UPDATE FOR YOUR NEEDS)
#' 
index.pred.fid <- "WSH"
index.ts  <- 390

#' ----- specify adult adult vs juv
#' 
# get age at maturity for plotting juv and ad diet 
prm.text.age.mat <- grep(paste(index.pred.fid,"_age_mat",sep = ""), readLines(con = path.bio.prm, warn = FALSE), value = TRUE) # get age_mat value
prm.text.age.mat <- gsub("^\\S+\\s+|\\s+\\S+$", "", prm.text.age.mat) # get age_mat value
start.age <- as.numeric(prm.text.age.mat) # get age_mat value (0-9)
remove(prm.text.age.mat) # purge 
age.juv <- c(0:(start.age - 1))
age.ad <- c(start.age:9)
remove(start.age) # purge 

#' ----- subset spatial biomass data frame by juv and ad cohorts
#' 
biomass.space.sub  <- biomass.space[which(biomass.space$FG == index.pred.fid &
                                           biomass.space$Time == index.ts),]
biomass.space.subj <- biomass.space[which(biomass.space$FG == index.pred.fid &
                                            biomass.space$age %in% age.juv &
                                            biomass.space$Time == index.ts),]
biomass.space.suba <- biomass.space[which(biomass.space$FG == index.pred.fid &
                                            biomass.space$age %in% age.ad &
                                            biomass.space$Time == index.ts),]
head(biomass.space.sub) 
head(biomass.space.subj) 
head(biomass.space.suba) 

#' ----- aggregate to get predator biomass by polygon
#' 
biomass.space.sub.sum <- aggregate(biomass.space.sub$Biomass, 
                                   by=list(box_id=biomass.space.sub$box_id), 
                                   FUN=sum)
names(biomass.space.sub.sum) <- c("box_id", "biomass")
biomass.space.subj.sum <- aggregate(biomass.space.subj$Biomass, 
                                    by=list(box_id=biomass.space.subj$box_id), 
                                    FUN=sum)
names(biomass.space.subj.sum) <- c("box_id", "biomass_juv")
biomass.space.suba.sum <- aggregate(biomass.space.suba$Biomass, 
                                    by=list(box_id=biomass.space.suba$box_id), 
                                    FUN=sum)
names(biomass.space.suba.sum) <- c("box_id", "biomass_ad")

#' ----- merge data for plot
#'
biomass.space.sub.forplot <- merge(biomass.space.sub.sum, 
                                   biomass.space.subj.sum, 
                                   by = 'box_id')
biomass.space.preddist <- merge(biomass.space.sub.forplot, 
                           biomass.space.suba.sum, 
                           by = 'box_id')
head(biomass.space.preddist)
remove(biomass.space.sub,biomass.space.subj,biomass.space.suba,
       biomass.space.sub.sum,biomass.space.subj.sum,biomass.space.suba.sum,
       biomass.space.sub.forplot)


###* ------------------------------ 
###* For the specified predator, determine the prey groups 
###* 
###* 

#' ----- extract pprey data from prm
#' 
data.bio.prm  <- readLines(path.bio.prm) #' get data in bio prm file 
pPREY <- text2num(data.bio.prm, 'pPREY', Vector=TRUE, pprey = TRUE)
pPREY <- as.data.frame(pPREY)
colnames(pPREY)  <- c(as.character(grp$Code), 'DLsed', 'DRsed', 'DCsed')
#Ava.mat2 <- text2num(data.bio.prm, 'pprey', Vector=TRUE, pprey = FALSE) # predation ON seagrass

#' ----- Get diet parameterization for specified predator
# all prey items for fid: 
pPREY.sub <- pPREY[str_detect(rownames(pPREY),index.pred.fid),]
pPREY.sub <- pPREY.sub[,which(colSums(pPREY.sub)>0)] # isolate prey groups 
#pPREY.sub <- (pPREY.sub / rowSums(pPREY.sub)*100)   # scale so all rows sum to one (1) 
grps.prey <- names(pPREY.sub)

# # make a barplot
# pPREY.sub$variable <- rownames(pPREY.sub)
# pPREY.sub.forplot <- gather(pPREY.sub, 
#                             prey, 
#                             value, 
#                             names(pPREY.sub)[1]:names(pPREY.sub)[length(names(pPREY.sub))-1], 
#                             factor_key=TRUE)


# ---- sub biomass.space data by prey groups 

# all prey

#grps.prey.sub  <- grps.prey[!((grps.prey %in% c('DLsed', 'DRsed', 'DCsed')))] 
biomass.space.sub.prey <- biomass.space[which(biomass.space$FG %in% grps.prey &
                                                biomass.space$Time == index.ts),]
biomass.space.sub.prey.sum <- aggregate(biomass.space.sub.prey$Biomass, 
                                        by=list(box_id = biomass.space.sub.prey$box_id, 
                                                FG     = biomass.space.sub.prey$FG), 
                                        FUN=sum)
names(biomass.space.sub.prey.sum) <- c("box_id", "FG", "Biomass")
head(biomass.space.sub.prey.sum)
# convert
biomass.space.preycomp <- spread(biomass.space.sub.prey.sum, FG, Biomass)
head(biomass.space.preycomp)

remove(biomass.space.sub.prey.sum,biomass.space.sub.prey)


###* ------------------------------ 
###* Prep data to make a spatial map of the predator distribution against prey
###* polygon specific compositions
###* 

#' ----- prep spatial data
bgm <- bgmfile(path.bgm)
bgm.boxSpatial.0 <- boxSpatial(bgm) # first, need to pass data through boxSpatial() funciton to isolate spaital data 
bgm.boxSpatial <- spTransform(bgm.boxSpatial.0, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) # transform bgm data to match the baseline map (to be made later)
attributes                 <- c("botz", "area", "box_id", "boundary") # key attributes for plotting 
bgm.boxSpatial.sub         <- bgm.boxSpatial[,attributes]             # subset spatial data by attributes 
bgm.boxSpatial.sub@data$id <- rownames(bgm.boxSpatial.sub@data)       # match row names 
bgm.boxSpatial.sub@data[bgm.boxSpatial.sub@data$boundary == TRUE,"botz"] <- NA # optional - I like to set depth of boundary boxes to NAs for cleaner plots 

remove(bgm,bgm.boxSpatial.0,bgm.boxSpatial,attributes)

#' ----- merge spatial data to PREDATOR biomass data - FOR geom_polygon PART OF PLOT
# add / merge data(s) for plotting: 
bgm.boxSpatial.sub.merge <- merge(bgm.boxSpatial.sub, 
                                  biomass.space.preddist, 
                                  by="box_id")
head(bgm.boxSpatial.sub.merge)
# use the tidy() function in broom to convert boxSpatial data
bgm.boxSpatial.sub.merge.tidy <- tidy(bgm.boxSpatial.sub.merge)
# use left_join in dplry to merge and create the data for ggplot
plotdata_preddist <- left_join(bgm.boxSpatial.sub.merge.tidy,
                       bgm.boxSpatial.sub.merge@data,
                       by="id")
head(plotdata_preddist)

remove(bgm.boxSpatial.sub.merge,bgm.boxSpatial.sub.merge.tidy)

#' ----- Find centroids for labels
centroids <- aggregate(cbind(long,lat) ~ id, data=plotdata_preddist, FUN=mean)

#' ----- merge spatial data to PREY data - FOR geom_scatterpie PART OF PLOT

# all prey

# # add / merge data(s) for plotting: 
bgm.boxSpatial.sub.merge <- merge(bgm.boxSpatial.sub, 
                                  biomass.space.preycomp, 
                                  by="box_id")
head(bgm.boxSpatial.sub.merge)
# next, use the tidy() function in broom to convert boxSpatial data
bgm.boxSpatial.sub.merge.tidy <- tidy(bgm.boxSpatial.sub.merge)
# next, use left_join in dplry to merge
plotdata <- left_join(bgm.boxSpatial.sub.merge.tidy,
                      bgm.boxSpatial.sub.merge@data,
                      by="id");head(plotdata)
# remove duplicates - cleaner image
plotdata <- plotdata[!duplicated(plotdata$id), ] 
# replace lat and long with centroids 
plotdata <- left_join(plotdata[,-c(1,2)],
                       centroids,
                       by="id");head(plotdata)
# remove levels
biomass.space.preddist$box_id <- as.numeric(as.character(biomass.space.preddist$box_id))
# merge
plotdata_preycomp <- left_join(plotdata,
                               biomass.space.preddist,
                               by="box_id");head(plotdata_preycomp)

remove(bgm.boxSpatial.sub.merge,bgm.boxSpatial.sub.merge.tidy,plotdata)

#' ----- palette 
getPalette = colorRampPalette(brewer.pal(9, "Set1"))



###* ------------------------------ 
###* 
###* MAKE BASE PLOT
###* 

# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
#install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
#install.packages("tools")
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library(rbgm)  # to read bgm file
library(broom) # need for tidy()
library(dplyr)      # need for left_join()

#'
#' World data
#'
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
#'
#' States data
#'
library("maps")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
head(states)

#states <- cbind(states, st_coordinates(st_centroid(states))) # this errored but I am not sure if I actually need it... 

library("tools")
states$ID <- toTitleCase(states$ID)
head(states)
#'
#' Make base plot
#'

map.baseline <- ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  geom_sf(data = states, fill = NA) + 
  #annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-99, -80), ylim = c(17.5, 31)) +
  xlab("") +
  ylab("")

map.baseline
#st_crs(world)



###* ------------------------------ 
###* 
###* MERGE BASE PLOT WITH PREDATOR DISTRIBUTION AND PREY BIOMASS COMPOSITION
###* 

map.baseline + 
  #ggplot() + 
  geom_polygon(data = plotdata_preddist,
               aes(long, lat, group=group,
                   #fill=as.factor(botz)), # this needs to be a as.factor() to force ggplot to make discrete colors
                   fill=(biomass/area)), 
               #fill = "grey",
               color="black") +
  #scale_fill_brewer(palette="PuBu",direction = -1,na.value="grey") +
  scale_fill_continuous(low = "white", high = "darkgreen",na.value="lightgrey") +
  labs(fill = paste0(index.pred.fid," mt/area")) +
  #theme_void()  # when plotting over the map.baseline, this needs to be dropped !! 
  #
  new_scale("fill") + # This is to split the filling so we have two legends - VERY IMPORTANT HERE!
  #
  geom_scatterpie(aes(x=long, 
                      y=lat, 
                      group=group), 
                  data=plotdata_preycomp[which((plotdata_preycomp$biomass)>1),],
                  #data=plotdata_preycomp,
                  cols=grps.prey.sub,
                  sorted_by_radius = TRUE) + 
  #coord_equal() + # when plotting over the map.baseline, you do not need this! 
  scale_fill_manual(values = getPalette(18)) +
  # misc
  #theme(legend.position = c(0.85, 0.23)) + 
  labs(title = paste0("Predator ",index.pred.fid," spatial distribution against prey composition per polygon"), 
       fill = paste0(index.pred.fid," Prey")) 



#'
#'
#' Lets say you want to plot juv predator biomass against juv predator prey items
#'
#'

# get prey for juv predator
pPREY.subj <- pPREY[str_detect(rownames(pPREY),paste0(index.pred.fid,"1")),]
pPREY.subj <- pPREY.subj[,which(colSums(pPREY.subj)>0)] # isolate prey groups 
grps.preyj <- names(pPREY.subj)
# re-select biomass data
grps.preyj  <- grps.preyj[!((grps.preyj %in% c('DLsed', 'DRsed', 'DCsed')))] # drop the sed ones
biomass.space.preycomp.sub <- biomass.space.preycomp[,c("box_id",grps.preyj)]
# re-merge spatial data to PREY data - FOR geom_scatterpie PART OF PLOT
bgm.boxSpatial.sub.merge <- merge(bgm.boxSpatial.sub, 
                                  biomass.space.preycomp.sub, 
                                  by="box_id")
head(bgm.boxSpatial.sub.merge)
bgm.boxSpatial.sub.merge.tidy <- tidy(bgm.boxSpatial.sub.merge)
plotdata <- left_join(bgm.boxSpatial.sub.merge.tidy,
                      bgm.boxSpatial.sub.merge@data,
                      by="id");head(plotdata)
plotdata <- plotdata[!duplicated(plotdata$id), ] 
plotdata <- left_join(plotdata[,-c(1,2)],
                      centroids,
                      by="id");head(plotdata)
biomass.space.preddist$box_id <- as.numeric(as.character(biomass.space.preddist$box_id))
plotdata_preycomp_sub <- left_join(plotdata,
                                   biomass.space.preddist,
                                   by="box_id");head(plotdata_preycomp)
remove(bgm.boxSpatial.sub.merge,bgm.boxSpatial.sub.merge.tidy,plotdata)
# re-plot
map.baseline + 
  #ggplot() + 
  geom_polygon(data = plotdata_preddist,
               aes(long, lat, group=group,
                   fill=(biomass_juv/area)), # ploting juv predator biomass now
               color="black") +
  scale_fill_continuous(low = "white", high = "darkgreen",na.value="lightgrey") +
  labs(fill = paste0(index.pred.fid," mt/area")) +
  new_scale("fill") + 
  geom_scatterpie(aes(x=long, 
                      y=lat, 
                      group=group), 
                  data=plotdata_preycomp_sub[which((plotdata_preycomp_sub$biomass_juv)>1 ),], # plotting plotdata_preycomp_sub 
                  cols=grps.prey.sub,
                  sorted_by_radius = TRUE) + 
  scale_fill_manual(values = getPalette(18)) +
  labs(title = paste0("Predator ",index.pred.fid," (juv) spatial distribution against prey composition per polygon"), 
       fill = paste0(index.pred.fid," Prey")) 























