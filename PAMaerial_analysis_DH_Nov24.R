#PAM-aerial analysis (NB: ALL CODE RUNS FROM START TO END, though various data input files are required)
#Written by Danielle Harris and Magda Chudzinska
#This version edited on 4th Nov 2024
#The aim of this code is to:
#Step 1: Read in required R libraries (for all steps)
#Step 2: Load the spatial data and format as needed
#Step 3: Spatial modelling of DAS data (using MRSea R package)
#Step 4: Extract densities from the spatial model for each PAM station
#Step 5: Read in PAM data and format as needed
#Step 6: Bayesian data integration
#Step 7: Estimate densities and uncertainties from PAM data
#Step 8: Refit spatial model to PAM-derived densities

###############################################################################
#STEP 1: Read in required R libraries (for all steps)
##############################################################################
require(fields)
require(dplyr)
require(ggplot2)
require(sf)
library(tidyverse)
library(gridExtra)
library(MRSea)
library(mrds)
library(dsm)
library(splancs)
library(patchwork)
library(tweedie)
library(statmod)
require(mgcv)
require(raster)
library(runjags)
library(nimble)
library(MASS)
library(SHELF)
library(ggplot2)
library(coda)
library(MuMIn)
library(colorspace)

###############################################################################
#STEP 2: Load the spatial data and format as needed
################################################################################

#read in the digital survey data, digital aerial survey coordinates are in UTM30N
dd <- readRDS("C:/.../Digital_Survey_Data.rds")
obsdatad<-dd$obsdata
segdatad<-dd$segdata
preddatad<-dd$preddata

#read in the PAM positions, original data are in lat long
PAM <- read.table("C:/.../PAM_locations.txt", sep=",", header=T)

#convert the PAM data to a spatial object and project to UTMs
PAM_SF <- PAM %>% 
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = "+proj=longlat +datum=WGS84")
PAM_UTM <- st_transform(PAM_SF, crs = "+proj=utm +zone=30 +ellps=GRS80 +datum=NAD83")

#read in the shape files of the two survey squares
control <- st_read("C:/.../control block UTM30.shp")
impact <- st_read("C:/.../impact block UTM30.shp")

#make a plot of the survey squares and the locations of the PAM instruments
ggplot() + 
  geom_sf(data=control, lty=2) +
  geom_sf(data=impact, lty=2)+
  geom_point(data=segdatad, aes(x=x,y=y))+
  geom_point(data=obsdatad, aes(x=x,y=y), col="red")+
  geom_sf(data=PAM_UTM, col="orange")

################################################################################
#STEP 3: Spatial modelling of DAS data (using MRSea R package)
################################################################################

# at the time of analysis, no information was available on when the observations were taken, so the observations
# could not be assigned to a specific survey day. Therefore, all data were combined and the
# following assumptions were made:
# 1. Data in obsdatad are ordered in time - this is necessary to estimate autocorrelation
# 2. It was assumed that 'Effort' in obsdatad and segdatad is a cumulative effort per grid combining all
# four days of survey

# The final spatial model is given below. 
# All the model diagnostics and model selection was performed in a separate analysis

# To set up the model, all DAS observations are combined from the same grid
obsdatad <- obsdatad %>% mutate(grid_id_obs = group_indices(., x, y))

obsdatadSum <- obsdatad %>%
  group_by(grid_id_obs) %>% 
  mutate(response = sum(size)) %>%
  distinct(grid_id_obs, .keep_all=TRUE)

# Now merge with segdatad to also have zeros (it was assumed that zeros
# are true zeros)
# Data are merged by Sample.Label
# Not all columns are required in the merged document so a temp document is made first
# seg_obs is then the main input to modelling
temp <- segdatad %>% left_join(obsdatadSum, by="Sample.Label")

seg_obs <- segdatad
seg_obs$count <- temp$response
seg_obs$count[is.na(seg_obs$count)] <- 0 
seg_obs$AreaOffset <- seg_obs$Effort * seg_obs$Strip_Width/1000 # here it is assumed 
#that strip width refers to the entire width, not half strip width, and unit is km2

# SALSA requires the following column names to be included: response, x.pos, y.pos
seg_obs <- seg_obs %>% 
  rename("response" = "count") %>%
  rename("x.pos" = "x") %>%
  rename("y.pos" = "y")

# set up the prediction grid. This is the same as the surveyed area but it has to be a separate file
predgrid <- seg_obs[,c("x.pos","y.pos","depth","sediment","slope","AreaOffset","response")]
predgrid$AreaOffset <- 16 # each grid is 16 km2
gridSize <- 16

# fitting the model
profoutFull <- tweedie.profile(response ~ x.pos + y.pos + depth + slope + offset(log(AreaOffset)), 
                               data=seg_obs,
                               xi.vec = seq(1.01, 1.99, by=0.05), do.plot=FALSE)

initialModelTw <- gamMRSea(response ~ 1 + offset(log(AreaOffset)), family = tweedie(var.power=profoutFull$xi.max, link.power = 0), 
                           data = seg_obs)

# 1d model i.e., no lat and lon fitted but needed for 2D model fitting
salsa1dlistTwNS<-list(fitnessMeasure = 'AICtweedie', 
                      minKnots_1d = c(1,1), 
                      maxKnots_1d = c(3,4), 
                      startKnots_1d = c(1,1), 
                      degree = c(NA,NA),
                      gaps = c(0,0),
                      splines = c("ns", "ns"))

salsa1dRunTwNS<-runSALSA1D(initialModelTw, 
                           salsa1dlistTwNS, 
                           varlist = c("depth","slope"), 
                           datain = seg_obs,
                           removal = TRUE,
                           predictionData = predgrid, 
                           suppress.printout = TRUE)

final1dTw <- salsa1dRunTwNS$bestModel

# 2d model, lat and lon included
knotgrid<- getKnotgrid(coordData = cbind(seg_obs$x.pos, seg_obs$y.pos), 
                       numKnots = 300,
                       plot = FALSE)

distMats <- makeDists(cbind(seg_obs$x.pos, seg_obs$y.pos), knotgrid)


salsa2dlist<-list(fitnessMeasure = 'AICtweedie',
                  knotgrid = knotgrid,
                  startKnots=5,
                  minKnots=2,
                  maxKnots=15,
                  gap=0)

salsa2dRunTw<-runSALSA2D(model = final1dTw, 
                         salsa2dlist = salsa2dlist, 
                         d2k=distMats$dataDist,
                         k2k=distMats$knotDist)

final2dTw <- salsa2dRunTw$bestModel

#predictions for the final model
#predictions with 'Km' show density per km2
#predictions with 'Grid' show density per grid so porpoises/16 km2
predgrid$final2dTwPredKm <- (predict(object = final2dTw, 
                                     newdata = predgrid)[,1])/gridSize
predgrid$final2dTwPredGrid <- predict(object = final2dTw, 
                                      newdata = predgrid)

# calculating CVs, CIs
# these are CVs and CIs related to model uncertainty only
preddist<-makeDists(cbind(predgrid$x.pos, predgrid$y.pos),
                    knotgrid, knotmat=FALSE)$dataDist

boots<-do.bootstrap.cress.robust(final2dTw, predgrid, g2k=preddist, B=300)
cis<-makeBootCIs(boots)

predgrid$final2dTw_lower95_grid<-cis[,1]
predgrid$final2dTw_upper95_grid<-cis[,2]
predgrid$final2dTw_lower95_km<-cis[,1]/gridSize
predgrid$final2dTw_upper95_km<-cis[,2]/gridSize

predgrid$final2dTw_sd <- apply(boots, 1, sd)
predgrid$final2dTw_cv <- predgrid$final2dTw_sd/predgrid$final2dTwPredGrid

# plot some figures to look at uncertainty

g1 <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dTwPredKm)) +
  scale_fill_distiller(palette = "BrBG", breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4), limits=c(0,4)) +
  xlab("x") + ylab("y") + theme_bw()+
  guides(fill=guide_legend(title="Density [km2]"))+
  ggtitle("Mean pred")

g2 <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dTw_cv)) +
  scale_fill_distiller(palette = "BrBG") +
  xlab("x") + ylab("y") + theme_bw()+
  ggtitle("CV")+
  guides(fill=guide_legend(title="CV"))

g3 <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dTw_lower95_km)) +
  scale_fill_distiller(palette = "BrBG", breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4), limits=c(0,4)) +
  xlab("x") + ylab("y") + theme_bw()+
  ggtitle("lower CI")+
  guides(fill=guide_legend(title="Density [km2]"))

g4 <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dTw_upper95_km)) +
  scale_fill_distiller(palette = "BrBG", breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4), limits=c(0,4)) +
  xlab("x") + ylab("y") + theme_bw()+
  ggtitle("upper CI")+
  guides(fill=guide_legend(title="Density [km2]"))

g1 + g3 + g4 + plot_layout(nrow = 2,ncol=2, byrow = TRUE ,guides = 'collect')
g1 + g3 + g4 + g2 + plot_layout(nrow = 2,ncol=2, byrow = TRUE ,guides = 'collect')

################################################################################
#STEP 4: Extract densities from the spatial model for each PAM station
################################################################################

# extracting from mean predictions (densities per km2)

# first converting model predictions into raster
modelr <- cbind.data.frame(predgrid$x.pos, predgrid$y.pos,predgrid$final2dTwPredKm)
modelr <- rasterFromXYZ(modelr)
crs(modelr) <- crs(PAM_UTM)

# making rasters for CVs
modelrCVmod <- cbind.data.frame(predgrid$x.pos, predgrid$y.pos,predgrid$final2dTw_cv)
modelrCVmod <- rasterFromXYZ(modelrCVmod)
crs(modelrCVmod) <- crs(PAM_UTM)

# making rasters for CIs
modelrCIupper <- cbind.data.frame(predgrid$x.pos, predgrid$y.pos, predgrid$final2dTw_upper95_km)
modelrCIupper <- rasterFromXYZ(modelrCIupper)
crs(modelrCIupper) <- crs(PAM_UTM)

modelrCIlower <- cbind.data.frame(predgrid$x.pos, predgrid$y.pos, predgrid$final2dTw_lower95_km)
modelrCIlower <- rasterFromXYZ(modelrCIlower)
crs(modelrCIlower) <- crs(PAM_UTM)

#extract the information at the location of the PAM instruments
PAM_dens <- raster::extract(modelr,PAM_UTM, sp=T) # the NAs are PAM locations outside the spatial model
PAM_CVmod <- raster::extract(modelrCVmod,PAM_UTM, sp=T) # the NAs are PAM locations outside the spatial model
PAM_CIupper<-raster::extract(modelrCIupper,PAM_UTM, sp=T)
PAM_CIlower<-raster::extract(modelrCIlower,PAM_UTM, sp=T)

#create required data frames for further analysis
PAM_CIupper_df<-as.data.frame(raster::extract(modelrCIupper,PAM_UTM,sp=T))
PAM_CIlower_df<-as.data.frame(raster::extract(modelrCIlower,PAM_UTM,sp=T))
PAM_densDF <- as.data.frame(PAM_dens)
colnames(PAM_densDF) <- c("PAM_id","Block","Density_km2","x","y")
PAM_CVmodDF <- as.data.frame(PAM_CVmod)
colnames(PAM_CVmodDF) <- c("PAM_id","Block","CVmod","x","y")
PAM_densDF$CVmod <- PAM_CVmodDF$CVmod

#plot to check underlying density and densities assigned to the PAM locations.
ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=predgrid$final2dTwPredKm)) +
  scale_fill_distiller(palette = "BrBG") +
  xlab("x") + ylab("y") + theme_bw() +
  guides(fill=guide_legend(title="Density [km2]"))+
  geom_point(data=PAM_densDF,aes(x=x,y=y, size = Density_km2), col="black")

################################################################################
#STEP 5: Read in PAM data and format as needed
#This step could be simplified by providing a file with the number of PPS per CPOD per day across the required time period
################################################################################

#use the CPOD_PPM_all.df object - this has all the relevant data.
CPOD_PPM_all.df<-read.csv("C:/.../File4_CPOD_PPM_01Aug10_01Oct10.df.csv")

#Extract all days from the CPOD data
limit1<-as.Date("2010-08-01")
limit2<-as.Date("2010-10-01")
date.index<-which(unique(CPOD_PPM_all.df$CPOD_Date)>=limit1 & unique(CPOD_PPM_all.df$CPOD_Date)<=limit2)
#now have the unique dates for between the limits as defined above
#use this index to pull out the CPOD counts for each day for each CPOD
#for each CPOD ID, find the dates and record the summed PPS for that day. If there are multiple records, print to check

#create a results matrix with dimensions #CPOD x #days
PAM_densDF_edit<-PAM_densDF
PPS_alldates<-matrix(NA,dim(PAM_densDF_edit)[1],length(date.index))

#for each CPOD ID
for (i in 1:dim(PAM_densDF_edit)[1]){
  #for each date
  for (j in 1:length(date.index)){
    temp.index<-which(CPOD_PPM_all.df$CPOD_locID==PAM_densDF_edit$PAM_id[i] & CPOD_PPM_all.df$CPOD_Date==unique(CPOD_PPM_all.df$CPOD_Date)[date.index][j],arr.ind=TRUE)
    temp.index.length<-length(temp.index)
    if (temp.index.length==1){
      PPS_alldates[i,j]<-CPOD_PPM_all.df$CPOD_summedPPS_edit[temp.index]} 
    if (temp.index.length>1){
      PPS_alldates[i,j]<-sum(CPOD_PPM_all.df$CPOD_summedPPS_edit[temp.index])
      print(temp.index)
      print(CPOD_PPM_all.df$CPOD_summedPPS_edit[temp.index])
      print(CPOD_PPM_all.df$CPOD_MinsOn[temp.index])} 
  }} 

#remove NAs caused by CPODs with no associated DAS density (NAs still remain because of no CPOD data on a given date)
PPS_alldates_noNA<-PPS_alldates[!is.na(PAM_densDF_edit$Density_km2),]

#there are some multiple records here - just manually check that the second records are zeros, and for now can be ignored. The second CPODS either have zero effort, or only on for a few minutes.

#now have a matrix with number of PPS per CPOD per day across the entire time period.
#this can be used for estimating densities on a larger timescale

###################
#sum the counts across the 4 days for the Bayesian model

#first remove the NAs from the case study data
PAM_densDF_edit_noNA<-PAM_densDF_edit[!is.na(PAM_densDF_edit$Density_km2),]

#find the date entries using the date index and take data from the PPS_alldates_noNA
PPS_280810_index<-which(date.index==which(unique(CPOD_PPM_all.df$CPOD_Date)=="2010-08-28"))
PPS_190910_index<-which(date.index==which(unique(CPOD_PPM_all.df$CPOD_Date)=="2010-09-19"))
PPS_260910_index<-which(date.index==which(unique(CPOD_PPM_all.df$CPOD_Date)=="2010-09-26"))
PPS_270910_index<-which(date.index==which(unique(CPOD_PPM_all.df$CPOD_Date)=="2010-09-27"))
PPS_010810_index<-which(date.index==which(unique(CPOD_PPM_all.df$CPOD_Date)=="2010-08-01"))
PPS_010910_index<-which(date.index==which(unique(CPOD_PPM_all.df$CPOD_Date)=="2010-09-01"))
PPS_011010_index<-which(date.index==which(unique(CPOD_PPM_all.df$CPOD_Date)=="2010-10-01"))

PAM_densDF_edit_noNA$PPS_280810<-PPS_alldates_noNA[,PPS_280810_index]
PAM_densDF_edit_noNA$PPS_190910<-PPS_alldates_noNA[,PPS_190910_index]
PAM_densDF_edit_noNA$PPS_260910<-PPS_alldates_noNA[,PPS_260910_index]
PAM_densDF_edit_noNA$PPS_270910<-PPS_alldates_noNA[,PPS_270910_index]
PAM_densDF_edit_noNA$PPS_010810<-PPS_alldates_noNA[,PPS_010810_index]
PAM_densDF_edit_noNA$PPS_010910<-PPS_alldates_noNA[,PPS_010910_index]
PAM_densDF_edit_noNA$PPS_011010<-PPS_alldates_noNA[,PPS_011010_index]

#now sum PPS for the 4 days of the DAS survey
#some rows will have NAs due to having no CPOD data available (reason unknown). Need to remove these rows
PAM_densDF_edit_noNA$PPS_DAS_summed<-rowSums(cbind(PAM_densDF_edit_noNA$PPS_280810,PAM_densDF_edit_noNA$PPS_190910,PAM_densDF_edit_noNA$PPS_260910,PAM_densDF_edit_noNA$PPS_270910))

#use the PPS_DAS_summed column to remove the NAs - these had all 4 DAS days as NAs so indicate the columns where all data were NAs. 
#This will not work generally so will need to be generalised for further use. These results were manually checked.
PAM_densDF_edit_noNAs<-PAM_densDF_edit_noNA[!is.na(PAM_densDF_edit_noNA$PPS_DAS_summed),]

#produce colour accessible plot for the report
ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=predgrid$final2dTwPredKm)) +
  scale_fill_continuous_sequential(palette = "Teal", limits=c(0,max(predgrid$final2dTwPredKm))) +
  xlab("x") + ylab("y") + theme_bw() +
  guides(fill=guide_legend(title="Density [km2]"))+
  geom_point(data=PAM_densDF_edit_noNAs,aes(x=x,y=y), col="black")
################################################################################
#STEP 6: Bayesian data integration
################################################################################

#set number of observations (length of vector containing the porpoise positive second data)
nobs <- length(as.vector(t(PAM_densDF_edit_noNAs$PPS_DAS_summed)))

#Enter the number of PPS at each PAM instrument location as a vector
PPS_DAS_summed<-as.vector(t(PAM_densDF_edit_noNAs$PPS_DAS_summed))

#Set the time monitored, specifically between dawn and dusk (for 4 days in our case study)
t <- 12*60*60*4

#Define the densities from the DAS surveys
Dhat<-PAM_densDF_edit_noNAs$Density_km2

#Set up the Nimble model to estimate g0 and vp from the DAS + CPOD data

#Need this additional analysis to produce assumed prior g(0) parameters
#################################################################################
# From Teilmann et al. 2013, we expect the % time available to be between
# 42.5% and 61.5% with a median of 50.8% (these are monthly predictions). 
# We set the upper and lower values to be the 0.05 and 0.95 quantiles. 
# The lowest and highest % time per *individual* ranged between 41.3% and 73.3%, 
# so these were set as the 0.01 and 0.99 quantiles.
# Then fit various ditributions to these data, including Normal and Beta

fg0 <- fitdist(vals = c(0.413, .425, .508, .615, 0.733), 
               probs = c(0.01, 0.05, .50, .95, 0.99),
               lower = 0, upper = 1)

# check that this looks reasonable
hist(rbeta(1000, fg0$Beta$shape1, fg0$Beta$shape2), xlim = c(0,1))
################################################################################

code <- nimbleCode({
  
  # Priors for parameters to be estimated by the model
  g0 ~ dbeta(47.61495, 45.94263) # from a separate analysis to define a beta distribution for g0 
  #(see immediately below for details)
  vp ~ dunif(0, 0.03) # upper bound selected based on Jacobson et al. (2017)
  
  # Model
  for (i in 1:nobs){
    
    # assume density estimates are lognormally distributed
    D[i] ~ dlnorm(meanlog = log(Dhat[i]), sd = D.lsigma[i]) 
    D.lsigma[i] <- sqrt(log(1+(CV[i]^2))) # convert var to log scale
    
    # assume PPS are poisson
    PPS[i] ~ dpois(mu[i]) 
    mu[i] <- (D[i]/g0) * vp * time # linear predictor
    
  } # end for i
  
}) # end nimbleCode

#list of data for Nimble model: DAS densities (with CV) and PPS at each instrument location 
nimbleData <- list(Dhat = Dhat, 
                   PPS = PPS_DAS_summed, 
                   CV = PAM_densDF_edit_noNA$CVmod)

#time and number of observations listed as constants
nimbleConstants <- list(nobs = nobs, time = t)

#list of initial starting values
nimbleInits <- list(D = rep(1, nobs))

#list of parameters to estimate
nimbleParams <- list("D", "vp", "g0")

#setting up the model
model <- nimbleModel(code = code,
                     constants = nimbleConstants,
                     data = nimbleData,
                     inits = nimbleInits)

thin = 10
niter = 250000
nburnin = 200000
nchains = 4

#running the model
nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        constants = nimbleConstants, data = nimbleData,
                        thin = thin, niter = niter, nburnin = nburnin, 
                        nchains = nchains, samplesAsCodaMCMC = TRUE)

single.mcmc <- combine.mcmc(nimbleOut)

#Results

#check convergence (values much above 1 indicate lack of convergence, taken from ?gelman.diag)
gelman.diag(nimbleOut)

# check samples of g0 via plotting
plot(single.mcmc[,nobs+1])
hist(single.mcmc[,nobs+1], xlim = c(0, 1), main="", xlab="Estimate of g(0)", ylim=c(0,3500))

# check samples of vp via plotting
plot(single.mcmc[,nobs+2])
hist(single.mcmc[,nobs+2], xlim = c(0, 0.003), main="", xlab="Estimate of vp", ylim=c(0,3000))

# check estimates of D
qdf <- summary(single.mcmc)$quantiles
plot(Dhat, qdf[1:nobs, 3], xlab="Estimated relative density from DAS (animals.km2)", ylab="Updated relative DAS estimates (animals.km2)") # obs v pred
abline(c(0,1))

# check how the point estimate of vp translates to a relative density given PPS (this is a quick check that the values look reasonable compared to the initial DAS relative densities)
vp_Est <- qdf[nobs+2,3]
g0_Est <- qdf[nobs+1,3]
D_PAMEst <- (PPS_DAS_summed * g0_Est)/(vp_Est * t)
plot(Dhat, D_PAMEst, xlab="Estimated relative density from DAS (animals.km2)", ylab="Estimated relative density from PPS data using vp and g(0) (animals.km2)")
abline(c(0,1))

################################################################################
#STEP 7: Estimate densities and uncertainties from PAM data
################################################################################

#Several absolute density estimates are needed for comparison
#1) A design-based estimate from the DAS surveys
#2) A model-based estimate from the DAS surveys 
#3) A design-based estimate from the CPOD data, using the DAS-derived vp parameters

########################################################################
#1) Design based DAS survey estimate
#Using the transect lines for the encounter rate variance
#Also need an estimate of g(0)
#Use the information in seg_obs for transect and effort data
#Assume that transect labels are unique transects
#Assume that effort, and therefore the offset, is summed over 4 days 

dim(seg_obs)
n_DASabs<-sum(seg_obs$response)
effort_DASabs<-sum(seg_obs$AreaOffset)
g0_DASabs<-as.numeric(fg0$Normal[1]) #this is the same distribution used in Step 6
D_DASabs<-(n_DASabs/effort_DASabs)/g0_DASabs

#variance
#encounter rate variance - use the transects as samplers
seg_obs_Sum<-seg_obs %>%
  group_by(Transect.Label) %>%
  mutate(response_sum = sum(response)) %>%
  mutate(AreaOffset_sum = sum(AreaOffset)) %>%
  mutate(ERpertransect = sum(response)/sum(AreaOffset)) %>%
  distinct(Transect.Label, .keep_all=TRUE)

#variance of ER
varER_DAS<-var(seg_obs_Sum$response_sum)
CV_ER_DAS<-sqrt(varER_DAS/(n_DASabs^2))

#variance of g(0)
#use the sd from the normal distribution
varg0_DAS<-as.numeric(fg0$Normal[2]^2)
CV_g0_DAS<-sqrt(varg0_DAS/(g0_DASabs^2))

#now combine the CVs and estimate var of the density
varD_DAS<-(D_DASabs^2)*((CV_ER_DAS^2)+(CV_g0_DAS^2))

#now produce the CIs, assumed to be log-normal as in Buckland et al. (2001)
za<-1.96
varlnD<-log(1+(varD_DAS/(D_DASabs^2)))
C<-exp(za*sqrt(varlnD))
LCL_D_DAS<-D_DASabs/C
UCL_D_DAS<-D_DASabs*C

#################################################################
#2) A model-based estimate from the DAS surveys
#Two estimates will be derived (1) using the mean density across the whole prediction grid 
#and (2) using the mean density based on the CPOD locations only.

#mean density estimate for the whole prediction area, corrected for g(0)
D.vis.model<-mean(predgrid$final2dTwPredKm/g0_DASabs)

#mean density from the predictions at the CPOD locations, corrected for g(0)
D.vis.model.CPOD<-mean(PAM_densDF_edit_noNAs$Density_km2/g0_DASabs)

#variance

#uncertainty of g(0) needs to be built into the CIs.
#to achieve this, draw values of g(0) from the assumed distribution and combine these with the bootstrap results from the confidence interval values

#first, create 500 draws from the g(0) normal distribution
g0_500draws<-rnorm(500, g0_DASabs, sqrt(varg0_DAS)) 

#second, now adjust the boots matrix (no. of prediction grid cells x #bootstrap iterations for model variance) with each value of g(0)
for (i in 1:length(g0_500draws)){
  temp_boots<-(boots/gridSize)/g0_500draws[i]
  if (i==1) CIs_bootstrap_predgrid<-temp_boots
  if (i>1){CIs_bootstrap_predgrid<-cbind(CIs_bootstrap_predgrid,temp_boots)} 
  print(i)
}

#now have 300*500 values for each prediction grid. This combines both the model uncertainty and the g0 uncertainty
#now need to extract the 2.5 and 97.5 percentile for each grid cell

#for each surface, calculate the mean density
mean_bootstrap_predgrid_persurface<-colMeans(CIs_bootstrap_predgrid)

#now take the LCL and UCL from this
LCL_bootstrap_predgrid_surface<-quantile(mean_bootstrap_predgrid_persurface, probs=0.025)
UCL_bootstrap_predgrid_surface<-quantile(mean_bootstrap_predgrid_persurface, probs=0.975)
D.UCL.vis.model<-UCL_bootstrap_predgrid_surface
D.LCL.vis.model<-LCL_bootstrap_predgrid_surface

#need to now extract CPOD locations from this larger grid
CIs_bootstrap_predgrid_df <- cbind.data.frame(predgrid$x.pos, predgrid$y.pos, CIs_bootstrap_predgrid)
CIs_bootstrap_predgrid_ras <- rasterFromXYZ(CIs_bootstrap_predgrid_df)
crs(CIs_bootstrap_predgrid_ras) <- crs(PAM_UTM)

#extract the information at the location of the PAM instruments
PAM_CIs_bootstrap_predgrid_ras<- raster::extract(CIs_bootstrap_predgrid_ras,PAM_UTM, sp=T) # the NAs are PAM locations outside the spatial model
PAM_CIs_bootstrap_predgrid_df<-as.data.frame(PAM_CIs_bootstrap_predgrid_ras)
PAM_CIs_bootstrap_predgrid_df_noNA<-PAM_CIs_bootstrap_predgrid_df[!is.na(PAM_densDF_edit$Density_km2),]

#take only the surfaces
PAM_CIs_bootstrap_predgrid_df_noNA_edit<-PAM_CIs_bootstrap_predgrid_df_noNA[,3:150002]

#now take the means of all surfaces
mean_PAM_CIs_bootstrap_predgrid_df<-colMeans(PAM_CIs_bootstrap_predgrid_df_noNA_edit)

#now take the LCL and UCL from this
mean_PAM_CIs_bootstrap_predgrid_df_LCL<-quantile(mean_PAM_CIs_bootstrap_predgrid_df, probs=0.025)
mean_PAM_CIs_bootstrap_predgrid_df_UCL<-quantile(mean_PAM_CIs_bootstrap_predgrid_df, probs=0.975)
D.UCL.vis.model.CPOD<-mean_PAM_CIs_bootstrap_predgrid_df_UCL
D.LCL.vis.model.CPOD<-mean_PAM_CIs_bootstrap_predgrid_df_LCL

#################################################################
#3) Design-based estimate from the CPOD data using the estimate of vp 
#Need to estimate density using each estimate of vp to maintain the uncertainty in estimation.
#Then take the credible intervals from these.
#I have 20000 estimates of vp to use in the density estimator

#Relative density using the CPOD data
g0_quantiles<-qdf[nobs+1,]
vp_quantiles<-qdf[nobs+2,]

g0_allestimates<-single.mcmc[,nobs+1]
vp_allestimates<-single.mcmc[,nobs+2]

#create results vectors with length #days
Dcpod_rel_alldates<-c(NA,length(date.index))
Dcpod_abs_alldates<-c(NA,length(date.index))
sumPPS_cpod_alldates<-c(NA,length(date.index))
varPPS_cpod_alldates<-c(NA,length(date.index))

#for each CPOD ID
for (i in 1:dim(PPS_alldates_noNA)[1]){
  #for each date
  for (j in 1:length(date.index)){
    temp.Dcpod_rel.vector<-(sum(PPS_alldates_noNA[,j], na.rm=TRUE)/((12*60*60)*vp_allestimates*nobs))*g0_allestimates
    Dcpod_rel_alldates[j]<-summary(temp.Dcpod_rel.vector)$quantiles[3]
    temp.Dcpod_abs.vector<-(sum(PPS_alldates_noNA[,j], na.rm=TRUE)/((12*60*60)*vp_allestimates*nobs))
    Dcpod_abs_alldates[j]<-summary(temp.Dcpod_abs.vector)$quantiles[3]
    sumPPS_cpod_alldates[j]<-sum(PPS_alldates_noNA[,j], na.rm=TRUE)
    varPPS_cpod_alldates[j]<-var(PPS_alldates_noNA[,j], na.rm=TRUE)
    #the below is commented out as it gives the same CV as using the counts only, as the effort is the same for all instruments in this case
    #the code is kept for cases where encounter rates would differ from the count data.
    #ERcpod_alldates[j]<-(sum(PPS_alldates_noNA[,j], na.rm=TRUE))/(12*60*60*nobs) 
    #varER_cpod_alldates[j]<-var(PPS_alldates_noNA[,j]/(12*60*60*nobs), na.rm=TRUE)
  }
} 

#variance
#encounter rate variance - use the CPODS as samplers per day

#variance of encounter rates
#NB: using the counts only gives the same CV as calculating the encounter rate, as the effort is the same for all in this case
CV_PPS_CPOD_alldates<-sqrt(varPPS_cpod_alldates/(sumPPS_cpod_alldates^2))
#CV_ER_CPOD_alldates<-sqrt(varER_cpod_alldates/(ERcpod_alldates^2))

#variance of vp
#take the variance of the 20000 estimates
var_vp_CPOD_alldates<-var(single.mcmc[,nobs+2])
mean_vp_CPOD_alldates<-mean(single.mcmc[,nobs+2])

#now estimate the CV of vp
CV_vp_CPOD_alldates<-sqrt(var_vp_CPOD_alldates/(mean_vp_CPOD_alldates^2))

#now combine the CVs and estimate variance of the density
varD_CPOD_alldates<-(Dcpod_abs_alldates^2)*((CV_PPS_CPOD_alldates^2)+(CV_vp_CPOD_alldates^2))

#now produce the CIs, assuming lognormal confidence intervals
za<-1.96
varlnD_CPOD_alldates<-log(1+(varD_CPOD_alldates/(Dcpod_abs_alldates^2)))
C<-exp(za*sqrt(varlnD_CPOD_alldates))
LCL_D_CPOD_alldates<-Dcpod_abs_alldates/C
UCL_D_CPOD_alldates<-Dcpod_abs_alldates*C

#create a plot for the whole time series
Dates.num.all<-c(1:dim(PPS_alldates_noNA)[2])
Dates.all<-seq(limit1, limit2, by="days")
data.range<-range(LCL_D_CPOD_alldates,UCL_D_CPOD_alldates,LCL_D_DAS,UCL_D_DAS,D.UCL.vis.model,D.LCL.vis.model,D.UCL.vis.model.CPOD,D.LCL.vis.model.CPOD)
plot(c(1:(length(Dates.num.all)+3)),c(D_DASabs,D.vis.model.CPOD,D.vis.model,Dcpod_abs_alldates),xlab="Date", ylab="Density (animals/km2)", xaxt = 'n',pch=16, ylim=c(0,ceiling(data.range)[2]))
points(c(1:3),c(D_DASabs,D.vis.model.CPOD,D.vis.model),col="blue", pch=16)
Dates.subset<-c("01 Aug", "15 Aug","01 Sept", "15 Sept", "01 Oct")
Dates.subset.num<-c(4,18,35,49,65)
axis(1, at=Dates.subset.num, labels=Dates.subset, las=1, cex.axis=1)
#axis(1, at=c(1:(length(Dates.num.all)+3)), labels=c(as.character(Dates.all),"DAS_des","DAS_mod1","DAS_mod2"), las=2, cex.axis=0.5)
## draw arrows up and down of data points to show the error bars.
arrows(c(1:(length(Dates.num.all)+3)), c(UCL_D_DAS,D.UCL.vis.model,D.UCL.vis.model.CPOD,UCL_D_CPOD_alldates), c(1:(length(Dates.num.all)+3)), c(LCL_D_DAS,D.LCL.vis.model,D.LCL.vis.model.CPOD,LCL_D_CPOD_alldates), angle=90, code=3, length=0.06)
text(1,1.5,"DAS", col="blue")
text(35,1.5,"PAM with DAS-derived parameters")
abline(v=3.5,lty="dashed")

################################################################################
#STEP 8: Refit spatial model to PAM-derived densities
################################################################################

# Below are two steps:
# 8a. A simple rasterising of the densities estimated using the PAM data and estimate of vp. No model is fitted here but a density map is produced.
#    This approach can also be used to produce density maps for periods outside the DAS data collection.
#    The same prediction grid is used as for MRSea. If more than one PAM instrument was deployed in the same grid cell, a mean
#    density is used. This approach will leave grids where there is no PAM instrument as NAs. 
# 8b. The same MRSea model as fitted above is fitted to the PAM estimates.
#    A separate model selection should be done for each predicted month if any seasonal covariates were included.
#    But this is not the case here. 
#    To run the code in step 8b, step 8a has to be run first.

#### both these steps will be run for two dates: 
# one day when a DAS also took place (19 Sept 2010)
# one day outside DAS data collection (1 Oct 2010)

################################################
# Step 8a - rasterising PAM estimates

#create same objects for 19/09 and 1/10 and add to PAM_densDF_edit_noNAs
PAM_densDF_edit_noNAs$PPS_190910_dens<-PAM_densDF_edit_noNAs$PPS_190910/((12*60*60)*mean_vp_CPOD_alldates)
PAM_densDF_edit_noNAs$PPS_011010_dens<-PAM_densDF_edit_noNAs$PPS_011010/((12*60*60)*mean_vp_CPOD_alldates)

#str(PAM_densDF_edit_noNAs)
#PAM_densDF_edit_noNA$PPS_DAS_summed_dens <- PAM_densDF_edit_noNA$PPS_DAS_summed_dens$PPS_DAS_summed

# converting PAMs into spatial points
pamsSP <- PAM_densDF_edit_noNAs
coordinates(pamsSP) <- ~x+y

# creating blank raster
r <- raster(ext=extent(modelr),res=res(modelr))

# rasterizing
PAM_19sept <- rasterize(x=pamsSP, y=r, field=pamsSP$PPS_190910_dens, fun=function(x, ...) mean(x, na.rm=T))
PAM_19sept <- as.data.frame(PAM_19sept, xy=T)

PAM_1oct <- rasterize(x=pamsSP, y=r, field=pamsSP$PPS_011010_dens, fun=function(x, ...) mean(x, na.rm=T))
PAM_1oct <- as.data.frame(PAM_1oct, xy=T)

# plotting
ggplot(PAM_19sept) + geom_tile(aes(x=x, y=y, fill=layer)) +
  scale_fill_distiller(palette = "BrBG", breaks = seq(0,3,0.5), limits=c(0,3),na.value = 'white') +
  xlab("x") + ylab("y") + theme_bw()+
  guides(fill=guide_legend(title="Density [km2]"))+
  ggtitle("19/09/2010")#+

ggplot(PAM_1oct) + geom_tile(aes(x=x, y=y, fill=layer)) +
  scale_fill_distiller(palette = "BrBG", breaks = seq(0,3,0.5), limits=c(0,3),na.value = 'white') +
  xlab("x") + ylab("y") + theme_bw()+
  guides(fill=guide_legend(title="Density [km2]"))+
  ggtitle("01/10/2010")#+

################################################
# Step 8b - fitting the MRSea model to the new density estimates from the PAM data and estimated vp
# To keep the method consistent with fitting MRSea to DAS data, each PAM instrument was assigned to respective grid cells
# If more than one PAM instrument was placed in a grid cell, the mean of the densities was used (as in Step 8a)
# Similar object names have been used as for MRSea modelling for DAS (in Step 3)
# An offset is no longer explicitly needed as estimates are already standardised (in per km2)
# In the modelling below, the densities are transformed into abundances per grid cell. 
# A quasipoisson distribution was used (zeros are much reduced compared to DAS data, so quasipoisson was judged to be sufficient here)

#### 19 September
PAM_BAY_19sept <- PAM_19sept
PAM_BAY_19sept$layer <- as.integer(PAM_BAY_19sept$layer * gridSize)

# using only the grids which are in predgrid file
temp <- as.data.frame(modelr, xy=T)
toKeep <- which(!is.na(temp$predgrid.final2dTwPredKm))
PAM_BAY_19sept <- PAM_BAY_19sept[c(toKeep),]

# SALSA requires the following column names to be included: response, x.pos, y.pos
seg_obs_PAM_19sept <- PAM_BAY_19sept %>% 
  rename("response" = "layer") %>%
  rename("x.pos" = "x") %>%
  rename("y.pos" = "y")

# adding environmental variables to the data
seg_obs_PAM_19sept$depth <- seg_obs$depth
seg_obs_PAM_19sept$sediment <- seg_obs$sediment
seg_obs_PAM_19sept$slope <- seg_obs$slope

# removing grids with NAs i.e., grids where there were no PAM instruments
seg_obs_PAM_19sept <- seg_obs_PAM_19sept[!is.na(seg_obs_PAM_19sept$response),]
hist(seg_obs_PAM_19sept$response)

# fitting the model
initialModel_19sept <- gamMRSea(response ~ 1, family = quasipoisson(link = "log"),
                                data = seg_obs_PAM_19sept)

# 1d model
salsa1dlist<-list(fitnessMeasure = 'QBIC', 
                  minKnots_1d = c(1,1), 
                  maxKnots_1d = c(4,4), 
                  startKnots_1d = c(1,1), 
                  degree = c(NA,NA),
                  gaps = c(0,0),
                  splines = c("ns", "ns"))

salsa1dRun_19sept<-runSALSA1D(initialModel_19sept, 
                              salsa1dlist, 
                              varlist = c("depth","slope"), 
                              datain = seg_obs_PAM_19sept,
                              removal = TRUE,
                              #predictionData = predgrid, 
                              suppress.printout = TRUE)

summary(salsa1dRun_19sept$bestModel) 
salsa1dRun_19sept$keptvarlist

final1d_19sept <- salsa1dRun_19sept$bestModel

# 2d model
knotgrid_19sept<- getKnotgrid(coordData = cbind(seg_obs_PAM_19sept$x.pos, seg_obs_PAM_19sept$y.pos), 
                              numKnots = 300,
                              plot = FALSE)

distMats_19sept <- makeDists(cbind(seg_obs_PAM_19sept$x.pos, seg_obs_PAM_19sept$y.pos), knotgrid_19sept)

salsa2dlist_19sept<-list(fitnessMeasure = 'QBIC',#AICtweedie
                         knotgrid = knotgrid_19sept,
                         startKnots=5,
                         minKnots=2,
                         maxKnots=15,
                         gap=0)

salsa2dRun_19sept<-runSALSA2D(model = final1d_19sept, 
                              salsa2dlist = salsa2dlist_19sept, 
                              d2k=distMats_19sept$dataDist,
                              k2k=distMats_19sept$knotDist)

final2d_19sept <- salsa2dRun_19sept$bestModel
summary(final2d_19sept)

# predictions for the final model
# prediction grid is the same as the surveyed area but it has to be a separate file
# make sure the observation grid from DAS is used.
predgrid <- seg_obs[,c("x.pos","y.pos","depth","sediment","slope","response")]

# predictions with 'Km' show density per km2 (as in the paper)
# predictions with 'Grid' show density per grid so # porpoises/16 km2
preddist_19sept<-makeDists(cbind(predgrid$x.pos, predgrid$y.pos),
                           knotgrid_19sept, knotmat=FALSE)$dataDist

predgrid$final2dPredKm_19sept <- (predict(object = final2d_19sept, 
                                          newdata = predgrid,
                                          g2k=preddist_19sept)[,1])/gridSize
predgrid$final2dPredGrid_19sept <- predict(object = final2d_19sept, 
                                           newdata = predgrid,
                                           g2k=preddist_19sept)

## calculating CVs, CIs
boots_19sept<-do.bootstrap.cress.robust(final2d_19sept, predgrid, g2k=preddist_19sept, B=300)
cis_19sept<-makeBootCIs(boots_19sept)
predgrid$final2d_lower95_grid_19sept<-cis_19sept[,1]
predgrid$final2d_upper95_grid_19sept<-cis_19sept[,2]
predgrid$final2d_lower95_km_19sept<-cis_19sept[,1]/gridSize
predgrid$final2d_upper95_km_19sept<-cis_19sept[,2]/gridSize

predgrid$final2d_sd_19sept <- apply(boots_19sept, 1, sd)
predgrid$final2d_cv_19sept <- predgrid$final2d_sd_19sept/predgrid$final2dPredGrid_19sept

# plotting
maxD_19sept <- round(max(predgrid$final2dPredKm_19sept),2)

g1_19sept <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dPredKm_19sept)) +
  scale_fill_distiller(palette = "BrBG", breaks = seq(0,3,0.5), limits=c(0,round(maxD_19sept)),na.value = 'white') +
  xlab("x") + ylab("y") + theme_bw()+
  guides(fill=guide_legend(title="Density [km2]"))+
  ggtitle("19/09/2010")

g2_19sept <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2d_cv_19sept)) +
  scale_fill_distiller(palette = "BrBG") +
  xlab("x") + ylab("y") + theme_bw()+
  ggtitle("19/09/2010: uncertainty")+
  guides(fill=guide_legend(title="CV"))

#plot predicted densities again using more accessible colour palette

ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dPredKm_19sept)) +
  scale_fill_continuous_sequential(palette = "Teal", limits=c(0,round(maxD_19sept)),na.value = 'white') +
  xlab("x") + ylab("y") + theme_bw()+
  guides(fill=guide_legend(title="Density [km2]"))+
  ggtitle("19/09/2010")

#### 1 October
PAM_BAY_01oct <- PAM_1oct
PAM_BAY_01oct$layer <- as.integer(PAM_BAY_01oct$layer * gridSize)

# using only the grids which are in predgrid file
# already defined above in 19 Sept code
PAM_BAY_01oct <- PAM_BAY_01oct[c(toKeep),]

# SALSA requires the following column names to be included: response, x.pos, y.pos
seg_obs_PAM_01oct <- PAM_BAY_01oct %>% 
  rename("response" = "layer") %>%
  rename("x.pos" = "x") %>%
  rename("y.pos" = "y")

# adding environmental variables to the data
seg_obs_PAM_01oct$depth <- seg_obs$depth
seg_obs_PAM_01oct$sediment <- seg_obs$sediment
seg_obs_PAM_01oct$slope <- seg_obs$slope

# removing grids with NAs i.e., grids where there were no PAM instruments
seg_obs_PAM_01oct <- seg_obs_PAM_01oct[!is.na(seg_obs_PAM_01oct$response),]
hist(seg_obs_PAM_01oct$response)

# fitting the model
initialModel_01oct <- gamMRSea(response ~ 1, family = quasipoisson(link = "log"),
                               data = seg_obs_PAM_01oct)

# 1d model
# salsa1dlist already defined above for 19 Sept 
salsa1dRun_01oct<-runSALSA1D(initialModel_01oct, 
                             salsa1dlist, 
                             varlist = c("depth","slope"), 
                             datain = seg_obs_PAM_01oct,
                             removal = TRUE,
                             suppress.printout = TRUE)

summary(salsa1dRun_01oct$bestModel) 
salsa1dRun_01oct$keptvarlist
final1d_01oct <- salsa1dRun_01oct$bestModel

# 2d model
knotgrid_01oct<- getKnotgrid(coordData = cbind(seg_obs_PAM_01oct$x.pos, seg_obs_PAM_01oct$y.pos), 
                             numKnots = 300,
                             plot = FALSE)

distMats_01oct <- makeDists(cbind(seg_obs_PAM_01oct$x.pos, seg_obs_PAM_01oct$y.pos), knotgrid_01oct)


salsa2dlist_01oct<-list(fitnessMeasure = 'QBIC',
                        knotgrid = knotgrid_01oct,
                        startKnots=5,
                        minKnots=2,
                        maxKnots=15,
                        gap=0)

salsa2dRun_01oct<-runSALSA2D(model = final1d_01oct, 
                             salsa2dlist = salsa2dlist_01oct, 
                             d2k=distMats_01oct$dataDist,
                             k2k=distMats_01oct$knotDist)

final2d_01oct <- salsa2dRun_01oct$bestModel
summary(final2d_01oct)

# predictions for the final model
# predictions with 'Km' show density per km2 (as in the paper)
# predictions with 'Grid' show density per grid so # porpoises/16 km2

preddist_01oct<-makeDists(cbind(predgrid$x.pos, predgrid$y.pos),
                          knotgrid_01oct, knotmat=FALSE)$dataDist

predgrid$final2dPredKm_01oct <- (predict(object = final2d_01oct, 
                                         newdata = predgrid,
                                         g2k=preddist_01oct)[,1])/gridSize
predgrid$final2dPredGrid_01oct <- predict(object = final2d_01oct, 
                                          newdata = predgrid,
                                          g2k=preddist_01oct)

## calculating CVs, CIs
boots_01oct<-do.bootstrap.cress.robust(final2d_01oct, predgrid, g2k=preddist_01oct, B=300)
cis_01oct<-makeBootCIs(boots_01oct)
predgrid$final2d_lower95_grid_01oct<-cis_01oct[,1]
predgrid$final2d_upper95_grid_01oct<-cis_01oct[,2]
predgrid$final2d_lower95_km_01oct<-cis_01oct[,1]/gridSize
predgrid$final2d_upper95_km_01oct<-cis_01oct[,2]/gridSize

predgrid$final2d_sd_01oct <- apply(boots_01oct, 1, sd)
predgrid$final2d_cv_01oct <- predgrid$final2d_sd_01oct/predgrid$final2dPredGrid_01oct

# plotting
maxD_01oct <- round(max(predgrid$final2dPredKm_01oct),2)

g1_01oct <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dPredKm_01oct)) +
  scale_fill_distiller(palette = "BrBG", breaks = seq(0,round(maxD_01oct),1), limits=c(0,round(maxD_01oct)),na.value = 'white') +
  xlab("x") + ylab("y") + theme_bw()+
  guides(fill=guide_legend(title="Density [km2]"))+
  ggtitle("01/10/2010")

g2_01oct <- ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2d_cv_01oct)) +
  scale_fill_distiller(palette = "BrBG") +
  xlab("x") + ylab("y") + theme_bw()+
  ggtitle("19/09/2010: uncertainty")+
  guides(fill=guide_legend(title="CV"))

#plot predicted densities again using more accessible colour palette

ggplot(predgrid) + geom_tile(aes(x=x.pos, y=y.pos, fill=final2dPredKm_01oct)) +
  scale_fill_continuous_sequential(palette = "Teal", limits=c(0,round(maxD_01oct)),na.value = 'white') +
  xlab("x") + ylab("y") + theme_bw()+
  guides(fill=guide_legend(title="Density [km2]"))+
  ggtitle("01/10/2010")
