  library(unmarked)
  library(tidyr)
  library(sf)
  library(ggplot2)
  library(MuMIn)
  library(cowplot)
  
  # Read in the species data
  list.files("../Data/occupancy/", full = TRUE)
  
  species <- read.csv("../Data/occupancy//elephantSmallGridUnion.csv")
  species <- species[- which(species$SURVEY.NO %in% c(2012, 2019)),] # 2012 has missing GIS coordinates in raw data, 2019 not fully surveyed due to COVID19 so many NAs 
  head(species)
  
  # Load in the grid cell data and sample level variables
  nestedGrid <- read.csv("../Data/occupancy//nestedGrid_variables.csv")
  head(nestedGrid)
  
  largeGridVars <- read.csv("../Data/occupancy//largeGridVariables.csv")
  head(largeGridVars)
  
  sampledGrid <- read.csv("../Data/occupancy//sampledSmallGrid_correct.csv") 
  head(sampledGrid)
  
  head(species)
  
  # Check species small grid IDs are correct
  species$smallGridID %in% sampledGrid$smallGridID # yes
  sampledGrid$smallGridID %in% nestedGrid$smallGridI # yes
  
  # Create variable indicating if small grid cell was surveyed
  nestedGrid$sampledYN <- ifelse(nestedGrid$smallGridI %in% sampledGrid$smallGridID, 1, 0)
  
  # Note that in 2019 large grid cells were not surveyd as follows (RSTs J,K,L,M,N,O)
#notSurveyd2019 <- c(483, 508, 484, 509, 534, 559, 485, 510, 535, 560, 585, 366, 391, 416, 441, 367, 392, 417, 442, 368, 393, 418)
  #notSurveyd2019 #### 2019 removed
  
  speciesMerge_summarised <- unique(species[,c("smallGridID", "SURVEY.NO")])
  head(speciesMerge_summarised)
  
  # Create subset of only sampled cells but with the large grid var in place
  sampledGrid <- subset(nestedGrid, sampledYN == 1)
  
  # Create index of surveys
  indexList <- vector("list", length = length(unique(sampledGrid$largeGridI)))
  sampledGrid$index <- NA
  
  for(i in 1:length(unique(sampledGrid$largeGridI))){
    
    newDat <- subset(sampledGrid, largeGridI == unique(sampledGrid$largeGridI)[i])
    newDat$index <- 1:nrow(newDat)
    indexList[[i]] <- newDat
    
  }
  
  indexList
  
  sampledGridMerged <- do.call("rbind", indexList)
  sampledGridMerged # now have an observation index per small grid cell nested in the large grid cell#
  
  # Remove the environmental vars from nested grid as they relate to small grid cells, which is wrong
  nestedGrid <- nestedGrid[,-c(3:6)]
  head(nestedGrid)
  
  sampledGridMerged <- sampledGridMerged[,-c(3:6)]
  names(sampledGridMerged)[1] <- "smallGridID"
  
  # For each year get detections at each sampled small grid cell
  years <- as.character(unique(speciesMerge_summarised$SURVEY.NO))
  years <- years[order(years)]
  years <- as.factor(years)
  years <-droplevels(years)
  yearList <- vector("list", length = length(years))
  
  for(i in 1:length(years)){
    
    speciesPresence <- speciesMerge_summarised[speciesMerge_summarised$SURVEY.NO == years[i],]
    speciesPresence <- unique(speciesPresence) # Get presence only in each cell
    speciesPresence$presAbs <- 1
    speciesPresence <- merge(sampledGridMerged, speciesPresence)
    speciesPresence <- merge(sampledGridMerged, speciesPresence, all.x = TRUE)
    speciesPresence$SURVEY.NO <- years[i]
    speciesPresence$presAbs <- ifelse(is.na(speciesPresence$presAbs), 0, 1)
    
    speciesPresence <- data.frame(pivot_wider(speciesPresence, names_from = index, id_cols = largeGridI, values_from = presAbs))
    speciesPresence <- data.frame(speciesPresence[order(speciesPresence$largeGridI),], row.names = 1:nrow(speciesPresence))
    
    yearList[[i]] <- as.matrix(speciesPresence[,c(2:6)])
    
    
  }
  
  # Remove the 2019 cells that were not surveyed #### 2019 removed
  yearList[[length(yearList)]][which(largeGridVars$largeGridID %in% notSurveyd2019),] <- NA
  yearList <- do.call("rbind", yearList)
  
  ## Setup for unmarked
  M <- 96 # Number of sites
  J <- 5 # num secondary sample periods (small grid)
  T <- 4 # num primary sample periods (years)
  
  yy <- matrix(yearList, M, J*T)
  
  years <- c('08','16','17','18')

  year <- matrix(years,
                 nrow(yy), T, byrow=TRUE)
  
  siteCovs <-largeGridVars[,c(2,3,5)]
  siteCovs$meanHanFor <- scale(siteCovs$meanHanFor)
  siteCovs$meanDistR <- scale(siteCovs$meanDistR)
  siteCovs$meanElev <- scale(siteCovs$meanElev)
  
  # Check for multicollinearity
  pairs(siteCovs)
  cor(siteCovs$meanHanFor, siteCovs$meanElev) # keep elevation
  
  speciesFrame <- unmarkedMultFrame(
    y = yy,
    yearlySiteCovs = list(year = year),
    siteCovs = siteCovs,
    numPrimary=T)
  speciesFrame
  summary(speciesFrame)
  
  # A model with constant parameters
  fm0 <- colext(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, pformula = ~1, speciesFrame, control = list(maxit = 100000))
  
  # Like fm0, but with year-dependent detection (observer variation)
  fm1 <- colext(psiformula =  ~1, gammaformula = ~1, epsilonformula = ~1, pformula = ~year, speciesFrame, control = list(maxit = 100000))
  
  # Like fm0, but with year-dependent colonization and extinction
  fm2 <- colext(psiformula = ~1, gammaformula = ~year, epsilonformula = ~year, pformula = ~1, speciesFrame, control = list(maxit = 100000))
  
  # A fully time-dependent model
  fm3 <- colext(psiformula = ~1, gammaformula = ~year, epsilonformula = ~year, pformula = ~year, speciesFrame, control = list(maxit = 100000))
  
  # Like fm3 with mean distance to road-dependence of 1st-year occupancy
  fm4 <- colext(psiformula = ~meanDistR, gammaformula = ~year, epsilonformula = ~year, pformula = ~year, speciesFrame, control = list(maxit = 100000))
  
  # Like fm4 with mean elevation dependence of 1st-year occupancy
  fm5 <- colext(psiformula = ~meanElev, gammaformula = ~year, epsilonformula = ~year, pformula = ~year, speciesFrame, control = list(maxit = 100000))

  # Like fm5 with mean elevation and mean dist r dependence of 1st-year occupancy ### Doesn't work for elephants ####
  #fm6 <- colext(psiformula = ~meanElev + meanDistR, gammaformula = ~year, epsilonformula = ~year, pformula = ~year, speciesFrame, control = list(maxit = 100000))
  
  
  models <- fitList(
    'psi(.)gam(.)eps(.)p(.)' = fm0,
    'psi(.)gam(.)eps(.)p(Y)' = fm1,
    'psi(.)gam(Y)eps(Y)p(.)' = fm2,
    'psi(.)gam(Y)eps(Y)p(Y)' = fm3,
    'psi(R)gam(Y)eps(Y)p(Y)' = fm4,
    'psi(E)gam(Y)eps(Y)p(Y)' = fm5
    #'psi(ER)gam(Y)eps(Y)p(Y)' = fm6
    
  )
  
  coef(models)
  SE(models)
  ms <- modSel(models, null = 'psi(.)gam(.)eps(.)p(.)')
  ms
  ms <- modSel(models)
  ms
  
  modTabs <- data.frame(model = ms@Full$model,
                        nPars = ms@Full$nPars,
                        AIC = ms@Full$AIC,
                        delta = ms@Full$delta)
  
  write.csv(ms@Full, "../Results/elephantModCoefs.csv")
  write.csv(modTabs, "../Results/elephantModAIC.csv")
  
  
# Keep models with cumulative weight up to 90

    models <- fm5 # overwhelming support
    
  
coef(models)
models

# Plot covariate effects
par(mfrow = c(1,1))
  
# Elevation
  nd <- data.frame(meanElev = seq(min(siteCovs$meanElev), max(siteCovs$meanElev), length.out = nrow(siteCovs)), meanDistR = 0)
  predDistR <- unmarked::predict(models, type = "psi", newdata = nd, appendData = TRUE) # initial occupancy
  plot(plogis(predDistR$Predicted) ~ predDistR$meanElev, ylim = c(0.2,1), type = "l")
  lines(plogis(predDistR$upper) ~ predDistR$meanElev, lty = 2)
  lines(plogis(predDistR$lower) ~ predDistR$meanElev, lty = 2)

# What are the estimated occupancy probabilities over the four years (\(\Psi_t\)) and corresponding confidence intervals now?
projected(fm5)
fm5 <- nonparboot(fm5, B = 500)
SE(fm5)

tempChange <- data.frame(cbind(projected=projected(fm5)[2,], SE=fm5@projected.mean.bsse[2,]))
tempChange$upper <- tempChange$projected + 1.96 * tempChange$SE
tempChange$lower <- tempChange$projected - 1.96 * tempChange$SE
tempChange$year <- c(2008, 2016,2017,2018)

pdf("../Results/elephantAnnualOccu.pdf", width = 4, height = 4)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0,0))
  plot(projected ~ year, tempChange, ylim = c(0,1), ylab = expression(paste("Occupancy ", psi)),
       xlab = "Year", pch = 16)
  arrows(x0 = tempChange$year,
         x1 = tempChange$year,
         y0 = tempChange$lower,
         y1 = tempChange$upper, code =3, length = 0.05, angle = 90)
dev.off()

# Map
largeGridShp <- st_read("../Data/shapefiles/grid/largeGrid.shp")
ebo <- st_read("../Data/shapefiles/oldENP/Ebo National Park_old.shp")

newDat <- siteCovs
psiPred <- unmarked::predict(fm5, type = "psi", newdata = newDat, appendData = TRUE)
psiPred$id <- largeGridVars$largeGridID

# Check everything aligns with grid IDs
psiPred$id == largeGridShp$id
largeGridShp$Psi <- psiPred$Predicted

min(psiPred$Predicted)
max(psiPred$Predicted)

pdf("../Results/elephantOccuMap.pdf", width = 4, height = 4)

ggplot(data = largeGridShp) +
  geom_sf(aes(fill = Psi)) +
  geom_sf(data = ebo, col = "gray80", fill = NA) +
  scale_fill_viridis_c(option = "viridis", limits = c(0,1)) + # set from min max predicted psi
  theme_classic()

dev.off()

# Confidence intervals of top model
confintList <- list(
  confint(fm5, type = "psi"),
  confint(fm5, type = "col"),
  confint(fm5, type = "ext"),
  confint(fm5, type = "det"))

confintListDF <- do.call("rbind", confintList)

modCoefs <- data.frame(estimate = (coef(fm6)))

modCoefs <- cbind(modCoefs, confintListDF[,c(1,2)])

modCoefs$estTrans <- inv.logit(modCoefs[,c(1)])
modCoefs$lowerTrans <- inv.logit(modCoefs[,c(2)])
modCoefs$upperTrans <- inv.logit(modCoefs[,c(3)])
modCoefs <- round(modCoefs, digits = 2)

write.csv(modCoefs, "../Results/elephantTopOccuConfints.csv")

save.image(file="../Results/elephantOccupancy.RData")













