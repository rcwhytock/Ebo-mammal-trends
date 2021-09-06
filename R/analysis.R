# Large mammal data analysis Author: RCW Last
# updated 29/05/2019
library(MuMIn)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(mgcv)
library(visreg)
library(rgdal)
library(boot)

#### Read, tidy and summarise the raw data ####
LMdat <- read.csv("../Data/LM_data_March_2019.csv")

# Store vector of transect lengths (calculated
# externally)
transectLengths <- data.frame(
  RECCE.ID = LETTERS[1:23],
  RECCE.LENGTH = c(
    5.919,
    18.553,
    6.158,
    20.517,
    15.968,
    19.577,
    19.229,
    18.434,
    19.631,
    16.535,
    21.65,
    15.695,
    20.563,
    16.113,
    21.331,
    15.914,
    8.285,
    14.398,
    9.062,
    12.241,
    10.77,
    12.698,
    6.038
  )
)

str(LMdat)
transectLengths

LMdat <- merge(LMdat, transectLengths, all.x = T)
names(LMdat)

# Check all transects surveyed in all years
options(max.print = 2000)
table(LMdat$SURVEY.NO, LMdat$RECCE.ID)  # data from F was lost in 2017 due to rain, but otherwise fine

summary(LMdat)
unique(LMdat$YEAR)  # OK
unique(LMdat$SURVEY.NO)  # OK
unique(LMdat$MONTH)  # OK
unique(LMdat$SPECIES)
table(LMdat$SPECIES)  # most odd or uncertain records like 'H/BD' are very few and will exclude later

# Look at vegetation types in the different recces
table(LMdat$VEG.TYPE, LMdat$RECCE.ID)

# Check sign codes
barplot(table(LMdat$SIGN), las = 2)

# Look at species codes in the different recces
table(LMdat$SPECIES, LMdat$RECCE.ID)

# Look at species codes in the different recces
table(LMdat$SURVEY.NO, LMdat$YEAR)

# Look at summary of timing of surveys
table(LMdat$SURVEY.NO, LMdat$MONTH, LMdat$YEAR)

# Survey 1: January - June 2008 # change to 2008
# Survey 2: January - August 2012 ### Change to
# 2012 Survey 3: October - December 2016 and
# January - March 2017 ### Change to 2016 Survey 4:
# October - December 2017 and January - February
# 2018 ### Change to 2017 Survey 5: November -
# December 2018 and January to March 2019 ###
# Change to 2018

LMdat$SURVEY.NO <- ifelse(LMdat$SURVEY.NO == 6, 2019,
                          LMdat$SURVEY.NO)
LMdat$SURVEY.NO <- ifelse(LMdat$SURVEY.NO == 5, 2018,
                          LMdat$SURVEY.NO)
LMdat$SURVEY.NO <- ifelse(LMdat$SURVEY.NO == 4, 2017,
                          LMdat$SURVEY.NO)
LMdat$SURVEY.NO <- ifelse(LMdat$SURVEY.NO == 3, 2016,
                          LMdat$SURVEY.NO)
LMdat$SURVEY.NO <- ifelse(LMdat$SURVEY.NO == 2, 2012,
                          LMdat$SURVEY.NO)
LMdat$SURVEY.NO <- ifelse(LMdat$SURVEY.NO == 1, 2008,
                          LMdat$SURVEY.NO)

# Make bay duiker RD
LMdat$SPECIES <- ifelse(LMdat$SPECIES == "BD", "RD", LMdat$SPECIES )
# Need average canopy height for each transect (for
# later)
meanCanopy <- aggregate(as.numeric(LMdat$CANOPY) ~
                          LMdat$RECCE.ID + LMdat$SURVEY.NO,
                        FUN = mean)

head(meanCanopy)
colnames(meanCanopy) <- c("RECCE.ID", "SURVEY.NO",
                          "CANOPY")

# Convert VIS to numeric on 4 point scale
unique(LMdat$VIS)
summary(LMdat$VIS)
LMdat$VIS <- as.character(LMdat$VIS)
LMdat$VIS <- ifelse(LMdat$VIS == "VO", "4", LMdat$VIS)
LMdat$VIS <-
  ifelse(LMdat$VIS %in% c("0", "O", "o", "Vo"), "3", LMdat$VIS)
LMdat$VIS <- ifelse(LMdat$VIS == "C", "2", LMdat$VIS)
LMdat$VIS <- ifelse(LMdat$VIS == "VC", "1", LMdat$VIS)
LMdat$VIS <-
  ifelse(!LMdat$VIS %in% c("1", "2", "3", "4"), "2", LMdat$VIS)

summary(factor(LMdat$VIS))

# Get mean vis for each transect
meanVis <- aggregate(as.numeric(LMdat$VIS) ~
                       LMdat$RECCE.ID + LMdat$SURVEY.NO, FUN = mean)

colnames(meanVis) <- c("RECCE.ID", "SURVEY.NO",
                       "VIS")

#### for modeling (n > 100) are bay duiker (BA), blue
#### duiker (BD), chimpanzee (C), elephant (E), humans
#### (H), putty-nosed guenon (PU), red-river hog (RR)

LMdatReshapeAllSpecies <- droplevels(LMdat[LMdat$SIGN %in%
                                             c("F",
                                               "D",
                                               "VO",
                                               "O",
                                               #"T",
                                               "N",
                                               #"AC",
                                               #"CA",
                                               #"GS",
                                               #"SK",
                                               #"SN",
                                               "ROOTING"
                                             )
                                           &
                                             LMdat$AGE %in% c("F", "R", "O"),])


tableSummaries <-
  data.frame(table(
    LMdatReshapeAllSpecies$RECCE.ID,
    LMdatReshapeAllSpecies$SPECIES,
    LMdatReshapeAllSpecies$SURVEY.NO))

write.csv(pivot_wider(tableSummaries, values_from =  Freq, names_from = Var1), "../Data/tableSummaries.csv")

# Reshape the data for focal species
species <- c("BA","C", "E", "PU", "RR", "RD")

#

LMdatReshape <- droplevels(LMdat[LMdat$SPECIES %in%
                                   species & LMdat$SIGN %in%
                                   c("F",
                                     "D",
                                     "VO",
                                     "O",
                                     #"T",
                                     "N",
                                     #"AC",
                                     #"CA",
                                     #"GS",
                                     #"SK",
                                     #"SN",
                                     "ROOTING"
                                   )
                                 &
                                   LMdat$AGE %in% c("F", "R", "O"),])


# Simple summary tables of signs
simpleSummary <- table(LMdatReshape$SPECIES, LMdatReshape$RECCE.ID, LMdatReshape$SURVEY.NO)

summary(LMdatReshape)

LMdatReshape <- droplevels(LMdatReshape[, c(
  "RECCE.ID",
  "YEAR",
  "SURVEY.NO",
  "MONTH",
  "SPECIES",
  "SIGN",
  "AGE",
  "VEG.TYPE",
  "CANOPY",
  "RECCE.LENGTH",
  "VIS"
)])
LMdatReshape <- LMdatReshape[complete.cases(LMdatReshape),]


# Re-code veg type
LMdatReshape$VEG.TYPE <-
  as.factor(ifelse(
    as.character(LMdatReshape$VEG.TYPE) %in%
      c("FM", "OFS", "NFS", "FS"),
    as.character(LMdatReshape$VEG.TYPE),
    "OTHER"
  ))

LMdatReshape <- droplevels(LMdatReshape)
summary(LMdatReshape)
head(LMdatReshape)

# Get counts of species in each transect and each
# session
LMdatCountSp <- data.frame(table(
  LMdatReshape$RECCE.ID,
  LMdatReshape$SPECIES,
  LMdatReshape$SURVEY.NO
))
head(LMdatCountSp)

par(mfrow = c(1, 1))
plot(Freq ~ Var3, LMdatCountSp)
colnames(LMdatCountSp) <- c("RECCE.ID", "SPECIES",
                            "SURVEY.NO", "COUNT.SP")

# Join species,canopy data and vis data
LMdatCountSp$UniqueID <- paste(LMdatCountSp$RECCE.ID,
                               LMdatCountSp$SURVEY.NO, sep = "_")
meanCanopy$UniqueID <-
  paste(meanCanopy$RECCE.ID, meanCanopy$SURVEY.NO,
        sep = "_")
meanVis$UniqueID <- paste(meanVis$RECCE.ID, meanVis$SURVEY.NO,
                          sep = "_")

mergeDat <- merge(meanCanopy, LMdatCountSp, by = "UniqueID",
                  all.y = TRUE)

mergeDat <- merge(mergeDat, meanVis, by = "UniqueID", all.y = T)
head(mergeDat)

mergeDat <- mergeDat[, c(1, 4, 6, 8, 9, 10, 11)]
head(mergeDat)

mergeDat <- merge(mergeDat, transectLengths, by = "RECCE.ID")
head(mergeDat)

mamSummary <- mergeDat

# Make sure got 0s for all
# species/survey/month/recce combinations
mamSummary <- data.frame(mamSummary %>% complete(SPECIES,
                                                 nesting(SURVEY.NO,
                                                         RECCE.ID,
                                                         CANOPY,
                                                         RECCE.LENGTH)))
mamSummary[1:100,]

table(mamSummary$SURVEY.NO,
      mamSummary$SPECIES,
      mamSummary$RECCE.ID)


table(mamSummary$RECCE.ID, mamSummary$SURVEY.NO)


summary(mamSummary$SPECIES)

# Summary tables
sumStats <- aggregate(COUNT.SP ~ SPECIES + SURVEY.NO, FUN = sum, data = mamSummary)
sumStats$min <- aggregate(COUNT.SP ~ SPECIES + SURVEY.NO, FUN = min, data = mamSummary)[,3]
sumStats$max <- aggregate(COUNT.SP ~ SPECIES + SURVEY.NO, FUN = max, data = mamSummary)[,3]
sumStats$mean <- aggregate(COUNT.SP ~ SPECIES + SURVEY.NO, FUN = mean, data = mamSummary)[,3]
sumStats$sd <- aggregate(COUNT.SP ~ SPECIES + SURVEY.NO, FUN = sd, data = mamSummary)[,3]


write.csv(sumStats, "../Results/summaryStats.csv")

# Scale variables
mamSummaryScale <- data.frame(scale(mamSummary[, c(2, 4, 8)]))
mamSummaryScale$SPECIES <- mamSummary$SPECIES
mamSummaryScale$RECCE.ID <- as.factor(mamSummary$RECCE.ID)
mamSummaryScale$COUNT.SP <- mamSummary$COUNT.SP
mamSummaryScale$RECCE.LENGTH <- mamSummary$RECCE.LENGTH

head(mamSummaryScale)

##### Create models ####
modNull <- glmmTMB(
  COUNT.SP ~ 1 +
    (1 | RECCE.ID) +
    (SURVEY.NO | SPECIES),
  data = mamSummaryScale,
  family = nbinom2, REML = F
)

mod1 <- glmmTMB(
  COUNT.SP ~
    SURVEY.NO +
    VIS +
    offset(log(RECCE.LENGTH)) +
    (1 | RECCE.ID) +
    (SURVEY.NO | SPECIES),
  data = mamSummaryScale,
  family = nbinom2, REML = F
)

mod2 <- glmmTMB(
  COUNT.SP ~
    SURVEY.NO +
    I(SURVEY.NO ^ 2) +
    VIS +
    offset(log(RECCE.LENGTH)) +
    (1 | RECCE.ID) +
    (SURVEY.NO | SPECIES),
  data = mamSummaryScale,
  family = nbinom2, REML = F
)

mod3 <- glmmTMB(
  COUNT.SP ~
    SURVEY.NO +
    VIS +
    CANOPY +
    I(CANOPY ^ 2) +
    offset(log(RECCE.LENGTH)) +
    (1 | RECCE.ID) +
    (SURVEY.NO | SPECIES),
  data = mamSummaryScale,
  family = nbinom2, REML = F
)

mod4 <- glmmTMB(
  COUNT.SP ~
    SURVEY.NO +
    I(SURVEY.NO ^ 2) +
    VIS +
    CANOPY +
    I(CANOPY ^ 2) +
    offset(log(RECCE.LENGTH)) +
    (1 | RECCE.ID) +
    (SURVEY.NO | SPECIES),
  data = mamSummaryScale,
  family = nbinom2, REML = F
)


AICc(modNull, mod1, mod2, mod3, mod4)
summary(mod4) #

# Inference with REML = T
mod4 <- update(mod4, REML = T)
summary(mod4)


hist(resid(mod4))
plot(resid(mod4) ~ fitted(mod4))

# Calculate confidence intervals
b1 <- lme4::bootMer(mod4, FUN=function(x) fixef(x)$cond, nsim=500, .progress="txt")
if (requireNamespace("boot")) {
  boot.ci(b1, conf = 0.95, type="norm")
}

coefTab <- data.frame(coefTable(mod4))
coefTab$lowerCI <- NA
coefTab$upperCI <- NA

for(i in 1:nrow(coefTab)){
  
  coefTab[i, "lowerCI"] <- b1$t[order(b1$t[,i], decreasing = F),i][125]
  coefTab[i, "upperCI"] <- b1$t[order(b1$t[,i], decreasing = F),i][375]
  
}


# Create caterpillar plot
write.csv(coefTab, "../Results/mod4_coefs.csv")

# Plot species predictions
speciesCodes <- unique(mamSummary$SPECIES)
speciesCodes
speciesLabs <- as.character(speciesCodes)
speciesLabs <- c(
  "Blue duiker",
  "Chimpanzee",
  "Forest elephant",
  "Putty-nosed guenon",
  "Red duiker",
  "Red river hog"
)

plotPreds <- vector("list", length = length(speciesCodes))
for (i in 1:length(plotPreds)) {
  for (j in 1:23) {
    plotPreds[[i]] <- vector("list", length = 23)
  }
}


# Which transect is lcosest to 0 BLUP?
transectBLUPS <- exp(coef(mod4)$cond$RECCE.ID[1])
transectBLUPS$Name <- row.names(transectBLUPS)

RSTs <- read.csv("../Data/RST Coordinates.csv")
RSTs <- merge(RSTs, transectBLUPS)
names(RSTs)[7] <- "BLUP"
write.csv(RSTs, "../Data/transectBLUPS.csv", row.names = F)

for (i in 1:length(speciesCodes))
{
  for (j in 1:23) {
    newdat <-
      data.frame(
        SURVEY.NO = seq(
          min(mamSummaryScale$SURVEY.NO),
          max(mamSummaryScale$SURVEY.NO),
          length = 500
        ),
        RECCE.ID = factor(
          rep(LETTERS[1:23][j],
              length.out = 500),
          levels = levels(mamSummaryScale$RECCE.ID)
        ),
        SPECIES = factor(
          rep(speciesCodes[i], length.out = 500),
          levels = levels(mamSummaryScale$SPECIES)
        ),
        RECCE.LENGTH = mean(mamSummaryScale$RECCE.LENGTH),
        CANOPY = rep(mean(mamSummaryScale$CANOPY),
                     length = 500),
        VIS = rep(mean(mamSummaryScale$VIS),
                  length = 500)
        #HumanSign = rep(mean(mamSummaryScale$HumanSign),
        #length = 500)
      )
    
    plotPreds[[i]][[j]] <- predict(mod4, newdata = newdat)
    }
  
}

#### Plot results ####

# Create function to add alpha channel to colours
add.alpha <- function(col, alpha = 1)
{
  if (missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1],
        x[2], x[3], alpha = alpha))
}

newCol <- add.alpha("steelblue", alpha = 0.5)

# Setup device
par(mfrow = c(2, 3), oma = c(1,1,1,1), mar = c(8,5,1,1))

# Store the xlabs for unscale 
xlabs <- 2008:2019
xlabsScale <- (xlabs - mean(mamSummary$SURVEY.NO)) / sd(mamSummary$SURVEY.NO)


# Plot BLUPs
for (i in 1:length(plotPreds))
{
  newDatSpecies <- mamSummaryScale[mamSummaryScale$SPECIES ==
                                     speciesCodes[i],]
  
  
  plot(
    exp(plotPreds[[i]][[1]]) ~ newdat$SURVEY.NO,
    type = "l",
    lwd = 2,
    ylim = c(0,60),
    main = speciesLabs[i],
    col = newCol,
    ylab = "Relative encounter rate",
    xlab = "",
    cex.lab = 1.2,
    axes = F
  )
  
  mtext("Year", side = 1, line = 3.5, cex = 0.8)
  axis(2)
  axis(1, at = xlabsScale, labels = xlabs, las =2, cex.axis = 0.9)
  
  for (j in 1:23) {
    if (j == 1) {
    }
    else{
      lines(
        exp(plotPreds[[i]][[j]]) ~ newdat$SURVEY.NO,
        type = "l",
        lwd = 2,
        col = newCol
      )
    }
    
  }
  
  text(
    x = 2018.5,
    y = exp(plotPreds[[i]][[1]][500]),
    speciesLabs[i],
    cex = 0.9,
    pos = 4,
    col = rainbow(n = 7,
                  start = 0.3, end = 0.85)[i]
  )
  
  points(COUNT.SP ~ SURVEY.NO,
         data = newDatSpecies,
         col = newCol,
         pch = 16)
}

# Simple tables of species signs
table(mamSummaryScale$SPECIES, mamSummaryScale$RECCE.ID)

# Convert UTM to lat lon
LMdatReshapeAllSpecies$LAT <- as.numeric(as.character(LMdatReshapeAllSpecies$LAT))
LMdatReshapeAllSpecies$LONG <- as.numeric(as.character(LMdatReshapeAllSpecies$LONG))
LMdatReshapeAllSpecies_coords <- LMdatReshapeAllSpecies[complete.cases(LMdatReshapeAllSpecies$LAT),]

sputm <-
  SpatialPoints(
    data.frame(LAT = LMdatReshapeAllSpecies_coords$LAT, LON =  LMdatReshapeAllSpecies_coords$LONG),
    proj4string = CRS("+proj=utm +zone=32 +datum=WGS84")
  )
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
LMdatReshapeAllSpecies_coords <- cbind(LMdatReshapeAllSpecies_coords, data.frame(spgeo))
head(LMdatReshapeAllSpecies_coords)

names(LMdatReshapeAllSpecies_coords)[33] <- "LATITUDE"
names(LMdatReshapeAllSpecies_coords)[34] <- "LONGITUDE"

hist(LMdatReshapeAllSpecies_coords$LATITUDE) # remove anything greater than 20
hist(LMdatReshapeAllSpecies_coords$LONGITUDE) # remove anything greater than 10

LMdatReshapeAllSpecies_coords <- subset(LMdatReshapeAllSpecies_coords, LATITUDE < 20)
LMdatReshapeAllSpecies_coords <- subset(LMdatReshapeAllSpecies_coords, LONGITUDE < 5)

LMdatReshapeAllSpecies_coords <- droplevels(LMdatReshapeAllSpecies_coords)
unique(LMdatReshapeAllSpecies_coords$SPECIES)

table(LMdatReshapeAllSpecies_coords$SPECIES, LMdatReshapeAllSpecies_coords$SIGN)

# Get elephant GPS points
elephants <- subset(LMdatReshapeAllSpecies_coords, SPECIES == "E")
str(elephants)

elephants <- elephants[complete.cases(elephants$LAT),]
plot(elephants$LATITUDE, elephants$LONGITUDE)

head(elephants)
write.csv(elephants, "./Results/elephants.csv")

# Chimpanzees
chimpanzees <- subset(LMdatReshapeAllSpecies_coords, SPECIES == "C")
str(chimpanzees)
plot(chimpanzees$LAT ~ chimpanzees$LONG)

head(chimpanzees)
write.csv(chimpanzees, "./Results/chimpanzees.csv")
