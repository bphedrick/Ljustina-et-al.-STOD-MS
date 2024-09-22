#Code for:
#Eastern Musk Turtles exhibit multiple elevated stress responses to urbanization in Southern Louisiana
#Published in:
#DOI: 

#Load relevant packages
library(geomorph) #v.4.0.4
library(emmeans) #v. 1.8.9
library(data.table) #v. 1.14.2
library(ggstatsplot) #v. 0.12.1
library(ggplot2) #v. 3.4.4
library(RColorBrewer) #v. 1.1-3
library(ggpubr) #v. 0.6.0
library(lme4) #v. 1.1-33
library(dplyr) #v. 1.1.4
library(nlme) #v. 3.1-159

#Set a palette
palette(c("#0099FF", "red"))

#Set working directory
setwd("~/Desktop/Current Projects/Oliver-MuskTurtle/RCode/Ljustina et al., STOD CODE")


####Body Condition####

#Load data and log10 transform data and generate BMI (SCL/mass)
bodyCond <- read.csv("STOD_Field_Book_Fixed.csv")
    bodyCond$log10Mass <- log10(bodyCond$Mass_g)
    bodyCond$log10SCL <- log10(bodyCond$SCL_mm)
    bodyCond$BMI <- bodyCond$SCL_mm / bodyCond$Mass_g

#Body Condition Assessment 1: Body Mass
    
#Check whether body mass differs by site and by sex within sites
massANOVA <- aov(log10Mass ~ Site + Site:Sex, data = bodyCond) 
      summary(massANOVA)

TukeyHSD(massANOVA) #Tukey posthoc test to look at pairwise comparisons


#Make a plot for body mass by site
bodyMassPlot <- ggbetweenstats(
  data = bodyCond,
  x = Site,
  y = log10Mass) +
  labs(x = "", y = "log10(Body Mass)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 18)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.title = element_text(size =　18)) +
  theme(legend.position="none") + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

bodyMassPlot

#Reorganize data so you have sexes split
bodyCond$sexSite <- paste(bodyCond$Site, "_", bodyCond$Sex, sep = "")

sexSitePlot <- ggbetweenstats(
  data = bodyCond,
  x = sexSite,
  y = log10Mass) + 
  labs(x = "", y = "log10(Body Mass)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 18)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.title = element_text(size =　18)) +
  theme(legend.position="none") + 
  ggplot2::scale_color_manual(values = c("red", "red", "#0099FF", "#0099FF"))
  
sexSitePlot #This shows what the pairwise tests above showed


#Second look at the residuals of mass and SVL, which is robust to differences in body size

#Split turts by site into two separate dataframes in case there are different trajectories between sites
kennerTurts <- bodyCond[which(bodyCond$Site == "Kenner"), ]
joyceTurts <- bodyCond[which(bodyCond$Site == "Joyce"), ]

#lm of SCL and mass to extract residuals for Kenner
kennerLM <- lm(log10(SCL_mm) ~ log10(Mass_g), data = kennerTurts)
      plot(log10(SCL_mm) ~ log10(Mass_g), data = kennerTurts)

#lm of SCL and mass to extract residuals for Joyce
joyceLM <- lm(log10(SCL_mm) ~ log10(Mass_g), data = joyceTurts)
      plot(log10(SCL_mm) ~ log10(Mass_g), data = joyceTurts)

#Residuals dataframe
residDF <- data.frame(c(kennerLM$residuals, joyceLM$residuals), 
                         c(rep("Kenner", times = length(kennerLM$residuals)), rep("Joyce", times = length(joyceLM$residuals)))) 
      colnames(residDF) <- c("resids", "site")

#Run t-test comparing size residuals by site
t.test(resids ~ site, data = residDF)      

bmiResidsPlot <- ggbetweenstats(
  data = residDF,
  x = site,
  y = resids) + 
  labs(x = "", y = "BMI Residuals") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 18)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.title = element_text(size =　18)) +
  theme(legend.position="none") + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

bmiResidsPlot  



####Potential Stressors####

leeches <- na.omit(bodyCond[, 1:13]) #Get rid of turtles that didn't have leeches counted

sd(leeches[leeches$Site == "Joyce", ]$Leech.Count)
sd(leeches[leeches$Site == "Kenner", ]$Leech.Count)

#Run t-test between leech count and site
t.test(Leech.Count ~ Site, data = leeches)

#Make a plot for body mass by site
leechPlot <- ggbetweenstats(
  data = leeches,
  x = Site,
  y = Leech.Count) +
  labs(x = "", y = "Number of Leeches") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 18)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.title = element_text(size =　18)) +
  theme(legend.position="none") + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

leechPlot #There are more leeches on turtles at kenner than at Joyce 


#Compare water temp between sites based on our field measurements
t.test(WaterTemp_C ~ Site, data = bodyCond) #Not sig different


#Compare water quality data across sites
#Paired Mann-Whitney U (Wilcoxonrank-sum) tests on water quality data. Default 
#hypothesis test is two-tailed.

Water <- read.csv("Public_WaterQual_R.csv")

DO <- wilcox.test(DO ~ Site, data = Water, exact = FALSE)
Cond <- wilcox.test(Cond ~ Site, data = Water, exact = FALSE)
Temp <- wilcox.test(Temp ~ Site, data = Water, exact = FALSE)
Sal <- wilcox.test(Sal ~ Site, data = Water, exact = FALSE)
pH <- wilcox.test(pH ~ Site, data = Water, exact = FALSE)

#Put it in a matrix
waterQual <- matrix(nrow = 5, ncol = 2)
rownames(waterQual) <- c("DO", "Cond", "Temp", "Sal", "pH")
colnames(waterQual) <- c("stat", "p")

waterQual[1, 1] <- DO$statistic
waterQual[2, 1] <- Cond$statistic
waterQual[3, 1] <- Temp$statistic
waterQual[4, 1] <- Sal$statistic
waterQual[5, 1] <- pH$statistic
waterQual[1, 2] <- DO$p.value
waterQual[2, 2] <- Cond$p.value
waterQual[3, 2] <- Temp$p.value
waterQual[4, 2] <- Sal$p.value
waterQual[5, 2] <- pH$p.value

waterQual


####Shape Data and Prelim Analyses####
#Load TPS files
Dorsal <- readland.tps("STOD_Dorsal.TPS", specID = c("imageID"))
Ventral <- readland.tps("STOD_Ventral.TPS", specID = c("imageID"))

dorsalSiteRep <- read.csv("Site_Replicate_Dorsal.csv")
ventralSiteRep <- read.csv("Site_Replicate_Ventral.csv")

#Run Generalized Procrustes analysis
gpaDorsal <- gpagen(Dorsal, ProcD = FALSE, approxBE = TRUE)
gpaVentral <- gpagen(Ventral, ProcD = FALSE, approxBE = TRUE)

#Run a principal component analysis
pcaDorsal <- gm.prcomp(gpaDorsal$coords)
    dorsalPCASum <- summary(pcaDorsal)
pcaVentral <- gm.prcomp(gpaVentral$coords)
    ventralPCASum <- summary(pcaVentral)

#Creates geomorph dataframes
dorsalGDF <- geomorph.data.frame(shape = gpaDorsal$coords,
                                 size = log10(gpaDorsal$Csize),
                                 ind = dorsalSiteRep$individual,
                                 site = dorsalSiteRep$site,
                                 replicate = dorsalSiteRep$replicate)

ventralGDF <- geomorph.data.frame(shape = gpaVentral$coords,
                                 size = log10(gpaVentral$Csize),
                                 ind = ventralSiteRep$individual,
                                 site = ventralSiteRep$site,
                                 replicate = ventralSiteRep$replicate)

#PCA plots to explore shape space
par(mfrow = c(1, 2))

dorsalPCAPlot <- plot(pcaDorsal,
                      cex = 3,
                      pch = 21,
                      bg = as.factor(dorsalGDF$site),
                      col = c("black"),
                      xlab = paste("Principal Component 1 ", "(", sep = "",
                                   paste(round(dorsalPCASum$PC.summary$Comp1[2]*100, digits = 2),"%", ")", sep = "")),
                      ylab = paste("Principal Component 2 ", "(", sep = "",
                                   paste(round(dorsalPCASum$PC.summary$Comp2[2]*100, digits = 2),"%", ")", sep = "")),
                      cex.lab = 1,
                      cex.axis = 0.8)

ventralPCAPlot <- plot(pcaVentral,
                      cex = 3,
                      pch = 21,
                      bg = as.factor(ventralGDF$site),
                      col = c("black"),
                      xlab = paste("Principal Component 1 ", "(", sep = "",
                                   paste(round(ventralPCASum$PC.summary$Comp1[2]*100, digits = 2),"%", ")", sep = "")),
                      ylab = paste("Principal Component 2 ", "(", sep = "",
                                   paste(round(ventralPCASum$PC.summary$Comp2[2]*100, digits = 2),"%", ")", sep = "")),
                      cex.lab = 1,
                      cex.axis = 0.8)

#Look at size differences between groups for the carapace and plastron
violinDorsalDF <- data.frame(dorsalGDF$site, dorsalGDF$size)

violinDorsal <- ggplot(violinDorsalDF, aes(x = dorsalGDF$size, y = dorsalGDF$site)) +
  geom_violin(trim = FALSE) +
  coord_flip() +
  xlab("Carpace log10(Centroid Size)") +
  ylab("") +
  geom_boxplot(width = 0.2, fill = "white") +
  theme_classic() +
  scale_fill_manual(palette = "blue", "red")

violinDorsal


violinVentralDF <- data.frame(ventralGDF$site, ventralGDF$size)

violinVentral <- ggplot(violinVentralDF, aes(x = ventralGDF$size, y = ventralGDF$site)) +
  geom_violin(trim = FALSE) +
  coord_flip() +
  xlab("Plastron log10(Centroid Size)") +
  ylab("") +
  geom_boxplot(width = 0.2, fill = "white") +
  theme_classic() +
  scale_fill_manual(palette = "blue", "red")

violinVentral

            
####Asymmetry####

#Dorsal Asymmetry Analyses
#Import "LM Pairs csv file#
LMPairsDorsal <- read.csv("LMPairsSTODDorsal.csv", header = FALSE)

#Look at asymmetry generally for the carapace
dorsalSymm <- bilat.symmetry(A = shape, 
                             ind = ind, 
                             rep = rep, 
                             object.sym = TRUE, 
                             land.pairs = LMPairsDorsal, 
                             data = dorsalGDF, 
                             RRPP = TRUE, 
                             iter = 1000)

summary(dorsalSymm)

#Look at asymmetry for Joyce and Kenner separately for the carapace
#Kenner sites are 1:78 and Joyce sites are 79:162
#Kenner
dorsalKennerGDF <- geomorph.data.frame(shape = gpaDorsal$coords[, , 1:78],
                                       size = gpaDorsal$Csize[1:78],
                                       ind = dorsalSiteRep$individual[1:78],
                                       site = dorsalSiteRep$site[1:78],
                                       replicate = dorsalSiteRep$replicate[1:78])

#Joyce
dorsalJoyceGDF <- geomorph.data.frame(shape = gpaDorsal$coords[, , 79:162],
                                       size = gpaDorsal$Csize[79:162],
                                       ind = dorsalSiteRep$individual[79:162],
                                       site = dorsalSiteRep$site[79:162],
                                       replicate = dorsalSiteRep$replicate[79:162])

#Rerun within site asymmetry analyses
kennerDorsalSymm <- bilat.symmetry(A = shape, 
                                   ind = ind, 
                                   rep = rep, 
                                   object.sym = TRUE, 
                                   land.pairs = LMPairsDorsal, 
                                   data = dorsalKennerGDF, 
                                   RRPP = TRUE, 
                                   iter = 1000)
summary(kennerDorsalSymm)

#Repeating for Joyce#
joyceDorsalSymm <- bilat.symmetry(A = shape, 
                                   ind = ind, 
                                   rep = rep, 
                                   object.sym = TRUE, 
                                   land.pairs = LMPairsDorsal, 
                                   data = dorsalJoyceGDF, 
                                   RRPP = TRUE, 
                                   iter = 1000)
summary(joyceDorsalSymm)

#R square values suggest shell asymmetry explains greater amount of variation
#in Kenner turts than Joyce turts


##Now, we can go about determining whether magnitudes of shell asymmetry differ between sites
#First, determine whether landmark scheme is correctly formatted

#Flip the y coordinates for half of the landmarks to calculate Procrustes distances between sides
#1 is midline (red)
#2 is left side (blue)
#3 is right side (gold)

lmSideDorsal <- read.csv("lmSideDorsal.csv", header = FALSE)
    lmSideDorsal

#Set a palette and look that it worked
palette(c("red", "lightblue", "gold"))

plot(gpaDorsal$coords[, , 50], 
     pch = 21,
     cex = 3,
     bg = lmSideDorsal$V1)
#Looks good#

#Loop to make the second column (y axis) positive
dorsalGPAFlipped <- gpaDorsal$coords

#Take total length of array and divide by number of landmarks and dimensions
for(i in 1:dim(dorsalGPAFlipped)[3]){
  dorsalGPAFlipped[, 2, i] <- abs(dorsalGPAFlipped[, 2, i])
}

#This will show the x and y coordinates to make sure that it worked properly and 
#paired landmarks are reflected on top of one another
plot(dorsalGPAFlipped[, 2, 50] ~ dorsalGPAFlipped[, 1, 50],
     pch = 21,
     cex = 3,
     bg = lmSideDorsal$V1)

text(dorsalGPAFlipped[, 2, 50] ~ dorsalGPAFlipped[, 1, 50],
     labels = c(1:dim(dorsalGPAFlipped)[1]), cex = 1, col = c("black"))

#Calculate Procrustes distances between sides
#Procrustes distance is the sum of distances between corresponding landmarks of two shapes
#The distance between landmarks is sqrt((x2 - x1)^2 + (y2 - y1)^2)

#Make a set of lefts and a set of rights
colnames(LMPairsDorsal) <- c("left", "right") #These are lefts

left <- dorsalGPAFlipped[LMPairsDorsal$left, , ]
right <- dorsalGPAFlipped[LMPairsDorsal$right, , ]

lmPDist <- matrix(nrow = dim(left)[1], ncol = 1)
specDist <- matrix(nrow = dim(left)[3], ncol = 1)
allLMDists <- matrix(nrow= dim(left)[1], ncol = dim(left)[3])

#j is the individual. i is the landmark set of x, y coordinates
for(j in 1:dim(left)[3]){
  
  for(i in 1:nrow(lmPDist)){
    lmPDist[i] <- sqrt((left[i, 1, j] - right[i, 1, j])^2 + (left[i, 2, j] - right[i, 2, j])^2)
  }
  
  specDist[j] <- sum(lmPDist)
  
  allLMDists[, j] <- cbind(lmPDist)
}

dorsalSiteRep$specDist <- specDist

#Which landmark pair was the most asymmetric?
meanLMDist <- matrix(nrow = dim(left)[1], ncol = 1)

for(i in 1:dim(left)[1]){
  meanLMDist[i] <- mean(allLMDists[i, ])
  
}

LMPairsDorsal$mean <- meanLMDist
    LMPairsDorsal$mean

#Plot asymmetry magnitudes for Joyce and Kenner dorsal view as a violin plot    
    
#Average the magnitudes of asymmetry for the three replicates for the plot
    dorsalSiteRep
    
dorsalMeans <- dorsalSiteRep %>% group_by(individual) %>% summarise(mean = mean(specDist))
    
dorsalMeans$site <- c(dorsalKennerGDF$site[1:(length(dorsalKennerGDF$site)/3)], dorsalJoyceGDF$site[1:(length(dorsalJoyceGDF$site)/3)])

#plot it with individual means    
dorsalAsymPlot <- ggbetweenstats(
  data = dorsalMeans,
  x = site,
  y = mean) + 
  labs(x = "", y = "Asymmetry Magnitude") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 18)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.title = element_text(size =　18)) +
  theme(legend.position="none") + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

dorsalAsymPlot 

#method = ML fits the model using maximum likelihood, which is needed for model reduction 
#(not necessary with only one factor; in that case, or when you have arrived at the minimum adequate model, 
#you just leave out the method argument and it fits the model using REML).
DorsalAsymmANOVA <- lme(specDist ~ site, 
                        random = ~1|individual, 
                        method = "ML", 
                        data = dorsalSiteRep)

summary(DorsalAsymmANOVA)

PostSiteComp <- emmeans(DorsalAsymmANOVA, specs = "site")
    PostSiteComp

    
#Ventral Asymmetry Analyses
#Import "LM Pairs csv file#
LMPairsVentral <- read.csv("LMPairsSTODVentral.csv", header = FALSE)

ventralSymm <- bilat.symmetry(A = shape, 
                             ind = ind, 
                             replicate = replicate, 
                             object.sym = TRUE, 
                             land.pairs = LMPairsVentral, 
                             data = ventralGDF, 
                             RRPP = TRUE, 
                             iter = 1000)

summary(ventralSymm)

#Kenner
ventralKennerGDF <- geomorph.data.frame(shape = gpaVentral$coords[, , 1:102],
                                       size = gpaVentral$Csize[1:102],
                                       ind = ventralSiteRep$individual[1:102],
                                       site = ventralSiteRep$site[1:102],
                                       replicate = ventralSiteRep$replicate[1:102])

#Joyce
ventralJoyceGDF <- geomorph.data.frame(shape = gpaVentral$coords[, , 103:192],
                                      size = gpaVentral$Csize[103:192],
                                      ind = ventralSiteRep$individual[103:192],
                                      site = ventralSiteRep$site[103:192],
                                      replicate = ventralSiteRep$replicate[103:192])

#Rerun within site asymmetry analyses
kennerVentralSymm <- bilat.symmetry(A = shape, 
                                   ind = ind, 
                                   rep = replicate, 
                                   object.sym = TRUE, 
                                   land.pairs = LMPairsVentral, 
                                   data = ventralKennerGDF, 
                                   RRPP = TRUE, 
                                   iter = 1000)
summary(kennerVentralSymm)

#Repeating for Joyce#
joyceVentralSymm <- bilat.symmetry(A = shape, 
                                  ind = ind, 
                                  rep = replicate, 
                                  object.sym = TRUE, 
                                  land.pairs = LMPairsVentral, 
                                  data = ventralJoyceGDF, 
                                  RRPP = TRUE, 
                                  iter = 1000)
summary(joyceVentralSymm)

#R square values suggest shell asymmetry explains greater amount of variation
#in Kenner turts than Joyce turts#


##Now, we can go about determining whether magnitudes of shell asymmetry
##differ between sites##

#First, determine whether landmark scheme is correctly formatted#
#Flip the y coordinates for half of the landmarks to calculate PD 
#1 is midline (red)
#2 is left side (blue)
#3 is right side (gold)

lmSideVentral <- read.csv("lmSideVentral.csv", header = FALSE)
lmSideVentral

#Set a palette and look that it worked
palette(c("red", "lightblue", "gold"))

plot(gpaVentral$coords[, , 50], 
     pch = 21,
     cex = 3,
     bg = lmSideVentral$V1)
#Looks good#

#Loop to make the second column (y axis) positive
ventralGPAFlipped <- gpaVentral$coords

#Take total length of array and divide by number of landmarks and dimensions
for(i in 1:dim(ventralGPAFlipped)[3]){
  ventralGPAFlipped[, 2, i] <- abs(ventralGPAFlipped[, 2, i])
}

#This will show the x and y coordinates to make sure that it worked properly and 
#paired landmarks are reflected on top of one another
plot(ventralGPAFlipped[, 2, 50] ~ ventralGPAFlipped[, 1, 50],
     pch = 21,
     cex = 3,
     bg = lmSideVentral$V1)

text(ventralGPAFlipped[, 2, 50] ~ ventralGPAFlipped[, 1, 50],
     labels = c(1:dim(ventralGPAFlipped)[1]), cex = 1, col = c("black"))

#Calculate Procrustes distances between sides
#Procrustes distance is the sum of distances between corresponding landmarks of two shapes
#The distance between landmarks is sqrt((x2 - x1)^2 + (y2 - y1)^2)

#Make a set of lefts and a set of rights
colnames(LMPairsVentral) <- c("left", "right") #These are lefts

leftVentral <- ventralGPAFlipped[LMPairsVentral$left, , ]
rightVentral <- ventralGPAFlipped[LMPairsVentral$right, , ]

lmPDistVentral <- matrix(nrow = dim(leftVentral)[1], ncol = 1)
specDistVentral <- matrix(nrow = dim(leftVentral)[3], ncol = 1)
allLMDistsVentral <- matrix(nrow= dim(leftVentral)[1], ncol = dim(leftVentral)[3])

#j is the individual. i is the landmark set of x, y coordinates
for(j in 1:dim(leftVentral)[3]){
  
  for(i in 1:nrow(lmPDistVentral)){
    lmPDistVentral[i] <- sqrt((leftVentral[i, 1, j] - rightVentral[i, 1, j])^2 + (leftVentral[i, 2, j] - rightVentral[i, 2, j])^2)
  }
  
  specDistVentral[j] <- sum(lmPDistVentral)
  
  allLMDistsVentral[, j] <- cbind(lmPDistVentral)
}

ventralSiteRep$specDist <- specDistVentral

#Which landmark pair was the most asymmetric?

meanLMDistVentral <- matrix(nrow = dim(leftVentral)[1], ncol = 1)

for(i in 1:dim(leftVentral)[1]){
  meanLMDistVentral[i] <- mean(allLMDistsVentral[i, ])
  
}

LMPairsVentral$mean <- meanLMDistVentral
LMPairsVentral$mean

#Plot asymmetry magnitudes for Joyce and Kenner ventral view as a violin plot    

#Average the three replicates for the plot
ventralSiteRep

ventralMeans <- ventralSiteRep %>% group_by(individual) %>% summarise(mean = mean(specDist))

ventralMeans$site <- c(ventralKennerGDF$site[1:(length(ventralKennerGDF$site)/3)], ventralJoyceGDF$site[1:(length(ventralJoyceGDF$site)/3)])

#plot it with individual means    
ventralAsymPlot <- ggbetweenstats(
  data = ventralMeans,
  x = site,
  y = mean) + 
  labs(x = "", y = "Asymmetry Magnitude") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 18)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.title = element_text(size =　18)) +
  theme(legend.position="none") + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

ventralAsymPlot 

#method = ML fits the model using maximum likelihood, which is needed for model reduction 
#(not necessary with only one factor; in that case, or when you have arrived at the minimum adequate model, 
#you just leave out the method argument and it fits the model using REML).
VentralAsymmANOVA <- lme(specDist ~ site, 
                        random = ~1|individual, 
                        method = "ML", 
                        data = ventralSiteRep)

summary(VentralAsymmANOVA)

PostSiteCompVentral <- emmeans(VentralAsymmANOVA, specs = "site")
PostSiteCompVentral

    
    
####Cell Counts####

#Load cell count data
cellRatios <- read.csv("AllCellCountsFixed.csv")
  cellRatios$turtles <- paste(cellRatios$Turtle, "_", cellRatios$Site, sep = "")
  cellRatios$Ratio <- cellRatios$Heterophil/cellRatios$Lymphocyte

#plot it with individual means    
HLPlot <- ggbetweenstats(
  data = cellRatios,
  x = Site,
  y = Ratio) + 
  labs(x = "", y = "H:L Ratio") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 18)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.title = element_text(size =　18)) +
  theme(legend.position="none") + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

HLPlot 

#Check for significant difference in ratios by location
allCellsANOVA <- aov(Ratio ~ Site, data = cellRatios)

summary(allCellsANOVA)

emmeans(allCellsANOVA, specs = "Site")


