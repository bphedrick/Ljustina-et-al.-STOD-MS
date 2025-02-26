#Code for:
#Eastern Musk Turtles exhibit multiple elevated stress responses to urbanization in Southern Louisiana
#Published in:
#DOI: 

#Revised 2/25/25

#Load relevant packages
library(data.table) #v. 1.14.2
library(dplyr) #v. 1.1.4
library(emmeans) #v. 1.8.9
library(geomorph) #v.4.0.4
library(ggplot2) #v. 3.4.4
library(ggpubr) #v. 0.6.0
library(ggstatsplot) #v. 0.12.1
library(gridExtra)
library(lme4) #v. 1.1-33
library(nlme) #v. 3.1-159
library(RColorBrewer) #v. 1.1-3

#Set a palette
palette(c("#0099FF", "red","#D4AF37")) #Kenner, Joyce, St. Marks (museum specimens, only in asymmetry analyses)#

#Set working directory
setwd("~/Desktop/Current Projects/Oliver-MuskTurtle/RCode/STODFinalOL")

####Cell Counts####

#Load cell count data
cellRatios <- read.csv("AllCellCountsFixed.csv")
    cellRatios$turtles <- paste(cellRatios$Turtle, "_", cellRatios$Site, sep = "")
    cellRatios$Ratio <- cellRatios$Heterophil/cellRatios$Lymphocyte

#Check for significant difference in ratios by location with size and sex as additional factors
allCellsSexANOVA <- aov(Ratio ~ Site + Sex + SCL, data = cellRatios)
    summary(allCellsSexANOVA)
emmeans(aov(Ratio ~ Site, data = cellRatios), specs = "Site") #Only consider site since other variables were not significant
    
#plot it with individual means    
HLPlot <- ggbetweenstats(
  data = cellRatios,
  x = Site,
  y = Ratio) + 
  labs(x = "", y = "H:L Ratio") + 
  theme_classic() + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

HLPlot 


#Looking at just heterophils#
HeteroANOVA <- aov(Heterophil ~ Site + Sex + SCL, data = cellRatios)
    summary(HeteroANOVA)
emmeans(aov(Heterophil ~ Site, data = cellRatios), specs = "Site") #Just considers site

#plot it with individual means    
HPlot <- ggbetweenstats(
  data = cellRatios,
  x = Site,
  y = Heterophil) + 
  labs(x = "", y = "Heterophil count") + 
  theme_classic() + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

HPlot #Heterophil counts alone are significantly elevated in Kenner#

#Looking at just lymphocytes#
LymphoANOVA <- aov(Lymphocyte ~ Site + Sex + SCL, data = cellRatios)
    summary(LymphoANOVA)
emmeans(aov(Lymphocyte ~ Site, data = cellRatios), specs = "Site") #Just considers site
    
#plot it with individual means    
LPlot <- ggbetweenstats(
  data = cellRatios,
  x = Site,
  y = Lymphocyte) + 
  labs(x = "", y = "Lmyphocyte count") + 
  theme_classic() + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

LPlot #Lymphocyte counts alone are significantly elevated in Joyce#

grid.arrange(HPlot, LPlot, ncol=2)

####Body Condition####

#Load data and log10 transform data and generate BMI (SCL/mass)
bodyCond <- read.csv("STOD_Field_Book_Fixed.csv")

bodyCondJK <- bodyCond[(bodyCond$Site == "JWMA" | bodyCond$Site == "Kenner" ), ]
bodyCondJK$log10Mass <- log10(bodyCondJK$Mass_g)
bodyCondJK$log10SCL <- log10(bodyCondJK$SCL_mm)

#Check whether body mass differs by site and by sex within sites (does not include St. Marks museum specimens)
massANOVA <- aov(log10Mass ~ Site + Site:Sex, data = bodyCondJK)
summary(massANOVA)

TukeyHSD(massANOVA) #Tukey posthoc test to look at pairwise comparisons

#Reorganize data so you have sexes split for mass
bodyCondJK$sexSite <- paste(bodyCondJK$Site, "_", bodyCondJK$Sex, sep = "")

sexSiteMassPlot <- ggbetweenstats(
  data = bodyCondJK,
  x = sexSite,
  y = log10Mass) + 
  labs(x = "", y = "log10(Body Mass)") + 
  theme_classic() + 
  ggplot2::scale_color_manual(values = c("red", "red", "#0099FF", "#0099FF"))

sexSiteMassPlot #This shows what the pairwise tests above showed


#Second look at the residuals of mass and SCL, which is robust to differences in body size (St. Marks specimens not included)
turtLM <- lm(log10(SCL_mm) ~ log10(Mass_g), data = bodyCondJK)
      summary(turtLM)
      plot(log10(SCL_mm) ~ log10(Mass_g), data = bodyCondJK, pch = 21, bg = as.factor(bodyCondJK$Site), col = "black")

      t.test(turtLM$residuals ~ bodyCondJK$Site)
                            
      residsDF <- data.frame(resids = c(turtLM$residuals), 
                             site = c(bodyCondJK$Site)) 
      
bmiResidsPlot <- ggbetweenstats(
        data = residsDF, 
        x = site, 
        y = resids) + 
        labs(x = "", y = "BMI Residuals") + 
        theme_classic() + 
        ggplot2::scale_color_manual(values = c("red", "#0099FF"))
      
      bmiResidsPlot  #Printed at 7 x 5



####Shape Data and Prelim Analyses####
#Load TPS files
Dorsal2 <- readland.tps("STOD_Dorsal2.txt", specID = c("imageID"))
Ventral2 <- readland.tps("STOD_Ventral2.txt", specID = c("imageID"))

#Reorder site factors
dorsalSiteRep2 <- read.csv("Site_Replicate_Dorsal2.csv")
      dorsalSiteRep2$site <- factor(dorsalSiteRep2$site, levels=c("K", "J", "S"))

ventralSiteRep2 <- read.csv("Site_Replicate_Ventral2.csv")
      ventralSiteRep2$site <- factor(ventralSiteRep2$site, levels=c("K", "J", "S"))


#Run Generalized Procrustes analysis

gpaDorsal2 <- gpagen(Dorsal2, ProcD = FALSE, approxBE = TRUE)
gpaVentral2 <- gpagen(Ventral2, ProcD = FALSE, approxBE = TRUE)

plot(gpaDorsal2)
plot(gpaVentral2)

#Run a principal component analysis
pcaDorsal2 <- gm.prcomp(gpaDorsal2$coords)
dorsalPCASum2 <- summary(pcaDorsal2)
pcaVentral2 <- gm.prcomp(gpaVentral2$coords)
ventralPCASum2 <- summary(pcaVentral2)


#Creates geomorph dataframes
dorsalGDF2 <- geomorph.data.frame(shape = gpaDorsal2$coords,
                                  size = log10(gpaDorsal2$Csize),
                                  ind = dorsalSiteRep2$individual,
                                  site = dorsalSiteRep2$site,
                                  replicate = dorsalSiteRep2$replicate,
                                  habitat = dorsalSiteRep2$habitat,
                                  sex = dorsalSiteRep2$sex)

ventralGDF2 <- geomorph.data.frame(shape = gpaVentral2$coords,
                                   size = log10(gpaVentral2$Csize),
                                   ind = ventralSiteRep2$individual,
                                   site = ventralSiteRep2$site,
                                   replicate = ventralSiteRep2$replicate,
                                   habitat = ventralSiteRep2$habitat,
                                   sex = ventralSiteRep2$sex)

#PCA plots to explore shape space
par(mfrow = c(1, 2))

dorsalPCAPlot2 <- plot(pcaDorsal2,
                       cex = 3,
                       pch = 21,
                       bg = as.factor(dorsalGDF2$site),
                       col = c("black"),
                       xlab = paste("Principal Component 1 ", "(", sep = "",
                                    paste(round(dorsalPCASum2$PC.summary$Comp1[2]*100, digits = 2),"%", ")", sep = "")),
                       ylab = paste("Principal Component 2 ", "(", sep = "",
                                    paste(round(dorsalPCASum2$PC.summary$Comp2[2]*100, digits = 2),"%", ")", sep = "")),
                       cex.lab = 1,
                       cex.axis = 0.8)

#text(pcaDorsal2$x[, 2] ~ pcaDorsal2$x[, 1])
     
summary(procD.lm(pcaDorsal2$x ~ dorsalGDF2$site))
summary(pairwise(procD.lm(pcaDorsal2$x ~ dorsalGDF2$site), groups = dorsalGDF2$site))

ventralPCAPlot2 <- plot(pcaVentral2,
                        cex = 3,
                        pch = 21,
                        bg = as.factor(ventralGDF2$site),
                        col = c("black"),
                        xlab = paste("Principal Component 1 ", "(", sep = "",
                                     paste(round(ventralPCASum2$PC.summary$Comp1[2]*100, digits = 2),"%", ")", sep = "")),
                        ylab = paste("Principal Component 2 ", "(", sep = "",
                                     paste(round(ventralPCASum2$PC.summary$Comp2[2]*100, digits = 2),"%", ")", sep = "")),
                        cex.lab = 1,
                        cex.axis = 0.8) #Removing outliers in plastron shape had no impact on asymmetry magnitudes, so they are present in subsequent analyses

summary(procD.lm(pcaVentral2$x ~ ventralGDF2$site))
summary(pairwise(procD.lm(pcaVentral2$x ~ ventralGDF2$site), groups = ventralGDF2$site))


####Asymmetry####

#Dorsal Asymmetry Analyses
#Import "LM Pairs csv file#
LMPairsDorsal <- read.csv("LMPairsSTODDorsal.csv", header = FALSE)

#Look at asymmetry generally for the carapace

dorsalSymm2 <- bilat.symmetry(A = shape, 
                              ind = ind, 
                              rep = rep, 
                              object.sym = TRUE, 
                              land.pairs = LMPairsDorsal, 
                              data = dorsalGDF2, 
                              RRPP = TRUE, 
                              iter = 1000)

summary(dorsalSymm2)

#Look at asymmetry for Joyce, Kenner, and St. Marks separately for the carapace
#Kenner sites are 1:78, Joyce sites are 79:162, St. Marks museum specimens are 163:249
#Kenner
dorsalKennerGDF2 <- geomorph.data.frame(shape = gpaDorsal2$coords[, , 1:78],
                                        size = gpaDorsal2$Csize[1:78],
                                        ind = dorsalSiteRep2$individual[1:78],
                                        site = dorsalSiteRep2$site[1:78],
                                        replicate = dorsalSiteRep2$replicate[1:78],
                                        sex = dorsalSiteRep2$sex[1:78])

#Joyce
dorsalJoyceGDF2 <- geomorph.data.frame(shape = gpaDorsal2$coords[, , 79:162],
                                       size = gpaDorsal2$Csize[79:162],
                                       ind = dorsalSiteRep2$individual[79:162],
                                       site = dorsalSiteRep2$site[79:162],
                                       replicate = dorsalSiteRep2$replicate[79:162],
                                       sex = dorsalSiteRep2$sex[79:162])

#St. Marks
dorsalMarksGDF2 <- geomorph.data.frame(shape = gpaDorsal2$coords[, , 163:249],
                                       size = gpaDorsal2$Csize[163:249],
                                       ind = dorsalSiteRep2$individual[163:249],
                                       site = dorsalSiteRep2$site[163:249],
                                       replicate = dorsalSiteRep2$replicate[163:249],
                                       sex = dorsalSiteRep2$sex[163:249])


#Rerun within site asymmetry analyses
kennerDorsalSymm2 <- bilat.symmetry(A = shape, 
                                    ind = ind, 
                                    rep = rep, 
                                    object.sym = TRUE, 
                                    land.pairs = LMPairsDorsal, 
                                    data = dorsalKennerGDF2, 
                                    RRPP = TRUE, 
                                    iter = 1000)
summary(kennerDorsalSymm2)

#Repeating for Joyce#
joyceDorsalSymm2 <- bilat.symmetry(A = shape, 
                                   ind = ind, 
                                   rep = rep, 
                                   object.sym = TRUE, 
                                   land.pairs = LMPairsDorsal, 
                                   data = dorsalJoyceGDF2, 
                                   RRPP = TRUE, 
                                   iter = 1000)
summary(joyceDorsalSymm2)

#Repeating for St. Marks#
marksDorsalSymm2 <- bilat.symmetry(A = shape, 
                                   ind = ind, 
                                   rep = rep, 
                                   object.sym = TRUE, 
                                   land.pairs = LMPairsDorsal, 
                                   data = dorsalMarksGDF2, 
                                   RRPP = TRUE, 
                                   iter = 1000)
summary(marksDorsalSymm2)

#R square values suggest shell asymmetry explains greater amount of variation
#in Kenner turts than Joyce and St. Marks turts


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

plot(gpaDorsal2$coords[, , 50], 
     pch = 21,
     cex = 3,
     bg = lmSideDorsal$V1)
#Looks good#

#Loop to make the second column (y axis) positive
dorsalGPAFlipped2 <- gpaDorsal2$coords

#Take total length of array and divide by number of landmarks and dimensions
for(i in 1:dim(dorsalGPAFlipped2)[3]){
  dorsalGPAFlipped2[, 2, i] <- abs(dorsalGPAFlipped2[, 2, i])
}

#This will show the x and y coordinates to make sure that it worked properly and 
#paired landmarks are reflected on top of one another
plot(dorsalGPAFlipped2[, 2, 50] ~ dorsalGPAFlipped2[, 1, 50],
     pch = 21,
     cex = 3,
     bg = lmSideDorsal$V1)

text(dorsalGPAFlipped2[, 2, 50] ~ dorsalGPAFlipped2[, 1, 50],
     labels = c(1:dim(dorsalGPAFlipped2)[1]), cex = 1, col = c("black"))

#Calculate Procrustes distances between sides
#Procrustes distance is the sum of distances between corresponding landmarks of two shapes
#The distance between landmarks is sqrt((x2 - x1)^2 + (y2 - y1)^2)

#Make a set of lefts and a set of rights
colnames(LMPairsDorsal) <- c("left", "right") #These are lefts

left <- dorsalGPAFlipped2[LMPairsDorsal$left, , ]
right <- dorsalGPAFlipped2[LMPairsDorsal$right, , ]

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

dorsalSiteRep2$specDist <- specDist

#Which landmark pair was the most asymmetric?
meanLMDist <- matrix(nrow = dim(left)[1], ncol = 1)

for(i in 1:dim(left)[1]){
  meanLMDist[i] <- mean(allLMDists[i, ])
  
}

LMPairsDorsal$mean <- meanLMDist
LMPairsDorsal$mean

#Plot asymmetry magnitudes for Joyce, Kenner and St. Marks dorsal view as a violin plot    

#Average the magnitudes of asymmetry for the three replicates for the plot

dorsalMeans <- dorsalSiteRep2 %>% group_by(individual) %>% summarise(mean = mean(specDist))

dorsalMeans$site <- c(dorsalKennerGDF2$site[1:(length(dorsalKennerGDF2$site)/3)], 
                      dorsalJoyceGDF2$site[1:(length(dorsalJoyceGDF2$site)/3)],
                      dorsalMarksGDF2$site[1:(length(dorsalMarksGDF2$site)/3)])

#plot by site with individual means    
dorsalAsymPlotSite <- ggbetweenstats(
  data = dorsalMeans,
  x = site,
  y = mean) + 
  labs(x = "", y = "Asymmetry Magnitude") + 
  theme_classic()+ 
  ggplot2::scale_color_manual(values = c("#0099FF", "red", "#D4AF37"))

dorsalAsymPlotSite

#Plot asymmetry magnitudes for male and female dorsal view as a violin plot    

#Average the magnitudes of asymmetry for the three replicates for the plot

dorsalMeans$sex <- c(dorsalKennerGDF2$sex[seq(1, length(dorsalKennerGDF2$sex), 3)],
                     dorsalJoyceGDF2$sex[seq(1, length(dorsalJoyceGDF2$sex), 3)],
                     dorsalMarksGDF2$sex[seq(1, length(dorsalMarksGDF2$sex), 3)])


#plot by sex with individual means    
dorsalAsymPlotSex <- ggbetweenstats(
  data = dorsalMeans,
  x = sex,
  y = mean) + 
  labs(x = "", y = "Asymmetry Magnitude") + 
  theme_classic()+ 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

dorsalAsymPlotSex


#method = ML fits the model using maximum likelihood, which is needed for model reduction 
#(not necessary with only one factor; in that case, or when you have arrived at the minimum adequate model, 
#you just leave out the method argument and it fits the model using REML).

#Comparing by site
DorsalAsymmSiteANOVA <- lme(specDist ~ site, 
                            random = ~1|individual, 
                            method = "ML", 
                            data = dorsalSiteRep2)

summary(DorsalAsymmSiteANOVA)

PostSiteComp <- emmeans(DorsalAsymmSiteANOVA, specs = "site")
PostSiteComp
pairs(PostSiteComp, simple = "each") #Pairwise comparisons

#Compare by habitat

DorsalAsymmHabANOVA <- lme(specDist ~ habitat, 
                           random = ~1|individual, 
                           method = "ML", 
                           data = dorsalSiteRep2)

summary(DorsalAsymmHabANOVA)

PostHabComp <- emmeans(DorsalAsymmHabANOVA, specs = "habitat")
PostHabComp


#Ventral Asymmetry Analyses
#Import "LM Pairs csv file#
LMPairsVentral <- read.csv("LMPairsSTODVentral.csv", header = FALSE)

ventralSymm2 <- bilat.symmetry(A = shape, 
                               ind = ind, 
                               replicate = replicate, 
                               object.sym = TRUE, 
                               land.pairs = LMPairsVentral, 
                               data = ventralGDF2, 
                               RRPP = TRUE, 
                               iter = 1000)

summary(ventralSymm2)

#Kenner
ventralKennerGDF2 <- geomorph.data.frame(shape = gpaVentral2$coords[, , 1:102],
                                         size = gpaVentral2$Csize[1:102],
                                         ind = ventralSiteRep2$individual[1:102],
                                         site = as.factor(ventralSiteRep2$site[1:102]),
                                         replicate = ventralSiteRep2$replicate[1:102],
                                         sex = ventralSiteRep2$sex[1:102])

#Joyce
ventralJoyceGDF2 <- geomorph.data.frame(shape = gpaVentral2$coords[, , 103:192],
                                        size = gpaVentral2$Csize[103:192],
                                        ind = ventralSiteRep2$individual[103:192],
                                        site = ventralSiteRep2$site[103:192],
                                        replicate = ventralSiteRep2$replicate[103:192],
                                        sex = ventralSiteRep2$sex[103:192])

#St. Marks
ventralMarksGDF2 <- geomorph.data.frame(shape = gpaVentral2$coords[, , 193:276],
                                        size = gpaVentral2$Csize[193:276],
                                        ind = ventralSiteRep2$individual[193:276],
                                        site = ventralSiteRep2$site[193:276],
                                        replicate = ventralSiteRep2$replicate[193:276],
                                        sex = ventralSiteRep2$sex[193:276])

#Rerun within site asymmetry analyses
kennerVentralSymm2 <- bilat.symmetry(A = shape, 
                                     ind = ind, 
                                     rep = replicate, 
                                     object.sym = TRUE, 
                                     land.pairs = LMPairsVentral, 
                                     data = ventralKennerGDF2, 
                                     RRPP = TRUE, 
                                     iter = 1000)
summary(kennerVentralSymm2)

#Repeating for Joyce#
joyceVentralSymm2 <- bilat.symmetry(A = shape, 
                                    ind = ind, 
                                    rep = replicate, 
                                    object.sym = TRUE, 
                                    land.pairs = LMPairsVentral, 
                                    data = ventralJoyceGDF2, 
                                    RRPP = TRUE, 
                                    iter = 1000)
summary(joyceVentralSymm2)

#Repeating for St. Marks
marksVentralSymm2 <- bilat.symmetry(A = shape, 
                                    ind = ind, 
                                    rep = replicate, 
                                    object.sym = TRUE, 
                                    land.pairs = LMPairsVentral, 
                                    data = ventralMarksGDF2, 
                                    RRPP = TRUE, 
                                    iter = 1000)
summary(marksVentralSymm2)



#R square values suggest shell asymmetry explains greater amount of variation
#in Kenner turts than Joyce turts and St Marks turts#


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

plot(gpaVentral2$coords[, , 50], 
     pch = 21,
     cex = 3,
     bg = lmSideVentral$V1)
#Looks good#

#Loop to make the second column (y axis) positive
ventralGPAFlipped <- gpaVentral2$coords

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

ventralSiteRep2$specDist <- specDistVentral

#Which landmark pair was the most asymmetric?

meanLMDistVentral <- matrix(nrow = dim(leftVentral)[1], ncol = 1)

for(i in 1:dim(leftVentral)[1]){
  meanLMDistVentral[i] <- mean(allLMDistsVentral[i, ])
  
}

LMPairsVentral$mean <- meanLMDistVentral
LMPairsVentral$mean

#Plot asymmetry magnitudes for Joyce and Kenner ventral view as a violin plot    

#Average the three replicates for the plot
ventralSiteRep2

ventralMeans <- ventralSiteRep2 %>% group_by(individual) %>% summarise(mean = mean(specDist))

ventralMeans$site <- c(ventralKennerGDF2$site[1:(length(ventralKennerGDF2$site)/3)], 
                       ventralJoyceGDF2$site[1:(length(ventralJoyceGDF2$site)/3)],
                       ventralMarksGDF2$site[1:(length(ventralMarksGDF2$site)/3)])

ventralMeans$sex <- c(ventralKennerGDF2$sex[seq(1, length(ventralKennerGDF2$sex), 3)],
                      ventralJoyceGDF2$sex[seq(1, length(ventralJoyceGDF2$sex), 3)],
                      ventralMarksGDF2$sex[seq(1, length(ventralMarksGDF2$sex), 3)])


#plot by site with individual means    
ventralAsymPlotSite <- ggbetweenstats(
  data = ventralMeans,
  x = site,
  y = mean) + 
  labs(x = "", y = "Asymmetry Magnitude") + 
  theme_classic()+ 
  ggplot2::scale_color_manual(values = c("#0099FF", "red", "#D4AF37"))

ventralAsymPlotSite

#plot by sex with individual means    
ventralAsymPlotSex <- ggbetweenstats(
  data = ventralMeans,
  x = sex,
  y = mean) + 
  labs(x = "", y = "Plastron Asymmetry Magnitude") + 
  theme_classic() + 
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

ventralAsymPlotSex 

#method = ML fits the model using maximum likelihood, which is needed for model reduction 
#(not necessary with only one factor; in that case, or when you have arrived at the minimum adequate model, 
#you just leave out the method argument and it fits the model using REML).

#Compare by site

#Reorder to compare Kenner as the first variable

VentralAsymmSiteANOVA <- lme(specDist ~ site, 
                             random = ~1|individual, 
                             method = "ML", 
                             data = ventralSiteRep2)

summary(VentralAsymmSiteANOVA)


PostSiteCompVentral <- emmeans(VentralAsymmSiteANOVA, specs = "site")
PostSiteCompVentral
pairs(PostSiteCompVentral, simple = "each")


#Compare by habitat
VentralAsymmANOVA <- lme(specDist ~ habitat, 
                         random = ~1|individual, 
                         method = "ML", 
                         data = ventralSiteRep2)

summary(VentralAsymmANOVA)

PostSiteCompVentral <- emmeans(VentralAsymmANOVA, specs = "habitat")
PostSiteCompVentral


####Potential Stressors####

bodyCond <- read.csv("STOD_Field_Book_Fixed.csv")

leeches <- na.omit(bodyCond[, 1:13]) #Get rid of turtles that didn't have leeches counted

#Standard deviation of leech counts by site
sd(leeches[leeches$Site == "JWMA", ]$Leech.Count)
sd(leeches[leeches$Site == "Kenner", ]$Leech.Count)

#Run t-test between leech count and site
LeechANOVA <- aov(Leech.Count ~ Site, data = leeches)
summary(LeechANOVA)

t.test(Leech.Count ~ Site, data = leeches)

#Make a plot for leech count by site
leechPlot <- ggbetweenstats(
  data = leeches,
  x = Site,
  y = Leech.Count) +
  labs(x = "", y = "Number of Leeches") + 
  theme_classic() +
  ggplot2::scale_color_manual(values = c("red", "#0099FF"))

leechPlot #There are more leeches on turtles at kenner than at Joyce 


#Compare water temp between sites based on our field measurements
t.test(WaterTemp_C ~ Site, data = bodyCond) #Not sig different

#Posthoc regression of leech count and eosinophil counts. 
cellRatios <- read.csv("AllCellCountsFixed.csv")
cellRatios$turtles <- paste(cellRatios$Turtle, "_", cellRatios$Site, sep = "")

leeches$ID <- paste(leeches$Turtle, "_", leeches$Site, sep = "") 
cellRatios$ID <- paste(cellRatios$Turtle, "_", cellRatios$Site, sep = "") 

mergedData <- merge(cellRatios, leeches, by = "ID")

mergedDataJoyce <- mergedData[mergedData$Site.x == "JWMA", ]
mergedDataKenner <- mergedData[mergedData$Site.x == "Kenner", ]

plot(Leech.Count ~ Eosinophil, pch = 21, bg = as.factor(Site.x), data = mergedData)
summary(lm(Leech.Count ~ Eosinophil, data = mergedData))

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


