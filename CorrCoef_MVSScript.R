#---------------
#all packages
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(broom)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(PCAtools)
#--------------
#Importing
DEC <- read_excel("C:/Users/ericr/Desktop/SUNY ESF/PFAS/In progress/DEC2023 PFAS Dataset Project/April 2023/NY PFAS in fish tissue_20230405.xlsx", 
                  sheet = "NY PFAS analytical results", na = "NA")

#--------------
#subsetting
#sticking with whole fish censored data for now, need to determine
#best way to do this after discussing then come back to this 
Data <- DEC[c(7,24,28,31,34,37,40,47,51)]
#if negative, make positive and add column for nondetects
Cen <- ifelse(Data$LAB_RESULT<0, T, F)
Data <- cbind(Data, Cen)
Data$LAB_RESULT <- abs(Data$LAB_RESULT)
#if values are -9 (ie not sampled) they are removed
Data1 <- subset(Data, Data$LAB_RESULT != 9)
#subsetting to just whole fish
Whole <- subset(Data1, FISH_PREP_TYPE_NAME %in% c("Whole", 
                                                  "Synthetic Whole Fish"))
#splitting
PFAS <- split(Whole, Whole$SUBSTANCE_NAME_ABBREVIATION)

#rotating
#remove censored T/F
Data2 <- Data[-10]
#rotate and keep names attached
PFASmatrix <- Data2 %>%
  group_by(SUBSTANCE_NAME_ABBREVIATION) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = SUBSTANCE_NAME_ABBREVIATION, values_from = LAB_RESULT) %>%
  select(-row)

PFASwhole <- subset(PFASmatrix, FISH_PREP_TYPE_NAME %in% c("Whole", "Synthetic Whole Fish"))
PFASwhole <- PFASwhole[c(1:7,11,12,13,14,15,16,19,21,23,30)]
PFASwhole <- na.omit(PFASwhole)

#split by location
PFASwhole1 <- PFASwhole[c(1,4,6:17)]
PFASwholeSplit <- split(PFASwhole1, PFASwhole1$WATERBODY_NAME)
#split by spp
PFASwholeSplitSpp <- split(PFASwhole1, PFASwhole1$COMMON_NAME)
#subset to only numbers
#removes columns 1:2 for all dataframes in the list by magic (?)
PFAStest<- lapply(PFASwholeSplit, "[", 3:14)
PFASspp <- lapply(PFASwholeSplitSpp, "[", 3:14)
#-------------
#plotting
#correlation plots by waterbody
funCorr = function(data, title) {
scale <- scale(data)
cor <- cor(scale)
ggcorrplot(cor,
           type = "lower",
           lab = TRUE,
           ggtheme = ggpubr::theme_classic2(),
           title = title,
           lab_size = 2.5) +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8)) +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.6))
}

#usage - by waterbody
names <- names(PFAStest)
for(i in names){print(funCorr(get(i,PFAStest), i))}
funCorr(PFASwhole1[c(3:14)], "All Locations")

#by species
namesSpp <- names(PFASspp)
for(i in namesSpp){print(funCorr(get(i,PFASspp), i))}

#PCA code - is the data best suited for multivariate analysis?
funPCA = function(data, title){
pca <- prcomp(data)
fviz_pca_biplot(pca, habillage = PFASwhole1$COMMON_NAME,
                  label = "var") +
    theme_classic2() + 
    theme(plot.title=element_text(face="bold", 
                                  size=15, hjust=0.6))
}
#doesn't work, issue with habillage, fix if needed
for(i in names){print(funPCA(get(i,PFAStest), i))}
funPCA(PFASwhole1[c(5:14)], "All")
