#------------
#Packages
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
#-----------
#Importing
DEC <- read_excel("C:/Users/ericr/Desktop/SUNY ESF/PFAS/In progress/DEC2023 PFAS Dataset Project/April 2023/NY PFAS in fish tissue_20230405.xlsx", 
                  sheet = "NY PFAS analytical results", 
                  na = "NA")
#------------
#subsetting
Data <- DEC[c(7,24,28,31,34,37,40,47,51)]
#if negative now NA to remove nondetects
Data$LAB_RESULT <- replace(Data$LAB_RESULT, which(Data$LAB_RESULT < 0), NA)
Data <- Data[!(is.na(Data$LAB_RESULT)), ]

#subset to just whole and synthetic whole, makes sense for comparing to length, 
#but need to confirm this is appropriate
Whole <- subset(Data, FISH_PREP_TYPE_NAME %in% c("Synthetic Whole Fish", "Whole"))

#subset to individual PFAS groups that are large enough to sample
table(Whole$SUBSTANCE_NAME_ABBREVIATION)
table(Whole$COMMON_NAME)
split <- split(Whole, Whole$SUBSTANCE_NAME_ABBREVIATION)
PFOS <- split$PFOS
PFOS1 <- split(PFOS, PFOS$COMMON_NAME)
#---------------
#plotting

funlogplot2 = function(data, title) {
  ggplot(data, aes(x = log(`LENGTH_MEAN (MM)`), y = log(LAB_RESULT))) +
    geom_point(aes(color = WATERBODY_NAME)) +
    theme_classic2() +
    labs(x = "Log - Length (mm)",
         y = "Log - Whole Fish Concentration (ppb)",
         title = title,
         color = "Waterbody Name") +
    geom_smooth(method = "lm", color = "grey") +
    theme(axis.text.y=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size=7),
          legend.title = element_text(size=10),
          plot.title=element_text(face="bold", 
                                  size=15, 
                                  hjust=0.6)) +
    guides(color = guide_legend(ncol = 1)) +
    stat_poly_eq(use_label(c("eq", "adj.R2", "p", "n")))
}
#xxx#2 is by waterbody instead of by species
#each individual PFAS group run for length to PFAS concentration

names <- names(split)
for(i in names){print(funlogPlot(get(i,split), i))}
#PFDA, PFDoA, PFHxS, PFNA, PFOA, PFOS, PFOSA, PFUnA

namesPFOS <- names(PFOS1)
for(i in namesPFOS){print(funlogplot2(get(i,PFOS1), i))}
