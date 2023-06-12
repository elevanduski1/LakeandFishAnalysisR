#---------------
#all packages
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(broom)

#--------------
#Importing
DEC <- read_excel("C:/Users/ericr/Desktop/SUNY ESF/PFAS/In progress/DEC2023 PFAS Dataset Project/April 2023/NY PFAS in fish tissue_20230405.xlsx", 
                  sheet = "NY PFAS analytical results", na = "NA")

#--------------
#subsetting
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
unique(Whole$COMMON_NAME) #46 total
unique(Data$COMMON_NAME) #56 total, lose 10 spp doing only Whole fish

#splitting
PFAS <- split(Whole, Whole$SUBSTANCE_NAME_ABBREVIATION)
PFOS <- PFAS$PFOS
PFOSspp <- split(PFOS, PFOS$COMMON_NAME)
PFOSwater <- split(PFOS, PFOS$WATERBODY_NAME)

#making the means and sd shown for pfos plot
#use mean for pfos only then use median for the rest?
MeansSDPFOS <- PFOS %>%
  group_by(COMMON_NAME) %>%
  summarize(mean=mean(LAB_RESULT, na.rm=TRUE), 
            sd=sd(LAB_RESULT, na.rm=TRUE))
MeansSDPFOS

#--------
#plotting
#plots for each pfas by species across all locations
funBar = function(data, title) {
  ggplot(data, aes(x = reorder(COMMON_NAME, -mean), y = mean)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(min = mean-sd, max = mean+sd, width=0.05)) +
    labs(x = "Species",
         y = "Concentration (ppb)",
         title = title) + 
    theme_classic2() +
    theme(axis.text.y=element_text(size=15),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, 
                                     color="Black", face="bold", size = 15),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          plot.title=element_text(face="bold", 
                                  size=15, 
                                  hjust=0.5),
          legend.position = "none")
}
#usage
funBar(MeansSDPFOS, "PFOS")

#plots each pfas group by species in boxplots
funBoxplot = function(data, title) {
  ggplot(data, aes(x = reorder(COMMON_NAME, LAB_RESULT, 
                               FUN=median, decreasing = TRUE), 
                   y = LAB_RESULT, fill = COMMON_NAME)) +
    geom_boxplot() +
    theme_classic2() +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, 
                                     color="Black", face="bold", size = 15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          plot.title=element_text(face="bold", 
                                  size=15, 
                                  hjust=0.5)) +
    labs(title = title,
         x="Species", 
         y= "Whole Fish Concentration (ppb)")
}

#plots each waterbody for 1 PFAS group by species
funWaterBox = function(data, title) {
  ggplot(data, aes(x = reorder(COMMON_NAME, LAB_RESULT, 
                               FUN=median, decreasing = TRUE), 
                   y = LAB_RESULT, fill = COMMON_NAME)) +
    geom_boxplot() +
    theme_classic2() +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, 
                                     color="Black", face="bold", size = 15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          plot.title=element_text(face="bold", 
                                  size=15, 
                                  hjust=0.5)) +
    labs(title = title,
         x="Species", 
         y=  "Whole Fish Concentration (ppb)")
}

#plots each pfas group by family
funFamBox = function(data, title) {
  ggplot(data, aes(x = reorder(FAMILY, LAB_RESULT, 
                               FUN=median, decreasing = TRUE), 
                   y = LAB_RESULT, fill = FAMILY)) +
    geom_boxplot() +
    #stat_summary(fun.data = give.n, geom = "text") +
    theme_classic2() +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, 
                                     color="Black", face="bold", size = 15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          plot.title=element_text(face="bold", 
                                  size=15, 
                                  hjust=0.5)) +
    labs(title = title,
         x="Family", 
         y= "Whole Fish Concentration (ppb)")
}

#determines n value for each
#give.n <- function(x){
#  return(c(y = mean(x), label = length(x)))
#}

#usage
#boxplots by pfas group
names <- names(PFAS)
for(i in names){print(funBoxplot(get(i,PFAS), i))}

funBoxplot(PFOS, "PFOS")
#boxplots by waterbody for one PFAS group 
namesPFOS <- names(PFOSwater)
for(i in namesPFOS){print(funWaterBox(get(i,PFOSwater), i))}

#boxplots by pfas group for taxa groups
for(i in names){print(funFamBox(get(i,PFAS), i))}
funFamBox(PFAS$`11Cl-PF3OUdS`, "11Cl-PF3OUdS")
