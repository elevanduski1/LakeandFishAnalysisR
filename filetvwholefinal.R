#---------------
#all packages
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(broom)
library(car)
library(emmeans)
 #Packages
#--------------
#Importing
DEC <- read_excel("C:/Users/ericr/Desktop/SUNY ESF/PFAS/In progress/DEC2023 PFAS Dataset Project/April 2023/NY PFAS in fish tissue_20230405.xlsx", 
                  +     sheet = "NY PFAS analytical results", 
                  +     na = "NA")
 #Importing
#------------- 
#subsetting
Data <- DEC[c(7,24,28,31,34,37,40,47,51)]
#if negative now NA to remove nondetects
Data$LAB_RESULT <- replace(Data$LAB_RESULT, which(Data$LAB_RESULT < 0), NA)
Data <- Data[!(is.na(Data$LAB_RESULT)), ]

#Subsetting by only synthetic whole fish and their filets, then combining
Whole <- subset(Data, FISH_PREP_TYPE_NAME == "Synthetic Whole Fish")
table(Whole$SUBSTANCE_NAME_ABBREVIATION)
Filet <- subset(Data, FISH_PREP_TYPE_NAME %in% c("Skin-Off Fillets", 
                                                       "Skin-On Fillets/Scaleless"))
SHFvFilet <- merge(Whole, Filet, by = c("SUBSTANCE_NAME_ABBREVIATION","FIELD_SAMPLE_NO"))
table(SHFvFilet$SUBSTANCE_NAME_ABBREVIATION)

#------------
#plotting
#Function for input of each species as the dataset and creates scatterplot
#can add or remove log around the aes of x and y to create log-log plots or not
funlogSpp = function(data, title) {
  ggplot(data, aes(x = log(LAB_RESULT.y), y = log(LAB_RESULT.x))) +
    geom_point(aes(color = WATERBODY_NAME.x)) +
    theme_classic2() +
    labs(x = "Filet Concentration (ppb)",
         y = "Whole Fish Concentration (ppb)",
         title = title,
         color = "Waterbody Name") +
    geom_smooth(method = "lm", color = "grey") +
    theme(axis.text.y=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size=10),
          legend.title = element_text(size=10),
          plot.title=element_text(face="bold", 
                                  size=15, 
                                  hjust=0.6)) +
    guides(color = guide_legend(ncol = 1)) +
    stat_poly_eq(use_label(c("eq", "adj.R2", "p", "n")))
}

#Split data from all filet to whole fish comparisons by common name

splitPFAS <- split(SHFvFilet, SHFvFilet$SUBSTANCE_NAME_ABBREVIATION)
PFOS <- splitPFAS$PFOS
splitPFOS <- split(PFOS, PFOS$COMMON_NAME.x)
PFDA <- splitPFAS$PFDA
splitPFDA <- split(PFDA, PFDA$COMMON_NAME.x)
PFUnA <- splitPFAS$PFUnA
splitPFUnA <- split(PFUnA, PFUnA$COMMON_NAME.x)
PFDoA <- splitPFAS$PFDoA
splitPFDoA <- split(PFDoA, PFDoA$COMMON_NAME.x)

#now that it is split can use individual split data frames made to create all the graphs
#creates all the species graphs at once
#PFOS
names <- names(splitPFOS)
for(i in names){print(funSpp(get(i,splitPFOS), i))}
#PFDA
namesPFDA <- names(splitPFDA)
for(i in namesPFDA){print(funSpp(get(i,splitPFDA), i))}
#PFUnA
namesPFUnA <- names(splitPFUnA)
for(i in namesPFUnA){print(funSpp(get(i,splitPFUnA), i))}
#PFDoA
namesPFDoA <- names(splitPFDoA)
for(i in namesPFDoA){print(funSpp(get(i,splitPFDoA), i))}

#can use the function by itself for individual graphs based on different subsets
SHF1 <- subset(PFOS, LAB_RESULT.y < 500)
funSpp(PFOS, "All Species")
funSpp(PFDA, "Log-Log All Species PFDA")
funSpp(PFUnA, "Log-Log All Species PFUnA")
funSpp(PFDoA, "All Species PFDoA")
#this is when funlogSpp was made, too little to late
funlogSpp(PFDoA, "Log-Log All Species PFDoA")

#Same function as before but separation by water body and points are by species
funlogWater = function(data, title) {
  ggplot(data, aes(x = log(LAB_RESULT.y), y = log(LAB_RESULT.x))) +
    geom_point(aes(color = COMMON_NAME.x)) +
    theme_classic2() +
    labs(x = "Filet Concentration (ppb)",
         y = "Whole Fish Concentration (ppb)",
         title = title,
         color = "Species") +
    geom_smooth(method = "lm", color = "grey") +
    theme(axis.text.y=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size=10),
          legend.title = element_text(size=10),
          plot.title=element_text(face="bold", 
                                  size=15, 
                                  hjust=0.6)) +
    guides(color = guide_legend(ncol = 1)) +
    stat_poly_eq(use_label(c("eq", "adj.R2", "p", "n")))
}

#Same as above but with water bodies representing points

#PFOS
splitPFOS.WB <- split(PFOS, PFOS$WATERBODY_NAME.x)
namesPFOS.WB <- names(splitPFOS.WB)
for(i in namesPFOS.WB){print(funWater(get(i,splitPFOS.WB), i))}

#PFDA
splitPFDA.WB <- split(PFDA, PFDA$WATERBODY_NAME.x)
namesPFDA.WB <- names(splitPFDA.WB)
for(i in namesPFDA.WB){print(funWater(get(i,splitPFDA.WB), i))}

#PFUnA
splitPFUnA.WB <- split(PFUnA, PFUnA$WATERBODY_NAME.x)
namesPFUnA.WB <- names(splitPFUnA.WB)
for(i in namesPFUnA.WB){print(funWater(get(i,splitPFUnA.WB), i))}

#PFDoA
splitPFDoA.WB <- split(PFDoA, PFDoA$WATERBODY_NAME.x)
namesPFDoA.WB <- names(splitPFDoA.WB)
for(i in namesPFDoA.WB){print(funWater(get(i,splitPFDoA.WB), i))}

#function by itself for all and logs

#PFOS
funWater(PFOS, "All PFOS")
funlogWater(PFOS, "Log-Log All PFOS")

#PFDA
funWater(PFDA, "All PFDA")
funlogWater(PFDA, "Log-Log All PFDA")

#PFUnA
funWater(PFUnA, "All PFUnA")
funlogWater(PFUnA, "Log-Log All PFUnA")

#PFDoA
funWater(PFDoA, "All PFDoA")
funlogWater(PFDoA, "Log-Log All PFDoA")

#----------
#Get all linear models, overwrote code for each lm of spp, water, and both

#PFOS
PFOSlmWaterbody <- PFOS %>% nest(data = c(-WATERBODY_NAME.x, -COMMON_NAME.x)) %>% 
  mutate(model = map(data, ~lm(LAB_RESULT.x~LAB_RESULT.y, data = .)), 
         tidied = map(model, tidy)) %>% 
  unnest(tidied)
PFOSlm3 <- PFOSlmWaterbody[-c(3,4)]
write.csv(PFOSlm3, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\PFOSlmBoth.csv")

#PFDA
PFDAlm1 <- PFDA %>% nest(data = c(-WATERBODY_NAME.x)) %>% 
  mutate(model = map(data, ~lm(LAB_RESULT.x~LAB_RESULT.y, data = .)), 
         tidied = map(model, tidy)) %>% 
  unnest(tidied)
PFDAlm3 <- PFDAlm1[-c(2,3)]
write.csv(PFDAlm3, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\PFDAlmwater.csv")

#PFUnA
PFUnAlm1 <- PFUnA %>% nest(data = c(-WATERBODY_NAME.x, -COMMON_NAME.x)) %>% 
  mutate(model = map(data, ~lm(LAB_RESULT.x~LAB_RESULT.y, data = .)), 
         tidied = map(model, tidy)) %>% 
  unnest(tidied)
PFUnAlm3 <- PFUnAlm1[-c(3,4)]
write.csv(PFUnAlm3, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\PFUnAlmboth.csv")

#PFDoA
PFDoAlm1 <- PFDoA %>% nest(data = c(-COMMON_NAME.x)) %>% 
  mutate(model = map(data, ~lm(LAB_RESULT.x~LAB_RESULT.y, data = .)), 
         tidied = map(model, tidy)) %>% 
  unnest(tidied)
PFDoAlm3 <- PFDoAlm1[-c(2,3)]
write.csv(PFDoAlm3, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\PFDoAlmspp.csv")

#-----------
#compare linear models with car::Anova

#by spp
Spplm1 <- lm(PFOSforlm$LAB_RESULT.x~PFOSforlm$FISH_PREP_TYPE_NAME.x + PFOSforlm$COMMON_NAME.x)
Sppaov1 <- car::Anova(Spplm1, type = 3)
Sppaov1
write.csv(Sppaov1, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\Sppaov1.csv")
summary(car::Anova(Spplm1, type = 3))

#by water body
waterlm1 <- lm(PFOSforlm$LAB_RESULT.x~PFOSforlm$FISH_PREP_TYPE_NAME.x + PFOSforlm$WATERBODY_NAME.x)
wateraov1 <- car::Anova(waterlm1, type = 3)
wateraov1
write.csv(wateraov1, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\wateraov1.csv")
summary(car::Anova(waterlm1, type = 3))

#there is a significant difference by species
#and by water body

#post-hoc
Spp.emm <- emmeans(Spplm1, ~ COMMON_NAME.x)
Spp.emm
spp.contrast <- contrast(Spp.emm, show.ests = TRUE)
write.csv(spp.contrast, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\sppContrast.csv")

#the ones with p< 0.05 are different lm, while those above are the same?
water.emm <- emmeans(waterlm1, ~ WATERBODY_NAME.x)
water.emm
water.contrast<- contrast(water.emm, show.ests = TRUE)
write.csv(water.contrast, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\WaterContrast.csv")

#most are significantly different, so unifying them wouldn't work?