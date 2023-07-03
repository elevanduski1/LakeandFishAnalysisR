#---------------
#all packages
library(writexl)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(broom)
library(car)
library(emmeans)
library(lme4)
library(lmerTest)
library(sjPlot)
#--------------
#Importing
DEC <- read_excel("C:/Users/ericr/Desktop/SUNY ESF/PFAS/In progress/DEC2023 PFAS Dataset Project/April 2023/NY PFAS in fish tissue_20230405.xlsx", 
                  +     sheet = "NY PFAS analytical results", 
                  +     na = "NA")
#------------- 
#subsetting
Data <- DEC[c(7,11,12,24,28,31,34,37,40,47,51, 57)]

#if negative now NA to remove nondetects
Data$LAB_RESULT <- replace(Data$LAB_RESULT, which(Data$LAB_RESULT < 0), NA)
Data <- Data[!(is.na(Data$LAB_RESULT)), ]

#determine sites purpose
Reference <- subset(Data, Data$SAMPLE_SELECTION_PURPOSE == "REF")
table(Reference$WATERBODY_NAME)

Potential <- subset(Data, Data$SAMPLE_SELECTION_PURPOSE == "POT")
table(Potential$WATERBODY_NAME)

#Subsetting by only synthetic whole fish and their filets, then combining
Whole <- subset(Data, FISH_PREP_TYPE_NAME %in% c("Synthetic Whole Fish", "Whole"))
table(Whole$SUBSTANCE_NAME_ABBREVIATION)
PFOS <- subset(Whole, Whole$SUBSTANCE_NAME_ABBREVIATION == "PFOS")
PFOS <- subset(PFOS, PFOS$LAT_DECIMAL_DEGREES != "NA")
write.csv(PFOS, "C:\\Users\\ericr\\PFOS.csv")
Filet <- subset(Data, FISH_PREP_TYPE_NAME %in% c("Skin-Off Fillets", 
                                                       "Skin-On Fillets/Scaleless"))
SHFvFilet <- merge(Whole, Filet, by = c("SUBSTANCE_NAME_ABBREVIATION","FIELD_SAMPLE_NO"))
table(SHFvFilet$SUBSTANCE_NAME_ABBREVIATION)

#------------
#plotting
#Function for input of each species as the dataset and creates scatterplot
#can add or remove log around the aes of x and y to create log-log plots or not
funSpp = function(data, title) {
  ggplot(data, aes(x = LAB_RESULT.y, y = LAB_RESULT.x)) +
    geom_point(aes(color = COMMON_NAME.x)) +
    theme_classic2() +
    labs(x = "Filet Concentration (ppb)",
         y = "Whole Fish Concentration (ppb)",
         title = title,
         color = "Waterbody Name") +
    geom_smooth(method = "lm", color = "grey", aes()) +
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
funPFASbySpp = function(data, title) {
  ggplot(data, aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                   color = COMMON_NAME.x)) +
    geom_point() +
    theme_classic2() +
    labs(x = "Filet Concentration (ppb)",
         y = "Whole Fish Concentration (ppb)",
         title = title,
         color = "Species",
         fill = "Species") +
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
    guides(color = guide_legend(ncol = 1))
}

funSppbyWB = function(data, title) {
  ggplot(data, aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = WATERBODY_NAME.x, 
                   color = WATERBODY_NAME.x)) +
    geom_point() +
    theme_classic2() +
    labs(x = "Filet Concentration (ppb)",
         y = "Whole Fish Concentration (ppb)",
         title = title,
         color = "Waterbody Name",
         fill = "Waterbody Name") +
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
#Each PFAS by species and lm individually by waterbody
funSppbyWB(splitPFOS$`American Eel`, "American Eel")
for(i in names){print(funSppbyWB(get(i,splitPFOS), i))}
for(i in namesPFDA){print(funSppbyWB(get(i,splitPFDA), i))}
for(i in namesPFDoA){print(funSppbyWB(get(i,splitPFDoA), i))}
for(i in namesPFUnA){print(funSppbyWB(get(i,splitPFUnA), i))}
for(i in namesPFOA){print(funSppbyWB(get(i,splitPFOA), i))}
for(i in namesPFNA){print(funSppbyWB(get(i,splitPFNA), i))}
for(i in namesPFHxS){print(funSppbyWB(get(i,splitPFHxS), i))}
for(i in namesPFOSA){print(funSppbyWB(get(i,splitPFOSA), i))}

#barplot code for individual species separated by waterbody
#will include code to integrate reference and known sites

#if ksi or ksr == known, elseif ref == reference, else == name
Data$Usage <- ifelse(Data$SAMPLE_SELECTION_PURPOSE == 'KSI', 'Known',
                    ifelse(Data$SAMPLE_SELECTION_PURPOSE == 'KSR', 'Known',
                           ifelse(Data$SAMPLE_SELECTION_PURPOSE == 'REF', 'Reference', Data$WATERBODY_NAME)))
PFOSBar <- subset(Data, Data$SUBSTANCE_NAME_ABBREVIATION == "PFOS")
#used only synthetic whole fish for this, as they will be the same samples as the filet
#we are comparing them to
PFOSBar <- subset(PFOSBar, PFOSBar$FISH_PREP_TYPE_NAME == "Whole")
PFOSBarSplit <- split(PFOSBar, PFOSBar$COMMON_NAME)
namesPFOS <- names(PFOSBarSplit)
funBar = function(data, title) {
  data1 <- data %>%
    group_by(Usage) %>%
    summarize(mean=mean(LAB_RESULT, na.rm=TRUE), 
              sd=sd(LAB_RESULT, na.rm=TRUE))
  
  plot <- ggplot(data1, aes(x = reorder(Usage, -mean), y = mean,
                            fill=factor(ifelse(Usage == "Known","Highlighted",
                                               ifelse(Usage == "Reference", "Highlighted", "Normal"))))) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(min = mean-sd, max = mean+sd, width=0.05)) +
    labs(x = "Waterbody Name",
         y = "Whole Fish Concentration (ppb)",
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
for(i in namesPFOS){print(funBar(get(i,PFOSBarSplit), i))}
#-------------
#all species colored points and lm all in one
funPFASbySpp(SHF1, "PFOS")
ForAllLmPFOS <- subset(SHFvFilet, SHFvFilet$SUBSTANCE_NAME_ABBREVIATION == "PFOS")
funPFASbySpp(ForAllLmPFOS, "PFOS")
ForAllLmPFDA <- subset(SHFvFilet, SHFvFilet$SUBSTANCE_NAME_ABBREVIATION == "PFDA")
funPFASbySpp(ForAllLmPFDA, "PFDA")
ForAllLmPFUnA <- subset(SHFvFilet, SHFvFilet$SUBSTANCE_NAME_ABBREVIATION == "PFUnA")
funPFASbySpp(ForAllLmPFUnA, "PFUnA")
ForAllLmPFDoA <- subset(SHFvFilet, SHFvFilet == "PFDoA")
funPFASbySpp(ForAllLmPFDoA, "PFDoA")
#---------------
#same as above but limiting to showing lm when n>= 5, but keeps showing all 
#points regardless of the n value
#PFOS
ForAllLmPFOSn5 <- ForAllLmPFOS %>% 
  group_by(COMMON_NAME.x) %>% 
  filter(n() >= 5)
#put count for each species 
CountPFOS <- ForAllLmPFOS %>%
  group_by(COMMON_NAME.x) %>%
  summarise(Count = n())

#plot
ggplot(ForAllLmPFOS, aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                 color = COMMON_NAME.x)) +
  geom_point() +
  theme_classic2() +
  labs(x = "Filet Concentration (ppb)",
       y = "Whole Fish Concentration (ppb)",
       title = "PFOS",
       fill = "Species and Count",
       caption = "Linear regression shown for Species with n ≥ 5") +
  geom_smooth(ForAllLmPFOSn5, mapping = aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                                  color = COMMON_NAME.x),
              method = "lm", color = "grey", inherit.aes = F) +
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
  guides(fill = guide_legend(ncol = 1), color = "none") +
  scale_fill_discrete(labels = paste(CountPFOS$COMMON_NAME.x, CountPFOS$Count))

#PFDA
ForAllLmPFDAn5 <- ForAllLmPFDA %>% 
  group_by(COMMON_NAME.x) %>% 
  filter(n() >= 5)
#put count for each species 
CountPFDA <- ForAllLmPFDA %>%
  group_by(COMMON_NAME.x) %>%
  summarise(Count = n())

#plot
ggplot(ForAllLmPFDA, aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                         color = COMMON_NAME.x)) +
  geom_point() +
  theme_classic2() +
  labs(x = "Filet Concentration (ppb)",
       y = "Whole Fish Concentration (ppb)",
       title = "PFDA",
       fill = "Species and Count",
       caption = "Linear regression shown for Species with n ≥ 5") +
  geom_smooth(ForAllLmPFDAn5, mapping = aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                                            color = COMMON_NAME.x),
              method = "lm", color = "grey", inherit.aes = F) +
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
  guides(fill = guide_legend(ncol = 1), color = "none") +
  scale_fill_discrete(labels = paste(CountPFDA$COMMON_NAME.x, CountPFDA$Count))

#PFDoA
ForAllLmPFDoAn5 <- ForAllLmPFDoA %>% 
  group_by(COMMON_NAME.x) %>% 
  filter(n() >= 5)
#put count for each species 
CountPFDoA <- ForAllLmPFDoA %>%
  group_by(COMMON_NAME.x) %>%
  summarise(Count = n())

#plot
ggplot(ForAllLmPFDoA, aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                         color = COMMON_NAME.x)) +
  geom_point() +
  theme_classic2() +
  labs(x = "Filet Concentration (ppb)",
       y = "Whole Fish Concentration (ppb)",
       title = "PFDoA",
       fill = "Species and Count",
       caption = "Linear regression shown for Species with n ≥ 5") +
  geom_smooth(ForAllLmPFDoAn5, mapping = aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                                            color = COMMON_NAME.x),
              method = "lm", color = "grey", inherit.aes = F) +
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
  guides(fill = guide_legend(ncol = 1), color = "none") +
  scale_fill_discrete(labels = paste(CountPFDoA$COMMON_NAME.x, CountPFDoA$Count))

#PFUnA
ForAllLmPFUnAn5 <- ForAllLmPFUnA %>% 
  group_by(COMMON_NAME.x) %>% 
  filter(n() >= 5)
#put count for each species 
CountPFUnA <- ForAllLmPFUnA %>%
  group_by(COMMON_NAME.x) %>%
  summarise(Count = n())

#plot
ggplot(ForAllLmPFUnA, aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                          color = COMMON_NAME.x)) +
  geom_point() +
  theme_classic2() +
  labs(x = "Filet Concentration (ppb)",
       y = "Whole Fish Concentration (ppb)",
       title = "PFUnA",
       fill = "Species and Count",
       caption = "Linear regression shown for Species with n ≥ 5") +
  geom_smooth(ForAllLmPFUnAn5, mapping = aes(x = LAB_RESULT.y, y = LAB_RESULT.x, fill = COMMON_NAME.x, 
                                             color = COMMON_NAME.x),
              method = "lm", color = "grey", inherit.aes = F) +
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
  guides(fill = guide_legend(ncol = 1), color = "none") +
  scale_fill_discrete(labels = paste(CountPFUnA$COMMON_NAME.x, CountPFUnA$Count))
#--------------
splitPFAS <- split(SHFvFilet, SHFvFilet$SUBSTANCE_NAME_ABBREVIATION)
PFOS <- splitPFAS$PFOS
splitPFOS <- split(PFOS, PFOS$COMMON_NAME.x)
PFDA <- splitPFAS$PFDA
splitPFDA <- split(PFDA, PFDA$COMMON_NAME.x)
PFUnA <- splitPFAS$PFUnA
splitPFUnA <- split(PFUnA, PFUnA$COMMON_NAME.x)
PFDoA <- splitPFAS$PFDoA
splitPFDoA <- split(PFDoA, PFDoA$COMMON_NAME.x)
PFOA <- splitPFAS$PFOA
splitPFOA <- split(PFOA, PFOA$COMMON_NAME.x)
PFNA <- splitPFAS$PFNA
splitPFNA <- split(PFNA, PFNA$COMMON_NAME.x)
PFHxS <- splitPFAS$PFHxS
splitPFHxS <- split(PFHxS, PFHxS$COMMON_NAME.x)
PFOSA <- splitPFAS$PFOSA
splitPFOSA <- split(PFOSA, PFOSA$COMMON_NAME.x)
PFPeA <- splitPFAS$PFPeA
splitPFPeA <- split(PFPeA, PFPeA$COMMON_NAME.x)
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
#PFOA
namesPFOA <- names(splitPFOA)
for(i in namesPFOA){print(funSpp(get(i,splitPFOA), i))}
#PFNA
namesPFNA <- names(splitPFNA)
for(i in namesPFNA){print(funSpp(get(i,splitPFNA), i))}
#PFHxS
namesPFHxS<- names(splitPFHxS)
for(i in namesPFHxS){print(funSpp(get(i,splitPFHxS), i))}
#PFOSA
namesPFOSA <- names(splitPFOSA)
for(i in namesPFOSA){print(funSpp(get(i,splitPFOSA), i))}
#PFPeA
namesPFPeA <- names(splitPFPeA)
for(i in namesPFPeA){print(funSpp(get(i,splitPFPeA), i))}

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
#NOT THE BEST WAY TO DO SO, KEEPING CODE FOR LEGACY IN CASE THIS COMES UP AGAIN
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
summary(spp.contrast)
write.csv(spp.contrast, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\sppContrast.csv")

#the ones with p< 0.05 are different lm, while those above are the same?
water.emm <- emmeans(waterlm1, ~ WATERBODY_NAME.x)
water.emm
water.contrast<- contrast(water.emm, show.ests = TRUE)
write.csv(water.contrast, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\WaterContrast.csv")

#most are significantly different, so unifying them wouldn't work?
#-------------
#lme4 linear mixed effect modelling 
#PFOS
#rise/run so whole (y axis) / filet (x axis)
#do emmeans on this model as post hoc testing
Whole.Filet <- SHF1$LAB_RESULT.x/SHF1$LAB_RESULT.y
PFOSlme4 <- cbind(SHF1, Whole.Filet)
PFOSlme4$COMMON_NAME.x <- as.factor(PFOSlme4$COMMON_NAME.x)
PFOSlme <- lmer(Whole.Filet ~ COMMON_NAME.x + (1|WATERBODY_NAME.x), 
                data = PFOSlme4)
sumPFOSlme <- summary(PFOSlme)
sumPFOSlme
#----------
#make column names just the names of the groups, not including the column
csumm <- coef(sumPFOSlme)
rownames(csumm) <- sub("^COMMON_NAME.x", "", rownames(csumm))
modelPFOS <-print(csumm[-1,], digits=4)
namesPFOSlme <- rownames(modelPFOS)
namesPFOSlme <- append(namesPFOSlme,'(Intercept)',after=0)
#output table model and export
tab_model(PFOSlme, show.stat = T, show.ngroups = T, show.obs = T, 
          pred.labels = namesPFOSlme, 
          file = "PFOSlme.pdf")
#----------
#post hoc
Spp.emm <- emmeans(PFOSlme, ~ COMMON_NAME.x)
Spp.emm
spp.contrast <- contrast(Spp.emm, show.ests = TRUE)
summary(spp.contrast)
write.csv(Spp.emm, "C:\\Users\\ericr\\emmeans.csv")
write.csv(summary(spp.contrast), "C:\\Users\\ericr\\contrast.csv")
knitr::kable(summary(spp.contrast), format = "html")
#lmertest is useful as it provides lmer output and p values, 
#but giving p values for mixed effects models is currently under discussion 
#as to whether it is appropriate, so useful but not the best?
summary(lmerTest::lmer(Whole.Filet ~ COMMON_NAME.x + (1|WATERBODY_NAME.x), 
               data = PFOSlme4))

#low variance by waterbody, means it explains little of the variance in the data?

#PFDA
Whole.FiletPFDA <- ForAllLmPFDA$LAB_RESULT.x/ForAllLmPFDA$LAB_RESULT.y
PFDAlme4 <- cbind(ForAllLmPFDA, Whole.FiletPFDA)
PFDAlme <- lmer(Whole.FiletPFDA ~ COMMON_NAME.x + (1|WATERBODY_NAME.x), 
                data = PFDAlme4)
summary(PFDAlme)
summary(lmerTest::lmer(Whole.FiletPFDA ~ COMMON_NAME.x + (1|WATERBODY_NAME.x), 
                       data = PFDAlme4))
PFDAlmerTest <- lmerTest::lmer(Whole.FiletPFDA ~ COMMON_NAME.x + (1|WATERBODY_NAME.x), 
                               data = PFDAlme4)
tab_model(PFDAlme, show.stat = T, show.ngroups = T, show.obs = T, 
          file = "PFOSlme.doc")
#post hoc
Spp.emmPFDA <- emmeans(PFDAlme, ~ COMMON_NAME.x)
Spp.emmPFDA
spp.contrastPFDA <- contrast(Spp.emmPFDA, show.ests = TRUE)
summary(spp.contrastPFDA)
write.csv(Spp.emmPFDA, "C:\\Users\\ericr\\emmeansPFDA.csv")
write.csv(summary(spp.contrastPFDA), "C:\\Users\\ericr\\contrastPFDA.csv")

#PFDoA
Whole.FiletPFDoA <- ForAllLmPFDoA$LAB_RESULT.x/ForAllLmPFDoA$LAB_RESULT.y
PFDoAlme4 <- cbind(ForAllLmPFDoA, Whole.FiletPFDoA)
PFDoAlme <- lmer(Whole.FiletPFDoA ~ COMMON_NAME.x + (1|WATERBODY_NAME.x), 
                data = PFDoAlme4)
summary(PFDoAlme)
PFDoAlme
tab_model(PFDoAlme, show.stat = T, show.ngroups = T, show.obs = T, 
          file = "PFDoAlme.doc")
#posthoc
Spp.emmPFDoA <- emmeans(PFDoAlme, ~ COMMON_NAME.x)
Spp.emmPFDoA
spp.contrastPFDoA <- contrast(Spp.emmPFDoA, show.ests = TRUE)
summary(spp.contrastPFDoA)
write.csv(Spp.emmPFDoA, "C:\\Users\\ericr\\emmeansPFDoA.csv")
write.csv(summary(spp.contrastPFDoA), "C:\\Users\\ericr\\contrastPFDoA.csv")

#PFUnA
Whole.FiletPFUnA <- ForAllLmPFUnA$LAB_RESULT.x/ForAllLmPFUnA$LAB_RESULT.y
PFUnAlme4 <- cbind(ForAllLmPFUnA, Whole.FiletPFUnA)
PFUnAlme <- lmer(Whole.FiletPFUnA ~ COMMON_NAME.x + (1|WATERBODY_NAME.x), 
                 data = PFUnAlme4)
tab_model(PFUnAlme, show.stat = T, show.ngroups = T, show.obs = T, 
          file = "PFUnAlme.doc")

#posthoc
Spp.emmPFUnA <- emmeans(PFUnAlme, ~ COMMON_NAME.x)
Spp.emmPFUnA
spp.contrastPFUnA <- contrast(Spp.emmPFUnA, show.ests = TRUE)
summary(spp.contrastPFUnA)
write.csv(Spp.emmPFUnA, "C:\\Users\\ericr\\emmeansPFUnA.csv")
write.csv(summary(spp.contrastPFUnA), "C:\\Users\\ericr\\contrastPFUnA.csv")
