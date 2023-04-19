#packages
library(tidyverse)
library(ggpubr)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library("factoextra")
library(readxl)
library(psych)
#Creating data subsets
Data <- Fish[c(7,24,31,34,37,40,47,51)]
PFAS <- Data[c(2,7,8)]
PFAS <- PFAS %>% 
  group_by(c(SUBSTANCE_NAME_ABBREVIATION)) %>% 
  mutate(row_id = row_number()) %>% 
  pivot_wider(id_cols = "SUBSTANCE_NAME_ABBREVIATION",  names_from = row_id,
              values_from = LAB_RESULT, 
              names_glue = "{.value}{row_id}") %>% 
  ungroup()
PFAS <- t(PFAS)
write.csv(PFAS, "C:\\Users\\ericr\\OneDrive\\Desktop\\SUNY ESF\\Spring 2023\\R&Repro\\Final\\PFAS.csv")
PFAS <- read_excel("PFAS.xlsx", na = "NA")

PFOS <- subset(Data, SUBSTANCE_NAME_ABBREVIATION == "PFOS")
head(PFAS)
head(Fish)

#plots
PFASscale <- scale(PFAS)
corr_matrix <- cor(PFAS)
corr_matrix
cor_test_mat <- corr.test(corr_matrix)$p
cor_test_mat

ggcorrplot(corr_matrix, 
           type = "lower",
           lab = TRUE,
           p.mat = cor_test_mat,
           insig = "pch",
           pch = 4,
           pch.cex = 10,
           ggtheme = ggpubr::theme_classic2(),
           title = "PFAS Correlation Matrix",
           lab_size = 2.5) +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(plot.title=element_text(face="bold", 
                             size=15, hjust=0.6))

ggsave("PFAScorrmatrix.jpg", dpi=300, device=NULL)

ggplot(PFOS, aes(x=PFOS$`WEIGHT_MEAN (G)`, y=PFOS$`LENGTH_MEAN (MM)`)) +
  geom_point() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  labs(title = "Length versus Weight",
       x = "Weight (g)",
       y = "Length (mm)") +
  geom_smooth(method = "lm", se = T)
ggsave("lengthtoweight.jpg", dpi=300, device=NULL)
lwlm <- lm(PFOS$`LENGTH_MEAN (MM)`~PFOS$`WEIGHT_MEAN (G)`)
summary(lwlm)

ggplot(PFOS, aes(x=PFOS$LAB_RESULT, y=PFOS$`LENGTH_MEAN (MM)`)) +
  geom_point() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  labs(title = "Length versus PFOS \n Concentration",
       x = "Concentration (ppb)",
       y = "Length (mm)") +
  geom_smooth(method = "lm", se = T)
ggsave("lengthtopfos.jpg", dpi=300, device=NULL)
pfosllm <- lm(PFOS$`LENGTH_MEAN (MM)`~PFOS$LAB_RESULT)
summary(pfosllm)

ggplot(PFOS, aes(x=PFOS$LAB_RESULT, y=PFOS$`WEIGHT_MEAN (G)`)) +
  geom_point() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  labs(title = "Weight versus PFOS \n Concentration",
       x = "Concentration (ppb)",
       y = "Weight (g)") +
  geom_smooth(method = "lm", se = T)
ggsave("weighttopfos.jpg", dpi=300, device=NULL)
pfoswlm <- lm(PFOS$`WEIGHT_MEAN (G)`~PFOS$LAB_RESULT)
summary(pfoswlm)

ggplot(PFOS, aes(x=PFOS$COMMON_NAME, y=PFOS$LAB_RESULT)) +
  geom_boxplot()

#looks ugly so lets subset
PFOSsub <- PFOS %>%
  count(COMMON_NAME) %>%
  top_n(10)

PFOSsubset <- PFOS %>% filter(COMMON_NAME == c("American Eel", "Bluegill", 
  "Brown Trout", "Common Carp", "Largemouth Bass", 
  "Pumpkinseed", "Walleye", "White Sucker", "Yellow Perch"))


ggplot(PFOS, aes(x=COMMON_NAME, y=LAB_RESULT, fill = COMMON_NAME)) +
  geom_boxplot() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PFOS Concentration \n by Species",
       y = "Concentration (ppb)",
       x = "Species",
       caption = "Cutoff at 2400, 1 value greater \n (7700 ppb BT @ Six Mile Creek)") +
  theme(legend.position = "none") +
  stat_compare_means(method = "anova", hjust=0, vjust = 2.2) +
scale_y_continuous(expand = c(0,0), limits = c(0,2400), breaks = seq(0, 2400, by=200))
ggsave("pfosspp.jpg", dpi=300, device=NULL)

PFOSsubwater <- PFOS %>%
  count(WATERBODY_NAME) %>%
  top_n(10)
PFOSsubsetwater <- subset(PFOS, WATERBODY_NAME == c("Hoosic River", "Little Hoosic River",
    "Beaverdam Lake", "Six Mile Creek", "Browns Pond", "Barge Canal", "Thayers Pond", 
    "Lockwood Basin", "Washington Lake", "Moodna Creek"))

ggplot(PFOSeels, aes(x=WATERBODY_NAME, y=LAB_RESULT, fill = WATERBODY_NAME)) +
  geom_boxplot() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PFOS Concentration in \n American Eels by Waterbody",
       y = "Concentration (ppb)",
       x = "Waterbody") +
       #caption = "Cutoff at 2400, 1 value greater \n (7700 ppb BT @ Six Mile Creek)") +
  theme(legend.position = "none")
  #stat_compare_means(method = "anova", hjust=0, vjust = 3.3)
  #scale_y_continuous(expand = c(0,0), limits = c(0,2400), breaks = seq(0, 2400, by=200))
ggsave("pfoswbeelsall.jpg", dpi=300, device=NULL)


fulllm <- lm(PFOS$LAB_RESULT~PFOS$`WEIGHT_MEAN (G)`+PFOS$`LENGTH_MEAN (MM)`)
summary(fulllm)

PFOSeel <- subset(PFOSsubsetwater, COMMON_NAME == "American Eel")

PFOSeels <- subset(PFOS, COMMON_NAME == "American Eel")
eellm <- lm(PFOSeels$LAB_RESULT~PFOSeels$`WEIGHT_MEAN (G)`+PFOSeels$`LENGTH_MEAN (MM)`)
summary(eellm)

PFOSeel %>% distinct(WATERBODY_NAME)
