library(tidyverse)
library(ggpubr)
summary(Derby)
Derby$SumPFAS <- rowSums(Derby[,c(7:45)])
Derby$SumPFASpergram <- Derby$SumPFAS/Derby$Weight

#multivariate, correlation, and pca
library(MVN)
mvn(Derby[,40:41], mvnTest = "mardia", univariateTest = "SW", 
                multivariatePlot = "qq", covariance = T)

library(corrr)
library(ggcorrplot)
library(FactoMineR)
DerbyNormal <- Derby[,c(4,5,6,8,27,28,29,31,33,34,35,36,37,39,40,43,44,45,46)]
DerbyNormal <- scale(DerbyNormal)
corr_matrix <- cor(DerbyNormal)
corr_matrix
library(psych)
cor_test_mat <- corr.test(corr_matrix)$p
cor_test_mat

ggcorrplot(corr_matrix,  p.mat = cor_test_mat, 
           hc.order = TRUE,
           type = "lower",
           lab = TRUE,
           ggtheme = ggpubr::theme_classic2(),
           title = "Seneca Lake PFAS Correlation Matrix",
           lab_size = 1.5) +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  theme(legend.key.size = unit(0.4, 'cm'))
?ggcorrplot
data.pca <- princomp(corr_matrix)
ggsave("PFAScorrmatrixwinsig.jpg", dpi=300, device=NULL)

pca <- prcomp(Derby[,c(8,27,28,29,31,33,34,35,36,37,39,40,43,44,45)], scale = T)

#PCA plot
library("factoextra")
fviz_eig(pca, addlabels = TRUE)
fviz_pca_ind(pca, habillage = Derby$Species, label = "none", 
             title = "PFAS PCA by Species") +
  theme_classic2() + 
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.6))
ggsave("PFASPCAbyspp.jpg", dpi=300, device=NULL)

#subsetting
Stivers<- subset(Derby, Marina=="Stivers")
Watkins <- subset(Derby, Marina=="Watkins")
LT <- subset(Derby, Species == "LT")
LLS <- subset(Derby, Species == "LLS")
BNT <-subset(Derby, Species == "BNT")
RBT <- subset(Derby, Species == "RBT")
#PFCA <- Derby[,c(1,2,3,4,5,27,29,30,33,34,36,38,40,42,44,45,46)]
#PFSA <- Derby[,c(1,2,3,4,5,28,31,32,35,37,39,41,43)]

#univariate testing - t.test
t.test(Stivers$Hg,Watkins$Hg)
t.test(Stivers$PFOS, Watkins$PFOS)
t.test(LT$Hg, BNT$Hg)
t.test(LT$Hg, RBT$Hg)
t.test(LT$PFOS, RBT$PFOS)

#ANOVA
Hgaov <- aov(Hg~Species, data=Derby)
summary(Hgaov)
Hglakeaov<- aov(Hg ~ Marina, data=Derby)
summary(Hglakeaov)
PFOSaov <- aov(PFOS~Species, data=Derby)
summary(PFOSaov)

PFOSaov <- aov(PFOS~Species, data=Derby)
summary(PFOSaov)
tukeyPFOS <- TukeyHSD(PFOSaov)
tukeyPFOS
cldPFOS <- multcompLetters4(PFOSaov, tukeyPFOS)
cldPFOS

PFUdAaov <- aov(PFUdA~Species, data=Derby)
summary(PFUdAaov)
tukeyPFUdA <-TukeyHSD(PFUdAaov)
tukeyPFUdA
cldPFUdA <- multcompLetters4(PFUdAaov, tukeyPFUdA)
print(cldPFNA)

summary(BNT)
summary(PFOSspp)
t.test(LT$PFOS, BNT$PFOS)
t.test(LT$PFOS, RBT$PFOS)
t.test(LT$PFOS, LLS$PFOS)
t.test(BNT$PFOS, RBT$PFOS)
t.test(BNT$PFOS, LLS$PFOS)
mean(LT$PFOS)
Derbymanova <- manova(cbind(8,27,28,29,34,35,36,37,38,44)~Species + Marina, data=Derby)
#variable lengths for each set, hard to do

#multiple comparison to add to plots if needed
library(multcompView)
cldPFNA <- multcompLetters4(PFNAaov, tukeyPFNA)
print(cldPFNA)


# plots
#PFXX by species general code, can do any ind chemical with this set by changing y
compare_means( data = Derby)


ggplot(Derby, aes(x=Species, y=Length, fill=Species)) +
  geom_boxplot() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.6)) +
  labs(title = "Length by Species",
       x = "Species",
       y = "Length (cm)") +
  stat_compare_means(method = "anova", hjust=0.5) #label = "p.format")
ggsave("lengthbyspp.jpg", dpi=300, device=NULL)

?stat_compare_means
#stacked barplot for all pfas by species
# y=concentration, fill=name, x=species
Derbysubset <- Derby[,c(1,2,3,7:45)]
derbymelt <-melt(Derbysubset)

ggplot(derbymelt, aes(x=Species, fill=variable, y=value)) +
  geom_bar(position = "stack", stat="identity") +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.6)) +
  labs(title = "PFAS Concentration by Species",
       x = "Species",
       y = "PFAS (ppb)",
       fill = "PFAS Name") +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.key.size = unit(0.3, 'cm'))
ggsave("PFASConcbysppbystack.jpg", dpi=300, device=NULL)

ggplot(derbymelt, aes(x=Marina, fill=variable, y=value)) +
  geom_bar(position = "stack", stat="identity") +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.6)) +
  labs(title = "PFAS Concentration by Marina",
       x = "Marina",
       y = "PFAS (ppb)",
       fill = "PFAS Name") +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.key.size = unit(0.3, 'cm'))
ggsave("PFASConcbymarinabystack.jpg", dpi=300, device=NULL)

#transposing and melting all pfas and all spp
PFAS <- t(Derby[7:45])
PFAS <-as.data.frame(PFAS)
PFAS[ "Name" ] <- rownames(PFAS)
PFAS <-melt(PFAS)

PFSAmelt <- t(PFSA[6:13])
PFSAmelt <-as.data.frame(PFSAmelt)
PFSAmelt[ "Name" ] <- rownames(PFSAmelt)
PFSAmelt <-melt(PFSAmelt)

#all species represented by chemical
ggplot(PFAS, aes(x=Name, y=value)) +
  geom_boxplot(varwidth=T) +
  theme_classic2() + 
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, color="Black", face="bold")) +
  labs(title = "PFAS Concentration by Compound",
       x = "PFAS Chemicals",
       y = "Concentration (ppb)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,50), breaks = seq(0, 50, by=10))
ggsave("PFASconcbycompound.jpg", dpi=300, device=NULL)

#ind spp represented, still by compound
#Lake Trout
LTmelt <- t(LT[7:45])
LTmelt <-as.data.frame(LTmelt)
LTmelt[ "Name" ] <- rownames(LTmelt)
library(reshape2)
LTmelt <-melt(LTmelt)

ggplot(LTmelt, aes(x=Name, y=value)) +
  geom_boxplot(varwidth=T) +
  theme_classic2() + 
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, color="Black", face="bold")) +
  labs(title = "PFAS Concentration by Compound in Lake Trout",
       x = "PFAS Chemicals",
       y = "Concentration (ppb)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,50), breaks = seq(0, 50, by=10))
ggsave("PFASconcinLT.jpg", dpi=300, device=NULL)

#Land Locked Salmon
LLSmelt <- t(LLS[7:45])
LLSmelt <-as.data.frame(LLSmelt)
LLSmelt[ "Name" ] <- rownames(LLSmelt)
library(reshape2)
LLSmelt <-melt(LLSmelt)

ggplot(LLSmelt, aes(x=Name, y=value)) +
  geom_boxplot(varwidth=T) +
  theme_classic2() + 
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, color="Black", face="bold")) +
  labs(title = "PFAS Concentration by Compound in \n Land Locked Salmon",
       x = "PFAS Chemicals",
       y = "Concentration (ppb)") +
  scale_y_continuous(expand = c(0,0))
ggsave("PFASconcinLLS.jpg", dpi=300, device=NULL)

#Brown Trout
BNTmelt <- t(BNT[7:45])
BNTmelt <-as.data.frame(BNTmelt)
BNTmelt[ "Name" ] <- rownames(BNTmelt)
library(reshape2)
BNTmelt <-melt(BNTmelt)

ggplot(BNTmelt, aes(x=Name, y=value)) +
  geom_boxplot(varwidth=T) +
  theme_classic2() + 
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, color="Black", face="bold")) +
  labs(title = "PFAS Concentration by Compound in \n Brown Trout",
       x = "PFAS Chemicals",
       y = "Concentration (ppb)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,5), breaks = seq(0, 5, by=1))
ggsave("PFASconcinBNT.jpg", dpi=300, device=NULL)

#Rainbow Trout
RBTmelt <- t(RBT[7:45])
RBTmelt <-as.data.frame(RBTmelt)
RBTmelt[ "Name" ] <- rownames(RBTmelt)
library(reshape2)
RBTmelt <-melt(RBTmelt)

ggplot(RBTmelt, aes(x=Name, y=value)) +
  geom_boxplot(varwidth=T) +
  theme_classic2() + 
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, color="Black", face="bold")) +
  labs(title = "PFAS Concentration by Compound in \n Rainbow Trout",
       x = "PFAS Chemicals",
       y = "Concentration (ppb)") +
  scale_y_continuous(expand = c(0,0)) #, limits = c(0,5), breaks = seq(0, 5, by=1))
ggsave("PFASconcinRBT.jpg", dpi=300, device=NULL)

#Boxplot of all species PFXX concentration by dependent variable (marina, spp, etc.)
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )


ggplot(Derby, aes(x=Species, y=PFOS, fill=Species)) +
  geom_boxplot() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  labs(title = "PFOS Concentration by Species",
       x = "Species",
       y = "PFOS (ppb)") +
  stat_compare_means(method = "anova", hjust=0.5)
ggsave("pfosconcbyspp.jpg", dpi=300, device=NULL)

#both marinas barplot by species
ggplot(Derby, aes(x=Marina, y=SumPFAS, fill=Marina)) +
  geom_boxplot() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  labs(title = "Summated PFAS Concentration by Marina",
       x = "Marina",
       y = "Summated PFAS (ppb)") +
  stat_compare_means(method = "t.test", hjust=-0.5)
ggsave("sumpfasconcbymarina.jpg", dpi=300, device=NULL)

#stivers barplot by species
ggplot(Stivers, aes(x=Species, y=SumPFAS, fill=Species)) +
  geom_boxplot() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  labs(title = "Summated PFAS Concentration \n at Stivers Marina",
       x = "Species",
       y = "PFAS (ppb)") +
  stat_compare_means(method = "anova", hjust=1.2)
ggsave("sumpfasconcstivers.pdf", dpi=300, device=NULL)

compare_means(SumPFAS ~ Species,  data = Stivers)

#watkins glen barplot by species
ggplot(Watkins, aes(x=Species, y=SumPFAS, fill=Species)) +
  geom_boxplot() +
  theme_classic2() +
  theme(plot.title=element_text(face="bold", 
                                size=15, hjust=0.5)) +
  labs(title = "Summated PFAS Concentration \n at Watkins Glen Marina",
       x = "Species",
       y = "PFAS (ppb)") +
  stat_compare_means(method = "anova", hjust=0.5)
ggsave("sumpfasconcwatkins.pdf", dpi=300, device=NULL)

compare_means(SumPFAS ~ Species,  data = Watkins)

#Linear model comparing 2 compounds to one another, then subsequent lm to add to plot
ggplot(Derby, aes(x=PFOS, y=PFHpS)) +
  geom_point() +
  theme_classic2() +
  geom_smooth(method="lm")

pfossumlm <- lm(Derby$PFOS~Derby$PFHpS)
summary(pfossumlm)
