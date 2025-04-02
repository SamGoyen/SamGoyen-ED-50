#setwd to source file location
library(lmerTest)
library(emmeans)
library(sjPlot)
library(drc)
library(Rmisc)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(reshape2)
library(tidyverse)
library(broom)
library(forcats)

PAM_data<-read.csv("AS-2022_CBASS_PAM.csv")
Allpam<-read.delim("/Users/mkwoo/Documents/Australia/Global search/All_pam_data.txt")


Allpam$Site<-as.factor(Allpam$Site)
Allpam$Geno<-as.factor(Allpam$Geno)
Allpam$Species<-as.factor(Allpam$Species)
Allpam$Temp<-as.factor(Allpam$Temp)

PAM_data<-subset(PAM_data, Timepoint == 'TP')

Allpam$Site_Species[Name=(Allpam$Site == "NthDirection") & (Allpam$Species=="A.hya")]<-"NthDirection_A.hya"
Allpam$Site_Species[Name=(Allpam$Site == "NthDirection") & (Allpam$Species=="S.pis")]<-"NthDirection_S.pis"
Allpam$Site_Species[Name=(Allpam$Site == "NthDirection") & (Allpam$Species=="P.lob")]<-"NthDirection_P.lob"
Allpam$Site_Species[Name=(Allpam$Site == "NthDirection") & (Allpam$Species=="P.ver")]<-"NthDirection_P.ver"

Allpam$Site_Species[Name=(Allpam$Site == "Moore") & (Allpam$Species=="A.hya")]<-"Moore_A.hya"
Allpam$Site_Species[Name=(Allpam$Site == "Moore") & (Allpam$Species=="S.pis")]<-"Moore_S.pis"
Allpam$Site_Species[Name=(Allpam$Site == "Moore") & (Allpam$Species=="P.lob")]<-"Moore_P.lob"
Allpam$Site_Species[Name=(Allpam$Site == "Moore") & (Allpam$Species=="P.ver")]<-"Moore_P.ver"

Allpam$Site_Species[Name=(Allpam$Site == "Davies") & (Allpam$Species=="A.hya")]<-"Davies_A.hya"
Allpam$Site_Species[Name=(Allpam$Site == "Davies") & (Allpam$Species=="S.pis")]<-"Davies_S.pis"
Allpam$Site_Species[Name=(Allpam$Site == "Davies") & (Allpam$Species=="P.lob")]<-"Davies_P.lob"
Allpam$Site_Species[Name=(Allpam$Site == "Davies") & (Allpam$Species=="P.ver")]<-"Davies_P.ver"

Allpam$Site_Species<-as.factor(Allpam$Site_Species)

NthDirection_data<-subset(Allpam, Site == 'NthDirection')
Moore_data<-subset(Allpam, Site == 'Moore')
Davies_data<-subset(Allpam, Site == 'Davies')

A.hya_data<-subset(Allpam, Species == 'A.hya')
S.pis_data<-subset(Allpam, Species == 'S.pis')
P.lob_data<-subset(Allpam, Species == 'P.lob')
P.ver_data<-subset(Allpam, Species == 'P.ver')

str(Allpam)

#Run individual genotype fits for each species/population

################################################
#################### NthDirection ####################
################################################

#### NthDirection_A.hya ####
NthDirection_A.hya <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="A.hya",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))
NthDirection_A.hya_pop <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="A.hya",], fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))

summary(NthDirection_A.hya_pop)
plot(NthDirection_A.hya_pop)

#extract coeffs by geno, then compute 95% CIs
NthDirection_A.hya_genocoeffs_50<-data.frame(ED(NthDirection_A.hya, c(50)))
NthDirection_A.hya_coeff_mean<-mean(NthDirection_A.hya_genocoeffs_50$Estimate)
NthDirection_A.hya_coeff_mean

NthDirection_A.hya_summary<-data.frame(CI(NthDirection_A.hya_genocoeffs_50$Estimate, ci=0.95))
NthDirection_A.hya_coeff_lower<-NthDirection_A.hya_summary[3,]
NthDirection_A.hya_coeff_upper<-NthDirection_A.hya_summary[1,]

#### NthDirection_S.pis ####
NthDirection_S.pis <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="S.pis",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))
NthDirection_S.pis_pop <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="S.pis",], fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))

summary(NthDirection_S.pis)
plot(NthDirection_S.pis)

#extract coeffs by geno, then compute 95% CIs
NthDirection_S.pis_genocoeffs_50<-data.frame(ED(NthDirection_S.pis, c(50)))
NthDirection_S.pis_coeff_mean<-mean(NthDirection_S.pis_genocoeffs_50$Estimate)
NthDirection_S.pis_coeff_mean

NthDirection_S.pis_summary<-data.frame(CI(NthDirection_S.pis_genocoeffs_50$Estimate, ci=0.95))
NthDirection_S.pis_coeff_lower<-NthDirection_S.pis_summary[3,]
NthDirection_S.pis_coeff_upper<-NthDirection_S.pis_summary[1,]

##### NthDirection_P.lob ##### 
NthDirection_P.lob <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="P.lob",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')))
NthDirection_P.lob_pop <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="P.lob",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(NthDirection_P.lob)
plot(NthDirection_P.lob)

#extract coeffs by geno, then compute 95% CIs
NthDirection_P.lob_genocoeffs_50<-data.frame(ED(NthDirection_P.lob, c(50)))
NthDirection_P.lob_coeff_mean<-mean(NthDirection_P.lob_genocoeffs_50$Estimate)
NthDirection_P.lob_coeff_mean

NthDirection_P.lob_summary<-data.frame(CI(NthDirection_P.lob_genocoeffs_50$Estimate, ci=0.95))
NthDirection_P.lob_coeff_lower<-NthDirection_P.lob_summary[3,]
NthDirection_P.lob_coeff_upper<-NthDirection_P.lob_summary[1,]

##### NthDirection_P.ver ##### 
NthDirection_P.ver <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="P.ver",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')))
NthDirection_P.ver_pop <- drm(PAM ~ Temp, data = NthDirection_data[NthDirection_data$Species=="P.ver",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(NthDirection_P.ver)
plot(NthDirection_P.ver)

#extract coeffs by geno, then compute 95% CIs
NthDirection_P.ver_genocoeffs_50<-data.frame(ED(NthDirection_P.ver, c(50)))
NthDirection_P.ver_coeff_mean<-mean(NthDirection_P.ver_genocoeffs_50$Estimate)
NthDirection_P.ver_coeff_mean

NthDirection_P.ver_summary<-data.frame(CI(NthDirection_P.ver_genocoeffs_50$Estimate, ci=0.95))
NthDirection_P.ver_coeff_lower<-NthDirection_P.ver_summary[3,]
NthDirection_P.ver_coeff_upper<-NthDirection_P.ver_summary[1,]

################################################
#################### Moore ####################
################################################

#### Moore_A.hya ####
Moore_A.hya <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="A.hya",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')))
Moore_A.hya_pop <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="A.hya",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Moore_A.hya)
plot(Moore_A.hya)

#extract coeffs by geno, then compute 95% CIs
Moore_A.hya_genocoeffs_50<-data.frame(ED(Moore_A.hya, c(50)))
Moore_A.hya_coeff_mean<-mean(Moore_A.hya_genocoeffs_50$Estimate)
Moore_A.hya_coeff_mean

Moore_A.hya_summary<-data.frame(CI(Moore_A.hya_genocoeffs_50$Estimate, ci=0.95))
Moore_A.hya_coeff_lower<-Moore_A.hya_summary[3,]
Moore_A.hya_coeff_upper<-Moore_A.hya_summary[1,]

##### Moore_S.pis  ##### 
Moore_S.pis <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="S.pis",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))
Moore_S.pis_pop <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="S.pis",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Moore_S.pis)
plot(Moore_S.pis)

#extract coeffs by geno, then compute 95% CIs
Moore_S.pis_genocoeffs_50<-data.frame(ED(Moore_S.pis, c(50)))
Moore_S.pis_coeff_mean<-mean(Moore_S.pis_genocoeffs_50$Estimate)
Moore_S.pis_coeff_mean

Moore_S.pis_summary<-data.frame(CI(Moore_S.pis_genocoeffs_50$Estimate, ci=0.95))
Moore_S.pis_coeff_lower<-Moore_S.pis_summary[3,]
Moore_S.pis_coeff_upper<-Moore_S.pis_summary[1,]

##### Moore_P.lob ##### 
Moore_P.lob <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="P.lob",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.7, 45))
Moore_P.lob_pop <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="P.lob",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Moore_P.lob)
plot(Moore_P.lob)

#extract coeffs by geno, then compute 95% CIs
Moore_P.lob_genocoeffs_50<-data.frame(ED(Moore_P.lob, c(50)))
Moore_P.lob_coeff_mean<-mean(Moore_P.lob_genocoeffs_50$Estimate)
Moore_P.lob_coeff_mean

Moore_P.lob_summary<-data.frame(CI(Moore_P.lob_genocoeffs_50$Estimate, ci=0.95))
Moore_P.lob_coeff_lower<-Moore_P.lob_summary[3,]
Moore_P.lob_coeff_upper<-Moore_P.lob_summary[1,]

##### Moore_P.ver ##### 
Moore_P.ver <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="P.ver",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))
Moore_P.ver_pop <- drm(PAM ~ Temp, data = Moore_data[Moore_data$Species=="P.ver",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Moore_P.ver)
plot(Moore_P.ver)

#extract coeffs by geno, then compute 95% CIs
Moore_P.ver_genocoeffs_50<-data.frame(ED(Moore_P.ver, c(50)))
Moore_P.ver_coeff_mean<-mean(Moore_P.ver_genocoeffs_50$Estimate)
Moore_P.ver_coeff_mean

Moore_P.ver_summary<-data.frame(CI(Moore_P.ver_genocoeffs_50$Estimate, ci=0.95))
Moore_P.ver_coeff_lower<-Moore_P.ver_summary[3,]
Moore_P.ver_coeff_upper<-Moore_P.ver_summary[1,]

##### ##### ##### ##### ##### 
##### ## Davies ## ##### 
##### ##### ##### ##### ##### 

#### Davies_A.hya ####
Davies_A.hya <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="A.hya",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.7, 45))
Davies_A.hya_pop <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="A.hya",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Davies_A.hya)
plot(Davies_A.hya)

#extract coeffs by geno, then compute 95% CIs
Davies_A.hya_genocoeffs_50<-data.frame(ED(Davies_A.hya, c(50)))
Davies_A.hya_coeff_mean<-mean(Davies_A.hya_genocoeffs_50$Estimate)
Davies_A.hya_coeff_mean

Davies_A.hya_summary<-data.frame(CI(Davies_A.hya_genocoeffs_50$Estimate, ci=0.95))
Davies_A.hya_coeff_lower<-Davies_A.hya_summary[3,]
Davies_A.hya_coeff_upper<-Davies_A.hya_summary[1,]

##### Davies_S.pis ##### 
Davies_S.pis <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="S.pis",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))
Davies_S.pis_pop <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="S.pis",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Davies_S.pis)
plot(Davies_S.pis)

#extract coeffs by geno, then compute 95% CIs
Davies_S.pis_genocoeffs_50<-data.frame(ED(Davies_S.pis, c(50)))
Davies_S.pis_coeff_mean<-mean(Davies_S.pis_genocoeffs_50$Estimate)
Davies_S.pis_coeff_mean

Davies_S.pis_summary<-data.frame(CI(Davies_S.pis_genocoeffs_50$Estimate, ci=0.95))
Davies_S.pis_coeff_lower<-Davies_S.pis_summary[3,]
Davies_S.pis_coeff_upper<-Davies_S.pis_summary[1,]

##### Davies_P.lob ##### 
Davies_P.lob <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="P.lob",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')))
Davies_P.lob_pop <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="P.lob",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Davies_P.lob)
plot(Davies_P.lob)

#extract coeffs by geno, then compute 95% CIs
Davies_P.lob_genocoeffs_50<-data.frame(ED(Davies_P.lob, c(50)))
Davies_P.lob_coeff_mean<-mean(Davies_P.lob_genocoeffs_50$Estimate)
Davies_P.lob_coeff_mean

Davies_P.lob_summary<-data.frame(CI(Davies_P.lob_genocoeffs_50$Estimate, ci=0.95))
Davies_P.lob_coeff_lower<-Davies_P.lob_summary[3,]
Davies_P.lob_coeff_upper<-Davies_P.lob_summary[1,]

##### Davies_P.ver ##### 
Davies_P.ver <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="P.ver",], curveid=Geno, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl = c(1000, 0.9, 45))
Davies_P.ver_pop <- drm(PAM ~ Temp, data = Davies_data[Davies_data$Species=="P.ver",], fct = LL.3(names = c('hill', 'max', 'ed50')))

summary(Davies_P.ver)
plot(Davies_P.ver)

#extract coeffs by geno, then compute 95% CIs
Davies_P.ver_genocoeffs_50<-data.frame(ED(Davies_P.ver, c(50)))
Davies_P.ver_coeff_mean<-mean(Davies_P.ver_genocoeffs_50$Estimate)
Davies_P.ver_coeff_mean

Davies_P.ver_summary<-data.frame(CI(Davies_P.ver_genocoeffs_50$Estimate, ci=0.95))
Davies_P.ver_coeff_lower<-Davies_P.ver_summary[3,]
Davies_P.ver_coeff_upper<-Davies_P.ver_summary[1,]

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Combine genotpye-ED50s into dataframe for statistical analysis ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

NthDirection_A.hya_ED50s<-data.frame(NthDirection_A.hya_genocoeffs_50[,1])
NthDirection_S.pis_ED50s<-data.frame(NthDirection_S.pis_genocoeffs_50[,1])
NthDirection_P.lob_ED50s<-data.frame(NthDirection_P.lob_genocoeffs_50[,1])
NthDirection_P.ver_ED50s<-data.frame(NthDirection_P.ver_genocoeffs_50[,1])

Moore_A.hya_ED50s<-data.frame(Moore_A.hya_genocoeffs_50[,1])
Moore_S.pis_ED50s<-data.frame(Moore_S.pis_genocoeffs_50[,1])
Moore_P.lob_ED50s<-data.frame(Moore_P.lob_genocoeffs_50[,1])
Moore_P.ver_ED50s<-data.frame(Moore_P.ver_genocoeffs_50[,1])

Davies_A.hya_ED50s<-data.frame(Davies_A.hya_genocoeffs_50[,1])
Davies_S.pis_ED50s<-data.frame(Davies_S.pis_genocoeffs_50[,1])
Davies_P.lob_ED50s<-data.frame(Davies_P.lob_genocoeffs_50[,1])
Davies_P.ver_ED50s<-data.frame(Davies_P.ver_genocoeffs_50[,1])

library(tidyverse)
Geno_ED50s<-list(NthDirection_A.hya_ED50s,NthDirection_S.pis_ED50s,NthDirection_P.lob_ED50s,NthDirection_P.ver_ED50s,Moore_A.hya_ED50s,
                 Moore_S.pis_ED50s,Moore_P.lob_ED50s,Moore_P.ver_ED50s,Davies_A.hya_ED50s,Davies_S.pis_ED50s,Davies_P.lob_ED50s,Davies_P.ver_ED50s) %>% 
  map(~ .x %>% 
        as.data.frame %>%
        rownames_to_column('rn')) %>% 
  reduce(full_join, by = 'rn') %>%
  column_to_rownames('rn')

Geno_ED50s<-Geno_ED50s %>% 
  dplyr::rename(NthDirection_A.hya=NthDirection_A.hya_genocoeffs_50...1.,
                NthDirection_S.pis=NthDirection_S.pis_genocoeffs_50...1.,
                NthDirection_P.lob=NthDirection_P.lob_genocoeffs_50...1.,
                NthDirection_P.ver=NthDirection_P.ver_genocoeffs_50...1.,
                Moore_A.hya=Moore_A.hya_genocoeffs_50...1.,
                Moore_S.pis=Moore_S.pis_genocoeffs_50...1.,
                Moore_P.lob=Moore_P.lob_genocoeffs_50...1.,
                Moore_P.ver=Moore_P.ver_genocoeffs_50...1.,
                Davies_A.hya=Davies_A.hya_genocoeffs_50...1.,
                Davies_S.pis=Davies_S.pis_genocoeffs_50...1.,
                Davies_P.lob=Davies_P.lob_genocoeffs_50...1.,
                Davies_P.ver=Davies_P.ver_genocoeffs_50...1.)

Geno_ED50s$Geno<-as.factor(1:nrow(Geno_ED50s))
str(Geno_ED50s)

Geno_ED50s_long<-melt(Geno_ED50s, id="Geno")

Geno_ED50s_long<-Geno_ED50s_long %>% 
  dplyr::rename(Site_Species= variable,
                ED50=value)

Geno_ED50s_long<-na.omit(Geno_ED50s_long)

#Compare groups statistically
ED50_mod<-aov(ED50 ~ Site_Species, Geno_ED50s_long)
summary(ED50_mod)
TukeyHSD(ED50_mod)

#### Create Rank table #### 
Rank_table<-Geno_ED50s_long %>%
  arrange(Geno, Site_Species) %>%
  group_by(Site_Species) %>% 
  mutate(rank_order = rank(desc(ED50))) %>%
  arrange(rank_order)

#write.csv(Rank_table,'Genotype_ranks_full.csv')

Top_performers<-Rank_table %>%
  group_by(Site_Species) %>%
  slice_min(rank_order, n=5)

Bottom_performers<-Rank_table %>%
  group_by(Site_Species) %>%
  slice_max(rank_order, n=5)

Top_bottom_performers<-rbind(Top_performers,Bottom_performers)

#write.csv(Top_bottom_performers,'Genotype_ranks_top-bottom-performers.csv')

#### Subset for individual genotype DRC plots of top/bottom performers ####
#NthDirection
NthDirection_A.hya_performers<-filter(Allpam,Site_Species == 'NthDirection_A.hya' & Geno %in% c(13,9,11,1,4,10,8,5,7,6))
NthDirection_S.pis_performers<-filter(Allpam,Site_Species == 'NthDirection_S.pis' & Geno %in% c(5,9,11,13,14,3,4,6,12,1))
NthDirection_P.lob_performers<-filter(Allpam,Site_Species == 'NthDirection_P.lob' & Geno %in% c(14,3,5,6,7,4,12,10,13,11))
NthDirection_P.ver_performers<-filter(Allpam,Site_Species == 'NthDirection_P.ver' & Geno %in% c(5,6,4,7,9,2,10,1,15,14))

NthDirection_performers<-rbind(NthDirection_A.hya_performers,NthDirection_S.pis_performers,NthDirection_P.lob_performers,NthDirection_P.ver_performers)

#Moore
Moore_A.hya_performers<-filter(Allpam,Site_Species == 'Moore_A.hya' & Geno %in% c(14,8,13,5,10,7,3,2,1,4))
Moore_S.pis_performers<-filter(Allpam,Site_Species == 'Moore_S.pis' & Geno %in% c(9,8,11,10,1,5,2,14,6,13))
Moore_P.lob_performers<-filter(Allpam,Site_Species == 'Moore_P.lob' & Geno %in% c(2,6,14,1,3,5,7,13,11,12))
Moore_P.ver_performers<-filter(Allpam,Site_Species == 'Moore_P.ver' & Geno %in% c(10,15,12,9,16,13,2,4,3,14))

Moore_performers<-rbind(Moore_A.hya_performers,Moore_S.pis_performers,Moore_P.lob_performers,Moore_P.ver_performers)

#Davies
Davies_A.hya_performers<-filter(Allpam,Site_Species == 'Davies_A.hya' & Geno %in% c(7,11,6,8,5,2,13,12,14,10))
Davies_S.pis_performers<-filter(Allpam,Site_Species == 'Davies_S.pis' & Geno %in% c(10,1,9,8,3,5,6,12,14,13))
Davies_P.lob_performers<-filter(Allpam,Site_Species == 'Davies_P.lob' & Geno %in% c(10,11,9,1,5,7,14,8,12,13))
Davies_P.ver_performers<-filter(Allpam,Site_Species == 'Davies_P.ver' & Geno %in% c(4,9,3,7,1,10,14,13,15,16))

Davies_performers<-rbind(Davies_A.hya_performers,Davies_S.pis_performers,Davies_P.lob_performers,Davies_P.ver_performers)

#For each genotype in each 'Site_Genus' grouping, run DRC
temp_x<- data.frame(temp = seq(28.5,39, length.out = 100))

NthDirection_Geno_DRCs <- NthDirection_performers %>%
  nest(dataset = c(-Site_Species, -Geno) ) %>% 
  mutate(model = map(dataset, ~drc::drm(formula = .x$PAM~.x$Temp, fct = LL.3(), upperl = c(1000, 0.8, 45),lowerl = c(0.1, 0.01, 26)) ),
         results = map(model,~ predict(.x,newdata = temp_x,conf.int = FALSE))) %>% 
  unnest(results)
NthDirection_Geno_DRCs$Temp<-rep(c(seq(28.5,39, length.out = 100)),times=40)

Moore_Geno_DRCs <- Moore_performers %>%
  nest(dataset = c(-Site_Species, -Geno) ) %>% 
  mutate(model = map(dataset, ~drc::drm(formula = .x$PAM~.x$Temp, fct = LL.3(), upperl = c(1000, 0.8, 45),lowerl = c(0.1, 0.01, 26)) ),
         results = map(model,~ predict(.x,newdata = temp_x,conf.int = FALSE))) %>% 
  unnest(results)
Moore_Geno_DRCs$Temp<-rep(c(seq(28.5,39, length.out = 100)),times=40)

Davies_Geno_DRCs <- Davies_performers %>%
  nest(dataset = c(-Site_Species, -Geno) ) %>% 
  mutate(model = map(dataset, ~drc::drm(formula = .x$PAM~.x$Temp, fct = LL.3(), upperl = c(1000, 0.8, 45),lowerl = c(0.1, 0.01, 26)) ),
         results = map(model,~ predict(.x,newdata = temp_x,conf.int = FALSE))) %>% 
  unnest(results)
Davies_Geno_DRCs$Temp<-rep(c(seq(28.5,39, length.out = 100)),times=40)

#NthDirection geno plots
library(RColorBrewer)
mypalette<-brewer.pal(10, "Spectral")
display.brewer.pal(10,"Spectral")

#low to high 10 total
mypalette<-mypalette[10:1]
mypalette

#### reorder Genotypes for geom_line, geom_vline, and geom_jitter data based on ranks not genotype numbers, then plot geno curves ####

#### NthDirection Ahya geno plot ####

#geom_jitter data
NthDirection_A.hya_performers %>%
  mutate(Geno = fct_relevel(Geno, '13','9','11','1','4','10','8','5','7','6'))

#geom_line data
NthDirection_A.hya_Geno_DRCs <- NthDirection_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='NthDirection_A.hya') 

levels(NthDirection_A.hya_Geno_DRCs$Geno)
NthDirection_A.hya_Geno_DRCs$Geno <- factor(NthDirection_A.hya_Geno_DRCs$Geno, levels=c('13','9','11','1','4','10','8','5','7','6'))

#geom_vline data
NthDirection_A.hya_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='NthDirection_A.hya')

levels(NthDirection_A.hya_performers_ED50s$Geno)
NthDirection_A.hya_performers_ED50s$Geno <- factor(NthDirection_A.hya_performers_ED50s$Geno, levels=c('13','9','11','1','4','10','8','5','7','6'))

#plot
ggplot() +
  geom_line(data = NthDirection_A.hya_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = NthDirection_A.hya_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = NthDirection_A.hya_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28.5,38.5), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("NthDirection Acropora hyacinthus") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### NthDirection Spis geno plot ####

#geom_jitter data
NthDirection_S.pis_performers %>%
  mutate(Geno = fct_relevel(Geno, '5','9','11','13','14','3','4','6','12','1'))

#geom_line data
NthDirection_S.pis_Geno_DRCs <- NthDirection_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='NthDirection_S.pis') 

levels(NthDirection_S.pis_Geno_DRCs$Geno)
NthDirection_S.pis_Geno_DRCs$Geno <- factor(NthDirection_S.pis_Geno_DRCs$Geno, levels=c('5','9','11','13','14','3','4','6','12','1'))

#geom_vline data
NthDirection_S.pis_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='NthDirection_S.pis')

levels(NthDirection_S.pis_performers_ED50s$Geno)
NthDirection_S.pis_performers_ED50s$Geno <- factor(NthDirection_S.pis_performers_ED50s$Geno, levels=c('5','9','11','13','14','3','4','6','12','1'))

#plot
ggplot() +
  geom_line(data = NthDirection_S.pis_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = NthDirection_S.pis_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = NthDirection_S.pis_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28.5,45), breaks=c(28,30,32,34,36,38,40,42,44)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("NthDirection_S.pis") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### NthDirection P.lob geno plot ####

#geom_jitter data
NthDirection_P.lob_performers %>%
  mutate(Geno = fct_relevel(Geno, '14','3','5','6','7','4','12','10','13','11'))

#geom_line data
NthDirection_P.lob_Geno_DRCs <- NthDirection_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='NthDirection_P.lob') 

levels(NthDirection_P.lob_Geno_DRCs$Geno)
NthDirection_S.pis_Geno_DRCs$Geno <- factor(NthDirection_P.lob_Geno_DRCs$Geno, levels=c('14','3','5','6','7','4','12','10','13','11'))

#geom_vline data
NthDirection_P.lob_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='NthDirection_P.lob')

levels(NthDirection_P.lob_performers_ED50s$Geno)
NthDirection_P.lob_performers_ED50s$Geno <- factor(NthDirection_P.lob_performers_ED50s$Geno, levels=c('14','3','5','6','7','4','12','10','13','11'))

#plot
ggplot() +
  geom_line(data = NthDirection_P.lob_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = NthDirection_P.lob_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = NthDirection_P.lob_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28.5,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("NthDirection_P.lob") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### NthDirection P.ver geno plot ####

#geom_jitter data
NthDirection_P.ver_performers %>%
  mutate(Geno = fct_relevel(Geno, '5','6','4','7','9','2','10','1','15','14'))

#geom_line data
NthDirection_P.ver_Geno_DRCs <- NthDirection_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='NthDirection_P.ver') 

levels(NthDirection_P.ver_Geno_DRCs$Geno)
NthDirection_P.ver_Geno_DRCs$Geno <- factor(NthDirection_P.ver_Geno_DRCs$Geno, levels=c('5','6','4','7','9','2','10','1','15','14'))

#geom_vline data
NthDirection_P.ver_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='NthDirection_P.ver')

levels(NthDirection_P.ver_performers_ED50s$Geno)
NthDirection_P.ver_performers_ED50s$Geno <- factor(NthDirection_P.ver_performers_ED50s$Geno, levels=c('5','6','4','7','9','2','10','1','15','14'))

#plot
ggplot() +
  geom_line(data = NthDirection_P.ver_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = NthDirection_P.ver_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = NthDirection_P.ver_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("NthDirection_P.ver") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Moore A.hya geno plot ####

#geom_jitter data
Moore_A.hya_performers %>%
  mutate(Geno = fct_relevel(Geno, '14','8','13','5','10','7','3','2','1','4'))

#geom_line data
Moore_A.hya_Geno_DRCs <- Moore_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Moore_A.hya') 

levels(Moore_A.hya_Geno_DRCs$Geno)
Moore_A.hya_Geno_DRCs$Geno <- factor(Moore_A.hya_Geno_DRCs$Geno, levels=c('14','8','13','5','10','7','3','2','1','4'))

#geom_vline data
Moore_A.hya_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Moore_A.hya')

levels(Moore_A.hya_performers_ED50s$Geno)
Moore_A.hya_performers_ED50s$Geno <- factor(Moore_A.hya_performers_ED50s$Geno, levels=c('14','8','13','5','10','7','3','2','1','4'))

#plot
ggplot() +
  geom_line(data = Moore_A.hya_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Moore_A.hya_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Moore_A.hya_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Moore_A.hya") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Moore S.pis geno plot ####

#geom_jitter data
Moore_S.pis_performers %>%
  mutate(Geno = fct_relevel(Geno, '9','8','11','10','1','5','2','14','6','13'))

#geom_line data
Moore_S.pis_Geno_DRCs <- Moore_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Moore_S.pis') 

levels(Moore_S.pis_Geno_DRCs$Geno)
Moore_S.pis_Geno_DRCs$Geno <- factor(Moore_S.pis_Geno_DRCs$Geno, levels=c('9','8','11','10','1','5','2','14','6','13'))

#geom_vline data
Moore_S.pis_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Moore_S.pis')

levels(Moore_S.pis_performers_ED50s$Geno)
Moore_S.pis_performers_ED50s$Geno <- factor(Moore_S.pis_performers_ED50s$Geno, levels=c('9','8','11','10','1','5','2','14','6','13'))

#plot
ggplot() +
  geom_line(data = Moore_S.pis_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Moore_S.pis_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Moore_S.pis_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Moore_S.pis") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Moore P.lob geno plot ####

#geom_jitter data
Moore_P.lob_performers %>%
  mutate(Geno = fct_relevel(Geno, '2','6','14','1','3','5','7','13','11','12'))

#geom_line data
Moore_P.lob_Geno_DRCs <- Moore_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Moore_P.lob') 

levels(Moore_P.lob_Geno_DRCs$Geno)
Moore_P.lob_Geno_DRCs$Geno <- factor(Moore_P.lob_Geno_DRCs$Geno, levels=c('2','6','14','1','3','5','7','13','11','12'))

#geom_vline data
Moore_P.lob_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Moore_P.lob')

levels(Moore_P.lob_performers_ED50s$Geno)
Moore_P.lob_performers_ED50s$Geno <- factor(Moore_P.lob_performers_ED50s$Geno, levels=c('2','6','14','1','3','5','7','13','11','12'))

#plot
ggplot() +
  geom_line(data = Moore_P.lob_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Moore_P.lob_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Moore_P.lob_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,41), breaks=c(28,30,32,34,36,38,40)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Moore_P.lob") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Moore P.ver geno plot ####

#geom_jitter data
Moore_P.ver_performers %>%
  mutate(Geno = fct_relevel(Geno, '10','15','12','9','16','13','2','4','3','14'))

#geom_line data
Moore_P.ver_Geno_DRCs <- Moore_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Moore_P.ver') 

levels(Moore_P.ver_Geno_DRCs$Geno)
Moore_P.ver_Geno_DRCs$Geno <- factor(Moore_P.ver_Geno_DRCs$Geno, levels=c('10','15','12','9','16','13','2','4','3','14'))

#geom_vline data
Moore_P.ver_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Moore_P.ver')

levels(Moore_P.ver_performers_ED50s$Geno)
Moore_P.ver_performers_ED50s$Geno <- factor(Moore_P.ver_performers_ED50s$Geno, levels=c('10','15','12','9','16','13','2','4','3','14'))

#plot
ggplot() +
  geom_line(data = Moore_P.ver_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Moore_P.ver_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Moore_P.ver_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Moore_P.ver") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Davies A.hya geno plot ####

#geom_jitter data
Davies_A.hya_performers %>%
  mutate(Geno = fct_relevel(Geno, '7','11','6','8','5','2','13','12','14','10'))

#geom_line data
Davies_A.hya_Geno_DRCs <- Davies_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Davies_A.hya') 

levels(Davies_A.hya_Geno_DRCs$Geno)
Davies_A.hya_Geno_DRCs$Geno <- factor(Davies_A.hya_Geno_DRCs$Geno, levels=c('7','11','6','8','5','2','13','12','14','10'))

#geom_vline data
Davies_A.hya_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Davies_A.hya')

levels(Davies_A.hya_performers_ED50s$Geno)
Davies_A.hya_performers_ED50s$Geno <- factor(Davies_A.hya_performers_ED50s$Geno, levels=c('7','11','6','8','5','2','13','12','14','10'))

#plot
ggplot() +
  geom_line(data = Davies_A.hya_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Davies_A.hya_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Davies_A.hya_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Davies_A.hya") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Davies S.pis geno plot ####

#geom_jitter data
Davies_S.pis_performers %>%
  mutate(Geno = fct_relevel(Geno, '10','1','9','8','3','5','6','12','14','13'))

#geom_line data
Davies_S.pis_Geno_DRCs <- Davies_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Davies_S.pis') 

levels(Davies_S.pis_Geno_DRCs$Geno)
Davies_S.pis_Geno_DRCs$Geno <- factor(Davies_S.pis_Geno_DRCs$Geno, levels=c('10','1','9','8','3','5','6','12','14','13'))

#geom_vline data
Davies_S.pis_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Davies_S.pis')

levels(Davies_S.pis_performers_ED50s$Geno)
Davies_S.pis_performers_ED50s$Geno <- factor(Davies_S.pis_performers_ED50s$Geno, levels=c('10','1','9','8','3','5','6','12','14','13'))

#plot
ggplot() +
  geom_line(data = Davies_S.pis_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Davies_S.pis_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Davies_S.pis_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Davies_S.pis") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Davies P.lob geno plot ####

#geom_jitter data
Davies_P.lob_performers %>%
  mutate(Geno = fct_relevel(Geno, '10','11','9','1','5','7','14','8','12','13'))

#geom_line data
Davies_P.lob_Geno_DRCs <- Davies_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Davies_P.lob') 

levels(Davies_P.lob_Geno_DRCs$Geno)
Davies_P.lob_Geno_DRCs$Geno <- factor(Davies_P.lob_Geno_DRCs$Geno, levels=c('10','11','9','1','5','7','14','8','12','13'))

#geom_vline data
Davies_P.lob_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Davies_P.lob')

levels(Davies_P.lob_performers_ED50s$Geno)
Davies_P.lob_performers_ED50s$Geno <- factor(Davies_P.lob_performers_ED50s$Geno, levels=c('10','11','9','1','5','7','14','8','12','13'))

#plot
ggplot() +
  geom_line(data = Davies_P.lob_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Davies_P.lob_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Davies_P.lob_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Davies_P.lob") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Davies P.ver geno plot ####

#geom_jitter data
Davies_P.ver_performers %>%
  mutate(Geno = fct_relevel(Geno, '4','9','3','7','1','10','14','13','15','16'))

#geom_line data
Davies_P.ver_Geno_DRCs <- Davies_Geno_DRCs %>%
  select(Geno,Site_Species,results,Temp)  %>%
  filter(Site_Species=='Davies_P.ver') 

levels(Davies_P.ver_Geno_DRCs$Geno)
Davies_P.ver_Geno_DRCs$Geno <- factor(Davies_P.ver_Geno_DRCs$Geno, levels=c('4','9','3','7','1','10','14','13','15','16'))

#geom_vline data
Davies_P.ver_performers_ED50s<- Top_bottom_performers %>%
  filter(Site_Species=='Davies_P.ver')

levels(Davies_P.ver_performers_ED50s$Geno)
Davies_P.ver_performers_ED50s$Geno <- factor(Davies_P.ver_performers_ED50s$Geno, levels=c('4','9','3','7','1','10','14','13','15','16'))

#plot
ggplot() +
  geom_line(data = Davies_P.ver_Geno_DRCs, aes(x = Temp, y = results, color = Geno), size = 1) +
  geom_jitter(data = Davies_P.ver_performers, aes(x = Temp, y = PAM, color = Geno), size = 1, width = 0.15) +
  geom_vline(data = Davies_P.ver_performers_ED50s, aes(xintercept = ED50, color = Geno), show.legend = FALSE, size = 1, alpha = 0.75) +
  scale_x_continuous(limits=c(28,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  ggtitle("Davies_P.ver") +
  scale_color_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

#### Combine ED50 data plus predict curves from models for population-level plotting ####

Coeff_means<-data.frame(NthDirection_A.hya_coeff_mean,NthDirection_S.pis_coeff_mean,NthDirection_P.lob_coeff_mean,NthDirection_P.ver_coeff_mean,Moore_A.hya_coeff_mean,Moore_S.pis_coeff_mean,Moore_P.lob_coeff_mean,
                        Moore_P.ver_coeff_mean,Davies_A.hya_coeff_mean,Davies_S.pis_coeff_mean,Davies_P.lob_coeff_mean,Davies_P.ver_coeff_mean)
Coeff_lowers<-data.frame(NthDirection_A.hya_coeff_lower,NthDirection_S.pis_coeff_lower,NthDirection_P.lob_coeff_lower,NthDirection_P.ver_coeff_lower,Moore_A.hya_coeff_lower,Moore_S.pis_coeff_lower,Moore_P.lob_coeff_lower,
                         Moore_P.ver_coeff_lower,Davies_A.hya_coeff_lower,Davies_S.pis_coeff_lower,Davies_P.lob_coeff_lower,Davies_P.ver_coeff_lower)
Coeff_uppers<-data.frame(NthDirection_A.hya_coeff_upper,NthDirection_S.pis_coeff_upper,NthDirection_P.lob_coeff_upper,NthDirection_P.ver_coeff_upper,Moore_A.hya_coeff_upper,Moore_S.pis_coeff_upper,Moore_P.lob_coeff_upper,
                         Moore_P.ver_coeff_upper,Davies_A.hya_coeff_upper,Davies_S.pis_coeff_upper,Davies_P.lob_coeff_upper,Davies_P.ver_coeff_upper)

NthDirection_A.hya_preddata = data.frame(temp = seq(29,38, length.out = 100))
NthDirection_A.hya_pred = as.data.frame(predict(NthDirection_A.hya_pop, newdata = NthDirection_A.hya_preddata, interval = 'confidence'))
NthDirection_A.hya_preddata = data.frame(NthDirection_A.hya_preddata, fvfm = NthDirection_A.hya_pred$Prediction, Lower = NthDirection_A.hya_pred$Lower, Upper = NthDirection_A.hya_pred$Upper)

NthDirection_S.pis_preddata = data.frame(temp = seq(29,42, length.out = 100))
NthDirection_S.pis_pred = as.data.frame(predict(NthDirection_S.pis_pop, newdata = NthDirection_S.pis_preddata, interval = 'confidence'))
NthDirection_S.pis_preddata = data.frame(NthDirection_S.pis_preddata, fvfm = NthDirection_S.pis_pred$Prediction, Lower = NthDirection_S.pis_pred$Lower, Upper = NthDirection_S.pis_pred$Upper)

NthDirection_P.lob_preddata = data.frame(temp = seq(29,38, length.out = 100))
NthDirection_P.lob_pred = as.data.frame(predict(NthDirection_P.lob_pop, newdata = NthDirection_P.lob_preddata, interval = 'confidence'))
NthDirection_P.lob_preddata = data.frame(NthDirection_P.lob_preddata, fvfm = NthDirection_P.lob_pred$Prediction, Lower = NthDirection_P.lob_pred$Lower, Upper = NthDirection_P.lob_pred$Upper)

NthDirection_P.ver_preddata = data.frame(temp = seq(29,38, length.out = 100))
NthDirection_P.ver_pred = as.data.frame(predict(NthDirection_P.ver_pop, newdata = NthDirection_P.ver_preddata, interval = 'confidence'))
NthDirection_P.ver_preddata = data.frame(NthDirection_P.ver_preddata, fvfm = NthDirection_P.ver_pred$Prediction, Lower = NthDirection_P.ver_pred$Lower, Upper = NthDirection_P.ver_pred$Upper)

Moore_A.hya_preddata = data.frame(temp = seq(29,38, length.out = 100))
Moore_A.hya_pred = as.data.frame(predict(Moore_A.hya_pop, newdata = Moore_A.hya_preddata, interval = 'confidence'))
Moore_A.hya_preddata = data.frame(Moore_A.hya_preddata, fvfm = Moore_A.hya_pred$Prediction, Lower = Moore_A.hya_pred$Lower, Upper = Moore_A.hya_pred$Upper)

Moore_S.pis_preddata = data.frame(temp = seq(29,38, length.out = 100))
Moore_S.pis_pred = as.data.frame(predict(Moore_S.pis_pop, newdata = Moore_S.pis_preddata, interval = 'confidence'))
Moore_S.pis_preddata = data.frame(Moore_S.pis_preddata, fvfm = Moore_S.pis_pred$Prediction, Lower = Moore_S.pis_pred$Lower, Upper = Moore_S.pis_pred$Upper)

Moore_P.lob_preddata = data.frame(temp = seq(29,38, length.out = 100))
Moore_P.lob_pred = as.data.frame(predict(Moore_P.lob_pop, newdata = Moore_P.lob_preddata, interval = 'confidence'))
Moore_P.lob_preddata = data.frame(Moore_P.lob_preddata, fvfm = Moore_P.lob_pred$Prediction, Lower = Moore_P.lob_pred$Lower, Upper = Moore_P.lob_pred$Upper)

Moore_P.ver_preddata = data.frame(temp = seq(29,38, length.out = 100))
Moore_P.ver_pred = as.data.frame(predict(Moore_P.ver_pop, newdata = Moore_P.ver_preddata, interval = 'confidence'))
Moore_P.ver_preddata = data.frame(Moore_P.ver_preddata, fvfm = Moore_P.ver_pred$Prediction, Lower = Moore_P.ver_pred$Lower, Upper = Moore_P.ver_pred$Upper)

Davies_A.hya_preddata = data.frame(temp = seq(29,38, length.out = 100))
Davies_A.hya_pred = as.data.frame(predict(Davies_A.hya_pop, newdata = Davies_A.hya_preddata, interval = 'confidence'))
Davies_A.hya_preddata = data.frame(Davies_A.hya_preddata, fvfm = Davies_A.hya_pred$Prediction, Lower = Davies_A.hya_pred$Lower, Upper = Davies_A.hya_pred$Upper)

Davies_S.pis_preddata = data.frame(temp = seq(29,38, length.out = 100))
Davies_S.pis_pred = as.data.frame(predict(Davies_S.pis_pop, newdata = Davies_S.pis_preddata, interval = 'confidence'))
Davies_S.pis_preddata = data.frame(Davies_S.pis_preddata, fvfm = Davies_S.pis_pred$Prediction, Lower = Davies_S.pis_pred$Lower, Upper = Davies_S.pis_pred$Upper)

Davies_P.lob_preddata = data.frame(temp = seq(29,38, length.out = 100))
Davies_P.lob_pred = as.data.frame(predict(Davies_P.lob_pop, newdata = Davies_P.lob_preddata, interval = 'confidence'))
Davies_P.lob_preddata = data.frame(Davies_P.lob_preddata, fvfm = Davies_P.lob_pred$Prediction, Lower = Davies_P.lob_pred$Lower, Upper = Davies_P.lob_pred$Upper)

Davies_P.ver_preddata = data.frame(temp = seq(29,38, length.out = 100))
Davies_P.ver_pred = as.data.frame(predict(Davies_P.ver_pop, newdata = Davies_P.ver_preddata, interval = 'confidence'))
Davies_P.ver_preddata = data.frame(Davies_P.ver_preddata, fvfm = Davies_P.ver_pred$Prediction, Lower = Davies_P.ver_pred$Lower, Upper = Davies_P.ver_pred$Upper)

##### Plot all populations/species together or separately, as desired ##### 

#### NthDirection plot ####
levels(NthDirection_data$Site_Species)
NthDirection_data$Site_Species = factor(NthDirection_data$Site_Species,levels(NthDirection_data$Site_Species)[c(4,2,1,3)]) 

NthDirection_plot<- ggplot() +
  geom_jitter(data = NthDirection_data, aes(x = Temp, y = PAM, color = Site_Species), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(28.5,42), breaks=c(28,30,32,34,36,38,40,42)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  geom_line(data = NthDirection_A.hya_preddata, aes(x = temp, y = fvfm), color = 'seagreen4', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_A.hya_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_A.hya_coeff_mean), color = 'seagreen4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_A.hya_coeff_lower, xmax=Coeff_uppers$NthDirection_A.hya_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_A.hya_coeff_mean, digits = 2)), x = 31, y = 0.35, show.legend = FALSE, color = 'seagreen4') +
  
  geom_line(data = NthDirection_S.pis_preddata, aes(x = temp, y = fvfm), color = '#B31E6F', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_S.pis_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#B31E6F', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_S.pis_coeff_mean), color = '#B31E6F', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_S.pis_coeff_lower, xmax=Coeff_uppers$NthDirection_S.pis_coeff_upper, ymin=-Inf, ymax=Inf, fill= '#B31E6F',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_S.pis_coeff_mean, digits = 2)), x = 31, y = 0.30, show.legend = FALSE, color = '#B31E6F') +
  
  geom_line(data = NthDirection_P.lob_preddata, aes(x = temp, y = fvfm), color = 'royalblue4', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_P.lob_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'royalblue4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_P.lob_coeff_mean), color = 'royalblue4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_P.lob_coeff_lower, xmax=Coeff_uppers$NthDirection_P.lob_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'royalblue4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_P.lob_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'royalblue4') +
  
  geom_line(data = NthDirection_P.ver_preddata, aes(x = temp, y = fvfm), color = 'darkorange4', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_P.ver_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'darkorange4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_P.ver_coeff_mean), color = 'darkorange4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_P.ver_coeff_lower, xmax=Coeff_uppers$NthDirection_P.ver_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'darkorange4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_P.ver_coeff_mean, digits = 2)), x = 31, y = 0.20, show.legend = FALSE, color = 'darkorange4') +
  
  ggtitle("NthDirection only") +
  scale_color_manual(values=c('seagreen4','#B31E6F','royalblue4','darkorange4')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()
NthDirection_plot

#### Moore plot ####
levels(Moore_data$Site_Species)
Moore_data$Site_Species = factor(Moore_data$Site_Species,levels(Moore_data$Site_Species)[c(4,2,1,3)]) 

Moore_plot<- ggplot() +
  geom_jitter(data = Moore_data, aes(x = Temp, y = PAM, color = Site_Species), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(28.5,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  geom_line(data = Moore_A.hya_preddata, aes(x = temp, y = fvfm), color = 'seagreen4', show.legend = FALSE) +
  geom_ribbon(data = Moore_A.hya_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_A.hya_coeff_mean), color = 'seagreen4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_A.hya_coeff_lower, xmax=Coeff_uppers$Moore_A.hya_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_A.hya_coeff_mean, digits = 2)), x = 31, y = 0.35, show.legend = FALSE, color = 'seagreen4') +
  
  geom_line(data = Moore_S.pis_preddata, aes(x = temp, y = fvfm), color = '#B31E6F', show.legend = FALSE) +
  geom_ribbon(data = Moore_S.pis_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#B31E6F', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_S.pis_coeff_mean), color = '#B31E6F', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_S.pis_coeff_lower, xmax=Coeff_uppers$Moore_S.pis_coeff_upper, ymin=-Inf, ymax=Inf, fill= '#B31E6F',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_S.pis_coeff_mean, digits = 2)), x = 31, y = 0.30, show.legend = FALSE, color = '#B31E6F') +
  
  geom_line(data = Moore_P.lob_preddata, aes(x = temp, y = fvfm), color = 'royalblue4', show.legend = FALSE) +
  geom_ribbon(data = Moore_P.lob_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'royalblue4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_P.lob_coeff_mean), color = 'royalblue4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_P.lob_coeff_lower, xmax=Coeff_uppers$Moore_P.lob_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'royalblue4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_P.lob_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'royalblue4') +
  
  geom_line(data = Moore_P.ver_preddata, aes(x = temp, y = fvfm), color = 'darkorange4', show.legend = FALSE) +
  geom_ribbon(data = Moore_P.ver_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'darkorange4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_P.ver_coeff_mean), color = 'darkorange4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_P.ver_coeff_lower, xmax=Coeff_uppers$Moore_P.ver_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'darkorange4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_P.ver_coeff_mean, digits = 2)), x = 31, y = 0.20, show.legend = FALSE, color = 'darkorange4') +
  
  ggtitle("Moore only") +
  scale_color_manual(values=c('seagreen4','#B31E6F','royalblue4','darkorange4')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()
Moore_plot

#### Davies plot ####
levels(Davies_data$Site_Species)
Davies_data$Site_Species = factor(Davies_data$Site_Species,levels(Davies_data$Site_Species)[c(4,2,1,3)]) 

Davies_plot<- ggplot() +
  geom_jitter(data = Davies_data, aes(x = Temp, y = PAM, color = Site_Species), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(28.5,39), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  geom_line(data = Davies_A.hya_preddata, aes(x = temp, y = fvfm), color = 'seagreen4', show.legend = FALSE) +
  geom_ribbon(data = Davies_A.hya_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_A.hya_coeff_mean), color = 'seagreen4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_A.hya_coeff_lower, xmax=Coeff_uppers$Davies_A.hya_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_A.hya_coeff_mean, digits = 2)), x = 31, y = 0.35, show.legend = FALSE, color = 'seagreen4') +
  
  geom_line(data = Davies_S.pis_preddata, aes(x = temp, y = fvfm), color = '#B31E6F', show.legend = FALSE) +
  geom_ribbon(data = Davies_S.pis_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#B31E6F', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_S.pis_coeff_mean), color = '#B31E6F', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_S.pis_coeff_lower, xmax=Coeff_uppers$Davies_S.pis_coeff_upper, ymin=-Inf, ymax=Inf, fill= '#B31E6F',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_S.pis_coeff_mean, digits = 2)), x = 31, y = 0.30, show.legend = FALSE, color = '#B31E6F') +
  
  geom_line(data = Davies_P.lob_preddata, aes(x = temp, y = fvfm), color = 'royalblue4', show.legend = FALSE) +
  geom_ribbon(data = Davies_P.lob_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'royalblue4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_P.lob_coeff_mean), color = 'royalblue4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_P.lob_coeff_lower, xmax=Coeff_uppers$Davies_P.lob_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'royalblue4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_P.lob_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'royalblue4') +
  
  geom_line(data = Davies_P.ver_preddata, aes(x = temp, y = fvfm), color = 'darkorange4', show.legend = FALSE) +
  geom_ribbon(data = Davies_P.ver_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'darkorange4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_P.ver_coeff_mean), color = 'darkorange4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_P.ver_coeff_lower, xmax=Coeff_uppers$Davies_P.ver_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'darkorange4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_P.ver_coeff_mean, digits = 2)), x = 31, y = 0.20, show.legend = FALSE, color = 'darkorange4') +
  
  ggtitle("Davies only") +
  scale_color_manual(values=c('seagreen4','#B31E6F','royalblue4','darkorange4')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()
Davies_plot

###A.hya among sites###

A.hya_sites_plot<- ggplot() +
  geom_jitter(data = A.hya_data, aes(x = Temp, y = PAM, color = Site_Species), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(28.5,38.5), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  geom_line(data = NthDirection_A.hya_preddata, aes(x = temp, y = fvfm), color = 'seagreen4', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_A.hya_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_A.hya_coeff_mean), color = 'seagreen4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_A.hya_coeff_lower, xmax=Coeff_uppers$NthDirection_A.hya_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_A.hya_coeff_mean, digits = 2)), x = 31, y = 0.35, show.legend = FALSE, color = 'seagreen4') +
  
  geom_line(data = Moore_A.hya_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Moore_A.hya_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_A.hya_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_A.hya_coeff_lower, xmax=Coeff_uppers$Moore_A.hya_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_A.hya_coeff_mean, digits = 2)), x = 31, y = 0.30, show.legend = FALSE, color = 'seagreen1') +
  
  geom_line(data = Davies_A.hya_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Davies_A.hya_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_A.hya_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_A.hya_coeff_lower, xmax=Coeff_uppers$Davies_A.hya_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_A.hya_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'seagreen1') +
  
  ggtitle("A.hya only") +
  scale_color_manual(values=c('seagreen4','#B31E6F','royalblue4')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()
A.hya_sites_plot

###S.pis among sites###

S.pis_sites_plot<- ggplot() +
  geom_jitter(data = S.pis_data, aes(x = Temp, y = PAM, color = Site_Species), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(28.5,38.5), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  geom_line(data = NthDirection_S.pis_preddata, aes(x = temp, y = fvfm), color = 'seagreen4', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_S.pis_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_S.pis_coeff_mean), color = 'seagreen4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_S.pis_coeff_lower, xmax=Coeff_uppers$NthDirection_S.pis_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_S.pis_coeff_mean, digits = 2)), x = 31, y = 0.35, show.legend = FALSE, color = 'seagreen4') +
  
  geom_line(data = Moore_S.pis_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Moore_S.pis_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_S.pis_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_S.pis_coeff_lower, xmax=Coeff_uppers$Moore_S.pis_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_S.pis_coeff_mean, digits = 2)), x = 31, y = 0.30, show.legend = FALSE, color = 'seagreen1') +
  
  geom_line(data = Davies_S.pis_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Davies_S.pis_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_S.pis_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_S.pis_coeff_lower, xmax=Coeff_uppers$Davies_S.pis_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_S.pis_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'seagreen1') +
  
  ggtitle("S.pis only") +
  scale_color_manual(values=c('seagreen4','#B31E6F','royalblue4')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()
S.pis_sites_plot

###P.lob among sites###

P.lob_sites_plot<- ggplot() +
  geom_jitter(data = P.lob_data, aes(x = Temp, y = PAM, color = Site_Species), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(28.5,38.5), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  geom_line(data = NthDirection_P.lob_preddata, aes(x = temp, y = fvfm), color = 'seagreen4', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_P.lob_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_P.lob_coeff_mean), color = 'seagreen4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_P.lob_coeff_lower, xmax=Coeff_uppers$NthDirection_P.lob_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_P.lob_coeff_mean, digits = 2)), x = 31, y = 0.35, show.legend = FALSE, color = 'seagreen4') +
  
  geom_line(data = Moore_P.lob_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Moore_P.lob_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_P.lob_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_P.lob_coeff_lower, xmax=Coeff_uppers$Moore_P.lob_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_P.lob_coeff_mean, digits = 2)), x = 31, y = 0.30, show.legend = FALSE, color = 'seagreen1') +
  
  geom_line(data = Davies_P.lob_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Davies_P.lob_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_P.lob_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_P.lob_coeff_lower, xmax=Coeff_uppers$Davies_P.lob_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_P.lob_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'seagreen1') +
  
  ggtitle("P.lob only") +
  scale_color_manual(values=c('seagreen4','#B31E6F','royalblue4')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()
P.lob_sites_plot

###P.ver among sites###

P.ver_sites_plot<- ggplot() +
  geom_jitter(data = P.ver_data, aes(x = Temp, y = PAM, color = Site_Species), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(28.5,38.5), breaks=c(28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  geom_line(data = NthDirection_P.ver_preddata, aes(x = temp, y = fvfm), color = 'seagreen4', show.legend = FALSE) +
  geom_ribbon(data = NthDirection_P.ver_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen4', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = NthDirection_P.ver_coeff_mean), color = 'seagreen4', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$NthDirection_P.ver_coeff_lower, xmax=Coeff_uppers$NthDirection_P.ver_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen4',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(NthDirection_P.ver_coeff_mean, digits = 2)), x = 31, y = 0.35, show.legend = FALSE, color = 'seagreen4') +
  
  geom_line(data = Moore_P.ver_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Moore_P.ver_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Moore_P.ver_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Moore_P.ver_coeff_lower, xmax=Coeff_uppers$Moore_P.ver_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Moore_P.ver_coeff_mean, digits = 2)), x = 31, y = 0.30, show.legend = FALSE, color = 'seagreen1') +
  
  geom_line(data = Davies_P.ver_preddata, aes(x = temp, y = fvfm), color = 'seagreen1', show.legend = FALSE) +
  geom_ribbon(data = Davies_P.ver_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen1', linetype=2, alpha = 0.2) +
  geom_vline(data = Coeff_means, aes(xintercept = Davies_P.ver_coeff_mean), color = 'seagreen1', show.legend = FALSE) +
  annotate("rect", xmin=Coeff_lowers$Davies_P.ver_coeff_lower, xmax=Coeff_uppers$Davies_P.ver_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'seagreen1',  alpha = 0.1) +
  geom_text(data = Coeff_means, aes(label=round(Davies_P.ver_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'seagreen1') +
  
  ggtitle("P.ver only") +
  scale_color_manual(values=c('seagreen4','#B31E6F','royalblue4')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()
P.ver_sites_plot