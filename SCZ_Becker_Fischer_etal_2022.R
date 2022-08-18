#### ANALYSIS OF DATA FOR MANUSCRIPT: Becker, Fischer et al., 2022 ########

## last update: 08.2022 (almaxi@gmail.com)

#### Library
library(lme4)
library(lmerTest)
library(corrplot)
library(sjPlot)
library(sjmisc)
library(knitr)
library(magrittr)
library(sjlabelled)      
library(sjmisc)                                                                                    
library(sjstats) 
library(ggeffects)
library(performance)
library(parameters)
library(betareg)
library(tidyverse)
library(glmmTMB)
library(MKinfer)
library(effectsize)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(emmeans)
library(ggstatsplot)
library(car)
library(lsmeans)
library(see)
library(glmmTMB)

rm(list= ls())

# read in data
setwd("C:/Users/CabezaLab-Berlin1/Google Drive/personal/Mx/SCZ/data/") #/CabezaLab-Berlin1
load('SCZ_Becker_Fischer_etal_2022_dataA.Rdata')
load('SCZ_Becker_Fischer_etal_2022_dataB.Rdata')
load('SCZ_Becker_Fischer_etal_2022_dataC.Rdata')

###################################-
data$group = factor(data$group, levels = c("E-book","2D", "3D" )) # damit die Werte nicht als Zahlen genommen werden
data$group3 = factor(data$group, levels = c("3D","E-book","2D"  )) # damit die Werte nicht als Zahlen genommen werden
data$group2 = ordered(data$group, levels = c("E-book","2D", "3D" )) #factor(data$group, levels = c("E-book","2D", "3D" )) # damit die Werte nicht als Zahlen genommen werden
data$pat_ctrl = as.factor(data$pat_ctrl)

################################################################################-

##### calculate average stars #################

# stars for 2D -> 10 subjects excluded due to stars >1 (4 subjects for >0)
 mean(data[data$Timepoint == 8 & data$group == "2D",]$Spielstand, na.rm = T)
 sd(data[data$Timepoint == 8 & data$group == "2D",]$Spielstand, na.rm = T)
 max(data[data$Timepoint == 8 & data$group == "2D",]$Spielstand, na.rm = T)
# stars for 3D -> 2 personen excluded due to stars >1 
 mean(data[data$Timepoint == 8 & data$group == "3D",]$Spielstand, na.rm = T)
 sd(data[data$Timepoint == 8 & data$group == "3D",]$Spielstand, na.rm = T)
 max(data[data$Timepoint == 8 & data$group == "3D",]$Spielstand, na.rm = T)

# correlate intraclass correlation for PANSS ratings
 raterdata = data.frame(data$Gesamt_score_R1, data$Gesamt_score_R2)
 cor.test(data$Gesamt_score_R1, data$Gesamt_score_R2, method = "spearman")

### fill - FUN variable with respective group mean (5 items)
data$SpielFB2_NK_m_fill = data$SpielFB2_NK_m

  Pat_all_Fun = mean(data[data$pat_ctrl== "patient" &  data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "patient" &  is.na(data$SpielFB2_NK_m)]= Pat_all_Fun
  HP_all_Fun = mean(data[data$pat_ctrl== "control" &  data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "control" &  is.na(data$SpielFB2_NK_m)]= HP_all_Fun
  
  Pat_3D_Fun = mean(data[data$pat_ctrl== "patient" & data$group == "3D" & data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "patient" & data$group == "3D" & is.na(data$SpielFB2_NK_m)]= Pat_3D_Fun
  
  Pat_2D_Fun = mean(data[data$pat_ctrl== "patient" & data$group == "2D" & data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "patient" & data$group == "2D" & is.na(data$SpielFB2_NK_m)]= Pat_2D_Fun
  
  Pat_EB_Fun = mean(data[data$pat_ctrl== "patient" & data$group == "E-book" & data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "patient" & data$group == "E-book" & is.na(data$SpielFB2_NK_m)]= Pat_EB_Fun
  
  HP_3D_Fun = mean(data[data$pat_ctrl== "control" & data$group == "3D" & data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "control" & data$group == "3D" & is.na(data$SpielFB2_NK_m)]= HP_3D_Fun
  
  HP_2D_Fun = mean(data[data$pat_ctrl== "control" & data$group == "2D" & data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "control" & data$group == "2D" & is.na(data$SpielFB2_NK_m)]= HP_2D_Fun
  
  HP_EB_Fun = mean(data[data$pat_ctrl== "control" & data$group == "E-book" & data$Timepoint == 8,]$SpielFB2_NK_m, na.rm =T)
  data$SpielFB2_NK_m_fill[data$pat_ctrl== "control" & data$group == "E-book" & is.na(data$SpielFB2_NK_m)]= HP_EB_Fun


#### check for group differences in terms of fun variable ####
ALL_LMER1 <- lm(SpielFB2_NK_m_fill ~ group3 , data=data[data$Spielstand_neu==0 & data$Timepoint ==8 ,], REML= FALSE)  # REML= FALSE
summary(ALL_LMER1)

  mean(data[data$Spielstand_neu==0 & data$Timepoint ==8  & data$group3 == "3D",]$SpielFB2_NK_m_fill, na.rm = T)
  sd(data[data$Spielstand_neu==0 & data$Timepoint ==8  & data$group3 == "3D",]$SpielFB2_NK_m_fill, na.rm = T)
  
  mean(data[data$Spielstand_neu==0 & data$Timepoint ==8  & data$group3 == "2D",]$SpielFB2_NK_m_fill, na.rm = T)
  sd(data[data$Spielstand_neu==0 & data$Timepoint ==8  & data$group3 == "2D",]$SpielFB2_NK_m_fill, na.rm = T)
  
  mean(data[data$Spielstand_neu==0 & data$Timepoint ==8  & data$group3 == "E-book",]$SpielFB2_NK_m_fill, na.rm = T)
  sd(data[data$Spielstand_neu==0 & data$Timepoint ==8  & data$group3 == "E-book",]$SpielFB2_NK_m_fill, na.rm = T)

#################################################################################-
####################### MAIN ANALYSES: Effects of training intervention on spatial orientation and cognition ####
#################################################################################-

############ 1) PANSS total score #############################################
  XX = data[data$Spielstand_neu!=1 & data$Timepoint != 16 & data$Mrdata != 2 & data$pat_ctrl == "patient" ,]
  qqp(data$Gesamt_score_Rm, "norm" ); 
  
  ALL_LMER1 <- lmer(Gesamt_score_Rm ~ group2+Timepoint + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  ALL_LMER2 <- lmer(Gesamt_score_Rm ~ group2*Timepoint + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  ALL_LMER3 <- lmer(Gesamt_score_Rm ~ group3*Timepoint + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  ALL_LMER4 <- lmer(Gesamt_score_Rm ~ group*Timepoint  + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  anova(ALL_LMER1, ALL_LMER2)
  summary(ALL_LMER2)
  summary(ALL_LMER3)
  summary(ALL_LMER4)
  
  hist(residuals(ALL_LMER2))
  qqnorm(residuals(ALL_LMER2)); qqline(residuals(ALL_LMER2)) #Ausrei_er!!!
  PANSS_plot2 <-  ggpredict(ALL_LMER2, c( "Timepoint","group2")) %>% plot(show.title=F, show.x.title=F)+ theme_classic()

   #plot results
      XX1 = data.frame(rep(0, length(predict(ALL_LMER2))))
      XX1$fit<- predict(ALL_LMER2)
      XX1$condition <- factor(ALL_LMER2@frame$group2, ordered = F, levels = c( "E-book","3D","2D" ))
      XX1$Timepoint  <- ALL_LMER2@frame$Timepoint 
      
    PANSS_plot <- ggplot(XX1, aes(x = Timepoint, y = fit,  color =  condition)) +
        geom_point() +stat_smooth( se= TRUE, na.rm = T) +
        labs(y="PANSS (total)", x= "Time (in months)" )+
        geom_jitter() + theme_classic(base_size = 15)
  
################### 1a) PANSS - GP Subscale #################################  
XX = data[data$Spielstand_neu ==0 & data$Timepoint != 16 & data$Mrdata != 2  ,]

  qqp(data$Sum_kog_Rm, "norm" ); 
  Kog_LMER1 <- lmer(Sum_kog_Rm ~ group2+Timepoint + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  Kog_LMER2 <- lmer(Sum_kog_Rm ~ group2*Timepoint + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  Kog_LMER3 <- lmer(Sum_kog_Rm ~ group3*Timepoint + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  Kog_LMER4 <- lmer(Sum_kog_Rm ~ group*Timepoint  + SpielFB2_NK_m_fill + Olanzapin_Equivalent + (1|ID), data=XX, REML= F)  
  anova(Kog_LMER1, Kog_LMER2)
  summary(Kog_LMER2)
  summary(Kog_LMER3)
  summary(Kog_LMER4)
  hist(residuals(Kog_LMER4))
  qqnorm(residuals(Kog_LMER2)); qqline(residuals(Kog_LMER2)) 
  tab_model(Kog_LMER2)
  PANSS_GS_plot2 <- ggpredict(Kog_LMER2, c( "Timepoint","group2")) %>% plot(show.title=F, show.x.title=F)+ theme_classic()

    #plot results
    XX1 = data.frame(rep(0, length(predict(Kog_LMER2))))
    XX1$fit<- predict(Kog_LMER2)
    XX1$condition <- factor(Kog_LMER2@frame$group2, ordered = F, levels = c( "E-book","3D","2D" ))
    XX1$Timepoint  <- Kog_LMER2@frame$Timepoint 
    
    PANSS_GS_plot <- ggplot(XX1, aes(x = Timepoint, y = fit,  color =  condition)) +
      geom_point() +geom_smooth( se= TRUE, na.rm = T) +
      labs(y="PANSS (GP)", x= "Time (in months)" )+
      geom_jitter() + theme_classic(base_size = 15)


#################  1b) PANSS - Negative Subscale #################################
XX = data[data$Spielstand_neu==0 & data$Timepoint != 16 & data$Mrdata != 2  ,]
  qqp(XX$Sum_neg_Rm, "norm" );  #  schlechter fit
  qqp(log(XX$Sum_neg_Rm), "norm"); hist(log(XX$Sum_neg_Rm))  #besser
  Neg_LMER1 <- lmer(Sum_neg_Rm ~ group2+Timepoint+ SpielFB2_NK_m_fill + Olanzapin_Equivalent +(1|ID),data=XX,  na.action= na.omit, REML= F)  
  Neg_LMER2 <- lmer(Sum_neg_Rm ~ group2*Timepoint+ SpielFB2_NK_m_fill + Olanzapin_Equivalent +(1|ID),data=XX,  na.action= na.omit, REML= F)  
  Neg_LMER3 <- lmer(Sum_neg_Rm ~ group3*Timepoint+ SpielFB2_NK_m_fill + Olanzapin_Equivalent +(1|ID),data=XX,  na.action= na.omit, REML= F)  
  Neg_LMER4 <- lmer(Sum_neg_Rm ~ group*Timepoint + SpielFB2_NK_m_fill + Olanzapin_Equivalent +(1|ID),data=XX,  na.action= na.omit, REML= F)  
  
  anova(Neg_LMER1, Neg_LMER2)
  summary(Neg_LMER2)
  summary(Neg_LMER3)
  summary(Neg_LMER4)

  hist(residuals(Neg_LMER2))
  tab_model(Neg_LMER2)
  qqnorm(residuals(Neg_LMER2)); qqline(residuals(Neg_LMER2))
  PANSS_NEG_plot2 <-  ggpredict(Neg_LMER2, c( "Timepoint","group2")) %>% plot(show.title=F, show.x.title=F)+ theme_classic()

  #plot results (use marginal means from Neg_LMER2 )
  XX1 = data.frame(rep(0, length(predict(Neg_LMER2))))
  XX1$fit<- predict(Neg_LMER2)
  XX1$group <- factor(Neg_LMER2@frame$group2, ordered = F, levels = c( "E-book","3D","2D" ))
  XX1$condition <- factor(XX1$group)
  XX1$Timepoint  <- Neg_LMER2@frame$Timepoint 
  
  PANSS_NEG_plot <- ggplot(XX1, aes(x = Timepoint, y = fit,  color =  condition)) +
    geom_point() +stat_smooth(se= TRUE, na.rm = T) +
    labs(y="PANSS (Neg)", x= "Time (in months)" )+
    geom_jitter() + theme_classic(base_size = 15)

#################  1c) PANSS - Positive Subscale  ################################# 
PX = data[data$Spielstand_neu==0 & data$Timepoint != 16 & data$Mrdata != 2  ,]
  qqp(log(PX$Sum_pos_Rm), "norm"); hist(log(PX$Sum_pos_Rm))  #besser
  Pos_LMER1 <- lmer(Sum_pos_Rm ~ group2+Timepoint + SpielFB2_NK_m+Olanzapin_Equivalent+(1|ID),data=PX, na.action= na.omit)   # REML= FALSE
  Pos_LMER2 <- lmer(Sum_pos_Rm ~ group2*Timepoint + SpielFB2_NK_m+Olanzapin_Equivalent+(1|ID),data=PX, na.action= na.omit)   # REML= FALSE
  anova(Pos_LMER1, Pos_LMER2)
  summary(Pos_LMER2)
  tab_model(Pos_LMER2)
  hist(residuals(Pos_LMER2))
  qqnorm(residuals(Pos_LMER2)); qqline(residuals(Pos_LMER2))


#################  2) RAS #################################
#### RAS total 
  XX = data[data$Spielstand_neu == 0 & data$Timepoint != 17 & data$pat_ctrl == "patient" & data$Mrdata != 2 ,]
  cc= mean(XX[!is.na(XX$RAS_total),]$RAS_total, na.rm = T) - 3*sd(XX[!is.na(XX$RAS_total),]$RAS_total, na.rm = T) 
  dd= mean(XX[!is.na(XX$RAS_total),]$RAS_total, na.rm = T) + 3*sd(XX[!is.na(XX$RAS_total),]$RAS_total, na.rm = T) 
  XX[(XX$RAS_total<= cc |  XX$RAS_total>= dd) &  !is.na(XX$RAS_total),]$RAS_total = NaN
  qqp((XX$RAS_total), "norm"); hist((XX$RAS_total))  #besser

  RAS_LMER1 <- lmer(RAS_total ~ group2+Timepoint+ SpielFB2_NK_m_fill + Olanzapin_Equivalent+(1|ID), data=XX, na.action= na.omit, REML= F)
  RAS_LMER2 <- lmer(RAS_total ~ group2*Timepoint+ SpielFB2_NK_m_fill + Olanzapin_Equivalent+(1|ID), data=XX, na.action= na.omit, REML= F)
  RAS_LMER3 <- lmer(RAS_total ~ group3*Timepoint+ SpielFB2_NK_m_fill + Olanzapin_Equivalent+(1|ID), data=XX, na.action= na.omit, REML= F)
  RAS_LMER4 <- lmer(RAS_total ~ group* Timepoint+ SpielFB2_NK_m_fill + Olanzapin_Equivalent+(1|ID), data=XX, na.action= na.omit, REML= F)
  
  anova(RAS_LMER1, RAS_LMER2)
  summary(RAS_LMER2)
  summary(RAS_LMER3)
  summary(RAS_LMER4)

  tab_model(RAS_LMER2)
  hist(residuals(RAS_LMER2))
  qqnorm(residuals(RAS_LMER2)); qqline(residuals(RAS_LMER2))
  RAS_plot2 <-ggpredict(RAS_LMER2, c( 'Timepoint',"group2" ))%>%plot(show.title=F, show.x.title=F)+theme_classic()
  ggpredict(RAS_LMER2, c("SpielFB2_NK_m_fill", "group2"))%>%plot()+theme_classic()

      #plot results (marginal means from RAS_LMER2)
      XX1 = data.frame(rep(0, length(predict(RAS_LMER2))))
      XX1$fit<- predict(RAS_LMER2)
      XX1$group <- factor(RAS_LMER2@frame$group2, ordered = F)
      XX1$Timepoint  <- RAS_LMER2@frame$Timepoint 
      
      RAS_plot <- ggplot(XX1, aes(x = Timepoint, y = fit,  color =  group)) +
        geom_point() +stat_smooth( se= TRUE, na.rm = T) +
        labs(y="RAS (total)", x= "Time (in months)" )+
        geom_jitter() + theme_classic(base_size = 15)

#### 2a) RAS subscale: Goal & Success Orientation #########
# RAS_GoalSuccessOrient, 
# RAS_WillingnessAskHelp, RAS_NoDominationBySymptoms, RAS_RelianceOnOthers, RAS_persConfidenceHope
  XX = data[data$Spielstand_neu == 0 & data$Timepoint != 17 & data$pat_ctrl == "patient" & data$Mrdata != 2  ,]
  XX$RAS_var = XX$RAS_GoalSuccessOrient 
  qqp((XX$RAS_var), "norm");
  hist((XX$RAS_var))  #besser
  cc= mean(XX[!is.na(XX$RAS_var),]$RAS_var, na.rm = T) - 3*sd(XX[!is.na(XX$RAS_var),]$RAS_var, na.rm = T) 
  dd= mean(XX[!is.na(XX$RAS_var),]$RAS_var, na.rm = T) + 3*sd(XX[!is.na(XX$RAS_var),]$RAS_var, na.rm = T) 
  XX[(XX$RAS_var<= cc  | XX$RAS_var>= dd)  & !is.na(XX$RAS_var),]$RAS_var = NaN
  
  RASvar_LMER1 <- lmer(RAS_var ~ group2+Timepoint+SpielFB2_NK_m_fill +Olanzapin_Equivalent+(1|ID),data=XX, na.action= na.omit, REML= F)   # REML= FALSE
  RASvar_LMER2 <- lmer(RAS_var ~ group2*Timepoint+SpielFB2_NK_m_fill +Olanzapin_Equivalent+(1|ID),data=XX, na.action= na.omit, REML= F)   # REML= FALSE
  hist(residuals(RASvar_LMER2))
  qqnorm(residuals(RASvar_LMER2)); qqline(residuals(RASvar_LMER2))
  anova(RASvar_LMER1, RASvar_LMER2)
  summary( RASvar_LMER2)
  ggpredict(RASvar_LMER2, c( "Timepoint", "group2")) %>% plot(show.title=F, show.x.title=F)+ theme_classic()

    
    #(marginal means from RASvar_LMER2)
    XX1 = data.frame(rep(0, length(predict(RASvar_LMER2))))
    XX1$fit<- predict(RASvar_LMER2)
    XX1$group <- factor(RASvar_LMER2@frame$group2, ordered = F)
    XX1$Timepoint  <- RASvar_LMER2@frame$Timepoint 
    
    RAS_VAR_plot <- ggplot(XX1, aes(x = Timepoint, y = fit,  color =  group)) +
      geom_point() +stat_smooth(method = "lm", se= TRUE, na.rm = T) +
      labs(y="RAS (Success & Goal Orientation)", x= "Time (in months)" )+
      geom_jitter() + theme_classic(base_size = 15)


############### 3) MATRICS ########################################################
  XX = data[ data$Spielstand_neu == 0 & data$Timepoint != 16 & data$Mrdata != 2 ,]
  XX$subMCCB = XX$Neurocog.Comp.AGT
  qqp(XX$subMCCB, "norm" );
  cc= mean(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) - 3*sd(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) 
  dd= mean(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) + 3*sd(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) 
  XX[(XX$subMCCB <= cc  | XX$subMCCB>= dd ) & !is.na(XX$subMCCB),]$subMCCB = NaN
  MCCB_LMER1 <- lmer(subMCCB ~ group2+Timepoint+pat_ctrl + SpielFB2_NK_m_fill+ (1|ID), data=XX, REML= F) 
  MCCB_LMER2 <- lmer(subMCCB ~ group2*Timepoint+pat_ctrl + SpielFB2_NK_m_fill+ (1|ID), data=XX, REML= F) 
  MCCB_LMER3 <- lmer(subMCCB ~ group2*Timepoint*pat_ctrl + SpielFB2_NK_m_fill+ (1|ID), data=XX, REML= F) 
  anova(MCCB_LMER1, MCCB_LMER2,MCCB_LMER3)
  summary(MCCB_LMER2)
  summary(MCCB_LMER3)
  
  hist(residuals(MCCB_LMER2))
  qqnorm(residuals(MCCB_LMER2)); qqline(residuals(MCCB_LMER2)) #Ausrei_er!!!
  tab_model(MCCB_LMER2)

#### 3a) MATRICS subdomain: Sustained Attention ########################
##  AV.AGT,WM.AGT, RPS.AGT, Vrbl.Lrng.AGT, Vis.Lrng.AGT # 

  XX = data[data$Spielstand_neu==0 & data$Timepoint != 16 & data$Mrdata != 2  ,]
  XX$subMCCB = XX$AV.AGT
  qqp(XX$subMCCB, "norm" );
  cc= mean(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) - 3*sd(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) 
  dd= mean(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) + 3*sd(XX[!is.na(XX$subMCCB),]$subMCCB, na.rm = T) 
  XX[(XX$subMCCB<= cc  | XX$subMCCB>= dd  ) & !is.na(XX$subMCCB),]$subMCCB = NaN
  
  MCCBsub_LMER1 <- lmer(subMCCB ~  group2+Timepoint+ pat_ctrl  +SpielFB2_NK_m_fill +(1|ID), data=XX, REML= F) 
  MCCBsub_LMER2 <- lmer(subMCCB ~  group2*Timepoint+ pat_ctrl  +SpielFB2_NK_m_fill +(1|ID), data=XX, REML= F) 
  MCCBsub_LMER3 <- lmer(subMCCB ~  group2*Timepoint* pat_ctrl  +SpielFB2_NK_m_fill +(1|ID), data=XX, REML= F) 
  MCCBsub_LMER4 <- lmer(subMCCB ~  group3*Timepoint+ pat_ctrl  +SpielFB2_NK_m_fill +(1|ID), data=XX, REML= F) 
  MCCBsub_LMER5 <- lmer(subMCCB ~  group* Timepoint+ pat_ctrl  +SpielFB2_NK_m_fill +(1|ID), data=XX, REML= F) 
  
  anova(MCCBsub_LMER1, MCCBsub_LMER2,MCCBsub_LMER3)
  summary(MCCBsub_LMER2)
  summary(MCCBsub_LMER4)
  summary(MCCBsub_LMER5)

  hist(residuals(MCCBsub_LMER2))
  qqnorm(residuals(MCCBsub_LMER2)); qqline(residuals(MCCBsub_LMER2)) #Ausrei_er!!!
  tab_model(MCCBsub_LMER2)
  Attention_plot2 <-ggpredict(MCCBsub_LMER2, c( "Timepoint","group2")) %>% plot(show.title=F, show.x.title=F)+ theme_classic()
  plot_model(MCCB_LMER2, show.values = TRUE)


    #plot results (marginal means from MCCBsub_LMER2)
    XX1 = data.frame(rep(0, length(predict(MCCBsub_LMER2))))
    XX1$fit<- predict(MCCBsub_LMER2)
    XX1$condition <- factor(MCCBsub_LMER2@frame$group2, ordered = F, levels = c( "E-book","3D","2D" ))
    XX1$Timepoint  <- MCCBsub_LMER2@frame$Timepoint 
    
    Attention_plot <- ggplot(XX1, aes(x = Timepoint, y = fit,  color =  condition)) +
      geom_point() +stat_smooth(se= TRUE, na.rm = T) +
      labs(y="MCCB: Attention", x= "Time (in months)" )+
      geom_jitter() + theme_classic(base_size = 15)

    
### 4) Brain Behavior Link: FC & Sustained Attention #################

  XX = FCdata

  XX= FCdata[!is.na(FCdata$AV.AGT_post),]
  XX$FC_m = (XX$FC_rRPFC_rHC_postpre + XX$FC_rRPFC_lHC_postpre+ XX$FC_lFEF_rHC_postpre + XX$FC_lFEF_lHC_postpre)/4
  XX$FC_best = (XX$FC_rRPFC_rHC_postpre + XX$FC_rRPFC_lHC_postpre+ XX$FC_lFEF_rHC_postpre )/3
  
  XX$Attention = scale(XX$AV.AGT_post - XX$AV.AGT_pre)
  XX$FC = XX$FC_best
  XX$Condition[XX$Condition == 1] = "E-book"
  XX$Condition[XX$Condition == 2] = "2D"
  XX$Condition[XX$Condition == 3] = "3D"
  XX$Condition = factor(XX$Condition, levels = c( "E-book","3D","2D"  )) #factor(data$group, levels = c("E-book","2D", "3D" )) # damit die Werte nicht als Zahlen genommen werden
  
  #check relationship with outliers
  cor.test(XX$Attention , XX$FC, method = c("spearman"))

  # extract outliers
  qqp(XX$Attention, "norm" );
  cc= mean(XX[!is.na(XX$Attention ),]$Attention, na.rm = T) - 3*sd(XX[!is.na(XX$Attention ),]$Attention, na.rm = T) 
  dd= mean(XX[!is.na(XX$Attention ),]$Attention, na.rm = T) + 3*sd(XX[!is.na(XX$Attention ),]$Attention, na.rm = T) 
  XX[(XX$Attention<= cc  | XX$Attention>= dd ) & !is.na(XX$Attention ),]$Attention = NaN
  qqp(XX$Attention, "norm" );
  
  qqp(XX$FC, "norm" );
  cc= mean(XX[!is.na(XX$FC),]$FC, na.rm = T) - 3*sd(XX[!is.na(XX$FC),]$FC, na.rm = T) 
  dd= mean(XX[!is.na(XX$FC),]$FC, na.rm = T) + 3*sd(XX[!is.na(XX$FC),]$FC, na.rm = T) 
  XX[(XX$FC<= cc  | XX$FC>= dd)  & !is.na(XX$FC),]$FC = NaN
  qqp(XX$FC, "norm" );

  summary(lm(Attention~FC_best, data = XX))

  # final plot
  YY = XX[!(is.na(XX$FC)| is.na(XX$Attention)), ] 
  AttFCplot=ggplot(YY[YY$Timepoint != 16, ], aes(y = FC, x = Attention,  color =  Condition)) +
    geom_point()  +geom_jitter() + theme(text = element_text(size = 21)) +   
     geom_smooth(method = 'lm')+ theme_classic(base_size = 21) +
    labs(x="MCCB:attention (post-pre)", y= "FC in HC-PFC (post-pre)" )
  
  
############## 5) Tunnel task: Spatial Orientation #####################################################

### create Table X in manuscript: proportion ego to allo in categorization ######-
data[data$anzahl_ego == 99 & !is.na(data$anzahl_ego),]$anzahl_ego = NA
data$anzahl_ego = data$anzahl_ego*10

  hist((data$anzahl_ego))  #besser
  PANSS4 = data[!(is.na(data$anzahl_ego)), ] 
  library(dplyr); library(hrbrthemes)
  
fig_S5 <- data[data$Timepoint == 8,] %>%
    ggplot( aes(x=anzahl_ego, fill=pat_ctrl)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_ipsum() + theme_classic()+
    labs(fill="")+ xlab("Tunnel task: % egocentric responses") 

# rename condition groups
tunneldata$group[tunneldata$group == "1"] <- "E-book"
tunneldata$group[tunneldata$group == "2"] <- "2D"
tunneldata$group[tunneldata$group == "3"] <- "3D"
tunneldata$group_fac = ordered(tunneldata$group, levels = c("E-book", "2D", "3D"))
tunneldata$group = factor(tunneldata$group, levels = c("E-book","2D", "3D" )) # damit die Werte nicht als Zahlen genommen werden
tunneldata$group3 = factor(tunneldata$group, levels = c("3D","E-book","2D" )) # damit die Werte nicht als Zahlen genommen werden

tunneldata$pat_ctrl = as.factor(tunneldata$pat_ctrl)

# redefine timepoints
tunneldata[tunneldata$Timepoint==1 & !is.na(tunneldata$Timepoint),]$Timepoint = 0
tunneldata[tunneldata$Timepoint==3 & !is.na(tunneldata$Timepoint),]$Timepoint = 16
tunneldata[tunneldata$Timepoint==2 & !is.na(tunneldata$Timepoint),]$Timepoint = 8
tunneldata[tunneldata$Timepoint==1.5 & !is.na(tunneldata$Timepoint),]$Timepoint = 4

tunneldata$Difficulty_ord = ordered(tunneldata$Difficulty, levels = c(1,2,3))
tunneldata$Difficulty_fac= as.factor(tunneldata$Difficulty)

tunneldata$ID = as.factor(tunneldata$ID)

#### create a variable
  tunneldata$SpielFB2_NK_m_fill = tunneldata$SpielFB2_NK_m
  
  Pat_3D_Fun = mean(tunneldata[tunneldata$pat_ctrl== "patient" & tunneldata$group == "3D" & tunneldata$Timepoint == 8 & tunneldata$Difficulty == 1,]$SpielFB2_NK_m, na.rm =T)
  tunneldata$SpielFB2_NK_m_fill[tunneldata$pat_ctrl== "patient" & tunneldata$group == "3D" & is.na(tunneldata$SpielFB2_NK_m)]= Pat_3D_Fun
  
  Pat_2D_Fun = mean(tunneldata[tunneldata$pat_ctrl== "patient" & tunneldata$group == "2D" & tunneldata$Timepoint == 8 & tunneldata$Difficulty == 1,]$SpielFB2_NK_m, na.rm =T)
  tunneldata$SpielFB2_NK_m_fill[tunneldata$pat_ctrl== "patient" & tunneldata$group == "2D" & is.na(tunneldata$SpielFB2_NK_m)]= Pat_2D_Fun
  
  Pat_EB_Fun = mean(tunneldata[tunneldata$pat_ctrl== "patient" & tunneldata$group == "E-book" & tunneldata$Timepoint == 8 & tunneldata$Difficulty == 1,]$SpielFB2_NK_m, na.rm =T)
  tunneldata$SpielFB2_NK_m_fill[tunneldata$pat_ctrl== "patient" & tunneldata$group == "E-book" & is.na(tunneldata$SpielFB2_NK_m)]= Pat_EB_Fun
  
  HP_3D_Fun = mean(tunneldata[tunneldata$pat_ctrl== "ctrl" & tunneldata$group == "3D" & tunneldata$Timepoint == 8 & tunneldata$Difficulty == 1,]$SpielFB2_NK_m, na.rm =T)
  tunneldata$SpielFB2_NK_m_fill[tunneldata$pat_ctrl== "ctrl" & tunneldata$group == "3D" & is.na(tunneldata$SpielFB2_NK_m)]= HP_3D_Fun
  
  HP_2D_Fun = mean(tunneldata[tunneldata$pat_ctrl== "ctrl" & tunneldata$group == "2D" & tunneldata$Timepoint == 8 & tunneldata$Difficulty == 1,]$SpielFB2_NK_m, na.rm =T)
  tunneldata$SpielFB2_NK_m_fill[tunneldata$pat_ctrl== "ctrl" & tunneldata$group == "2D" & is.na(tunneldata$SpielFB2_NK_m)]= HP_2D_Fun
  
  HP_EB_Fun = mean(tunneldata[tunneldata$pat_ctrl== "ctrl" & tunneldata$group == "E-book" & tunneldata$Timepoint == 8 & tunneldata$Difficulty == 1,]$SpielFB2_NK_m, na.rm =T)
  tunneldata$SpielFB2_NK_m_fill[tunneldata$pat_ctrl== "ctrl" & tunneldata$group == "E-book" & is.na(tunneldata$SpielFB2_NK_m)]= HP_EB_Fun

################################## calculate condition*time interaction #############-

XX = tunneldata[ tunneldata$Spielstand_neu==0 & tunneldata$Timepoint != 16 & tunneldata$Mrdata != 2  ,]

XX$tunnel = XX$diff_min_e_or_a_turn
qqp(log(XX$tunnel), "norm" );
cc1= mean(XX[!is.na(XX$tunnel)  ,]$tunnel, na.rm = T) - 3*sd(XX[!is.na(XX$tunnel) ,]$tunnel, na.rm = T) 
dd1= mean(XX[!is.na(XX$tunnel)  ,]$tunnel, na.rm = T) + 3*sd(XX[!is.na(XX$tunnel) ,]$tunnel, na.rm = T) 
XX[(XX$tunnel<= cc1  | XX$tunnel>= dd1 ) & !is.na(XX$tunnel ) ,]$tunnel = NaN

diff_turn123_LMER1 <- lmer(log(tunnel) ~ Timepoint+group_fac+ pat_ctrl +Difficulty_ord  +SpielFB2_NK_m_fill + (1|ID), data=XX, REML= F) 
  hist(residuals(diff_turn123_LMER1))
diff_turn123_LMER2 <- lmer(log(tunnel) ~ Timepoint*group_fac+ pat_ctrl +Difficulty_ord  +SpielFB2_NK_m_fill + (1|ID), data=XX, REML= F) 
  hist(residuals(diff_turn123_LMER2))
diff_turn123_LMER3 <- lmer(log(tunnel) ~ Timepoint*group_fac  *Difficulty_ord +pat_ctrl +SpielFB2_NK_m_fill + (1|ID), data=XX, REML= F) 
  hist(residuals(diff_turn123_LMER3))
diff_turn123_LMER4 <- lmer(log(tunnel) ~ Timepoint*group_fac *Difficulty_ord * pat_ctrl +SpielFB2_NK_m_fill + (1|ID), data=XX, REML= F) 
diff_turn123_LMER5 <- lmer(log(tunnel) ~ Timepoint*group3+ pat_ctrl+ Difficulty_ord  +SpielFB2_NK_m_fill + (1|ID), data=XX, REML= F) 
diff_turn123_LMER6 <- lmer(log(tunnel) ~ Timepoint*group+ pat_ctrl + Difficulty_ord +SpielFB2_NK_m_fill + (1|ID),  data=XX, REML= F) 

anova(diff_turn123_LMER1, diff_turn123_LMER2,diff_turn123_LMER3,diff_turn123_LMER4)
summary(diff_turn123_LMER2)
summary(diff_turn123_LMER5)
summary(diff_turn123_LMER6)
Tunnel_plot2 <- ggpredict(diff_turn123_LMER2, c("Timepoint","group_fac" )) %>% plot(show.title=F, show.x.title=F)+theme_classic()


  #plot results
  XX1 = data.frame(rep(0, length(predict(diff_turn123_LMER2))))
  XX1$fit<- predict(diff_turn123_LMER2)
  XX1$condition <- factor(diff_turn123_LMER2@frame$group_fac, ordered = F, levels = c( "E-book","3D","2D" ))
  XX1$Timepoint  <- diff_turn123_LMER2@frame$Timepoint 
  XX1$Difficulty_ord  <- diff_turn123_LMER2@frame$Difficulty_ord
  XX1 = XX1[XX1$Difficulty_ord == "3",]
  
  Tunnel_plot1 <- ggplot(XX1, aes(x = Timepoint, y = fit,  color =  condition)) + #, shape = Difficulty_ord
    geom_point() +stat_smooth( se= TRUE, na.rm = T) +
    labs(y="Spatial Orientation", x= "Time (in months)" )+
    geom_jitter() + theme_classic(base_size = 15)


############# plot all behavioral analyses together #####################
#### plot of pictures together ########-
library(ggpubr)

# for only 8 weeks
fig_1 <- ggarrange(Tunnel_plot2, Attention_plot2, PANSS_plot2, PANSS_GS_plot2, PANSS_NEG_plot2, RAS_plot2,
                  common.legend = T, legend = "bottom",
                  #labels = c("A", "B", "C", ),
                  ncol = 3, nrow = 2)

#for all 16 weeks
fig_S6 <- ggarrange(Tunnel_plot, Attention_plot, PANSS_plot, PANSS_GS_plot, PANSS_NEG_plot,RAS_plot,
                  common.legend = T, legend = "bottom",
                  ncol = 3, nrow = 2)

# tables for behavioral analyses (MCCB, Spatial Orientation, PANSS, RAS)
table_S3a <- tab_model( MCCBsub_LMER2,show.intercept = T, show.est =F, show.std = "std") 
table_S3b <- tab_model( diff_turn123_LMER2,show.intercept = T, show.est =F, show.std = "std") 
table_S4u5 <- tab_model(ALL_LMER2,Kog_LMER2, Neg_LMER2,  show.intercept = T, show.est =F, show.std = "std") # Neg_LMER2
table_S36 <- tab_model( RAS_LMER2, RASvar_LMER2,show.intercept = T, show.est =F, show.std = "std")

###############################################################################-
############# SUPPLEMENTARY ANALYSES ###########################################

############################ S1) Blood biomarkers (BDNF) ###############################
#BDNF_Conc_ng_per_mL

  XX = data[data$Spielstand_neu==0 & data$Timepoint != 16 & data$Mrdata != 2,]
  XX$X = XX$BDNF_Conc_ng_per_mL
  qqp(log(XX$X), "norm" );
  
  cc= median(XX[!is.na(XX$X),]$X, na.rm = T) - 3*sd(XX[!is.na(XX$X),]$X, na.rm = T) 
  dd= median(XX[!is.na(XX$X),]$X, na.rm = T) + 3*sd(XX[!is.na(XX$X),]$X, na.rm = T)   
  XX[(XX$X<= cc  | XX$X>= dd ) & !is.na(XX$X),]$X = NaN
  
  qqp(log(XX$X), "norm"); 
  hist(log(data$X))  
  BDNF_LMER1 <- glmmTMB(log(X) ~ group2+Timepoint+ pat_ctrl  +SpielFB2_NK_m_fill +(1|ID)  , data=XX, REML= F) 
  BDNF_LMER2 <- glmmTMB(log(X) ~ group2*Timepoint+ pat_ctrl  +SpielFB2_NK_m_fill +(1|ID)  , data=XX, REML= F) 
  BDNF_LMER3 <- glmmTMB(log(X) ~ group2*Timepoint* pat_ctrl  +SpielFB2_NK_m_fill +(1|ID)  , data=XX, REML= F) 
  anova(BDNF_LMER1, BDNF_LMER2,BDNF_LMER3)
  summary(BDNF_LMER1)
  plot(check_distribution(BDNF_LMER1))
  tab_model(BDNF_LMER1)
  hist(residuals(BDNF_LMER3))
  qqnorm(residuals(BDNF_LMER3)); qqline(residuals(BDNF_LMER3))


################################################################################-
#### S2) Magnetic Resonance Spectroscopy (MRS) ??? Glutamate ##########################
#	Gln_conc	
XX = data[data$Glu_SD<40  & data$Spielstand_neu==0 & data$Timepoint != 16 & data$Mrdata != 2 ,]

  # calculate fit deviation
  mean(XX$Glu_SD, na.rm = T )
  sd(XX$Glu_SD, na.rm = T )
  
  # clean extreme values
  XX[is.na(XX$Olanzapin_Equivalent),]$Olanzapin_Equivalent = 0
  XX$Neurotransmitter = XX$Glu_conc
  XX$NeurotransmitterSD = XX$Glu_SD
  qqp(XX$Neurotransmitter, "norm" );
  cc= median(XX[!is.na(XX$Neurotransmitter),]$Neurotransmitter, na.rm = T) - 3*sd(XX[!is.na(XX$Neurotransmitter),]$Neurotransmitter, na.rm = T) 
  dd= median(XX[!is.na(XX$Neurotransmitter),]$Neurotransmitter, na.rm = T) + 3*sd(XX[!is.na(XX$Neurotransmitter),]$Neurotransmitter, na.rm = T) 
  XX[(XX$Neurotransmitter<= cc  | XX$Neurotransmitter>= dd ) & !is.na(XX$Neurotransmitter),]$Neurotransmitter = NaN
  qqp(log(XX$Neurotransmitter), "norm" );  
  hist(log(XX$Neurotransmitter))  

# condition*time IA -> ns
  NT_LMER1 <- lmer(log(Neurotransmitter) ~ Timepoint+group2+ pat_ctrl +scale(SpielFB2_NK_m_fill) +(1|ID), data=XX, REML= F) 
  NT_LMER2 <- lmer(log(Neurotransmitter) ~ Timepoint*group2+ pat_ctrl +scale(SpielFB2_NK_m_fill) +(1|ID), data=XX, REML= F) 
  NT_LMER3 <- lmer(log(Neurotransmitter) ~ Timepoint*group2* pat_ctrl +scale(SpielFB2_NK_m_fill) +(1|ID), data=XX, REML= F) 
  plot(check_distribution(NT_LMER1))
  
  anova(NT_LMER1, NT_LMER2, NT_LMER3)
  summary(NT_LMER2)
  tab_model(NT_LMER1,  string.std.stat = "std. Statistic")
  hist(residuals(NT_LMER2))
  qqnorm(residuals(NT_LMER2)); qqline(residuals(NT_LMER2))
  ggpredict(NT_LMER2 , c( 'group2', 'Timepoint', 'pat_ctrl' )) %>% plot() + ggplot2::theme_classic()
  
      PANSS4 = XX[!(is.na(XX$group)| is.na(XX$Neurotransmitter)), ] 
      ggplot(PANSS4, aes(x = Timepoint, y = Neurotransmitter,  color =  group)) +
        geom_point() +theme_classic()+
        geom_jitter() + #geom_line()+
        stat_smooth( se= TRUE) +
        labs(y="Neurotransmitter", x= "Time (in weeks)" )
    
# correlate Glutamate with brain/behavior independent of condition
    
    # put data in short format
    XX = data[data$Spielstand_neu==0 ,] 
    XX1= reshape(XX, idvar = "ID", timevar = "Timepoint", direction = "wide")
    
    # calculate difference values
    XX1$Glu_SD_m = (XX1$Glu_SD.0+XX1$Glu_SD.8)/2 
    XX1$Glu_conc_diff = XX1$Glu_conc.8 - XX1$Glu_conc.0
    XX1$AV.AGT_diff = XX1$AV.AGT.8 - XX1$AV.AGT.0

    XX1$FC_rHC_lRPFC_diff=  XX1$X.connectivity.between.atlas.Hippocampus.r.and.networks.Salience.RPFC..L....32.45.27..at.post.0  -   XX1$X.connectivity.between.atlas.Hippocampus.r.and.networks.Salience.RPFC..L....32.45.27..at.pre.0                          
    XX1$FC_lHC_lRPFC_diff=XX1$X.connectivity.between.atlas.Hippocampus.l.and.networks.Salience.RPFC..L....32.45.27..at.post.0  - XX1$X.connectivity.between.atlas.Hippocampus.l.and.networks.Salience.RPFC..L....32.45.27..at.pre.0                         
    XX1$FC_rRPFC_lFEF_diff=XX1$X.connectivity.between.networks.DorsalAttention.FEF..L....27..9.64..and.networks.Salience.RPFC..R...32.46.27..at.pre.0  - XX1$X.connectivity.between.networks.DorsalAttention.FEF..L....27..9.64..and.networks.Salience.RPFC..R...32.46.27..at.pre.0 
    XX1$FC_lHC_rRPFC_diff=XX1$X.connectivity.between.atlas.Hippocampus.l.and.networks.Salience.RPFC..R...32.46.27..at.post.0  - XX1$X.connectivity.between.atlas.Hippocampus.l.and.networks.Salience.RPFC..R...32.46.27..at.pre.0                           
    XX1$FC_lHC_ACC_diff=XX1$X.connectivity.between.atlas.Hippocampus.l.and.networks.Salience.ACC..0.22.35..at.post.0 - XX1$X.connectivity.between.atlas.Hippocampus.l.and.networks.Salience.ACC..0.22.35..at.pre.0                                  
    XX1$FC_lLPFC_lRPFC_diff=XX1$X.connectivity.between.networks.FrontoParietal.LPFC..L....43.33.28..and.networks.Salience.RPFC..L....32.45.27..at.post.0  - XX1$X.connectivity.between.networks.FrontoParietal.LPFC..L....43.33.28..and.networks.Salience.RPFC..L....32.45.27..at.pre.0  
    XX1$FC_rHC_rRPFC_diff=XX1$X.connectivity.between.atlas.Hippocampus.r.and.networks.Salience.RPFC..R...32.46.27..at.post.0  - XX1$X.connectivity.between.atlas.Hippocampus.r.and.networks.Salience.RPFC..R...32.46.27..at.pre.0                           
    XX1$FC_MPFC_lRPFC_diff=XX1$X.connectivity.between.networks.DefaultMode.MPFC..1.55..3..and.networks.Salience.RPFC..L....32.45.27..at.post.0  - XX1$X.connectivity.between.networks.DefaultMode.MPFC..1.55..3..and.networks.Salience.RPFC..L....32.45.27..at.pre.0           
    
    # calculate mean values of HC-PFC FC
    XX1$FC_m = (XX1$FC_rHC_lRPFC_diff + XX1$FC_lHC_lRPFC_diff+ XX1$FC_rRPFC_lFEF_diff + XX1$FC_lHC_rRPFC_diff + XX1$FC_lHC_ACC_diff +  XX1$FC_lLPFC_lRPFC_diff +XX1$FC_rHC_rRPFC_diff + XX1$FC_MPFC_lRPFC_diff )/8

    XX1$FC = XX1$FC_m
    qqp((XX1$FC), "norm" );
    dd= median(XX1[!is.na(XX1$FC),]$FC, na.rm = T) + 3*sd(XX1[!is.na(XX1$FC),]$FC, na.rm = T)   #  schlechter fit
    XX1[ XX1$FC>= dd  & !is.na(XX1$FC),]$FC = NaN
    qqp((XX1$FC), "norm");
    
# change in glutamate predicting change in HC-PFC FC  
    
    # first calculate residuals
    FC_LMER3a <- glmmTMB(FC ~ pat_ctrl.0   +scale( SpielFB2_NK_m_fill.8 ), data=XX1[XX1$Glu_SD_m <40,], REML= FALSE)  # REML= FALSE  # REML= FALSE
    FC_LMER3b <- glmmTMB(FC ~ Glu_conc_diff+pat_ctrl.0  +scale( SpielFB2_NK_m_fill.8 ), data=XX1[XX1$Glu_SD_m <40,], REML= FALSE)  # REML= FALSE  # REML= FALSE
    FC_LMER3c <- glmmTMB(FC ~ Glu_conc_diff*pat_ctrl.0  +scale( SpielFB2_NK_m_fill.8 ), data=XX1[XX1$Glu_SD_m <40,], REML= FALSE)  # REML= FALSE  # REML= FALSE
    anova(FC_LMER3a, FC_LMER3b,FC_LMER3c) ;  hist(residuals(FC_LMER3b)); summary(FC_LMER3b)#p=.011

    YY = data.frame(rep('x',length( FC_LMER3b$frame$FC)))
    YY$Y= FC_LMER3b$frame$FC #
    YY$X  =  FC_LMER3b$frame$Glu_conc_diff
    YY$group = FC_LMER3b$frame$pat_ctrl.0 

    Glu_plot1 <-ggplot(YY, aes(x =X,y=Y, color = group)) + geom_point() +  stat_smooth(method = "lm", se= TRUE) +
                         labs(y="FC in HC-PFC (post-pre)",x= "adjusted glutamate concentration (post-pre)")+ geom_jitter() +
                           theme_classic() +  theme(text = element_text(size=12))
    
    # do final correlation between Glutamate & FC
    require(rcompanion)
    spearmanRho(x=YY$Y, y=YY$X, ci = TRUE, method = "spearman", conf = 0.95, histogram = T, R=10000) #-0.215   -0.412 -0.00235 
    cor.test(YY$Y, YY$X, method = c('spearman'))
    
# change in glutamate predicting change in attention
    qqp((XX1$AV.AGT_diff), "norm" );
    cc= median(XX1[!is.na(XX1$AV.AGT_diff),]$AV.AGT_diff, na.rm = T) - 3*sd(XX1[!is.na(XX1$AV.AGT_diff),]$AV.AGT_diff, na.rm = T) 
    dd= median(XX1[!is.na(XX1$AV.AGT_diff),]$AV.AGT_diff, na.rm = T) + 3*sd(XX1[!is.na(XX1$AV.AGT_diff),]$AV.AGT_diff, na.rm = T)   #  schlechter fit
    XX1[ XX1$AV.AGT_diff>= dd  & !is.na(XX1$AV.AGT_diff),]$AV.AGT_diff = NaN
    qqp((XX1$AV.AGT_diff), "norm");
    
    # first calculate residuals
        FC_LMER4a <- glmmTMB(AV.AGT_diff ~ pat_ctrl.0  +scale( SpielFB2_NK_m_fill.8 ), data=XX1[XX1$Glu_SD_m <40 ,], REML= FALSE)  
    FC_LMER4b <- glmmTMB(AV.AGT_diff  ~ Glu_conc_diff+pat_ctrl.0  +scale( SpielFB2_NK_m_fill.8 ), data=XX1[XX1$Glu_SD_m <40  ,], REML= FALSE)  
    FC_LMER4c <- glmmTMB(AV.AGT_diff  ~ Glu_conc_diff*pat_ctrl.0  +scale( SpielFB2_NK_m_fill.8 ), data=XX1[XX1$Glu_SD_m <40  ,], REML= FALSE)  
    
    anova(FC_LMER4a,FC_LMER4b,FC_LMER4c) ;  hist(residuals(FC_LMER4b)); summary(FC_LMER4b)#p=.011
    summary(FC_LMER4c)
    ZZ = data.frame(rep('x',length(FC_LMER4c$frame$AV.AGT_diff)))
    ZZ$Y= FC_LMER4c$frame$AV.AGT_diff 
    ZZ$X  =  FC_LMER4c$frame$Glu_conc_diff
    ZZ$group = FC_LMER4c$frame$pat_ctrl.0 
    Glu_plot2 <-ggplot(ZZ, aes(x =X,y=Y, color = group)) + geom_point() +  stat_smooth(method = "lm", se= TRUE) +labs(y="MCCB: attention (post-pre)", 
                                                                                                                      x= "adjusted glutamate concentration (post-pre)")+ geom_jitter() +theme_classic() +  theme(text = element_text(size=12))
    # do final correlation between Glutamate & FC
    spearmanRho(x=ZZ$Y, y=ZZ$X, ci = TRUE, method = "spearman", conf = 0.95, histogram = T, R=10000) #-0.215   -0.412 -0.00235 
    cor.test(ZZ$Y, ZZ$X, method = c('spearman'))
    
    # plotting both figures together
    library(ggpubr)
    figX <- ggarrange(Glu_plot1, Glu_plot2,
                      common.legend = T, legend = "bottom",
                      #labels = c("A", "B", "C", ),
                      ncol = 2, nrow = 1)
  

################################################################################-
########################## S3) T1 MRI - brain volume  ##########################

XX = data[data$Spielstand_neu==0 & data$Timepoint != 16 ,]


# Hippocampus
  cc= mean(XX[!is.na(XX$lHip),]$lHip, na.rm = T) - 3*sd(XX[!is.na(XX$lHip),]$lHip, na.rm = T) 
  dd= mean(XX[!is.na(XX$lHip),]$lHip, na.rm = T) + 3*sd(XX[!is.na(XX$lHip),]$lHip, na.rm = T) 
  XX[(XX$lHip<= cc  | XX$lHip>= dd)  & !is.na(XX$lHip),]$lHip = NaN
lHip_LMER <- glmmTMB(lHip ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE
  
  cc= mean(XX[!is.na(XX$rHip),]$rHip, na.rm = T) - 3*sd(XX[!is.na(XX$rHip),]$rHip, na.rm = T) 
  dd= mean(XX[!is.na(XX$rHip),]$rHip, na.rm = T) + 3*sd(XX[!is.na(XX$rHip),]$rHip, na.rm = T) 
  XX[(XX$rHip<= cc  | XX$rHip>= dd)  & !is.na(XX$rHip),]$rHip = NaN
rHip_LMER <- glmmTMB(rHip ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl   +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Anterior Cingulate
  cc= mean(XX[!is.na(XX$lAntCinGy),]$lAntCinGy, na.rm = T) - 3*sd(XX[!is.na(XX$lAntCinGy),]$lAntCinGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$lAntCinGy),]$lAntCinGy, na.rm = T) + 3*sd(XX[!is.na(XX$lAntCinGy),]$lAntCinGy, na.rm = T) 
  XX[(XX$lAntCinGy<= cc  | XX$lAntCinGy>= dd)  & !is.na(XX$lAntCinGy),]$lAntCinGy = NaN
lAntCinGy_LMER <- glmmTMB(lAntCinGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rAntCinGy),]$rAntCinGy, na.rm = T) - 3*sd(XX[!is.na(XX$rAntCinGy),]$rAntCinGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$rAntCinGy),]$rAntCinGy, na.rm = T) + 3*sd(XX[!is.na(XX$rAntCinGy),]$rAntCinGy, na.rm = T) 
  XX[(XX$rAntCinGy<= cc  | XX$rAntCinGy>= dd)  & !is.na(XX$rAntCinGy),]$rAntCinGy = NaN
rAntCinGy_LMER <- glmmTMB(rAntCinGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Frontal Operculum
  cc= mean(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) - 3*sd(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) 
  dd= mean(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) + 3*sd(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) 
  XX[(XX$lFroOpe<= cc  | XX$lFroOpe>= dd)  & !is.na(XX$lFroOpe),]$lFroOpe = NaN
lFroOpe_LMER <- glmmTMB(lFroOpe ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) - 3*sd(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) 
  dd= mean(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) + 3*sd(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) 
  XX[(XX$rFroOpe<= cc  | XX$rFroOpe>= dd)  & !is.na(XX$rFroOpe),]$rFroOpe = NaN
rFroOpe_LMER <- glmmTMB(rFroOpe ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Anterior Insula
  cc= mean(XX[!is.na(XX$lAntIns),]$lAntIns, na.rm = T) - 3*sd(XX[!is.na(XX$lAntIns),]$lAntIns, na.rm = T) 
  dd= mean(XX[!is.na(XX$lAntIns),]$lAntIns, na.rm = T) + 3*sd(XX[!is.na(XX$lAntIns),]$lAntIns, na.rm = T) 
  XX[(XX$lAntIns<= cc  | XX$lAntIns>= dd)  & !is.na(XX$lAntIns),]$lAntIns = NaN
lAntIns_LMER <- glmmTMB(lAntIns ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rAntIns),]$rAntIns, na.rm = T) - 3*sd(XX[!is.na(XX$rAntIns),]$rAntIns, na.rm = T) 
  dd= mean(XX[!is.na(XX$rAntIns),]$rAntIns, na.rm = T) + 3*sd(XX[!is.na(XX$rAntIns),]$rAntIns, na.rm = T) 
  XX[(XX$rAntIns<= cc  | XX$rAntIns>= dd)  & !is.na(XX$rAntIns),]$rAntIns = NaN
rAntIns_LMER <- glmmTMB(rAntIns ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Middle Frontal Gyrus  
  cc= mean(XX[!is.na(XX$lMidFroGy),]$lMidFroGy, na.rm = T) - 3*sd(XX[!is.na(XX$lMidFroGy),]$lMidFroGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$lMidFroGy),]$lMidFroGy, na.rm = T) + 3*sd(XX[!is.na(XX$lMidFroGy),]$lMidFroGy, na.rm = T) 
  XX[(XX$lMidFroGy<= cc  | XX$lMidFroGy>= dd)  & !is.na(XX$lMidFroGy),]$lMidFroGy = NaN
lMidFroGy_LMER <- glmmTMB(lMidFroGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rMidFroGy),]$rMidFroGy, na.rm = T) - 3*sd(XX[!is.na(XX$rMidFroGy),]$rMidFroGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$rMidFroGy),]$rMidFroGy, na.rm = T) + 3*sd(XX[!is.na(XX$rMidFroGy),]$rMidFroGy, na.rm = T) 
  XX[(XX$rMidFroGy<= cc  | XX$rMidFroGy>= dd)  & !is.na(XX$rMidFroGy),]$rMidFroGy = NaN
rMidFroGy_LMER <- glmmTMB(rMidFroGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Inferior Frontal Gyrus
  cc= mean(XX[!is.na(XX$lInfFroOrbGy),]$lInfFroOrbGy, na.rm = T) - 3*sd(XX[!is.na(XX$lInfFroOrbGy),]$lInfFroOrbGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$lInfFroOrbGy),]$lInfFroOrbGy, na.rm = T) + 3*sd(XX[!is.na(XX$lInfFroOrbGy),]$lInfFroOrbGy, na.rm = T) 
  XX[(XX$lInfFroOrbGy<= cc  | XX$lInfFroOrbGy>= dd)  & !is.na(XX$lInfFroOrbGy),]$lInfFroOrbGy = NaN
lInfFroOrbGy_LMER <- glmmTMB(lInfFroOrbGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rInfFroOrbGy),]$rInfFroOrbGy, na.rm = T) - 3*sd(XX[!is.na(XX$rInfFroOrbGy),]$rInfFroOrbGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$rInfFroOrbGy),]$rInfFroOrbGy, na.rm = T) + 3*sd(XX[!is.na(XX$rInfFroOrbGy),]$rInfFroOrbGy, na.rm = T) 
  XX[(XX$rInfFroOrbGy<= cc  | XX$rInfFroOrbGy>= dd)  & !is.na(XX$rInfFroOrbGy),]$rInfFroOrbGy = NaN
rInfFroOrbGy_LMER <- glmmTMB(rInfFroOrbGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Anterior Oribital Gyrus
  cc= mean(XX[!is.na(XX$lAntOrbGy),]$lAntOrbGy, na.rm = T) - 3*sd(XX[!is.na(XX$lAntOrbGy),]$lAntOrbGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$lAntOrbGy),]$lAntOrbGy, na.rm = T) + 3*sd(XX[!is.na(XX$lAntOrbGy),]$lAntOrbGy, na.rm = T) 
  XX[(XX$lAntOrbGy<= cc  | XX$lAntOrbGy>= dd)  & !is.na(XX$lAntOrbGy),]$lAntOrbGy = NaN
lAntOrbGy_LMER <- glmmTMB(lAntOrbGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rAntOrbGy),]$rAntOrbGy, na.rm = T) - 3*sd(XX[!is.na(XX$rAntOrbGy),]$rAntOrbGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$rAntOrbGy),]$rAntOrbGy, na.rm = T) + 3*sd(XX[!is.na(XX$rAntOrbGy),]$rAntOrbGy, na.rm = T) 
  XX[(XX$rAntOrbGy<= cc  | XX$rAntOrbGy>= dd)  & !is.na(XX$rAntOrbGy),]$rAntOrbGy = NaN
rAntOrbGy_LMER <- glmmTMB(rAntOrbGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Frontal Operculum
  cc= mean(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) - 3*sd(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) 
  dd= mean(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) + 3*sd(XX[!is.na(XX$lFroOpe),]$lFroOpe, na.rm = T) 
  XX[(XX$lFroOpe<= cc  | XX$lFroOpe>= dd)  & !is.na(XX$lFroOpe),]$lFroOpe = NaN
lFroOpe_LMER <- glmmTMB(lFroOpe ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) - 3*sd(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) 
  dd= mean(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) + 3*sd(XX[!is.na(XX$rFroOpe),]$rFroOpe, na.rm = T) 
  XX[(XX$rFroOpe<= cc  | XX$rFroOpe>= dd)  & !is.na(XX$rFroOpe),]$rFroOpe = NaN
rFroOpe_LMER <- glmmTMB(rFroOpe ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Medial Frontal Cbr
  cc= mean(XX[!is.na(XX$lMedFroCbr),]$lMedFroCbr, na.rm = T) - 3*sd(XX[!is.na(XX$lMedFroCbr),]$lMedFroCbr, na.rm = T) 
  dd= mean(XX[!is.na(XX$lMedFroCbr),]$lMedFroCbr, na.rm = T) + 3*sd(XX[!is.na(XX$lMedFroCbr),]$lMedFroCbr, na.rm = T) 
  XX[(XX$lMedFroCbr<= cc  | XX$lMedFroCbr>= dd)  & !is.na(XX$lMedFroCbr),]$lMedFroCbr = NaN
lMedFroCbr_LMER <- glmmTMB(lMedFroCbr ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rMedFroCbr),]$rMedFroCbr, na.rm = T) - 3*sd(XX[!is.na(XX$rMedFroCbr),]$rMedFroCbr, na.rm = T) 
  dd= mean(XX[!is.na(XX$rMedFroCbr),]$rMedFroCbr, na.rm = T) + 3*sd(XX[!is.na(XX$rMedFroCbr),]$rMedFroCbr, na.rm = T) 
  XX[(XX$rMedFroCbr<= cc  | XX$rMedFroCbr>= dd)  & !is.na(XX$rMedFroCbr),]$rMedFroCbr = NaN
rMedFroCbr_LMER <- glmmTMB(rMedFroCbr ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

#Superior Medial Frontal Gyrus
  cc= mean(XX[!is.na(XX$lSupMedFroGy),]$lSupMedFroGy, na.rm = T) - 3*sd(XX[!is.na(XX$lSupMedFroGy),]$lSupMedFroGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$lSupMedFroGy),]$lSupMedFroGy, na.rm = T) + 3*sd(XX[!is.na(XX$lSupMedFroGy),]$lSupMedFroGy, na.rm = T) 
  XX[(XX$lSupMedFroGy<= cc  | XX$lSupMedFroGy>= dd)  & !is.na(XX$lSupMedFroGy),]$lSupMedFroGy = NaN
lSupMedFroGy_LMER <- glmmTMB(lSupMedFroGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE
  
  cc= mean(XX[!is.na(XX$rSupMedFroGy),]$rSupMedFroGy, na.rm = T) - 3*sd(XX[!is.na(XX$rSupMedFroGy),]$rSupMedFroGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$rSupMedFroGy),]$rSupMedFroGy, na.rm = T) + 3*sd(XX[!is.na(XX$rSupMedFroGy),]$rSupMedFroGy, na.rm = T) 
  XX[(XX$rSupMedFroGy<= cc  | XX$rSupMedFroGy>= dd)  & !is.na(XX$rSupMedFroGy),]$rSupMedFroGy = NaN
rSupMedFroGy_LMER <- glmmTMB(rSupMedFroGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

#Posterior Insula
  cc= mean(XX[!is.na(XX$lPosIns),]$lPosIns, na.rm = T) - 3*sd(XX[!is.na(XX$lPosIns),]$lPosIns, na.rm = T) 
  dd= mean(XX[!is.na(XX$lPosIns),]$lPosIns, na.rm = T) + 3*sd(XX[!is.na(XX$lPosIns),]$lPosIns, na.rm = T) 
  XX[(XX$lPosIns<= cc  | XX$lPosIns>= dd)  & !is.na(XX$lPosIns),]$lPosIns = NaN
lPosIns_LMER <- glmmTMB(lPosIns ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX )  # REML= FALSE
 
  cc= mean(XX[!is.na(XX$rPosIns),]$rPosIns, na.rm = T) - 3*sd(XX[!is.na(XX$rPosIns),]$rPosIns, na.rm = T) 
  dd= mean(XX[!is.na(XX$rPosIns),]$rPosIns, na.rm = T) + 3*sd(XX[!is.na(XX$rPosIns),]$rPosIns, na.rm = T) 
  XX[(XX$rPosIns<= cc  | XX$rPosIns>= dd)  & !is.na(XX$rPosIns),]$rPosIns = NaN
rPosIns_LMER <- glmmTMB(rPosIns ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Middle Cingulate Gyrus  
  cc= mean(XX[!is.na(XX$lMidCinGy),]$lMidCinGy, na.rm = T) - 3*sd(XX[!is.na(XX$lMidCinGy),]$lMidCinGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$lMidCinGy),]$lMidCinGy, na.rm = T) + 3*sd(XX[!is.na(XX$lMidCinGy),]$lMidCinGy, na.rm = T) 
  XX[(XX$lMidCinGy<= cc  | XX$lMidCinGy>= dd)  & !is.na(XX$lMidCinGy),]$lMidCinGy = NaN
lMidCinGy_LMER <- glmmTMB(lMidCinGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rMidCinGy),]$rMidCinGy, na.rm = T) - 3*sd(XX[!is.na(XX$rMidCinGy),]$rMidCinGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$rMidCinGy),]$rMidCinGy, na.rm = T) + 3*sd(XX[!is.na(XX$rMidCinGy),]$rMidCinGy, na.rm = T) 
  XX[(XX$rMidCinGy<= cc  | XX$rMidCinGy>= dd)  & !is.na(XX$rMidCinGy),]$rMidCinGy = NaN
rMidCinGy_LMER <- glmmTMB(rMidCinGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Frontal Pole 
  cc= mean(XX[!is.na(XX$lFroPo),]$lFroPo, na.rm = T) - 3*sd(XX[!is.na(XX$lFroPo),]$lFroPo, na.rm = T) 
  dd= mean(XX[!is.na(XX$lFroPo),]$lFroPo, na.rm = T) + 3*sd(XX[!is.na(XX$lFroPo),]$lFroPo, na.rm = T) 
  XX[(XX$lFroPo<= cc  | XX$lFroPo>= dd)  & !is.na(XX$lFroPo),]$lFroPo = NaN
lFroPo_LMER <- glmmTMB(lFroPo ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rFroPo),]$rFroPo, na.rm = T) - 3*sd(XX[!is.na(XX$rFroPo),]$rFroPo, na.rm = T) 
  dd= mean(XX[!is.na(XX$rFroPo),]$rFroPo, na.rm = T) + 3*sd(XX[!is.na(XX$rFroPo),]$rFroPo, na.rm = T) 
  XX[(XX$rFroPo<= cc  | XX$rFroPo>= dd)  & !is.na(XX$rFroPo),]$rFroPo = NaN
rFroPo_LMER <- glmmTMB(rFroPo ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl  +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

# Superior Frontal Gyrus 
  cc= mean(XX[!is.na(XX$lSupFroGy),]$lSupFroGy, na.rm = T) - 3*sd(XX[!is.na(XX$lSupFroGy),]$lSupFroGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$lSupFroGy),]$lSupFroGy, na.rm = T) + 3*sd(XX[!is.na(XX$lSupFroGy),]$lSupFroGy, na.rm = T) 
  XX[(XX$lSupFroGy<= cc  | XX$lSupFroGy>= dd)  & !is.na(XX$lSupFroGy),]$lSupFroGy = NaN
lSupFroGy_LMER <- glmmTMB(lSupFroGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE

  cc= mean(XX[!is.na(XX$rSupFroGy),]$rSupFroGy, na.rm = T) - 3*sd(XX[!is.na(XX$rSupFroGy),]$rSupFroGy, na.rm = T) 
  dd= mean(XX[!is.na(XX$rSupFroGy),]$rSupFroGy, na.rm = T) + 3*sd(XX[!is.na(XX$rSupFroGy),]$rSupFroGy, na.rm = T) 
  XX[(XX$rSupFroGy<= cc  | XX$rSupFroGy>= dd)  & !is.na(XX$rSupFroGy),]$rSupFroGy = NaN
rSupFroGy_LMER <- glmmTMB(rSupFroGy ~ scale(GM) + group2*scale(Timepoint)+ pat_ctrl +scale(SpielFB2_NK_m_fill)+(1|ID), data=XX)  # REML= FALSE


# export tables
table_S8a = tab_model(lHip_LMER, rHip_LMER, lInfFroOrbGy_LMER,rInfFroOrbGy_LMER, lAntOrbGy_LMER,rAntOrbGy_LMER,
          lAntCinGy_LMER, rAntCinGy_LMER, lMidCinGy_LMER, rMidCinGy_LMER,lSupMedFroGy_LMER,rSupMedFroGy_LMER, lSupFroGy_LMER,rSupFroGy_LMER, 
          show.intercept = F, show.est =T, show.obs = F, robust = T, show.re.var= F)#p.style = "stars", 

table_S8b  = tab_model(lFroOpe_LMER, rFroOpe_LMER,lFroPo_LMER, rFroPo_LMER, lMidFroGy_LMER, rMidFroGy_LMER, lMedFroCbr_LMER,rMedFroCbr_LMER, 
          lAntIns_LMER, rAntIns_LMER,lPosIns_LMER, rPosIns_LMER,
          show.intercept = F,  show.est =T,show.obs = F, robust = T, show.re.var= F)#p.style = "stars",

# FDR adjust p value for GM group difference 
p=c (.450, .004, .022, .001, .993, .864, .923, .235, .001, .001, .088, .493, .847, .161, .833, .990, .001, .030, .066, .003, .059, .859, .830, .190, .898, .274)	
p.adjust(p, method = c('fdr'), n = length(p))    
