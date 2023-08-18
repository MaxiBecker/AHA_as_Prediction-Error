#### R-Script FOR MANUSCRIPT Becker, Wang & Cabeza, 2023 ############################################
#
# Assumptions: AHA-Dimensions correspond to (meta-cognitive) PE
#
########################################################################-

rm(list = ls())

## libraries #################################################################

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
  library(tidyverse)

## load / prepare data ########################################################

  setwd('C:/Users/Maxi/Google Drive/studies/Maxi/AHAPE/GitHub_4_publication')
  load("BeckerWangCabeza_2023_data.Rdata")
  CRA$Rate_Response[CRA$Rate_Response==7] = NA  #7 = was already solved during rating
  
  CRA1 = CRA[!is.na(CRA$subj_solution_likelihood),]
  min_SSL=min(CRA1$subj_solution_likelihood, na.rm = T)
  max_SSL = max(CRA1$subj_solution_likelihood, na.rm = T)
  CRA1$subj_solution_likelihood_norm = (CRA1$subj_solution_likelihood - min_SSL) /( max_SSL - min_SSL) 
  CRA1$PE = 1-CRA1$subj_solution_likelihood_norm


##############################################################################-
# calculating composite AHA score

 library(lavaan)
 #CRA_all = CRA[CRA$cor2 != 3,] %>% group_by(subject) %>% summarise_all(funs(mean), na.rm =T)   
 
 insightfac <-   'AHA  =~ certain + sudden + surprising + pleasing   
                   AHA ~~ 1*AHA
                  pleasing ~~ surprising
                  '
     #fit
     CRA_AHA = CRA[!is.na(CRA$sudden) & !is.na(CRA$pleasing) ,]
     fit <- lavaan(insightfac, data=CRA, 
                   auto.var=T, auto.fix.first=F ,
                   auto.cov.lv.x=T, estimator = "MLR" )
     summary(fit, fit.measures=TRUE, standardized = T, rsquare = T)
     modindices(fit, sort = TRUE)
     #standardizedSolution(fit)
     CRA_AHA$AHA_facload = predict(fit)
     inspect(fit,what="std")$lambda
     lavInspect(fit, what = "est")$theta
     
############################################################################-
### Descriptive Statistics #######
     
hist(CRA$Rate_Response)
CRA_all = CRA %>% group_by(subject) %>% summarise_all(funs(mean), na.rm =T)   

 CRA_sd= CRA %>%
    group_by(subject) %>%
    summarise_at(c( "sudden", "pleasing", "certain", "surprising", "Rate_Response"), sd, na.rm = TRUE) 

  summary(CRA_all)
  sum(CRA_all$sex)
  sd(CRA_all$age, na.rm = T)

    mean(CRA_all$cor1, na.rm = T)
    sd(CRA_all$cor1, na.rm = T)
    
    mean(CRA_all$cor2, na.rm = T)
    sd(CRA_all$cor2, na.rm = T)
    
    median(CRA_all$RT, na.rm = T)
    sd(CRA_all$RT, na.rm = T)
    
    median(CRA_all$RT_correct, na.rm = T)
    sd(CRA_all$RT_correct, na.rm = T)
    
    mean(CRA_all$surprising, na.rm = T)
    sd(CRA_all$surprising, na.rm = T)
    
    mean(CRA_all$sudden, na.rm = T)
    sd(CRA_all$sudden, na.rm = T)
    
    mean(CRA_all$certain, na.rm = T)
    sd(CRA_all$certain, na.rm = T)
    
    mean(CRA_all$pleasing, na.rm = T)
    sd(CRA_all$pleasing, na.rm = T)


######################################################################## ##
## PE - Actual analysis #############################################   ##
   

## Solution rating predicting accuracy ####################################################################################################################################
  M0_cor2 <- glmmTMB(cor2 ~ blocknum + scale(RT) + (1 | subject) + (1 | cues), data= CRA1[CRA1$RT>2  & !is.na(CRA1$subj_solution_likelihood),], family = binomial(link = "logit"), na.action = na.omit)
  M1_cor2 <- glmmTMB(cor2 ~ blocknum + scale(RT)+subj_solution_likelihood + (1 | subject) + (1 | cues), data= CRA1[CRA1$RT>2,], family = binomial(link = "logit"), na.action = na.omit)

  plot(check_distribution(M1_cor2))
  check_collinearity(M1_cor2)
  anova(M0_cor2,M1_cor2)
  #hist(residuals(M1))
  tab_model(M1_cor2, show.std = T)
  ggpredict(M2_cor2 , c("subj_solution_likelihood", "RT" ),main = "") %>% plot() + ggplot2::theme_classic()

## AHA: PLEASURE ~ PE ###########################################################################################################################
  M0a_pleasure <- glmmTMB(pleasing ~ RT       +blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  M0b_pleasure <- glmmTMB(pleasing ~ RT  +cor2+blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  M1_pleasure <- glmmTMB(pleasing ~ RT+PE+cor2+blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  M2_pleasure <- glmmTMB(pleasing ~ RT+PE*cor2+blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  
  plot(check_distribution(M2_pleasure))
  check_collinearity(M2_pleasure)
  anova(M0a_pleasure,M0b_pleasure,M1_pleasure, M2_pleasure)
  hist(residuals(M1_pleasure))

  summary(M0b_pleasure)
  tab_model(M0b_pleasure)
  tab_model(M1_pleasure)
  tab_model(M2_pleasure)
  
  Pleasure_ggpredict<- ggpredict(M1_pleasure , c( 'PE', "cor2" ))
  Pleasure_finalPlot2= ggplot(Pleasure_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.1,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Pleasure",y = "Pleasure (probability)", fill = "Accuracy")+
    xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))

## AHA: CERTAINTY ~ PE ############################################################################################################
  M0a_certain <-glmmTMB(certain ~ RT         +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M0b_certain <-glmmTMB(certain ~ RT   +cor2 +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M1_certain <- glmmTMB(certain ~ RT+PE+cor2 +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M2_certain <- glmmTMB(certain ~ RT+PE*cor2 +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  
  plot(check_distribution(M1_certain))#residuals are normal
  anova(M0a_certain,M0b_certain,M1_certain, M2_certain)
  summary(M1_certain)
  tab_model(M0b_certain,show.std  = T)
  tab_model(M1_certain,show.std  = T)
  tab_model(M2_certain,show.std  = T)
  
  Certain_ggpredict<- ggpredict(M2_certain , c( 'PE','cor2' ))
  Certain_finalPlot2= ggplot(Certain_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.1,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Certainty",y = "Certainty", fill = "Accuracy")+
    xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)") )+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))

## AHA: SUDDENNESS ~ PE ###################################################################################################################
  M0a_sudden<- glmmTMB(sudden ~RT+                  +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M0b_sudden<- glmmTMB(sudden ~RT+    cor2          +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M1_sudden <- glmmTMB(sudden ~RT+ PE+cor2           +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M2_sudden <- glmmTMB(sudden ~RT+ PE*cor2           +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  # Exploratory!
  M3_sudden <- glmmTMB(sudden ~RT*PE+ PE*cor2        +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M4_sudden <- glmmTMB(sudden ~RT*PE+ PE*cor2+cor2*RT+blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M5_sudden <- glmmTMB(sudden ~RT*PE*cor2            +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  plot(check_distribution(M2_sudden))
  check_collinearity(M1_sudden)
  anova(M0a_sudden,M0b_sudden,M1_sudden, M2_sudden,M3_sudden,M4_sudden,M5_sudden)

  #hist(residuals(M1_sudden))
  summary(M1_sudden)
  tab_model(M0b_sudden, show.std  = T)
  tab_model(M1_sudden, show.std  = T)
  tab_model(M2_sudden, show.std  = T)
  tab_model(M5_sudden, show.std  = T)

  Sudden_ggpredict<- ggpredict(M2_sudden , c( 'PE',  "cor2"))
  Sudden_finalPlot2= ggplot(Sudden_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.10,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Suddenness",y = "Suddenness", fill = "Accuracy")+
                  xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))+
                 
  ##Fig.S1
  Sudden_ggpredict<- ggpredict(M5_sudden , c( 'PE',  "cor2", 'RT'))
  Fig_S1= ggplot(Sudden_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.10,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Suddenness",y = "Suddenness", fill = "Accuracy")+
    xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))+
    facet_wrap(~ facet, nrow=1)  
  
## AHA: SURPRISE ~ PE #####################################
 M0a_surprise <- glmmTMB(surprising ~RT        +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2& !is.na(CRA1$cor2),],na.action=na.omit)
 M0b_surprise <- glmmTMB(surprising ~RT   +cor2+blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2& !is.na(CRA1$cor2),],na.action=na.omit)
  M1_surprise <- glmmTMB(surprising ~RT+PE+cor2+blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2& !is.na(CRA1$cor2),],na.action=na.omit)
  M2_surprise <- glmmTMB(surprising ~RT+PE*cor2+blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2& !is.na(CRA1$cor2),],na.action=na.omit)

  plot(check_distribution(M3_surprise))
  anova(M0a_surprise,M0b_surprise,M1_surprise, M2_surprise)
  check_collinearity(M3_surprise)
  hist(residuals(M3_surprise))
  summary(M2_surprise)
  tab_model(M0b_surprise,show.std  = T)

  Surprise_ggpredict<- ggpredict(M2_surprise , c( 'PE',   'cor2'))
  Surprise_finalPlot2= ggplot(Surprise_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.1,size = .8,
    position=position_dodge(.18)) + theme_classic() +labs(title = "Surprise",y = "Surprise", fill = "Accuracy")+
    xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))
               

  
#### GLMM results - output tables ###########################################

  tab_model(M1_certain,M1_pleasure, M1_sudden,M2_surprise,  show.std = T  )

### GLMM results - output plots #####################################################
  
  library(ggpubr)
  Fig3= ggarrange(  Certain_finalPlot2,Pleasure_finalPlot2,Sudden_finalPlot2, Surprise_finalPlot2,common.legend = T, legend = "bottom",ncol = 2, nrow = 2)
  

