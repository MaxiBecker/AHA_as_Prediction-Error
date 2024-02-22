#### R-Script FOR MANUSCRIPT Becker, Wang & Cabeza, 2024 ############################################

# Assumptions: AHA-Dimensions correspond to (meta-cognitive) PE
# last update: 22.02.2024
# (c) almaxi@gmail.com
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

## load / prepare data ############################################################

  setwd('C:/Users/Maxi/Google Drive/studies/Maxi/IN_REVISION/AHAPE/GitHub_4_publication')
  load("BeckerWangCabeza_2023_data.Rdata")
  
  CRA$subj_solution_likelihood[CRA$subj_solution_likelihood==6] = NA  #6 = was already solved during rating
  
  CRA1 = CRA[!is.na(CRA$subj_solution_likelihood),]
  min_SSL=min(CRA1$subj_solution_likelihood, na.rm = T)
  max_SSL = max(CRA1$subj_solution_likelihood, na.rm = T)
  CRA1$subj_solution_likelihood_norm = (CRA1$subj_solution_likelihood - min_SSL) /( max_SSL - min_SSL) 
  CRA1$PE = 1-CRA1$subj_solution_likelihood_norm

  unique(CRA1$subject)

##############################################################################-
#### AHA Measurement model: calculating composite AHA score ####

  #partial correlation with random effect subject
  CRA1_clean <- CRA[CRA$cor2 ==1 & CRA$RT>2,]
  sudden1_fit <- lme4::lmer(sudden ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA1_clean,na.action=na.omit )
  pleasing1_fit <- lme4::glmer(pleasing ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA1_clean,family = binomial (link = "logit"), na.action=na.omit )
  certain1_fit <- lme4::lmer(certain ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA1_clean,na.action=na.omit )
  surprising1_fit <- lme4::lmer(surprising ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA1_clean,na.action=na.omit )
  
  CRA_partial<- data.frame(residuals(sudden1_fit), residuals(pleasing1_fit),residuals(certain1_fit), residuals(surprising1_fit))
  names(CRA_partial) <- c("sudden", "pleasing", "certain", "surprising")
  
 library(lavaan)
 insightfac <-   'AHA  =~ certain + sudden + surprising + pleasing   
                   AHA ~~ 1*AHA
                   pleasing ~~ surprising'
                  
     #fit
     fit <- lavaan(insightfac, data=CRA_partial, 
                   auto.var=T, auto.fix.first=F ,
                   auto.cov.lv.x=T, estimator = "MLR" )
     summary(fit, fit.measures=TRUE, standardized = T)
     modindices(fit, sort = TRUE)
     #standardizedSolution(fit)
     CRA$AHA_facload = predict(fit)
     inspect(fit,what="std")$lambda
     lavInspect(fit, what = "est")$theta
     
     #check correlation plot
     CRA1_clean_dim = CRA1_clean[,c("certain" , "sudden" , "surprising" , "pleasing" )]
     cor.test(CRA1_clean_dim$surprising,CRA1_clean_dim$sudden, method = "spearman")
     cor.test(CRA1_clean_dim$surprising,CRA1_clean_dim$pleasing, method = "kendall")
     cor.test(CRA1_clean_dim$surprising,CRA1_clean_dim$certain, method = "spearman")
    
############################################################################-
### Descriptive Statistics #######
     
hist(CRA$Rate_Response)
CRA_all = CRA %>% group_by(subject) %>% summarise_all(funs(mean), na.rm =T)   

 CRA_sd= CRA %>%
    group_by(subject) %>%
    summarise_at(c( "sudden", "pleasing", "certain", "surprising", "Rate_Response"), sd, na.rm = TRUE) #"insight",

  summary(CRA_all)
  sum(CRA_all$sex)
  CRA_age = CRA[,c("age","subject")] %>% group_by(subject) %>% summarise_all(funs(mean), na.rm =T)   
  mean(CRA_all$age, na.rm = T)
  median(CRA_all$age, na.rm = T)
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
  M2_cor2 <- glmmTMB(cor2 ~ blocknum + scale(RT)*subj_solution_likelihood + (1 | subject) + (1 | cues), data= CRA1[CRA1$RT>2,], family = binomial(link = "logit"), na.action = na.omit)
  
  plot(check_distribution(M1_cor2))
  check_collinearity(M2_cor2)
  anova(M0_cor2,M1_cor2,M2_cor2)
  #hist(residuals(M1))
  tab_model(M1_cor2, show.std = T)
  ggpredict(M2_cor2 , c("subj_solution_likelihood", "RT" ),main = "") %>% plot() + ggplot2::theme_classic()

## AHA: PLEASURE ~ PE ###########################################################################################################################
  M0a_pleasure<- glmer(pleasing ~ RT       +blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  M0b_pleasure<- glmer(pleasing ~ RT   +cor2+blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  M1_pleasure <- glmer(pleasing ~ RT+PE+cor2+blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  M2_pleasure <- glmer(pleasing ~ RT+PE*cor2+blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  M3_pleasure <- glmer(pleasing ~ RT*PE*cor2+blocknum+(1|subject)+(1|cues),data= CRA1[CRA1$RT& !is.na(CRA1$cor2),],family = binomial(link = "logit"),na.action = na.omit)
  
  plot(check_distribution(M2_pleasure))
  check_collinearity(M2_pleasure)
  anova(M0a_pleasure,M0b_pleasure,M1_pleasure, M2_pleasure,M3_pleasure)
  hist(residuals(M1_pleasure))
  parameters::equivalence_test(M1_pleasure) #undecided
  
  summary(M0b_pleasure)
  tab_model(M0b_pleasure)
  tab_model(M1_pleasure)
  tab_model(M2_pleasure)
  tab_model(M3_pleasure)
  
  Pleasure_ggpredict<- ggpredict(M1_pleasure , c( 'PE', "cor2" ))
  Pleasure_finalPlot2= ggplot(Pleasure_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.1,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Pleasure",y = "Pleasure (probability)", fill = "Accuracy")+
    xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))

  #Pleasure_ggpredict<- ggpredict(M1_pleasure , c( 'PE', "cor2" ))
  Pleasure_dot = plot(Pleasure_ggpredict, show_ci = T, show_data = TRUE, jitter = T)+ggplot2::theme_classic()
  
## AHA: CERTAINTY ~ PE ############################################################################################################
  M0a_certain <-lmer(certain ~ RT        +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M0b_certain <-lmer(certain ~ RT+cor2+  +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M1_certain <- lmer(certain ~ RT+cor2+PE+blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M2_certain <- lmer(certain ~ RT+cor2*PE+blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M3_certain <- lmer(certain ~ RT*cor2*PE+blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  
  plot(check_distribution(M1_certain))#residuals are normal
  anova(M0a_certain,M0b_certain,M1_certain, M2_certain,M3_certain)
  summary(M1_certain)
  tab_model(M0b_certain,show.std  = T)
  tab_model(M1_certain,show.std  = T)
  tab_model(M2_certain,show.std  = T)
  tab_model(M3_certain,show.std  = T)
  
  Certain_ggpredict<- ggpredict(M2_certain , c( 'PE','cor2' ))
  Certain_finalPlot2= ggplot(Certain_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.1,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Certainty",y = "Certainty", fill = "Accuracy")+
    xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)") )+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))

  Certain_dot = plot(Certain_ggpredict, show_ci = T, show_data = TRUE, jitter = T)+ggplot2::theme_classic()

## AHA: SUDDENNESS ~ PE ###################################################################################################################
  M0a_sudden<- lmer(sudden ~RT+        +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M0b_sudden<- lmer(sudden ~RT+    cor2+blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M1_sudden <- lmer(sudden ~RT+ PE+cor2+blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M2_sudden <- lmer(sudden ~RT+ PE*cor2+blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  M3_sudden <- lmer(sudden ~RT*PE*cor2 +blocknum + (1|subject) +(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),], na.action = na.omit)
  plot(check_distribution(M2_sudden))
  check_collinearity(M1_sudden)
  anova(M0a_sudden,M0b_sudden,M1_sudden, M2_sudden,M3_sudden,M4_sudden,M5_sudden)

  #hist(residuals(M1_sudden))
  summary(M1_sudden)
  tab_model(M0b_sudden, show.std  = T)
  tab_model(M1_sudden, show.std  = T)
  tab_model(M2_sudden, show.std  = T)
  tab_model(M5_sudden, show.std  = T)
  Sudden_finalPlot=ggpredict (M5_sudden , c("PE", 'cor2', "RT"),main = "") %>% plot() + ggplot2::theme_classic()
  
  Sudden_ggpredict<- ggpredict(M5_sudden , c( 'PE',  "cor2"))
  Sudden_finalPlot2= ggplot(Sudden_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.10,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Suddenness",y = "Suddenness", fill = "Accuracy")+
                  xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))

      Sudden_dot = plot(Sudden_ggpredict, show_ci = T, show_data = TRUE, jitter = T)+ggplot2::theme_classic()
  

  ##Fig.S1
  Sudden_ggpredict<- ggpredict(M5_sudden , c( 'PE',  "cor2", 'RT'))
  Fig_S1= ggplot(Sudden_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.10,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Suddenness",y = "Suddenness", fill = "Accuracy")+
    xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))+
    facet_wrap(~ facet, nrow=1)  
  
## AHA: SURPRISE ~ PE #####################################
  
 M0a_surprise <- lmer(surprising ~RT         +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
 M0b_surprise <- lmer(surprising ~RT   +cor2 +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M1_surprise <- lmer(surprising ~RT+PE+cor2 +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M2_surprise <- lmer(surprising ~RT+PE*cor2 +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)
  M3_surprise <- lmer(surprising ~RT*PE*cor2 +blocknum+(1|subject)+(1|cues),data=CRA1[CRA1$RT>2 & !is.na(CRA1$cor2),],na.action=na.omit)

  plot(check_distribution(M3_surprise))
  anova(M0a_surprise,M0b_surprise,M1_surprise, M2_surprise,M3_surprise)
  check_collinearity(M3_surprise)
  hist(residuals(M3_surprise))
  summary(M2_surprise)
  tab_model(M0b_surprise,show.std  = T)
  tab_model(M3_surprise,show.std  = T)
  ggpredict(M3_surprise , c("PE", 'cor2', 'RT'),main = "") %>% plot() + ggplot2::theme_classic()
  Surprise_finalPlot=ggpredict(M2_surprise , c("PE", 'cor2'),main = "")# %>% plot() + ggplot2::theme_classic()
  #emmeans(M3_surprise, list(pairwise ~ diff_expected_outcome | cor2), adjust = "tukey")
  
  Surprise_ggpredict<- ggpredict(M2_surprise , c( 'PE',   'cor2'))
  Surprise_finalPlot2= ggplot(Surprise_ggpredict, aes(x= x, y = predicted , fill= group)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T) + #, color = "black"
    geom_errorbar(aes(ymin= conf.low, ymax=conf.high), width=.1,size = .8,
                  position=position_dodge(.18)) + theme_classic() +labs(title = "Surprise",y = "Surprise", fill = "Accuracy")+
                  xlab(expression(PE[meta]~": Δ(actual – expected solution likel.)"))+ scale_fill_manual(breaks = c(0,1),values=c( "#ecc19c", "#1e847f"))
                  #facet_wrap(~ facet, nrow=1)   #scale_fill_manual(breaks = c(0,1),values=c( "#f98578", "#37c0c6"))+

  Surprise_dot = plot(Surprise_ggpredict, show_ci = T, show_data = TRUE, jitter = T)+ggplot2::theme_classic()
  
#### GLMM results - output tables ###########################################

  tab_model(M1_certain,M1_pleasure, M1_sudden,M2_surprise,  show.std = T  )

### GLMM results - output plots #####################################################
  
  library(ggpubr)
  Fig3= ggarrange(  Certain_finalPlot2,Pleasure_finalPlot2,Sudden_finalPlot2, Surprise_finalPlot2,common.legend = T, legend = "bottom",ncol = 2, nrow = 2)
  REvision2_4= ggarrange(  Certain_dot,Pleasure_dot,Sudden_dot, Surprise_dot,common.legend = T, legend = "bottom",ncol = 2, nrow = 2)
  

  ###############################################################################-
  ########### CONTROL EXPERIMENT 1 - CRA without Rating task #########################
  rm(list=ls())
  CRA_ctrl = read.table("ahape_CRAs_control_norating.csv", sep = ";", dec= ".", header = T, na.strings=c("", " ", "NA", "NAN"))
  load("BeckerWangCabeza_2023_data_ctrl.Rdata")
  
  #calculate accuracy over subjects-> to identify outliers (<25%)
  CRAdata_ID <- CRA_ctrl %>% group_by(subject) %>% dplyr::summarise(across(c("cor2", "sudden", "certain", "surprising", "pleasing"), ~mean(., na.rm=T))) 
  
  hist(CRAdata_ID$cor2)
  outliers_cor2 = CRAdata_ID[CRAdata_ID$cor2 <.25,]$subject 
  participants_cor2 = CRAdata_ID[CRAdata_ID$cor2 >=.25,]$subject
  
  #calculate variance in AHA experience over subjects-> to identify outliers (=0)
  CRAdata_ID_sd <- CRA_ctrl %>% group_by(subject) %>% dplyr::summarise(across(c("sudden", "certain", "surprising", "pleasing"), ~sd(., na.rm=T))) 
  outliers_sd = CRAdata_ID_sd[CRAdata_ID_sd$sudden ==0 | CRAdata_ID_sd$pleasing ==0 |
                              CRAdata_ID_sd$certain ==0 |CRAdata_ID_sd$surprising ==0 ,]$subject 
  
  n= unique(CRAdata_ID[!CRAdata_ID$subject %in% outliers_total,])
  
  outliers_total = as.vector(na.omit(c(outliers_sd, outliers_cor2)))
  
  # compute demographics without outliers
  CRA_ctrl_ID <- CRA_ctrl[!CRA_ctrl$subject %in% outliers_total,] %>% group_by(subject) %>% dplyr::summarise(across(c("sex", 'age'), ~mean(., na.rm=T))) 
  mean(CRA_ctrl_ID$age, na.rm = T);  sd(CRA_ctrl_ID$age, na.rm = T)
  median(CRA_ctrl_ID$age, na.rm = T); 
  max(CRA_ctrl_ID$age, na.rm = T); min(CRA_ctrl_ID$age, na.rm = T)
  mean(CRA_ctrl_ID$sex, na.rm = T);
  
  #regress out RT, blocknum & random effects from AHA dimensions
  CRA_clean= CRA_ctrl[CRA_ctrl$cor2 == 1 & CRA_ctrl$RT >2 & !CRA_ctrl$subject %in% outliers_total  ,] 
  
  sudden_fit <- lme4::lmer(sudden ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA_clean,na.action=na.omit )
  pleasing_fit <- lme4::lmer(pleasing ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA_clean,na.action=na.omit )
  certain_fit <- lme4::lmer(certain ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA_clean,na.action=na.omit )
  surprising_fit <- lme4::lmer(surprising ~ scale(RT) + blocknum+(1|subject) +(1|cues),data=CRA_clean,na.action=na.omit )
  
  CRA_ctrl_partial<- data.frame(residuals(sudden_fit), residuals(pleasing_fit),residuals(certain_fit), residuals(surprising_fit))
  names(CRA_ctrl_partial) <- c("sudden", "pleasing", "certain", "surprising")

  ###############-
  #### AHA Measurement model: calculating composite AHA score ####
  library(lavaan)

  insightfac <-   'AHA  =~ sudden + surprising + pleasing + certain  
                   pleasing ~~ surprising
                   AHA ~~ 1*AHA                 '
                 
  #fit
  fit_ctrl <- lavaan(insightfac, data=CRA_ctrl_partial,
                auto.var=T, auto.fix.first=F ,
                auto.cov.lv.x=T, estimator = "MLR" )
  summary(fit_ctrl, fit.measures=TRUE, standardized = T, rsquare = T)
  modindices(fit_ctrl, sort = TRUE)
  
  #standardizedSolution(fit)
  CRA_ctrl_partial$AHA_facload = predict(fit)
  inspect(fit_ctrl,what="std")$lambda
  lavInspect(fit_ctrl, what = "est")$theta
  
  #check intercorrelations
  cor.test(CRA_ctrl_partial$surprising, CRA_ctrl_partial$sudden, method="spearman")
  cor.test(CRA_clean$surprising, CRA_clean$pleasing, method="kendall")
  cor.test(CRA_clean$surprising, CRA_clean$certain, method="spearman")
  

        
        