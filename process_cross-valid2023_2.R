#Looped 1000x (or whatever) XV of various models
#Cross validation on fd4 and other datasets
#Testing various model forms
#Need cfi_process3 to be run first, or empirical_cfi_update5.Rmd (preferred)
#Updated July 2023

require(dplyr)
require(stringr)
require(lubridate)

setwd("C:/Dan/_Remote_projects/ccp_git/ccp-cfi/")

#Use fd.xv for cross-validation; keep fd.recall just in case
fd.recall <- fd4

sharp <- filter(fd4, str_detect(ExpProject, "Sharpsand IM")) %>% 
  pull(num)  #original Sharpsand Immature stands only
sharp.sm <- filter(fd4, str_detect(ExpProject, "Sharpsand SM")) %>% 
  pull(num)  #Sharpsand - immature
icfme <- 77:87 #order is A, 1, 3, 4, 5, 6, 7, 8a, 8b, 9, 2; fixed in icfme process

#Cross-validation 
#
#kWayCrossValidation - let's see how this compares

require(vtreat)
require(stringr)

#Number of loops for cross-validation (~20 for testing; 1000 for production)
Times=1000
nSplits=4

#chg names so that models can work
fd4$FSG <- fd4$FSG1
fd4$CFI <- fd4$CFI1

fd4 <- filter(fd4, num <=113)

#fmlaFF is ffmc model
fComplex <- formula(CFI ~ ws + I(FSG^1.5) * I(log(SFC)) + ws:MC.SA + I(FSG^1.5):MC.SA + I(log(SFC)):MC.SA)

#Set up formulas for competing model forms
f1 <- formula(CFI ~ ws + FSG + SFC + MC.FFMC)
f2 <- formula(CFI ~ ws + FSG + SFC + MC.SA)
f3 <- formula(CFI ~ ws + I(FSG^1.5) + SFC + MC.FFMC)
f4 <- formula(CFI ~ ws + I(FSG^1.5) + SFC + MC.SA)
f5 <- formula(CFI ~ ws + I(FSG^1.5) + SFC + ws:MC.FFMC)
f6 <- formula(CFI ~ ws + I(FSG^1.5) + SFC + ws:MC.SA)
f7 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC)) + ws:MC.FFMC)
f8 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC)) + ws:MC.SA)

#for f9:f11 - FSG[sharp] = FSG=2, and FSG-1 for sharp-SM
#E.g. 4.3 m to 2.3 m for sharp-IM (1970s)
#5.3 to 4.3 for sharp-SM

fd5 <- fd4
fd5$FSG[sharp] <- fd4$FSG[sharp]-2
fd5$FSG[sharp.sm] <- fd4$FSG[sharp.sm]-1
fd5$FSG[sharp.th] <- fd4$FSG[sharp.th]-0.25  #not always used; also tried 0.3
fd5$FSG[icfme] <- filter(fd5, num %in% icfme) %>% 
  mutate(FSG=case_when(
    FSG > 2 ~ FSG-1, 
    TRUE ~ FSG)) %>% pull(FSG) #icfme FSG > 2 m or so, reduce by 1m, else leave alone

#Use fd5
f9 <- formula(CFI ~ ws + I(FSG^1.5) + SFC.CLS + ws:MC.FFMC)  #maybe use MC.SA as well?
f10 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC)) + ws:MC.FFMC)
f11 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC)) + ws:MC.SA)

#new - still with fd5  f12 or f11.b
f11.b <- formula(CFI ~ ws + I(FSG^2) + I(log(SFC)) + ws:MC.SA)  #FSG^2
f12 <- formula(CFI ~ ws + I(FSG^2) + I(log(SFC)) + ws:MC.SA)  #used to be f11.b

#new - fd4
f0 <- formula(CFI ~ ws + FSG)  #model A
f0.5 <- formula(CFI ~ ws + FSG + MC.FFMC)  #model B

#start with fd4, then calculate SFC2=SFC + SFC.LF
#Compare with base models (e.g. M7, M8)

fd6 <- fd4 %>%
  mutate(SFC2=SFC + SFC.LF)

f13 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC2)) + ws:MC.FFMC)
f14 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC2)) + ws:MC.SA)

#LF with multiplier ##Changed Nov 2023
#Assume -2m for Sharpsand-IM stands is correct due to snags
#FSG already calculated to include snag influence from base rescale equation
#Multiplier:  

#Crap getting myself confused!! Incomplete Nov. 11 2023

#(fd5$FSG[sharp.th] %>% mean() / (fd5$FSG[sharp.th] %>% mean() - fd5$FSG[sharp] %>% mean()))^1.5
#C & A simplest assumption: LCBH=4, FSG=2 (z=4, z-zL=2)
#(4/(4-2))^1.5
#2.8284

#Second calc, with actual stand values from Stocks measurements:
#(z/(z-zL))^1.5
#z=4.38, zL=3.45/2
#zL=ladder fuel centroid; at Sharpsand IM=1.725
#(4.38/(4.38-1.725))^1.5
#2.118919

#fd5$FSG[sharp.sm] <- fd4$FSG[sharp.sm]-1
#Equivalent doesn't work as well
#(fd4$FSG[sharp.sm] %>% mean() / (fd5$FSG[sharp.sm] %>% mean()))^1.5
#1.3689 
M = 4 #or so? See below

fd7 <- fd4 %>% 
  mutate(SFC3 = SFC + SFC.LF*M)

f15 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC3)) + ws:MC.FFMC)
f16 <- formula(CFI ~ ws + I(FSG^1.5) + I(log(SFC3)) + ws:MC.SA)


#forms <- c(f1, f2, f3, f4, f5, f6, f7, f8)

#New Mod Table (copied and edited from empirical_cfi_update5):
#Table 3 - fitted coefficients model table
#Models 7, 8, 9, 10, 11

#m7 and m8: fd4
model7 <- glm(f7, data=fd4, family=binomial)
preds.m7 <-predict(model7, type="response")
fd.m7 <- mutate(fd4, Pr_CFI=round(preds.m7, digits=4), 
                Corr=round(Pr_CFI, digits=0)==CFI)
acc.m7 <- nrow(filter(fd.m7, Corr==TRUE))/nrow(fd4)

model8 <- glm(f8, data=fd4, family=binomial)
preds.m8 <-predict(model8, type="response")
fd.m8 <- mutate(fd4, Pr_CFI=round(preds.m8, digits=4), 
                Corr=round(Pr_CFI, digits=0)==CFI)
acc.m8 <- nrow(filter(fd.m8, Corr==TRUE))/nrow(fd4)

#m9, m10, m11: fd5 with reduced FSG
model9 <- glm(f9, data=fd5, family=binomial)
preds.m9 <-predict(model9, type="response")
fd.m9 <- mutate(fd5, Pr_CFI=round(preds.m9, digits=4), 
                Corr=round(Pr_CFI, digits=0)==CFI)
acc.m9 <- nrow(filter(fd.m9, Corr==TRUE))/nrow(fd4)

model10 <- glm(f10, data=fd5, family=binomial)
preds.m10 <-predict(model10, type="response")
fd.m10 <- mutate(fd5, Pr_CFI=round(preds.m10, digits=4), 
                 Corr=round(Pr_CFI, digits=0)==CFI)
acc.m10 <- nrow(filter(fd.m10, Corr==TRUE))/nrow(fd4)

model11 <- glm(f11, data=fd5, family=binomial)
preds.m11 <-predict(model11, type="response")
fd.m11 <- mutate(fd5, Pr_CFI=round(preds.m11, digits=4), 
                 Corr=round(Pr_CFI, digits=0)==CFI)
acc.m11 <- nrow(filter(fd.m11, Corr==TRUE))/nrow(fd4)

model12 <- glm(f12, data=fd5, family=binomial)
preds.m12 <-predict(model12, type="response")
fd.m12 <- mutate(fd5, Pr_CFI=round(preds.m12, digits=4), 
                  Corr=round(Pr_CFI, digits=0)==CFI)
acc.m12 <- nrow(filter(fd.m12, Corr==TRUE))/nrow(fd4)

model13 <- glm(f13, data=fd6, family=binomial)
preds.m13 <-predict(model13, type="response")
fd.m13 <- mutate(fd6, Pr_CFI=round(preds.m13, digits=4), 
                 Corr=round(Pr_CFI, digits=0)==CFI)
acc.m13 <- nrow(filter(fd.m13, Corr==TRUE))/nrow(fd4)

model14 <- glm(f14, data=fd6, family=binomial)
preds.m14 <-predict(model14, type="response")
fd.m14 <- mutate(fd6, Pr_CFI=round(preds.m14, digits=4), 
                 Corr=round(Pr_CFI, digits=0)==CFI)
acc.m14 <- nrow(filter(fd.m14, Corr==TRUE))/nrow(fd4)

model15 <- glm(f15, data=fd7, family=binomial)
preds.m15 <-predict(model15, type="response")
fd.m15 <- mutate(fd7, Pr_CFI=round(preds.m15, digits=4), 
                 Corr=round(Pr_CFI, digits=0)==CFI)
acc.m15 <- nrow(filter(fd.m15, Corr==TRUE))/nrow(fd4)

model16 <- glm(f16, data=fd7, family=binomial)
preds.m16 <-predict(model16, type="response")
fd.m16 <- mutate(fd7, Pr_CFI=round(preds.m16, digits=4), 
                 Corr=round(Pr_CFI, digits=0)==CFI)
acc.m16 <- nrow(filter(fd.m16, Corr==TRUE))/nrow(fd4)

#Trying to fit M using variations on old models 5, 6
#get previously fitted coefficients (from Perrakis et al. 2023, SuppMat)

f5c.intercept<- -4.8496
fd4.m5 <- mutate(fd4, 
                 linCombo=f5c.intercept + 1.2658*ws + -0.397*FSG^1.5 + 
                   (-0.0662*MC.FFMC*ws) + 1.7128*SFC)

#All assigned coefficients and intercept become lumped into a linear combination variable, which is
#like the new intercept

f5con <- formula(CFI ~ linCombo -1 + SFC.LF)

#Model is fitted to last remaining variable, SFC.LF
model5con <- glm(formula=f5con, data=fd4.m5, family=binomial)
#4.2227/1.7128=2.465

#model6:
f6c.intercept<- -3.8869
fd4.m6 <- mutate(fd4, 
                 linCombo=f6c.intercept + 1.0597*ws + -0.38*FSG^1.5 + 
                   (-0.0493*MC.SA*ws) + 1.3516*SFC)

f6con <- formula(CFI ~ linCombo -1 + SFC.LF)

#Model is fitted to last remaining variable, SFC.LF
model6con <- glm(formula=f6con, data=fd4.m6, family=binomial)
#5.8251/1.3516 = 4.31  #fitted coefficient for LF divided by previously fitted coef. for SFC

#ok, around 3 should be optimal

preds.m5con <-predict(model5con, type="response")
fd.m5con <- mutate(fd4.m5, Pr_CFI=round(preds.m5con, digits=4), 
                   Corr=round(Pr_CFI, digits=0)==CFI)
acc.m5con <- nrow(filter(fd.m5con, Corr==TRUE))/nrow(fd4)


preds.m6con <-predict(model6con, type="response")
fd.m6con <- mutate(fd4.m6, Pr_CFI=round(preds.m6con, digits=4), 
                 Corr=round(Pr_CFI, digits=0)==CFI)
acc.m6con <- nrow(filter(fd.m6con, Corr==TRUE))/nrow(fd4)





#Get coefs
m7 <- model7$coefficients 
m8 <- model8$coefficients 
m9 <- model9$coefficients 
m10 <- model10$coefficients
m11 <- model11$coefficients
m13 <- model13$coefficients
m14 <- model14$coefficients
m15 <- model15$coefficients
m16 <- model16$coefficients

#accuracy for sf/cf per model
mods <- list(fd.m7, fd.m8, fd.m9, fd.m10, fd.m11, fd.m13, fd.m14, fd.m15, fd.m16)
true.cf <- NULL
true.sf <- NULL

for(i in 1:length(mods)) {
  true.cf[i] <- filter(mods[[i]], Corr==TRUE & CFI==1) %>% nrow()
}

for(i in 1:length(mods)) {
  true.sf[i] <- filter(mods[[i]], Corr==TRUE & CFI==0) %>% nrow()
}

true.ft <- data.frame(true.cf, true.sf) %>%
  mutate(model=c('m7', 'm8', 'm9', 'm10', 'm11', 'm13', 'm14', 'm15', 'm16')) #expand to m15 as needed

#Get model coefficients
model.coefs0 <- tibble(m7, m8, m10, m11, m13, m14, m15, m16) %>% t() #transpose for table
modNames <- model.coefs0[,1] %>% names   #clunky, but works; note - treats all vars like m7 (no FFMC)

vars <- names(m7)

model.coefs <- as.data.frame(model.coefs0)
colnames(model.coefs) <- vars
row.names(model.coefs) <- 1:8

#Good for most models - grab coefficients
model.coefs <- model.coefs %>% 
  select(1:3, 5, 4) %>% signif(5) %>%  #5 sig figs; chg order to put SFC at end
  mutate(model=modNames) %>%
  select(model, everything())  #resort columns for table

colnames(model.coefs) <- c('model', 'I', 'ws', 'FSG^1.5', 'ws X mc', 'log(SFC)')

m7sig <- summary(model7)$coefficients[,4]  #get p-values for each estimate
m8sig <- summary(model8)$coefficients[,4]
m10sig <- summary(model10)$coefficients[,4]
m11sig <- summary(model11)$coefficients[,4]
m13sig <- summary(model13)$coefficients[,4]
m14sig <- summary(model14)$coefficients[,4]
m15sig <- summary(model15)$coefficients[,4]
m16sig <- summary(model16)$coefficients[,4]

sig.vals0 <- data.frame(m7sig, m8sig, m10sig, m11sig, m13sig, m14sig, m15sig, m16sig) %>%
  t() #transpose for table

sig.vals1 <- as.data.frame(sig.vals0)
row.names(sig.vals1) <- 1:8  #normal row numbers
colnames(sig.vals1) <- c('sI', 'sw', 'sF', 'sC', 'sM')

sig.vals <- select(sig.vals1, sI, sw, sF, sM, sC) %>%
  signif(3) %>%
  mutate(model=modNames) 

aic.vals <- c(model7$aic, model8$aic, model10$aic, model11$aic, model13$aic, 
              model14$aic, model15$aic, model16$aic)
mod.list <- list(model7, model8, model10, model11, model13, model14, model15, model16)
nag.full <- lapply(mod.list, rcompanion::nagelkerke)

nagR2 <- NULL
for(i in 1:8) {
  nagR2=c(nagR2, nag.full[[i]][2][[1]][3])} %>% 
  unlist()  #nagR2 becomes vector of Nagelkerke R2 values for 8 models

aic.df <- data.frame(model=c(modNames), aic.vals=signif(aic.vals, 3), nagR2=signif(nagR2, 3)) 

sig.stars <- function(x) {
  case_when(
    x < 0.001 ~ '***',
    x < 0.01 ~ '**',
    x < 0.05 ~ '*',
    x < 0.1 ~ '.',
    TRUE ~ '-' )
}

mod.table <- left_join(model.coefs, sig.vals, by='model') %>%
  select(model, I, sI, ws, sw, 'FSG^1.5', sF, 'ws X mc', sM, 'log(SFC)', sC) %>%
  mutate(Acc = c(acc.m7, acc.m8, acc.m10, acc.m11, acc.m13, acc.m14, acc.m15, acc.m16) %>% 
           signif(3)) %>%
  left_join(aic.df, by='model') %>%
  left_join(true.ft, by='model') %>% #not sure if I need this; not kept in table 
  mutate(sI.s=sig.stars(sI),
         sw.s=sig.stars(sw),
         sF.s=sig.stars(sF),
         sM.s=sig.stars(sM),
         sC.s=sig.stars(sC)) %>% 
  select(model, I, sI.s, ws, sw.s, 'FSG^1.5', sF.s, 'ws X mc', sM.s, 'log(SFC)',
         sC.s, Acc, aic.vals, nagR2)  #keep Acc or true.cf; not both


#######################


#Monte Carlo Cross-Validation

##Redo as necessary

correct <- NULL
mcc <- NULL
aic <- NULL


#0: f0
form=f0
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.0 <- model
  correct.0 <- mean(correct)
  mcc.0 <- mean(mcc)
  aic.0 <- mean(aic)
  seed.0 <- xx
}

#0.5: f0.5
form=f0.5
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.05 <- model
  correct.05 <- mean(correct)
  mcc.05 <- mean(mcc)
  aic.05 <- mean(aic)
  seed.05 <- xx
}


#1: f1
form=f1
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.a <- model
  correct.a <- mean(correct)
  mcc.a <- mean(mcc)
  aic.a <- mean(aic)
  seed.a <- xx
}


#2: f2
form=f2
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)

  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.b <- model
  correct.b <- mean(correct)
  mcc.b <- mean(mcc)
  aic.b <- mean(aic)
  seed.b <- xx
}
#end big loop; 

#3: f3
form=f3
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)

    #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.c <- model
  correct.c <- mean(correct)
  mcc.c <- mean(mcc)
  aic.c <- mean(aic)
  seed.c <- xx
}
#end big loop; 


#4: f4
form=f4
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)

  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.d <- model
  correct.d <- mean(correct)
  mcc.d <- mean(mcc)
  aic.d <- mean(aic)
  seed.d <- xx
}
#end big loop; 


#5: f5
form=f5
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.e <- model
  correct.e <- mean(correct)
  mcc.e <- mean(mcc)
  aic.e <- mean(aic)
  seed.e <- xx
}
#end big loop; 


#6: f6
form=f6
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.f <- model
  correct.f <- mean(correct)
  mcc.f <- mean(mcc)
  aic.f <- mean(aic)
  seed.f <- xx
}
#end big loop; 


#7: f7
form=f7
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.g <- model
  correct.g <- mean(correct)
  mcc.g <- mean(mcc)
  aic.g <- mean(aic)
  seed.g <- xx
}
#end big loop; 


#8: f8
form=f8
fd.xv <- fd4 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.h <- model
  correct.h <- mean(correct)
  mcc.h <- mean(mcc)
  aic.h <- mean(aic)
  seed.h <- xx
}
#end big loop; 


#9: f9 - Note fd5
form=f9
fd.xv <- fd5 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }

  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.i <- model
  correct.i <- mean(correct)
  mcc.i <- mean(mcc)
  aic.i <- mean(aic)
  seed.i <- xx
}
#end big loop; 


#10: f10 - Note data=fd5
form=f10
fd.xv <- fd5 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.j <- model
  correct.j <- mean(correct)
  mcc.j <- mean(mcc)
  aic.j <- mean(aic)
  seed.j <- xx
}
#end big loop; 


#11: f11 - Note data=fd5
form=f11
fd.xv <- fd5 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.k <- model
  correct.k <- mean(correct)
  mcc.k <- mean(mcc)
  aic.k <- mean(aic)
  seed.k <- xx
}
#end big loop; 

#11b: f11b - Note data=fd5
form=f11.b
fd.xv <- fd5 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.kb <- model
  correct.kb <- mean(correct)
  mcc.kb <- mean(mcc)
  aic.kb <- mean(aic)
  seed.kb <- xx
}




#TO HERE FOR Empirical CFI Update####################################

#12: f12 - Not used?? Same as f11.b
form=f12
fd.xv <- fd6 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC2, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.l <- model
  correct.l <- mean(correct)
  mcc.l <- mean(mcc)
  aic.l <- mean(aic)
  seed.l <- xx
}
#end big loop; 


#13: f13 - Note data=fd6, SFC2
form=f13
fd.xv <- fd6 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC2, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.m <- model
  correct.m <- mean(correct)
  mcc.m <- mean(mcc)
  aic.m <- mean(aic)
  seed.m <- xx
}
#end big loop; 


#14: f14 - Note data=fd6, SFC2 
form=f14
fd.xv <- fd6 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC2, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.n <- model
  correct.n <- mean(correct)
  mcc.n <- mean(mcc)
  aic.n <- mean(aic)
  seed.n <- xx
}
#end big loop; 


#15: f15 - Note data=fd7, SFC3
form=f15
fd.xv <- fd7 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC3, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.o <- model
  correct.o <- mean(correct)
  mcc.o <- mean(mcc)
  aic.o <- mean(aic)
  seed.o <- xx
}
#end big loop; 

#16: f16 - Note data=fd7, SFC3
form=f16
fd.xv <- fd7 %>% select(num, fire, ws, MC.SA, MC.FFMC, SFC.CLS, FSG, SFC, SFC3, SFC.LF, CFI)

for(rep in 1:Times) {
  splitPlan <- kWayCrossValidation(nrow(fd.xv), nSplits, NULL, NULL)
  
  # initialize the prediction vector
  xx <- .Random.seed
  fd.xv$predxv <- 0  
  
  for(i in 1:4) {
    split <- splitPlan[[i]]
    model <- glm(formula=form, data = fd.xv[split$train, ], 
                 family=binomial(link="logit"))
    fd.xv$predxv[split$app] <- predict.glm(model, newdata = fd.xv[split$app, ],
                                           type="response")
  }
  
  fdxv<-mutate(fd.xv, predxv=round(predxv, digits=4)) %>%
    mutate(CorrXv=round(predxv, digits=0)==CFI)
  
  aic[rep]<-model$aic
  correct[rep]<-nrow(filter(fdxv, CorrXv==TRUE))/nrow(fdxv)
  ##
  #Confusion matrix and Matthews correlation coefficient (mcc)
  tp<-nrow(filter(fdxv, CorrXv==TRUE & CFI==1))  #TRUE crown fires
  tn<-nrow(filter(fdxv, CorrXv==TRUE & CFI==0))  #TRUE surface fires
  fp<-nrow(filter(fdxv, CorrXv==FALSE & CFI==1))  #FALSE missed crown fires
  fn<-nrow(filter(fdxv, CorrXv==FALSE & CFI==0))  #FALSE missed surface fires
  mcc[rep]<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  #rename for each
  model.p <- model
  correct.p <- mean(correct)
  mcc.p <- mean(mcc)
  aic.p <- mean(aic)
  seed.p <- xx
}
#end big loop; 


#End all

#Bring results to summary table - AllFigs script

require(stringr)

metrics <- c('correct.', 'mcc.', 'aic.')
model.order <- c(letters[1:18])
CFI.letters <- rep(model.order, each=length(metrics))


#CFI.models <- str_c(metrics, CFI.letters)

#latest version - no .l (duplicate of kb/11.b removed); jumps from .kb to .m
accuracy.vals <- c(correct.0, correct.05, correct.a, correct.b, correct.c, correct.d, correct.e, correct.f, correct.g, correct.h, 
                   correct.i, correct.j, correct.k, correct.kb, correct.m, correct.n, correct.o, correct.p) 
mcc.vals <- c(mcc.0, mcc.05, mcc.a, mcc.b, mcc.c, mcc.d, mcc.e, mcc.f, mcc.g, mcc.h, mcc.i, mcc.j, mcc.k, mcc.kb, mcc.m, 
              mcc.n, mcc.o, mcc.p)
aic.vals <- c(aic.0, aic.05, aic.a, aic.b, aic.c, aic.d, aic.e, aic.f, aic.g, aic.h, aic.i, aic.j, aic.k, aic.kb, aic.m,
              aic.n, aic.o, aic.p)
modelform <- c(f0, f0.5, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f11.b, f13, f14, f15, f16) 
runnum <- c('f0', 'f0.5', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f11.b', 
  'f13', 'f14', 'f15', 'f16')
Dataset <- c(rep('1.Base', times=10), rep('2.FSG Adj.', times=4), rep('3.LF to SFC', times=2),
             rep('3.LF x M to SFC', times=2))  


reorder <- 1:18 #reorder <- c(4, 3, 2, 1, 6, 8, 9, 10, 11, 5, 12, 13, 7)  #huge pain to do this correctly!

#XV.archive <- data.frame(as.character(modelform), accuracy.vals, mcc.vals, aic.vals, Dataset)
#XV.seed <- data.frame(seed.a, seed.b, seed.c, seed.d, seed.e, seed.f, seed.g, seed.h, seed.i, seed.j, 
#                      seed.k, seed.l, seed.r)

# write.csv(XV.archive, './XVarchive.csv' )
# write.csv(XV.seed, './XVseed.csv' )

#Full table summary, comparing all models
XV.summary <- data.frame(model.order, reorder, Model=as.character(modelform) %>% str_replace('CFI ~ ', ''), 
                         runnum, 
                         Accuracy=round(accuracy.vals, 4), MCC=round(mcc.vals, 4), AIC=round(aic.vals, 2), 
                         Dataset) %>%
  mutate(MCC.rank=rank(length(MCC) + 1 - rank(MCC)),  #isn't there a descending rank function?
         AIC.rank=rank(AIC)) %>% 
  arrange(reorder) %>% select(-model.order, -reorder) 


# Simplified table summary, comparing only models from base and FSG-adjust datasets 
# (no ladder fuel adjustments) for empirical CFI paper
XV.sum.simple <- XV.summary %>%
  slice(1:14) %>%
  select(-c(MCC.rank, AIC.rank)) %>%
  mutate(MCC.rank=rank(length(MCC) + 1 - rank(MCC)), 
         AIC.rank=rank(AIC),
         AIC=round(AIC,1))
  
XV.summary <- mutate(XV.summary, 
                     Accuracy=round(Accuracy, 3), MCC=round(MCC, 3), AIC=round(AIC, 1))

#new version, only good model forms, base and up for 2024
#9, 10, 12, 13, 15, 16, 17, 18
model.order2 <- letters[c(9, 10, 12, 13, 15, 16, 17, 18)]
reorder2 <- c(9, 10, 12, 13, 15, 16, 17, 18)
runnum2 <- runnum[reorder2]
XV.sum2 <- data.frame(model.order2, reorder2, 
                      Model=as.character(modelform[reorder2]) %>% str_replace('CFI ~ ', ''), 
                      runnum2, Accuracy=round(accuracy.vals[reorder2], 4), 
                      MCC=round(mcc.vals[reorder2], 4), 
                      AIC=round(aic.vals[reorder2], 2), 
                         Dataset[reorder2]) %>%
  mutate(MCC.rank=rank(length(MCC) + 1 - rank(MCC)),  #isn't there a descending rank function?
         AIC.rank=rank(AIC)) %>% 
  arrange(reorder2) %>% select(-model.order2, -reorder2) 

#changed to '_2' Apr 10, 2023 after adding 3 PNFI-RP fires
#Changed to '_3' on May 28, 2023 after changing sharp-SM SFC to estimate CFC portion
#Changed to '_4' on July 13, 2023 after fixing small errors in Sharp-SM SFC
#Added '_5' on Dec 6 2023 after solving for M using 'linearCombo' method
#Added '_2024' on Jan 25 2024 after cleaning up main M analyis in LF RMD file
#fixed LF calculations (fd_current _process2) and reran. fire_data_may2023b.csv has corrections

#write.csv(XV.summary, './tables/xv_summary5.csv') - not to be trusted; may have screwed up order
#write.csv(XV.sum.simple, './tables/xv_sum_simple4.csv') - same
#write.csv(XV.sum2, './tables/xv_sum2_2024.csv')


seeds <- list(seed.a, seed.b, seed.c, seed.d, seed.e, seed.f, seed.g, seed.h, seed.i, seed.j, seed.k, seed.kb, seed.l, 
              seed.m, seed.n, seed.o)

#write.csv(seeds, 'seeds4.csv')

#Final model output - WTF? why m16
model16 <- glm(f15, data=fd7, family=binomial)
final.mod <- model16

write.csv(final.mod$coefficients, './tables/final_mod.csv')
