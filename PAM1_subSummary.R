#### Data on participant for PAM1 (rebound) part ####
load(file='Z://PD_motor//subj_data///alldata.RData')

pam1data <- subset(alldata,PAM1=='yes')
pam1data <- subset(pam1data, MEG_ID != "345") # Did not complet task / discarded from all analysis
pam1data <- subset(pam1data, MEG_ID != "328") # Does not have PAM1 
pam1data <- subset(pam1data, MEG_ID != "393") # Bad MEG data. Excluded from analysis.
pam1data$MEG_ID <- factor(pam1data$MEG_ID)

# Analysis
library(BayesFactor)

pam1gender <- xtabs(~sex+Sub_type,pam1data)
fisher.test(pam1gender)
contingencyTableBF(pam1gender, sampleType = 'indepMulti',fixedMargin='cols')   # Testbetween colums (i.e group)

aggregate(pam1data$age, by=list(pam1data$Sub_type), mean)
aggregate(pam1data$age, by=list(pam1data$Sub_type), range)
aggregate(pam1data$age, by=list(pam1data$Sub_type), sd)
t.test(age~Sub_type, pam1data)
ttestBF(formula=age~Sub_type, data=pam1data, paired=F, rscale='wide')

aggregate(pam1data$MoCA, by=list(pam1data$Sub_type), mean)
aggregate(pam1data$MoCA, by=list(pam1data$Sub_type), sd)
t.test(MoCA~Sub_type, pam1data)
ttestBF(formula=MoCA~Sub_type, data=pam1data, paired=F, rscale='wide')

aggregate(pam1data$HADS_angst, by=list(pam1data$Sub_type), mean)
aggregate(pam1data$HADS_angst, by=list(pam1data$Sub_type), sd)
t.test(HADS_angst~Sub_type, pam1data)
ttestBF(formula=HADS_angst~Sub_type, data=pam1data, paired=F, rscale='wide')

aggregate(pam1data$HADS_depression, by=list(pam1data$Sub_type), mean)
aggregate(pam1data$HADS_depression, by=list(pam1data$Sub_type), sd)
t.test(HADS_depression~Sub_type, pam1data)
ttestBF(formula=HADS_depression~Sub_type, data=pam1data, paired=F, rscale='wide')

## UPDRS moved to "UPDRS_stats" script

mean(pam1patients$disease_dur, na.rm=T)
range(pam1patients$disease_dur, na.rm=T)
median(pam1patients$LEDD, na.rm=T)
mean(pam1patients$LEDD, na.rm=T)
sd(pam1patients$LEDD, na.rm=T)

sleepBefore <- c(pam1data$sleep_1_pam5,pam1data$sleep_2_pam5)
sleepAfter <- c(pam1data$sleep_1_pam1,pam1data$sleep_2_pam1)
session <- as.factor(c(rep(1,length(pam1data$sleep_1_pam5)),rep(2,length(pam1data$sleep_2_pam5))))

sleepData <- data.frame(id=pam1data$MEG_ID,group=pam1data$Sub_type, session, sleepBefore,sleepAfter)

mano <- manova(cbind(sleepBefore,sleepAfter) ~ session*group + Error(id), data=sleepData )
summary(mano)

library(lme4)
sleepLm <- lmer(sleepAfter~sleepBefore*session*group+(1|id),data=sleepData,REML = F)
summary(sleepLm)

aggregate(sleepData$sleepBefore, by=list(session=sleepData$session), mean, na.rm=T)
aggregate(sleepData$sleepBefore, by=list(session=sleepData$session), sd, na.rm=T)
aggregate(sleepData$sleepAfter, by=list(session=sleepData$session), mean, na.rm=T)
aggregate(sleepData$sleepAfter, by=list(session=sleepData$session), sd, na.rm=T)

library(ggplot2)
ggplot(sleepData, aes(x=sleepBefore, y=sleepAfter, color=session, shape=group))+
  geom_point()+theme_bw()+
  geom_smooth(method=lm)


