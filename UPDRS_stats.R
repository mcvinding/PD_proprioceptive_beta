setwd('Z://PD_motor//subj_data//beta_rebound//UPDRS_stats')

#Load data
data <- read.csv('Z://PD_motor//subj_data//beta_rebound//pow_summary//summary_pow.txt',sep=';',header=F)
names(data ) <- c('id','group','session','base','desync','rebound')
data$id <- as.factor(data$id)
data$group <- as.factor(data$group)
data$session <- as.factor(data$session)
levels(data$group) <- c('PD','ctrl')

library(xlsx)
raw.UPDRS <- read.xlsx('Z://PD_motor//subj_data//UPDRS_raw.xlsx',1,header=T)
raw.UPDRS$n.id = as.factor(raw.UPDRS$id)
raw.UPDRS$session <- as.factor(raw.UPDRS$session)
# levels(raw.UPDRS$session) <- c('off','on')

# raw.UPDRS$F1 <- as.numeric(as.character(raw.UPDRS$F1))
# raw.UPDRS$F2 <- as.numeric(as.character(raw.UPDRS$F2))
# raw.UPDRS$F3 <- as.numeric(as.character(raw.UPDRS$F3))
# raw.UPDRS$F4 <- as.numeric(as.character(raw.UPDRS$F4))
# raw.UPDRS$F5 <- as.numeric(as.character(raw.UPDRS$F5))
# raw.UPDRS$F6 <- as.numeric(as.character(raw.UPDRS$F6))
# raw.UPDRS$F7 <- as.numeric(as.character(raw.UPDRS$F7))

load(file='Z://PD_motor//subj_data//alldata.RData')

# Prepare data
updrs.long <- data.frame(id = rep(alldata$MEG_ID,2),
                         UPDRS.old = c(alldata$UPDRS_off,alldata$UPDRS_on),
                         group = rep(alldata$Sub_type,2),
                         session = rep(c(1,2),each=length(alldata$MEG_ID)),
                         n.id = rep(alldata$ID,2),
                         hand = rep(alldata$PAM_hand,2),
                         HY_stage = rep(alldata$HY_stage,2))

updrs.data <- merge(updrs.long,raw.UPDRS,by.x=c('n.id','session'), by.y=c('n.id','session'))

u.data <- merge(data,updrs.data, by.x=c('id','session'), by.y=c('id.x','session'))

# Flip lateralized factors for left hand subjects
F4.flip <- ifelse(u.data$hand=="left",u.data$F5,u.data$F4)
F5.flip <- ifelse(u.data$hand=="left",u.data$F4,u.data$F5)

u.data$F4 <- F4.flip
u.data$F5 <- F5.flip

save(u.data,file='uData.Rdata')

### New division ###
# Another way to divide the MDS-UPDRS-III is: 
# .	tremor (sum of items 15-18)
u.data$tremor = u.data$X3.15_postTremorHand_R+
  u.data$X3.15_postTremorHand_L+
  u.data$X3.16_kinTremorHands_L+
  u.data$X3.16_kinTremorHands_R+
  u.data$X3.17_restTremorAmp_L.J+
  u.data$X3.17_restTremorAmp_RUE+
  u.data$X3.17_restTremorAmp_LUE+
  u.data$X3.17_restTremorAmp_RLE+
  u.data$X3.17_restTremorAmp_LLE+
  u.data$X3.18_consRestTremor
# .	rigidity (item 3)
u.data$rigid = u.data$X3.3_rigidity_neck+
  u.data$X3.3_rigidity_RUE+
  u.data$X3.3_rigidity_LUE+
  u.data$X3.3_rigidity_RLE+
  u.data$X3.3_rigidity_LLE
# .	bradykinesia (sum of items 2, 4-9 and 14) 
u.data$brady = u.data$X3.2_facial_exp+
  u.data$X3.4_fingerTap_L+
  u.data$X3.4_fingerTap_R+
  u.data$X3.5_moveHands_L+
  u.data$X3.5_movemHands_R+
  u.data$X3.6_pro.supMoveHands_L+
  u.data$X3.6_pro.supMoveHands_R+
  u.data$X3.7_toeTap_R+
  u.data$X3.7_toeTap_L+
  u.data$X3.8_legAgil_L+
  u.data$X3.8_legAgil_R+
  u.data$X3.9_chair+
  u.data$X3.14_globalSpont
# .	axial (sum of items 1 and 9-13). 
u.data$axial = u.data$X3.11_freeze+
  u.data$X3.9_chair+
  u.data$X3.10_gait+
  u.data$X3.11_freeze+
  u.data$X3.12_postSablil+
  u.data$X3.13_posture

PD.data = subset(u.data, u.data$group.x=="PD")
PD.data$group.x <- factor(PD.data$group.x)
PD.data$id <- factor(PD.data$id)

# Super-long data
library(reshape2)
u.long <- melt(PD.data,
               # ID variables - all the variables to keep but not split apart on
               id.vars=c("id","rebound","session"),
               # The source columns
               measure.vars=c("F1", "F2", "F3", "F4", "F5", "F6", "F7" ),
                # Name of the destination column that will identify the original
                # column that the measurement came from
                variable.name="factor",
                value.name="score"
)

u.long2 <- melt(PD.data,
               # ID variables - all the variables to keep but not split apart on
               id.vars=c("id","rebound","session"),
               # The source columns
               measure.vars=c("tremor", "rigid", "brady", "axial" ),
               # Name of the destination column that will identify the original
               # column that the measurement came from
               variable.name="factor",
               value.name="score"
)

# Get summary of tremor
aggregate(u.long2$score, by=list(session=u.long2$session, factor=u.long2$factor),mean)
aggregate(u.long2$score, by=list(session=u.long2$session, factor=u.long2$factor),median)
aggregate(u.long2$score, by=list(session=u.long2$session, factor=u.long2$factor),range)


#################### Stats (mean) ##############################
load(file='uData.Rdata')

# Not used
library(lme4)

all.mod1 <- lmer(rebound~Total + (1|id), data=u.data, REML=F, subset = u.data$group.x=="PD")
summary(all.mod1)
all.mod0 <- lmer(rebound~1 + (1|id), data=u.data, REML=F, subset = u.data$group.x=="PD")
anova(all.mod0,all.mod1)
all.mod2 <- lmer(rebound~Total+session + (1|id), data=u.data, REML=F, subset = u.data$group.x=="PD")
anova(all.mod0,all.mod1,all.mod2)
all.mod3 <- lmer(rebound~Total*session + (1|id), data=u.data, REML=F, subset = u.data$group.x=="PD")
anova(all.mod0,all.mod1,all.mod2,all.mod3)


# Summary
# pam1patients <- subset(pam1data,Sub_type=='patient')
t.test(PD.data$Total~PD.data$session, paired=T)
ttestBF(x=PD.data$Total[PD.data$session=="1"],y=PD.data$Total[PD.data$session=="2"], data=pam1data, paired=T, rscale='wide')

t.test(x=pam1data$UPDRS_off[pam1data$Sub_type=='control'] , y=pam1patients$UPDRS_on, paired=F)
ttestBF(x=pam1data$UPDRS_off[pam1data$Sub_type=='control'], y=PD.data$Total[PD.data$session=="2"], paired=F, rscale='wide')
ttestBF(x=pam1data$UPDRS_off[pam1data$Sub_type=='control'], y=PD.data$Total[PD.data$session=="1"], paired=F, rscale='wide')

aggregate(u.data$Total, by=list(u.data$group.x,u.data$session), mean)
aggregate(u.data$Total, by=list(u.data$group.x,u.data$session), sd)

# Bayes factor (used)
load(file='uData.Rdata')
library(BayesFactor)

bf0 <- lmBF(rebound~id, data=PD.data, whichRandom="id")
bf1 <- lmBF(rebound~session+id, data=PD.data, whichRandom="id")
bf1/bf0

PD.data.x <- PD.data[complete.cases(PD.data),]
bf0 <- lmBF(rebound~id, data=PD.data, whichRandom="id")
bf1 <- lmBF(rebound~session*id, data=PD.data, whichRandom=c("id","session"))
bfT <- lmBF(rebound~Total+session+id, data=PD.data, whichRandom=c("id","session"))

bfT/bf1

bf.f1 <- lmBF(rebound~F1+session*id, whichRandom=c("id","session"),data=PD.data)
bf.f1/bf1
bf.f2 <- lmBF(rebound~F2+session*id, whichRandom=c("id","session"),data=PD.data)
bf.f2/bf1
bf.f3 <- lmBF(rebound~F3+session*id, whichRandom=c("id","session"),data=PD.data)
bf.f3/bf1
bf.f4 <- lmBF(rebound~F4+session*id, whichRandom=c("id","session"),data=PD.data)
bf.f4/bf1
bf.f5 <- lmBF(rebound~F5+session*id, whichRandom=c("id","session"),data=PD.data)
bf.f5/bf1
bf.f6 <- lmBF(rebound~F6+session*id, whichRandom=c("id","session"),data=PD.data)
bf.f6/bf1
bf.f7 <- lmBF(rebound~F7+session*id, whichRandom=c("id","session"),data=PD.data)
bf.f7/bf1

mF1.post <- posterior(bf.f1, iterations = 2000)
summary(mF1.post)
mF2.post <- posterior(bf.f2, iterations = 2000)
summary(mF2.post)
mF3.post <- posterior(bf.f3, iterations = 2000)
summary(mF3.post)
mF4.post <- posterior(bf.f4, iterations = 2000)
summary(mF4.post)
mF5.post <- posterior(bf.f5, iterations = 2000)
summary(mF5.post)
mF6.post <- posterior(bf.f6, iterations = 2000)
summary(mF6.post)
mF7.post <- posterior(bf.f7, iterations = 2000)
summary(mF7.post)

mT.post <- posterior(bfT, iterations = 2000)
summary(mT.post)

# Conventional correlation (just for show)
cor.test(PD.data.x$F1,PD.data.x$rebound)
cor.test(PD.data.x$F2,PD.data.x$rebound)
cor.test(PD.data.x$F3,PD.data.x$rebound)
cor.test(PD.data.x$F4,PD.data.x$rebound)
cor.test(PD.data.x$F5,PD.data.x$rebound)
cor.test(PD.data.x$F6,PD.data.x$rebound)
cor.test(PD.data.x$F7,PD.data.x$rebound)
cor.test(PD.data.x$Total,PD.data.x$rebound)

##
save.image(".RData")

### New division ###
bf0 <- lmBF(rebound~id, data=PD.data, whichRandom="id")
bf1 <- lmBF(rebound~session*id, data=PD.data, whichRandom=c("id","session"))

bf.tremor <- lmBF(rebound~tremor+session*id, whichRandom=c("id","session"),data=PD.data)
bf.tremor/bf1
bf.tremor.post <- posterior(bf.tremor, iterations = 2000)

bf.rigid <- lmBF(rebound~rigid+session*id, whichRandom=c("id","session"),data=PD.data)
bf.rigid/bf1
bf.rigid.post <- posterior(bf.rigid, iterations = 2000)

bf.axial <- lmBF(rebound~axial+session*id, whichRandom=c("id","session"),data=PD.data)
bf.axial/bf1
bf.axial.post <- posterior(bf.axial, iterations = 2000)

bf.brady <- lmBF(rebound~brady+session*id, whichRandom=c("id","session"),data=PD.data)
bf.brady/bf1
bf.brady.post <- posterior(bf.brady, iterations = 2000)

# Score
bf0 <- lmBF(score~id, data=u.long2, whichRandom="id")
bf1 <- lmBF(score ~ factor+session*id, data=u.long2, whichRandom=c("id","session"))
bf.F <- lmBF(score ~ factor+rebound+session*id, data=u.long2, whichRandom=c("id","session"))
bf.F/bf1
bf.Fx <- lmBF(score ~ factor*rebound+session*id, data=u.long2, whichRandom=c("id","session"))
bf.Fx/bf1

mT.post <- posterior(bf.Fx, iterations = 2000)
summary(mT.post)

# Conventional correlation (just for show)
cor.test(PD.data$tremor,PD.data$rebound)
cor.test(PD.data$rigid,PD.data$rebound)
cor.test(PD.data$axial,PD.data$rebound)
cor.test(PD.data$brady,PD.data$rebound)


plot(PD.data$tremor,PD.data$rebound)
plot(PD.data$rigid,PD.data$rebound)
plot(PD.data$axial,PD.data$rebound)
plot(PD.data$brady,PD.data$rebound)


## Each item
bf0 <- lmBF(rebound~session*id, data=PD.data, whichRandom=c("id","session"))

bf.x1 <- lmBF(rebound~X3.1+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x1/bf0
bf.x2 <- lmBF(rebound~X3.2+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x2/bf0
bf.x3 <- lmBF(rebound~X3.3+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x3/bf0
bf.x4 <- lmBF(rebound~X3.4+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x4/bf0
bf.x5 <- lmBF(rebound~X3.5+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x5/bf0
bf.x6 <- lmBF(rebound~X3.6+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x6/bf0
bf.x7 <- lmBF(rebound~X3.7+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x7/bf0
bf.x8 <- lmBF(rebound~X3.8+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x8/bf0
bf.x9 <- lmBF(rebound~X3.9+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x9/bf0
bf.x10 <- lmBF(rebound~X3.10+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x10/bf0
bf.x11 <- lmBF(rebound~X3.11+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x11/bf0
bf.x12 <- lmBF(rebound~X3.12+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x12/bf0
bf.x13 <- lmBF(rebound~X3.13+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x13/bf0
bf.x14 <- lmBF(rebound~X3.14+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x14/bf0
bf.x15 <- lmBF(rebound~X3.15+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x15/bf0
bf.x16 <- lmBF(rebound~X3.16+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x16/bf0
bf.x17 <- lmBF(rebound~X3.17+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x17/bf0
bf.x18 <- lmBF(rebound~X3.18+session*id, whichRandom=c("id","session"),data=PD.data)
bf.x18/bf0

############################ PLOTS #######################################
library(ggplot2)

### Plot F1
N=dim(mF1.post)[1]
x = 0:max(PD.data$F1)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF1.post[i,1]+mF1.post[i,2]*x}
df = data.frame(F1=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$F1)+1)

g1 <- ggplot(df, aes(F1, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F1, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F1))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('F1: Midline function')


### Plot F2
N=dim(mF2.post)[1]
x = 0:max(PD.data$F2)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF2.post[i,1]+mF2.post[i,2]*x}
df = data.frame(F2=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$F2)+1)

g2 <- ggplot(df, aes(F2, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F2, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F2))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('F2: Rest tremor')

### Plot F3
N=dim(mF3.post)[1]
x = 0:max(PD.data$F3)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF3.post[i,1]+mF3.post[i,2]*x}
df = data.frame(F3=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$F3)+1)

g3 <- ggplot(df, aes(F3, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F3, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F3))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('F3: Rigidity')


### Plot F4
N=dim(mF4.post)[1]
x = 0:max(PD.data$F4)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF4.post[i,1]+mF4.post[i,2]*x}
df = data.frame(F4=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$F4)+1)

g4 <- ggplot(df, aes(F4, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F4, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F4))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('F4: Bradykinesia right')

### Plot F5
N=dim(mF5.post)[1]
x = 0:max(PD.data$F5)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF5.post[i,1]+mF5.post[i,2]*x}
df = data.frame(F5=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$F5)+1)

g5 <- ggplot(df, aes(F5, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F5, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F5))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('F5: Bradykinesia left')

### Plot F6
N=dim(mF6.post)[1]
x = 0:max(PD.data$F6)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF6.post[i,1]+mF6.post[i,2]*x}
df = data.frame(F6=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$F6)+1)

g6 <- ggplot(df, aes(F6, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F6, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F6))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('F6: Postural and kinetic tremor')

### Plot F7
N=dim(mF7.post)[1]
x = 0:max(PD.data$F7)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF7.post[i,1]+mF7.post[i,2]*x}
df = data.frame(F7=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$F7)+1)

g7 <- ggplot(df, aes(F7, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F7, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F7))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('F7: Lower limb bradykinesia')

### Plot Total
N=dim(mT.post)[1]
x = 0:max(PD.data$Total)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mT.post[i,1]+mT.post[i,2]*x}
df = data.frame(Total=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$Total)+1)

gT <- ggplot(df, aes(Total, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=Total, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$Total))+ylim(-0.2,0.4)+
  ggtitle('Total UPDRS-III')+
  xlab('Total score') + ylab('Beta rebound %-change')


# Add layout
lay <- theme_bw() + theme(legend.position = "none",
                          text = element_text(size = 11, family="Helvetica"),
                          title = element_text(size = 12, vjust = 1.5, face="bold",lineheight = NULL),
                          # axis.text = element_text(size=10),
                          axis.title = element_text(size = 11, vjust = .5, face="plain"),
                          # axis.title.x = element_text(face="plain", size=4),
                          # axis.title.y = element_text(face="plain", size=4),
                          axis.text = element_text(face="bold", size=11),
                          panel.grid = element_blank())
  
  
g1 <- g1 + lay
g2 <- g2 + lay
g3 <- g3 + lay
g4 <- g4 + lay
g5 <- g5 + lay
g6 <- g6 + lay
g7 <- g7 + lay
gT <- gT + lay
  
ggsave("F1.jpeg", plot=g1, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F2.jpeg", plot=g2, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F3.jpeg", plot=g3, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F4.jpeg", plot=g4, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F5.jpeg", plot=g5, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F6.jpeg", plot=g6, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F7.jpeg", plot=g7, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("Total.jpeg", plot=gT, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)

save.image(".RData")


## Plot new division
### Plot F1
N=dim(bf.tremor.post)[1]
x = 0:max(PD.data$tremor)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = bf.tremor.post[i,1]+bf.tremor.post[i,2]*x}
df = data.frame(tremor=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$tremor)+1)

ggplot(df, aes(tremor, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=tremor, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$tremor))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('Tremor')

N=dim(bf.rigid.post)[1]
x = 0:max(PD.data$rigid)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = bf.rigid.post[i,1]+bf.rigid.post[i,2]*x}
df = data.frame(rigid=rep(x,N),rebound=unlist(y))
df$f = rep(1:N,each=max(PD.data$rigid)+1)

ggplot(df, aes(rigid, rebound)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=rigid, y=rebound, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$rigid))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta rebound %-change')+
  ggtitle('Rigid')




