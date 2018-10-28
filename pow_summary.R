setwd('Z://PD_motor//subj_data//beta_rebound//pow_summary') #NB! This is run from Windows partion

data <- read.csv('summary_pow.txt',sep=';',header=F)
names(data ) <- c('id','group','session','base','desync','rebound')
data$id <- as.factor(data$id)
data$group <- as.factor(data$group)
data$session <- as.factor(data$session)
levels(data$group) <- c('PD','ctrl')

library(lme4)
bslneDiff <- lmer(base~group*session+(1|id),data=data,REML=F)
summary(bslneDiff)

# desync ~ baseline
desyncAnalys <- lmer(desync ~ base*group*session+(1|id),data=data, REML=F)
summary(desyncAnalys)

desyncAnalys1 <- update(desyncAnalys, .~. -base:group:session)
anova(desyncAnalys,desyncAnalys1)
desyncAnalys2 <- update(desyncAnalys1, ~. -base:group)
anova(desyncAnalys1,desyncAnalys2)
desyncAnalys3 <- update(desyncAnalys2, ~. -base:session)
anova(desyncAnalys2,desyncAnalys3)
desyncAnalys4 <- update(desyncAnalys3, ~. -session:group)
anova(desyncAnalys3,desyncAnalys4)
desyncAnalys5 <- update(desyncAnalys4, ~. -session)
anova(desyncAnalys5,desyncAnalys4)
desyncAnalys6 <- update(desyncAnalys5, ~. -group)
anova(desyncAnalys6,desyncAnalys5)
desyncAnalys7 <- update(desyncAnalys6, ~. -base)
anova(desyncAnalys7,desyncAnalys6) #                            [!]

reboundAnalys <- lmer(rebound ~ base*group*session+(1|id),data=data,REML=F)
reboundAnalysC <- lmer(rebound ~ I(base-mean(base))*group*session+(1|id),data=data,REML=F)

summary(reboundAnalys)

reboundAnalys1 <- update(reboundAnalys, ~. -base:group:session)
anova(reboundAnalys,reboundAnalys1)
reboundAnalys2 <- update(reboundAnalys1, ~. -base:group)
anova(reboundAnalys2,reboundAnalys1)
reboundAnalys3 <- update(reboundAnalys2, ~. -base:session)
anova(reboundAnalys3,reboundAnalys2)
reboundAnalys4 <- update(reboundAnalys3, ~. -session:group)
anova(reboundAnalys4,reboundAnalys3) #                          [!]
reboundAnalys5 <- update(reboundAnalys3, ~. -base)
anova(reboundAnalys5,reboundAnalys3)

### Bayes ###
library(BayesFactor)

## Test difference in baseline
ttestBF(formula=base~group, data=data, paired=F, rscale='wide')

# Test effect of baseline on rebound
mod0 <- lmBF(rebound ~ id, data=data, whichRandom="id")
mod1a <- lmBF(rebound ~ id + group, data=data, whichRandom="id")
mod1b <- lmBF(rebound ~ id + base, data=data, whichRandom="id")
mod2 <- lmBF(rebound ~ id + base+group, data=data, whichRandom="id")
mod2b <- lmBF(rebound ~ id + group:session, data=data, whichRandom="id")
mod3 <- lmBF(rebound ~ id + base+group:session, data=data, whichRandom="id" )
mod4 <- lmBF(rebound ~ id + base:group:session, data=data, whichRandom="id")
mod4b <- lmBF(rebound ~ id + base:group+group:session, data=data, whichRandom="id")

bf0 <- mod1b/mod0  #Baseline without any group info
bf0
bf2 <- mod2/mod1a  # Test effect of baseline vs group
bf2
bf3 <- mod3/mod2b  # Test effect of baseline vs group:session
bf3
bf4 <- mod4/mod2b  # Test effect of full-interaction model vs group:session
bf4

posterior()


# Test effect of beta ERD on rebound
mod0 <- lmBF(rebound ~ id, data=data, whichRandom="id") 
mod1a <- lmBF(rebound ~ id + group, data=data, whichRandom="id")
mod1b <- lmBF(rebound ~ id + desync, data=data, whichRandom="id")
mod2 <- lmBF(rebound ~ id + desync+group, data=data, whichRandom="id")
mod2b <- lmBF(rebound ~ id + group:session, data=data, whichRandom="id")
mod3 <- lmBF(rebound ~ id + desync+group:session, data=data, whichRandom="id")
mod4 <- lmBF(rebound ~ id + desync:group:session, data=data, whichRandom="id")

bf0 <- mod1b/mod0  #ERD without any group info
bf0
bf2 <- mod2/mod1a  # Test effect of baseline pwr vs group
bf2
bf3 <- mod3/mod2b  # Test effect of baseline pwr vs group+session
bf3
bf4 <- mod4/mod2b  # Test effect and interaction of baeline vs group+session
bf4

# Test effect of baseline on beta ERD
mod0 <- lmBF(desync ~ id, data=data, whichRandom="id") 
mod1a <- lmBF(desync ~ id + group, data=data, whichRandom="id")
mod1b <- lmBF(desync ~ id + base, data=data, whichRandom="id")
mod2 <- lmBF(desync ~ id + base+group, data=data, whichRandom="id")
mod2b <- lmBF(desync ~ id + group+session, data=data, whichRandom="id")
mod3 <- lmBF(desync ~ id + base+group+session, data=data, whichRandom="id")

bf0 <- mod1b/mod0  #Baseline without any group info
bf0
bf2 <- mod2/mod1a  # Test effect of baseline vs group
bf2
bf3 <- mod3/mod2b  # Test effect of baeline vs group+session
bf3

## Plot for inspection
library(ggplot2)

ggplot(data, aes(base, desync, color=group))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  theme_bw()

ggplot(data, aes(base, rebound, color=group,shape=session))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  theme_bw()

ggplot(data, aes(group, base, color=session))+
  geom_point(position=position_dodge(width=.5))
ggplot(data, aes(group, desync, color=session))+
  geom_point(position=position_dodge(width=.5))

ggplot(data, aes(group, rebound, color=session))+
  geom_point(position=position_dodge(width=.5))

## Pretty plot
rebound.sum <- aggregate(data$rebound, by=list(group=data$group,session=data$session),mean)
rebound.sum$median <- aggregate(data$rebound, by=list(group=data$group,session=data$session),median)$x
rebound.sum$sd <- aggregate(data$rebound, by=list(group=data$group,session=data$session),sd)$x
rebound.sum$qt <- aggregate(data$rebound, by=list(group=data$group,session=data$session),quantile, probs=c(0.05,0.95))$x
rebound.sum$qt05 <- rebound.sum$qt[,1]
rebound.sum$qt95 <- rebound.sum$qt[,2]


ggplot(data,aes(x=group, y=rebound, color=id, group=session))+
  geom_point(position = position_jitterdodge(jitter.width = 1,dodge.width=0.5),size=2)+
  # stat_summary(fun.y="median",geom="crossbar",
  #              mapping=aes(ymin=..y.., ymax=..y..), width=0.5,
  #              position=position_dodge(width = 0.3),show.legend = FALSE)+
  stat_summary(fun.y="mean",geom="crossbar",
               mapping=aes(ymin=..y.., ymax=..y..), width=0.3,
               position=position_dodge(width = 0.3),show.legend = FALSE)
  # stat_summary(fun.y="quantile",fun.args = list(probs=0.05),geom="crossbar",
  #              mapping=aes(ymin=..y.., ymax=..y..), width=0.2,
  #              position=position_dodge(width = 0.3),show.legend = FALSE)+
  # stat_summary(fun.y="quantile",fun.args = list(probs=0.95),geom="crossbar",
  #              mapping=aes(ymin=..y.., ymax=..y..), width=0.2,
  #              position=position_dodge(width = 0.3),show.legend = FALSE)
  # 
  # geom_crossbar(data=rebound.sum,aes(x=group,y=median ,ymin=qt05,ymax=qt95), width = 0.5)+
  # geom_errorbar(data=p.summary,aes(x=Task,y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2)+
  ylab("Beta-rebound change")+
  xlab('Group')+ylim(c(.5, 1))+
  labs(title="Performance", tag="B") + guides(fill=F)+
  publish_theme + theme(legend.position = "none")
ggsave('pct_correct_all.png', pct_corrt,
       device='png',width=4,height=3, units='cm', dpi = 600, scale = 3.5)





