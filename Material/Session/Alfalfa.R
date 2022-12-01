###############################
## Alfalfa Experiment 
###############################

rm(list=ls())  # Removes all variables in memory
setwd("C:/Users/sgeza/OneDrive/Desktop/ASReml_UFV/Material/Session")

alfalfa<-read.table("ALFALFA.TXT",header=TRUE,na.string='*')
head(alfalfa)  
summary(alfalfa)
str(alfalfa)

# Creating factors
alfalfa$Variety<-as.factor(alfalfa$Variety)
alfalfa$Block<-as.factor(alfalfa$Block)
str(alfalfa)

# EDA - explatoraty data analysis
hist(alfalfa$Resp)
boxplot(alfalfa$Resp~alfalfa$Block)
table(alfalfa$Variety)
table(alfalfa$Variety,alfalfa$Block)

# Analysis using ASReml
library(asreml)

# asreml(fixed=y~<fixed effects>,
#        random=~<random effects>,
#        residual=<error structure>,
#        data=<mydata>)

# Variety Random
# y = mu + Block + Variety + e

mr <- asreml(fixed=Resp~Block,
             random=~Variety,
             residual=~idv(units),
             data=alfalfa)
plot(mr)
summary(mr)$varcomp

# H2 = s2g/(s2g+s2)
vpredict(mr, H2~V1/(V1+V2))

wald.asreml(mr, denDF='numeric')

# Solutions: BLUE, BLUP
BLUE <- summary(mr, coef=TRUE)$coef.fixed
BLUE
BLUP <- summary(mr, coef=TRUE)$coef.random
BLUP

# Predictions ('LSmeans')
pp <- predict.asreml(mr, classify='Variety', vcov=TRUE)
ls(pp)
pp$pvals # mu^ + Variety^
pp$vcov

# LRT
# H0: H2 = 0, H1: H2 > 0

# 2*(51.7370 - 44.878)

m0 <- asreml(fixed=Resp~Block,
             #random=~Variety,
             residual=~idv(units),
             data=alfalfa)
summary(m0)$varcomp

lrt.asreml(mr,m0, boundary=TRUE)


#######################
