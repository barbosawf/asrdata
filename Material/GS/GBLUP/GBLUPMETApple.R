################################
##  Genomic Selection: GBLUP  ##
################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/ASReml_UFV/Material/GS/GBLUP")
library(asreml)
library(ASRgenomics)

# Reading phenotypic data
datag <- read.table(file='phenotAppleMET.txt', header=TRUE)
head(datag)

# Reading AHAT inverse (sparse format) 
list.files()
load(file='Gmat.RData')

# Checking Ginverse
head(Ginv.bend)
head(attr(Ginv.bend, "rowNames"))
head(attr(Ginv.bend, "colNames"))
attr(Ginv.bend, "INVERSE")

# Matching data
check <- match.kinship2pheno(K=G_bend, pheno.data=datag,
                             indiv='INDIV', clean=FALSE, mism=TRUE)

# Defining factors
head(datag)
datag$INDIV<-as.factor(datag$INDIV)
datag$SITE<-as.factor(datag$SITE)
str(datag)

########################
# Single-SITE GBLUP Model in ASReml-R

# y = mu + INDIV + e
s1 <- datag[datag$SITE == 'HB',]
head(s1)

modelGBLUP<-asreml(fixed=FF~1,
                  random=~vm(INDIV,Ginv.bend),
                  residual=~idv(units),
                  workspace=128e06,
                  na.action=na.method(y="include"),
                  data=s1)
plot(modelGBLUP)
summary(modelGBLUP)$varcomp
vpredict(modelGBLUP,h2_1~V1/(V1+V2))

###################
# MET as Bivariate

# y = mu + SITE + INDIV:SITE + e
modelMET<-asreml(fixed=FF~SITE,
                 random=~vm(INDIV,Ginv.bend):corgh(SITE),
                 residual=~dsum(~units|SITE),
                 workspace=128e06,
                 na.action=na.method(y="include"),
                 data=datag)
modelMET<-update.asreml(modelMET)
summary(modelMET)$varcomp
wald(modelMET)
plot(modelMET)

vpredict(modelMET,h2_HB~V2/(V2+V4))
vpredict(modelMET,h2_MOT~V3/(V3+V5))

BLUP <- summary(modelMET,coef=TRUE)$coef.random
View(BLUP)

preds <- predict.asreml(modelMET, classify='INDIV:SITE')$pvals
head(preds)
predsA <- predict.asreml(modelMET, classify='INDIV')$pvals
head(predsA)

###### EXTRA ####
#################
##################
# Bivariate Model

# Reading Phenotypic data
data(pheno.apple)
datag <- pheno.apple
head(datag)

datag$INDIV <- as.factor(datag$INDIV)

# [y1 y2] = trait + INDIV:trait + e:trait
modelBV<-asreml(fixed=cbind(JUI_HB,FF_HB)~trait,
                random=~vm(INDIV,Ginv.bend):corgh(trait),
                #residual=~id(units):diag(trait),
                residual=~id(units):corgh(trait),
                workspace=128e06,
                na.action=na.method(y="include"),
                data=datag)
plot(modelBV)
summary(modelBV)$varcomp
