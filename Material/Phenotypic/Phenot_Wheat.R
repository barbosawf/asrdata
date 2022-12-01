###############################
## Wheat Data Experiment 
###############################

rm(list=ls())  # Removes all variables in memory
setwd("C:/Users/sgeza/OneDrive/Desktop/ASReml_UFV/Material/Phenotypic")

rm(list=ls()) 
library(asreml)
library(ASRdata)

pheno.wheat <- ASRdata::pheno.wheat
st.data <- pheno.wheat[pheno.wheat$location == "Alliance",]
head(st.data)

# Reference
m0 <- asreml(fixed=yield~rep+at(check,'1'):gen,
             random=~+at(check,'0'):gen+rep:ibk+rep:row+rep:col,
             residual=~idv(units),
             data=st.data)
summary(m0)$varcomp
vpredict(m0,H2~V4/(V1+V2+V3+V4+V5))
wald.asreml(m0, dendf='numeric')

preds0 <- predict(m0, classify='check:gen', levels=list(check='0'))$pvals
head(preds0,10)

# Asjusted Means Model
m1 <- asreml(fixed=yield~rep+gen,
             random=~rep:ibk+rep:row+rep:col,
             residual=~idv(units),
             data=st.data)
summary(m1)$varcomp
wald.asreml(m1, dendf='numeric')

preds <- predict(m1, classify='gen', vcov=TRUE)
head(preds$pvals)
View(preds$pvals)

