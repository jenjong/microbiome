# test
source("comp_logistic.R")
group.index <- c(1,1,2,2,2,3,4,4,4,4)
n = 100
p = length(group.index)  
set.seed(1)
x = matrix(rnorm(n*p),n,p)
y = sample(c(0,1),n,replace = T)
compLogistic.lasso(y, x, group.index,lambda = 1)
compLogistic.lasso(y, x, 1:p,lambda = 10)

library(glmnet)
# compare the result of logistic lasso
glmnet(x,y, family = "binomial", lambda = 5/n, standardize = F)$beta
compLogistic.lasso(y, x, 1:p,lambda = 5)

#  compLogistic
compLogistic.lasso(y, x, group.index,lambda = 5)
library(grpreg)
# group lasso
grpreg(x,y,group = group.index, penalty = "grLasso",
       family = "binomial", lambda = 5/n)$beta
