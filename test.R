# test # test code (by JJJ)
rm(list = ls()) ; gc()
source("./library/comp_logistic.R")
group.index <- c(1,1,2,2,2,3,4,4,4,4)
n = 100
p = length(group.index)  
set.seed(1)
x = matrix(rnorm(n*p),n,p)
y = sample(c(0,1),n,replace = T)
# com
complasso.Logistic(y, x, group.index,lambda = 1)

# compare the result of logistic lasso
library(glmnet)
glmnet(x,y, family = "binomial", lambda = 5/n, standardize = F)$beta
complasso.Logistic(y, x, 1:p,lambda = 5)

#  compLogistic
complasso.Logistic(y, x, group.index,lambda = 5)
# group lasso
library(grpreg)
drop(grpreg(x,y,group = group.index, penalty = "grLasso",
       family = "binomial", lambda = 5/n)$beta)


