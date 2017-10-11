# penalty function

# simulation code
# group information (should be sorted!)
group.index <- c(1,1,2,2,2,3,4,4,4,4)
n = 100
p = length(group.index)  
set.seed(1)
x = matrix(rnorm(n*p),n,p)
y = sample(c(0,1),n,replace = T)
lambda.var <- 1
max.iter = 10
tol.var = 1e-4
# check the group information and design matrix
gidx <- unique(group.index)
gsize <- as.integer( table(group.index) )
gloc.s <- cumsum(c(1,gsize[-length(gsize)]))
gloc.e <- cumsum(gsize)
group_info <- data.frame(group_Index = gidx, group_Size = gsize,
                         loc1 = gloc.s, loc2 = gloc.e)
# coefficients name
coefName <- paste("b",rep(gidx, gsize), 
                  unlist(lapply(gsize, FUN = seq, from = 1, by = 1)),
                  sep = '_')

# penalty matrix
Dmat <- diag(1,p)
for (i in 1:length(gidx))
{
  if( gsize[i] <= 2 ) next
  Dmat[gloc.e[i], (gloc.s[i]:gloc.e[i])] <- 1
}
Dmat <- Dmat[, - gloc.e[gsize>1]]
Dmat <- Dmat[-gloc.e[gsize == 2],]
# insert intercept term
Dmat <- cbind(0,Dmat)
# names of coefficients
colnames(Dmat) <- c("intercept", coefName[- gloc.e[gsize>1]])

# Desigm matrix
for (i in 1:length(gidx))
{
  if( gsize[i] == 1 ) next
  x[,(gloc.s[i]:gloc.e[i])] <- x[,(gloc.s[i]:gloc.e[i])] - x[, gloc.e[i]]
}
x <- x[, - gloc.e[gsize>1]]
# intercept
x<- cbind(1,x)

# genlasso algorithm
# initialize
beta.r <- beta.c <- rep(0, ncol(x))

  for( iter in 1:max.iter)
  {
    xb <- drop(x%*%beta.r)
    theta.vec<- 1/(1+exp(-xb))
    tilde.y <- 2*(y - theta.vec)
    fit = genlasso(tilde.y, x, Dmat)
    gamma.r <- drop(coef.genlasso(fit, lambda.var)$beta)  
    beta.r <- 2*gamma.r + beta.r
    if( iter>10 & iter%%10==0 )
    {
    sol.diff <- norm(beta.c - beta.r,"2")/norm(beta.c,"2") 
    if (sol.diff < tol.var) break else beta.c<-beta.r
    }
  }
# truncate solution
beta.r[abs(beta.r) < 1e-8] = 0
# intercept 
beta.0 <- beta.r[1]
beta.r <- beta.r[-1]

# recover coeff
beta.vec <- rep(0,p)
beta.vec[ -gloc.e[gsize>1] ] <- beta.r
i1 <- 1
for (i in 1:length(gidx))
{
  if( gsize[i] == 1 ) next
  beta.vec[gloc.e[i]] <- 1-sum( beta.r[i1:(i1+gsize[i]-2)] )
  i1 <- i1 + gsize[i]-1
}

beta.vec
beta.0
