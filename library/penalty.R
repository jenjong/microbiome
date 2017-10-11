# penalty function

# simulation code
# group information (should be sorted!)
group.index <- c(1,1,2,2,2,3,4,4,4,4)
# p = 10
n = 100
p = length(group.index)  
set.seed(1)
x = matrix(rnorm(n*p),n,p)
y = sample(c(0,1),n,replace = T)

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
colnames(Dmat) <- coefName[- gloc.e[gsize>1]]

# insert intercept term


# Desigm matrix
for (i in 1:length(gidx))
{
  if( gsize[i] == 1 ) next
  x[,(gloc.s[i]:gloc.e[i])] <- x[,(gloc.s[i]:gloc.e[i])] - x[, gloc.e[i]]
}
x <- x[, - gloc.e[gsize>1]]

# genlasso algorithm
# solve the problem for each lambda? inefficient!
beta.r <- rep(0, ncol(x))

  xb <- drop(x%*%beta.r)
  theta.vec<- 1/(1+exp(-xb))
  tilde.y <- 2*(y - theta.vec)
  
  fit = genlasso(y,x,Dmat[-2,])
  coef.genlasso(fit, lambda=1)
#fit = genlasso(y,x,D)



# check convergence

# intercept 

# recover coeff
  # beta.r <- rnorm(7)
  # beta.vec <- rep(0,p)
beta.vec[ -gloc.e[gsize>1] ] <- beta.r
i1 <- 1
for (i in 1:length(gidx))
{
  if( gsize[i] == 1 ) next
  beta.vec[gloc.e[i]] <- 1-sum( beta.r[i1:(i1+gsize[i]-2)] )
  i1 <- i1 + gsize[i]-1
}
