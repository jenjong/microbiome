# penalty function
makeD <- function(group.index)
{
  gidx <- unique(group.index)
  gsize <- as.integer( table(group.index) )
  gloc.s <- cumsum(c(1,gsize[-length(gsize)]))
  gloc.e <- cumsum(gsize)
  group_info <- data.frame(group_Index = gidx, group_Size = gsize,
                           loc1 = gloc.s, loc2 = gloc.e)
  coefName <- paste("b",rep(gidx, gsize), 
                    unlist(lapply(gsize, FUN = seq, from = 1, by = 1)),
                    sep = '_')
  
  # penalty matrix
  Dmat <- diag(1,length(group.index))
  for (i in 1:length(gidx))
  {
    if( gsize[i] <= 2 ) next
    Dmat[gloc.e[i], (gloc.s[i]:gloc.e[i])] <- 1
  }
  if (any(gsize>1)) Dmat <- Dmat[, - gloc.e[gsize>1], drop = F]
  if (any(gsize == 2)) Dmat <- Dmat[-gloc.e[gsize == 2],,drop = F]
  # insert intercept term
  Dmat <- cbind(0,Dmat)
  # names of coefficients
  if (any(gsize>1)) colnames(Dmat) <- c("intercept", coefName[- gloc.e[gsize>1]])
  return(Dmat)  
}

# Design
makeX <- function(x, group.index)
  
{
  gidx <- unique(group.index)
  gsize <- as.integer( table(group.index) )
  gloc.s <- cumsum(c(1,gsize[-length(gsize)]))
  gloc.e <- cumsum(gsize)
  
  for (i in 1:length(gidx))
  {
    if( gsize[i] == 1 ) next
    x[,(gloc.s[i]:gloc.e[i])] <- x[,(gloc.s[i]:gloc.e[i])] - x[, gloc.e[i]]
  }
  if ( any(gsize>1) ) x <- x[, - gloc.e[gsize>1]]
  # intercept
  x<- cbind(1,x)
  return(x)
}


# compositional logistic regression
compLogistic.lasso <- function(y,x,group.index,lambda=1,
                               max.iter = 100, tol.var = 1e-4)
{
  n = nrow(x)
  p = ncol(x)
  lambda.var = lambda
  if (length(group.index) != p) stop("the length of group.index should be equal to ncol(x)!")
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
  
  if (any(gsize == 1)) 
  {
    tmp <- coefName[gloc.s[which(gsize==1)]]
    coefName[gloc.s[which(gsize==1)]] <- substr(tmp,1, nchar(tmp)-2   )
  }
  
  
  X <- makeX(x, group.index)
  tilde.X <- X/2
  Dmat <-makeD(group.index)
  beta.r <- beta.c <- rep(0, ncol(X))
  for( iter in 1:max.iter)
  {
    xb <- drop(X%*%beta.r)
    theta.vec<- 1/(1+exp(-xb))
    tilde.y <- 2*(y - theta.vec) + drop(tilde.X%*%beta.r)
    fit = genlasso(tilde.y, tilde.X, Dmat)
    beta.r <- drop(coef.genlasso(fit, lambda.var)$beta)  
    if( iter>10 & iter%%10==0 )
    {
      sol.diff <- norm(beta.c - beta.r,"2")/norm(beta.c,"2") 
      if (sol.diff < tol.var) break else beta.c<-beta.r
    }
  }
  
  beta.r[abs(beta.r) < 1e-8] = 0
  # intercept 
  beta.0 <- beta.r[1]
  beta.r <- beta.r[-1]
  
  # recover coeff
  beta.vec <- rep(0,p)
  if (any(gsize>1)) beta.vec[ - gloc.e[gsize>1] ] <- beta.r  else beta.vec = beta.r
  
  i1 <- 1
  for (i in 1:length(gidx))
  {
    if( gsize[i] == 1 ) next
    beta.vec[gloc.e[i]] <-  -sum( beta.r[i1:(i1+gsize[i]-2)] )
    i1 <- i1 + gsize[i]-1
  }
  names(beta.vec) <- coefName
  return( list(Intercept = beta.0, beta.vec = beta.vec))
}
  