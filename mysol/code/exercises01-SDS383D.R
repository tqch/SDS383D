library(Matrix)
library(microbenchmark)
library(knitr)


lm_inv <- function(X,y,weights){
  W = diag(weights)
  A = t(X)%*%W%*%X
  b = t(X)%*%W%*%y
  A_inv = solve(A)
  return(A_inv%*%b)
}

lm_LU <- function(X, y, weights){
  W = diag(weights)
  A = t(X)%*%W%*%X
  p = nrow(A)
  b = t(X)%*%W%*%y
  eluA = expand(lu(A))
  b = b[sort(apply(eluA$P,1,function(x){which(x!=0)}),index.return=T)$ix]
  b = forwardsolve(eluA$L,b)
  b = backsolve(eluA$U,b)
  return(b)
}

lm_chol <- function(X, y, weights){
  W = diag(weights)
  A = t(X)%*%W%*%X
  p = nrow(A)
  b = t(X)%*%W%*%y
  U = base::chol(A)
  b = forwardsolve(t(U),b)
  b = backsolve(U,b)
  return(b)
}

bm_setup <- function(n, p){
  beta = rnorm(p)
  sigma = rexp(n)
  X = matrix(rnorm(n*(p-1)), nrow=n)
  X = cbind(rep(1,n), X)
  y = X%*%beta+sigma*rnorm(n)
  weights = 1/sigma^2
}

ps = c(50,100,500,1000)
bm_inv = list()
for(i in 1:4){
  p = ps[i]
  n = 2*p
  cat(sprintf("p=%d:\n",p))
  bm_inv[[i]] = microbenchmark(lm_inv(X,y,weights),times=10,setup={
    beta = rnorm(p)
    sigma = rexp(n)
    X = matrix(rnorm(n*(p-1)), nrow=n)
    X = cbind(rep(1,n), X)
    y = X%*%beta+sigma*rnorm(n)
    weights = 1/sigma^2
  })
}
df_inv = do.call(rbind,lapply(bm_inv,function(x){summary(x,unit="ms")[,c(-1,-8)]}))
df_inv = data.frame(p=ps, df_inv)
tab_inv = kable(df_inv, "latex", booktabs=T)

bm_lu = list()
for(i in 1:4){
  p = ps[i]
  n = 2*p
  cat(sprintf("p=%d:\n",p))
  bm_lu[[i]] = microbenchmark(lm_LU(X,y,weights),times=10,setup={
    beta = rnorm(p)
    sigma = rexp(n)
    X = matrix(rnorm(n*(p-1)), nrow=n)
    X = cbind(rep(1,n), X)
    y = X%*%beta+sigma*rnorm(n)
    weights = 1/sigma^2
  })
}
df_lu = do.call(rbind,lapply(bm_lu,function(x){summary(x,unit="ms")[,c(-1,-8)]}))
df_lu = data.frame(p=ps, df_lu)
tab_lu = kable(df_lu, "latex", booktabs=T)

bm_chol = list()
for(i in 1:4){
  p = ps[i]
  n = 2*p
  cat(sprintf("p=%d:\n",p))
  bm_chol[[i]] = microbenchmark(lm_inv(X,y,weights),times=10,setup={
    beta = rnorm(p)
    sigma = rexp(n)
    X = matrix(rnorm(n*(p-1)), nrow=n)
    X = cbind(rep(1,n), X)
    y = X%*%beta+sigma*rnorm(n)
    weights = 1/sigma^2
  })
}
df_chol = do.call(rbind,lapply(bm_chol,function(x){summary(x,unit="ms")[,c(-1,-8)]}))
df_chol = data.frame(p=ps, df_chol)
tab_chol = kable(df_chol, "latex", booktabs=T)