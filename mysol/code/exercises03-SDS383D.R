library(magrittr)
library(mvtnorm)

bldg = read.csv("./greenbuildings.csv")
bldg = bldg %>% dplyr::select(
  c("Rent","leasing_rate","green_rating","City_Market_Rent","age","class_a","class_b")
) %>% dplyr::mutate(revenue_per_sqft=Rent*leasing_rate/100)

n = nrow(bldg)
X = cbind(rep(1,n),as.matrix(bldg%>%dplyr::select(
  !c("Rent","leasing_rate","revenue_per_sqft"))))
p = ncol(X)
K = diag(rep(0.01,p))
d = 1
eta = d
y = as.matrix(bldg$revenue_per_sqft)
m = rep(0,p)
K_hat = K+t(X)%*%X
K_hat_inv = solve(K_hat)
m_hat = K%*%matrix(m)+t(X)%*%y # no K_hat_inv%*% in front
df = d+n
eta_star = t(m)%*%K%*%matrix(m)+t(y)%*%matrix(y)-t(m_hat)%*%K_hat_inv%*%m_hat
K_star = df*K_hat/drop(eta_star)
Sigma = solve(K_star)

beta_posterior_mean = K_hat_inv%*%m_hat
fitted = X%*%beta_posterior_mean
residuals = y-fitted

hist(residuals, breaks=50, freq=F)
plot(fitted,residuals)

mod = lm(revenue_per_sqft~green_rating+City_Market_Rent+age+class_a+class_b, data=bldg)
summary(mod)

gibbs_sampler <- function(steps,burnin,init_params,X,y){
  # 1:n lambda
  # (n+1):(n+p) beta
  # n+p+1 omega
  h = 1
  d = 1
  eta = d
  n = nrow(X)
  p = ncol(X)
  beta = init_params[(n+1):(n+p)]
  omega = init_params[n+p+1]
  chain = matrix(0, nrow=steps-burnin, ncol=n+p+1)
  m = rep(0, p)
  K = rep(0.01, p)
  for(i in 1:steps){
    sqerr = drop((y-X%*%matrix(beta))^2)
    lambda = rgamma(n,shape=(h+1)/2,rate=1)/(h+omega*sqerr)*2
    omega = rgamma(
      1,
      shape=d+n+p,
      rate=(eta+sum((beta-m)^2*K)+sqerr*lambda)/2
    )
    K_star = diag(K)+t(X)%*%(X*lambda)
    K_star_inv = solve(K_star)
    m_star = K_star_inv%*%(diag(K)%*%matrix(m)+t(X*lambda)%*%matrix(y)) # no K_star_inv in front
    beta = rmvnorm(1, mean=m_star, sigma=K_star_inv/omega, method="chol")
    if(i > burnin){
      chain[i-burnin,] = c(lambda,beta,omega)
    }
  }
  return(chain)
}


init_lambda = rep(NaN, n)
init_beta = rep(0, p)
init_omega = rep(1, 1)
init_params = c(init_lambda,init_beta,init_omega)
chain = gibbs_sampler(
  1000,
  100,
  init_params,
  X,
  y
)

myconfint <- function(x,ci=0.95){
  func <- function(z){
    return(c(
      quantile(z,(1-ci)/2),
      quantile(z,(1+ci)/2)
    ))
  }
  return(apply(x,2,func))
}

myconfint(chain[,(n+1):(n+p)], 0.95)

beta_posterior_mean_ = colMeans(chain[,(n+1):(n+p)])
fitted_ = X%*%beta_posterior_mean_
residuals_ = y-fitted_
hist(residuals_, breaks=50, freq=F)
plot(fitted_,residuals_)


par(mfrow=c(2,1),mai=c(0.5,0.5,0.5,0.5))
hist(residuals, breaks=50, freq=F)
hist(residuals_, breaks=50, freq=F)

par(mfrow=c(2,1),mai=c(0.5,0.5,0.5,0.5))
plot(fitted,residuals)
plot(fitted_,residuals_)