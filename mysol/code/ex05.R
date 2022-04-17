cheese = read.csv("./cheese.csv")
n_stores = length(unique(cheese$store))
data = data.frame(
  x = log(cheese$price),
  y = log(cheese$vol),
  d = cheese$disp,
  idx = factor(cheese$store, labels=c(1:n_stores))
)

n = nrow(data)
ni = sapply(1:n_stores, function(x){sum(data$idx==x)})


init_params = list(
  mu = mean(data$y),
  s = var(data$y),
  a = sapply(1:n_stores,
             function(x){mean(data$y[(data$d==0)&(data$idx==x)])}),
  # theta = mean(data$y[data$d==1])-mean(data$y[data$d==0]),
  theta = sapply(
    1:n_stores,
    function(x){mean(
      data$y[(data$d==1)&(data$idx==x)])-mean(
      data$y[(data$d==0)&(data$idx==x)])}
  ),
  m = mean(data$y/data$x),
  v = var(data$y/data$x),
  b = sapply(1:n_stores,
             function(x){mean((data$y/data$x)[data$d==0&data$idx==x])}),
  # lambda = mean((data$y/data$x)[data$d==1])-mean((data$y/data$x)[data$d==0]),
  lambda = sapply(
    1:n_stores,
    function(x){mean(
      data$y/data$x[(data$d==1)&(data$idx==x)])-mean(
        data$y/data$x[(data$d==0)&(data$idx==x)])}
  ),
  tau = 1
)

init_params$a[is.na(init_params$a)] = 0
init_params$b[is.na(init_params$b)] = 0
init_params$theta[is.na(init_params$theta)] = 0
init_params$lambda[is.na(init_params$lambda)] = 0

pos = cumsum(c(0,unlist(lapply(init_params, length))))
p = tail(pos, 1)
k = n_stores
# n1 = sum(data$d)
nz = sapply(1:k, function(z){sum(data$d[data$idx==z])})
sxi = sapply(1:k, function(z){sum(data$x[data$idx==z]^2)})
# sxst = sum(data$x[data$d==1]^2)
sxist = sapply(1:k, function(z){sum(data$x[(data$idx==z)&(data$d==1)]^2)})


m0 = 0
v0 = 100
a0 = 2
b0 = 1

# gibbs <- function(n_iter, init_params, data=data){
#   # unpack all the elements in the list
#   list2env(data, envir=environment())
#   list2env(init_params, envir=environment())
#   # parameter matrix
#   hist = matrix(NA, nrow=n_iter, ncol=p)
#   for(i in 1:n_iter){
#     mu = sqrt(s)*rnorm(1) + mean(a)
#     yhat = a[idx]+theta*d+(b[idx]+lambda*d)*x
#     sa = sum((a-mu)^2)
#     sy = sum((y-yhat)^2)
#     s = 1/rgamma(1, shape=(k+n+1)/2, rate=(sa+theta^2+sy/tau)/2)
#     ynoa = theta*d+(b[idx]+lambda*d)*x
#     ybara = sapply(1:k, function(z){mean((y-ynoa)[idx==z])})
#     temp1 = (mu+ni*ybara/tau)/s
#     temp2 = s/(1+ni/tau)
#     a = rnorm(k)*sqrt(temp2) + temp1*temp2
#     ynoth = a[idx]+(b[idx]+lambda*d)*x
#     ybarth = mean((y-ynoth)[d==1])
#     temp1 = n1*ybarth/tau
#     temp2 = 1/(1+n1/tau)
#     theta = rnorm(1)*sqrt(s*temp2) + temp1*temp2
#     temp1 = m0/v0+k*mean(b)/v
#     temp2 = 1/(1/v0+k/v)
#     m = rnorm(1)*sqrt(temp2) + temp1*temp2
#     sb = sum((b-m)^2)
#     v = 1/rgamma(1, shape=(k+n+1)/2, rate=(sb+lambda^2)/2)
#     ynob = a[idx]+theta*d+lambda*d*x
#     ybarb = sapply(1:k, function(z){sum((x*(y-ynob))[idx==z])})
#     temp1 = m/v+ybarb/tau/s
#     temp2 = 1/(1/v+sxi/tau/s)
#     b = rnorm(k)*sqrt(temp2) + temp1*temp2
#     ynolbd = a[idx]+theta*d+b[idx]*x
#     ybarlbd = sum(((y-ynolbd)*x)[d==1])
#     temp1 = ybarlbd/tau/s
#     temp2 = 1/(1/v+sxst/tau/s)
#     lambda = rnorm(1)*sqrt(temp2) + temp1*temp2
#     yhat = a[idx]+theta*d+(b[idx]+lambda*d)*x
#     sy = sum((y-yhat)^2)
#     tau = 1/rgamma(1, shape=a0+n/2, rate=b0+sy/2/s)
#     hist[i,] = c(mu, s, a, theta, m, v, b, lambda, tau)
#   }
#   return(hist)
# }

gibbs <- function(n_iter, init_params, data=data){
  # unpack all the elements in the list
  list2env(data, envir=environment())
  list2env(init_params, envir=environment())
  # parameter matrix
  hist = matrix(NA, nrow=n_iter, ncol=p)
  for(i in 1:n_iter){
    mu = rnorm(1)*sqrt(s) + mean(a)
    yhat = a[idx]+theta[idx]*d+(b[idx]+lambda[idx]*d)*x
    sa = sum((a-mu)^2)
    sy = sum((y-yhat)^2)
    s = 1/rgamma(1, shape=(2*k+n)/2, rate=(sa+sum(theta^2)+sy/tau)/2)
    ynoa = theta[idx]*d+(b[idx]+lambda[idx]*d)*x
    ybara = sapply(1:k, function(z){sum((y-ynoa)[idx==z])})
    temp1 = (mu+ybara/tau)/s
    temp2 = s/(1+ni/tau)
    a = rnorm(k)*sqrt(temp2) + temp1*temp2
    ynoth = a[idx]+(b[idx]+lambda[idx]*d)*x
    ybarth = sapply(1:k, function(z){sum((y-ynoth)[(idx==z)&(d==1)])})
    temp1 = ybarth/tau
    temp2 = 1/(1+nz/tau)
    theta = rnorm(k)*sqrt(s*temp2) + temp1*temp2
    temp1 = m0/v0+k*mean(b)/v
    temp2 = 1/(1/v0+k/v)
    m = rnorm(1)*sqrt(temp2) + temp1*temp2
    sb = sum((b-m)^2)
    v = 1/rgamma(1, shape=(2*k)/2, rate=(sb+sum(lambda^2))/2)
    ynob = a[idx]+theta[idx]*d+lambda[idx]*d*x
    ybarb = sapply(1:k, function(z){sum((x*(y-ynob))[idx==z])})
    temp1 = m/v+ybarb/tau/s
    temp2 = 1/(1/v+sxi/tau/s)
    b = rnorm(k)*sqrt(temp2) + temp1*temp2
    ynolbd = a[idx]+theta[idx]*d+b[idx]*x
    ybarlbd = sapply(1:k, function(z){sum(((y-ynolbd)*x)[(idx==z)&(d==1)])})
    temp1 = ybarlbd/tau/s
    temp2 = 1/(1/v+sxist/tau/s)
    lambda = rnorm(k)*sqrt(temp2) + temp1*temp2
    yhat = a[idx]+theta[idx]*d+(b[idx]+lambda[idx]*d)*x
    sy = sum((y-yhat)^2)
    tau = 1/rgamma(1, shape=a0+n/2, rate=b0+sy/2/s)
    hist[i,] = c(mu, s, a, theta, m, v, b, lambda, tau)
  }
  return(hist)
}

n_iter = 1000
hist = gibbs(n_iter, init_params, data)

burnin = 200
posterior_mean = colMeans(hist[burnin:n_iter, ])
posterior_mean = list(
  a = posterior_mean[(pos["a"]-k+1):pos["a"]],
  b = posterior_mean[(pos["b"]-k+1):pos["b"]],
  theta = posterior_mean[(pos["theta"]-k+1):pos["theta"]],
  lambda = posterior_mean[(pos["lambda"]-k+1):pos["lambda"]]
)

yhat = posterior_mean$a[data$idx] + posterior_mean$theta[data$idx]*data$d +
  (posterior_mean$b[data$idx]+posterior_mean$lambda[data$idx]*data$d)*data$x
mse = mean((data$y-yhat)^2)

library(latex2exp)
mycolors = rainbow(k)

png("./figs/logvol_vs_logprice.png", width=1200, height=800)
plot(
  data$x, data$y, col=mycolors[data$idx],
  xlab=TeX("$\\log P_{ij}$"), ylab=TeX("$\\log Q_{ij}$")
)


png("./figs/cheese_fitted.png", width=2200, height=1600)
par(mfrow=c(8,k/8))
colors = c("red", "green")

for(i in 1:k){
  datai = data[data$idx==i,]
  plot(datai$x, datai$y, col=colors[datai$d+1])
  abline(a=posterior_mean$a[i],b=posterior_mean$b[i],xlab=TeX("$\\log P_{ij}$"), ylab=TeX("$\\log Q_{ij}$"))
  abline(a=posterior_mean$a[i]+posterior_mean$theta[i],b=posterior_mean$b[i]+posterior_mean$lambda[i])
}
dev.off()