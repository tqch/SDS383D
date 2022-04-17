polls = read.csv("./polls.csv")
data = na.omit(polls)
data$org = NULL
data$year = NULL
data$survey = NULL
data$weight = NULL
data$state = factor(data$state)

data$edu = factor(data$edu)
data$age = factor(data$age)
X = model.matrix(bush~.-1-state, data)
k = length(unique(data$state))
p = ncol(X)
n = nrow(X)

state_idx = factor(data$state, labels=1:k)

states = levels(data$state)
mean_state = sapply(states, function(st){mean(data$bush[data$state==st])})
ni = sapply(1:k, function(st){sum(state_idx==st)})

mu = qnorm(mean_state)
mu[is.infinite(mu)] = 3
beta = matrix(0, nrow=k, ncol=p)
z = mu[state_idx]+rowSums(X*beta[state_idx,])
u = rep(NA, n)
u[y0] = -abs(rnorm(length(y0)))
u[y1] = abs(rnorm(length(y1)))

init_params = list(
  mu=mu,
  beta=beta,
  sigma2=rep(1, p),
  u=u
)

y0 = which(data$bush==0)
y1 = which(data$bush==1)

gibbs <- function(n_iter, init_params, data){
  list2env(init_params, envir=environment())
  list2env(data, envir=environment())
  hist = list()
  
  pb = txtProgressBar(min=0, max=n_iter, style=3, width=50, char="=")
  
  for(i in 1:n_iter){
    for(j in 1:k){
      sb = diag(1/sigma2) + t(X[state_idx==j,])%*%X[state_idx==j,]
      sbinv = solve(sb)
      b = colSums((X*(u-mu[state_idx]))[state_idx==j,])
      beta[j,] = t(chol(sbinv))%*%rnorm(p) + solve(sb,b)
    }
    xb = rowSums(X*beta[state_idx,])
    bmu = u-xb
    bmui = sapply(1:k, function(st){mean(bmu[state_idx==st])})
    mu = rnorm(k)*sqrt(1/ni) + bmui
    z = mu[state_idx] + xb
    inds = y0
    while(TRUE){
      u[inds] = rnorm(length(inds))+z[inds]
      if(all(u[inds]<=0)){
        break
      } else {
        inds = inds[which(u[inds]>0)]
      }
    }
    inds = y1
    while(TRUE){
      u[inds] = rnorm(length(inds))+z[inds]
      if(all(u[inds]>0)){
        break
      } else {
        inds = inds[which(u[inds]<=0)]
      }
    }
    # sbl = colSums(do.call(rbind, beta)^2)
    sbl = colSums(beta^2)
    sigma2 = 1/rgamma(p, shape=k/2)*sbl/2
    hist[[i]] = list(mu=mu, beta=beta, sigma2=sigma2, u=u)
    setTxtProgressBar(pb, i)
  }
  return(hist)
}

n_iter = 1000
hist = gibbs(n_iter, init_params, list(X=X, y=data$bush))

burnin = 200
beta_mean = 0
for(i in burnin:n_iter){
  beta_mean = beta_mean + hist[[i]]$beta
}
beta_mean = beta_mean/(n_iter-burnin)

beta_var = 0
for(i in burnin:n_iter){
  beta_var = beta_var + (hist[[i]]$beta-beta_mean)^2
}
beta_var = beta_var/(n_iter-burnin)

beta_std = sqrt(beta_var)

beta_effect = matrix("?", nrow=k, ncol=p)
beta_effect[beta_mean>1.96*beta_std] = "+"
beta_effect[beta_mean<-1.96*beta_std] = "-"

mu_mean = colMeans(do.call(rbind,lapply(hist[burnin:n_iter],function(x){x$mu})))
pred = (mu_mean[state_idx]+rowSums(X*beta_mean[state_idx,]))>0
acc0 = mean(pred==data$bush)
mod = glm(bush~.,data=data,family="binomial")
pred = predict(mod, type="response")>0.5
acc1 = mean(pred==data$bush)