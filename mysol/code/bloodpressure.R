library(latex2exp)

bloodpressure = read.csv("./bloodpressure.csv")

summary_level = function(df, level_ind, val_ind, func){
  c(by(df[,val_ind], df[,level_ind], func))
}

mycond = function(y,rel){
  func = function(x){
    do.call(rel,list(x,y))
  }
  return(func)
}

myfilter = function(df,val_col,cond_col,cond){
  df[val_col][cond(df[cond_col])]
}

treatment1 = myfilter(bloodpressure, "systolic", "treatment", mycond(1,"=="))
treatment2 = myfilter(bloodpressure, "systolic", "treatment", mycond(2,"=="))

t.test(treatment1,treatment2)

mean_by_subj = summary_level(bloodpressure, 3, 2, mean)

t.test(mean_by_subj[1:10],mean_by_subj[11:20])


sq_by_subj = summary_level(bloodpressure, 3, 2, function(x){sum(x^2)})
num_records = summary_level(bloodpressure, 3, 2, length)
group_indicator = summary_level(bloodpressure, 3, 4, mean) - 1

gibbs_sampler = function(
  steps,
  burnin,
  init_params,
  n_i=num_records,
  p=20,
  y_bar=mean_by_subj,
  y_sumsq=sq_by_subj,
  x=group_indicator
){
  # params[1]: s
  # params[2]: kappa
  # params[3]: mu
  # params[4]: beta
  # params[5:4+p]: theta
  
  n = sum(n_i)
  s = init_params$s
  kappa = init_params$kappa
  mu = init_params$mu
  beta = init_params$beta
  
  chain = matrix(0, nrow=steps-burnin, ncol=4+p)
  
  for(i in 1:steps){
    m = kappa*(mu+beta*x)+n_i*y_bar
    theta = rnorm(p)*sqrt(s/(kappa+n_i))+m/(kappa+n_i)
    s_theta = sum((theta-mu-beta*x)^2)
    s_y = sum(y_sumsq-2*n_i*y_bar*theta+n_i*theta^2)
    v_s = s_y+kappa*s_theta
    s = v_s/rgamma(1,shape=(n+p)/2)/2
    v_kappa = 1+s_theta/s
    kappa = rgamma(1,shape=(p+1)/2)/v_kappa*2
    theta_bar = mean(theta)
    mu = rnorm(1)*sqrt(s/kappa/p)+theta_bar-beta*mean(x)
    beta = rnorm(1)*sqrt(s/kappa/sum(x^2))+mean(x*(theta-mu))
    if(steps>burnin){
      chain[i-burnin,] = c(s,kappa,mu, beta, theta)
    }
  }
  return(chain)
}

chain = gibbs_sampler(1000, 100, init_params=list(
  s=var(bloodpressure$systolic),
  kappa=1,
  mu=mean(bloodpressure$systolic),
  beta=0
))

hist(
  chain[,4],
  xlim=c(-10,10),
  main=TeX("Posterior distribution of $\\beta$"),xlab=TeX("$\\beta$")
)

# 95% credible interval
c(quantile(chain[,4],0.025),quantile(chain[,4],0.975))