mathtest = read.csv("./mathtest.csv")
summary_level = function(df, level_ind, val_ind, func){
  c(by(df[,val_ind], df[,level_ind], func))
}

myvar = function(x){
  var(x)*(length(x)-1)/length(x)
}

school_size = summary_level(mathtest, 1, 2, length)
avg_score = summary_level(mathtest, 1, 2, mean)
var_score = summary_level(mathtest, 1, 2, myvar)


plot(
  school_size,
  avg_score,
  type="p",
  xlab="School size",
  ylab="Average test score"
)


sumsq_score = summary_level(mathtest, 1, 2, function(x){sum(x^2)})
gibbs_sampler = function(
  steps,
  burnin,
  init_params,
  n_i=school_size,
  p=100,
  y_bar=avg_score,
  y_sumsq=sumsq_score
){
  # params[1]: s
  # params[2]: kappa
  # params[3]: mu
  # params[4:3+p]: theta
  
  n = sum(n_i)
  s = init_params$s
  kappa = init_params$kappa
  mu = init_params$mu
  
  chain = matrix(0, nrow=steps-burnin, ncol=3+p)
  
  for(i in 1:steps){
    m = kappa*mu+n_i*y_bar
    theta = rnorm(p)*sqrt(s/(kappa+n_i))+m/(kappa+n_i)
    s_theta = sum((theta-mu)^2)
    s_y = sum(y_sumsq-2*n_i*y_bar*theta+n_i*theta^2)
    v_s = s_y+kappa*s_theta
    s = v_s/rgamma(1,shape=(n+p)/2)/2
    v_kappa = 1+s_theta
    kappa = rgamma(1,shape=(p+1)/2)/v_kappa*2
    theta_bar = mean(theta)
    mu = rnorm(1)*sqrt(s/kappa/p)+theta_bar
    if(steps>burnin){
      chain[i-burnin,] = c(s,kappa,mu,theta)
    }
  }
  return(chain)
}

chain = gibbs_sampler(1000, 100, init_params=list(
  s=var(mathtest$mathscore),
  kappa=1,
  mu=mean(mathtest$mathscore)
))

avg_kappa = colMeans(outer(chain[,2],school_size,function(x,y){1/(1+y/x)}))
plot(
  school_size,
  avg_kappa,
  type="p",
  xlab="School size",
  ylab="Average shrinkage coefficient"
)

