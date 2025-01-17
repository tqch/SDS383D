# helper function

sigmoid <- function(x){1/(1+exp(-x))}

get_model <- function(data, family){
  y = data$y
  X = data$X
  n_params = ncol(X)
  if(family=="binomial"){
    func <- function(weights=NULL){
      if(is.null(weights)){
        return(n_params)
      } else {
        p = c(sigmoid(rowSums(X*t(weights))))
        loglk = mean(y*log(p)+(1-y)*log(1-p))
        grad = colMeans(X*(y-p))
        return(list(loglk=loglk, grad=grad))
      }
    }
  }
}

gradient_ascent <- function(
  model,
  init_weights=NULL,
  max_iter=1000,
  step_size=0.001,
  tol=1e-6,
  adam=T,
  beta=c(0.9,0.999),
  eps=1e-8
){
  if(is.null(init_weights)){
    weights = rep(0, model())
  } else {
    weights = init_weights
  }
  loglks = rep(0,max_iter)
  stop_iter = max_iter
  convg = F
  for(i in 1:max_iter){
    eval = model(weights)
    loglks[i] = eval$loglk
    if(adam){
      if(i==1){
        m = (1-beta[1])*eval$grad
        v = (1-beta[2])*eval$grad^2
      } else {
        m = beta[1]*m + (1-beta[1])*eval$grad
        v = beta[2]*v + (1-beta[2])*eval$grad^2
      }
      m_bc = m/(1-beta[1]^i)
      v_bc = v/(1-beta[2]^i)
      weights = weights + step_size*m_bc/(sqrt(v_bc)+eps)
    } else {
      weights = weights + step_size*eval$grad  
    }
    if(i>1){
      print(loglks[i])
    }
    if((i>1) && (abs((loglks[i]-loglks[i-1])/loglks[i-1])<tol)){
      convg = T
      stop_iter = i
      break
    }
  }
  return(list(
    convergence=convg, 
    log_likelihood=loglks[1:stop_iter], 
    weights=weights
    ))
}

whitening <- function(x){
  x.means = colMeans(x)
  x.sds = apply(x, 2, sd)
  out = sweep(x, 2, x.means, "-")
  out = sweep(out, 2, x.sds, "/")
  return(out)
}

wdbc = read.csv("./WDBC.csv")
outcome = ifelse(wdbc[,2]=="B",0,1)
features = cbind(rep(1,nrow(wdbc)),whitening(wdbc[,3:12]))
data = list(y=outcome,X=features)
model = get_model(data, family="binomial")
trained = gradient_ascent(model)
fitted = sigmoid(rowSums(features*t(trained$weights)))