# gd

sigmoid <- function(x){1/(1+exp(-x))}
log_sigmoid <- function(x){-log(1+exp(-x))}

get_model <- function(data, family){
  y = data$y
  X = data$X
  n_params = ncol(X)
  if(family=="binomial"){
    func <- function(weights=NULL){
      if(is.null(weights)){
        return(n_params)
      } else {
        eta = rowSums(sweep(X, 2, weights, "*"))
        p = sigmoid(eta)
        loglk = mean(y*log_sigmoid(eta)+(1-y)*log_sigmoid(-eta))
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
  step_size=0.01,
  tol=1e-6,
  adam=F,
  beta=c(0.9,0.999),
  eps=1e-8,
  pho=0.9,
  c=0.1
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
    if(adam){
      eval = model(weights)
      loglks[i] = eval$loglk
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
      # weights = weights + step_size*eval$grad  
      alpha = step_size
      if(i==1){
        eval = model(weights)
      }
      loglks[i] = eval$loglk
      while(T){
        neweval = model(weights + alpha*eval$grad)
        if(is.nan(neweval$loglk)){
          alpha = alpha*pho
          next
        }
        if(neweval$loglk>=(eval$loglk+c*alpha*sum(eval$grad^2))){
          eval = neweval
          weights = weights + alpha*eval$grad
          break
        } else {
          alpha = alpha*pho
        }
      }
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
trained = gradient_ascent(model, adam=T, max_iter=10000)
plot(trained$log_likelihood, type="l", lwd=3)
fitted = sigmoid(rowSums(features*t(trained$weights)))
cat(sprintf("Accuracy:%f", mean((fitted>0.5)==outcome)))


# newton's method

sigmoid <- function(x){1/(1+exp(-x))}
log_sigmoid <- function(x){-log(1+exp(-x))}

get_model <- function(data, family){
  y = data$y
  X = data$X
  n_samples = nrow(X)
  n_params = ncol(X)
  if(family=="binomial"){
    func <- function(weights=NULL){
      if(is.null(weights)){
        return(n_params)
      } else {
        eta = as.matrix(X)%*%matrix(weights)
        p = drop(sigmoid(eta))
        loglk = mean(y*log_sigmoid(eta)+(1-y)*log_sigmoid(-eta))
        grad = colMeans(X*(y-p))
        hessian = -t(X)%*%as.matrix(X*p*(1-p))/n_samples
        return(list(loglk=loglk, grad=grad, hessian=hessian))
      }
    }
  }
}

newton <- function(
  model,
  init_weights=NULL,
  max_iter=1000,
  step_size=1,
  tol=1e-6,
  pho=0.9,
  c=0.1
){
  if(is.null(init_weights)){
    weights = rep(0.1, model())
  } else {
    weights = init_weights
  }
  loglks = rep(0, max_iter)
  stop_iter = max_iter
  convg = F
  for(i in 1:max_iter){
    alpha = step_size
    if(i==1){
      eval = model(weights)
    }
    # eval = model(weights)
    loglks[i] = eval$loglk
    g = -solve(eval$hessian, eval$grad)
    weights = weights + step_size*g
    while(T){
      neweval = model(weights + alpha*g)
      if(is.nan(neweval$loglk)){
        alpha = alpha*pho
        next
      }
      if(neweval$loglk>=(eval$loglk+c*alpha*sum(eval$grad*g))){
        eval = neweval
        weights = weights + alpha*g
        break
      } else {
        alpha = alpha*pho
      }
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
trained = newton(model)
fitted = sigmoid(rowSums(features*t(trained$weights)))
cat(sprintf("Accuracy:%f", mean((fitted>0.5)==outcome)))