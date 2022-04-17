# helper function

sigmoid <- function(x){1/(1+exp(-x))}

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
        p = c(sigmoid(rowSums(X*t(weights))))
        loglk = mean(y*log(p)+(1-y)*log(1-p))
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
    weights = rep(0, model())
  } else {
    weights = init_weights
  }
  loglks = rep(0,max_iter)
  stop_iter = max_iter
  convg = F
  for(i in 1:max_iter){
    alpha = step_size
    if(i==1){
      eval = model(weights)
    }
    loglks[i] = eval$loglk
    g = -solve(eval$hessian, eval$grad)
    print(sum(g*eval$grad))
    while(T){
      neweval = model(weights + alpha*g)
      if(neweval$loglk>=(eval$loglk+c*alpha*sum(eval$grad*g))){
        eval = neweval
        break
      } else {
        alpha = alpha*pho
      }
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
trained = newton(model)
fitted = sigmoid(rowSums(features*t(trained$weights)))
cat(sprintf("Accuracy:%f", mean((fitted>0.5)==outcome)))