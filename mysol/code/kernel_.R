library("scales")

get_kernel <- function(type, bandwidth){
  if(type=="gaussian"){
    kernel <- function(x, y){
      if(!is.null(dim(x))){
        if(nrow(y)==1){
          dif = abs(sweep(x, 1, drop(y), "-"))/bandwidth
        }} else {
          dif = abs(x-y)/bandwidth
        }
      if(is.null(dim(dif))){
        return(exp(-sum(dif^2)/2)/sqrt(2*pi)/bandwidth)
      } else {
        return(exp(-rowSums(dif^2)/2)/sqrt(2*pi)/bandwidth)
      }
    }
    return(kernel)
  }
}


local_linear_estimator <- function(data, kernel){
  X = data$X
  y = data$y
  n = nrow(X)
  A = outer(X, X, function(x,y){y-x})
  predict <- function(x){
    weights = sapply(X, function(y){kernel(y, x)})
    weights = (A%*%((X-x)*weights))*weights
    fitted = t(y)%*%weights
    fitted = fitted/sum(weights)
    return(fitted)
  }
  return(list(predict=predict))
}

pairwise_inner_product <- function(X, kernel){
  n = length(X)
  grid = expand.grid(1:n,1:n)
  inds = grid$Var1<=grid$Var2
  grid = grid[inds,]
  func <- function(inds){
    ind_1 = inds[1]
    ind_2 = inds[2]
    return(kernel(X[ind_1],X[ind_2]))
  }
  upper = apply(grid, 1, func)
  m = matrix(0, nrow=n, ncol=n)
  m[inds] = upper
  m = m + t(m*upper.tri(m))
  return(m)
}

loocv <- function(data, kernel_type, bandwidths){
  X = data$X
  y = data$y
  cvs = c()
  for(bwd in bandwidths){
    kernel = get_kernel(kernel_type, bwd)
    weights = pairwise_inner_product(X, kernel)
    A = outer(X, X, function(x,y){y-x})
    cv = 0
    for(i in length(X)){
      # weights_i = outer(weights[i,], weights[i,], "*")
      # weights_i = weights_i * outer(X, X, function(x,y){y-x})
      # weights_i = rowSums(sweep(weights_i, 2, X-X[i], "*")
      # fitted_i = sum(weights[-i,-i]%*%y[-i])/sum(weights[-i,-i])
      # cv = cv + (y[i]-fitted_i)^2
      weights_i = (A[-i,-i] %*% (weights[i,-i]*(X[-i]-X[i])))*weights[i,-i]
      fitted_i = t(y[-i])%*%weights_i
      fitted_i = fitted_i/sum(weights_i)
      cv = cv + (y[i]-fitted_i)^2
    }
    cvs = c(cvs, cv)
  }
  return(cvs)
}

bandwidths = c(1, 3, 10, 15, 30)

utilities = read.csv("./utilities.csv")
data = data.frame(X=utilities$temp, y=utilities$gasbill/utilities$billingdays)

cvs = loocv(data, "gaussian", bandwidths)

best_kernel = get_kernel("gaussian", bandwidths[which.min(cvs)])

model = local_linear_estimator(data, best_kernel)
fitted = sapply(data$X, model$predict)
residuals = data$y-fitted

# png("./figs/local_linear_regression_residuals_vs_fitted.png", width=1200, height=800)

plot(fitted, residuals,
     ylab="residuals", xlab="fitted values", cex = 3)

# dev.off()

# png("./figs/local_linear_regression_residuals_vs_x.png", width=1200, height=800)

plot(data$X, residuals, ylab="residuals", xlab="x", cex=3)

# dev.off()

# png("./figs/local_linear_regression_95_confint.png", width=1200, height=800)

xlim = c(min(data$X)-5, max(data$X)+5)
ylim = c(0, max(data$y)+1)
plot(data$X, data$y, type="p",xlim=xlim,ylim=ylim,xlab="x",ylab="y", pch=1)

par(new=TRUE)
sorted = sort(data$X, index.return=TRUE)
plot(sorted$x, fitted[sorted$ix],
     type="l",xlim=xlim,ylim=ylim,xlab="x",ylab="y",xaxt="n",yaxt="n",lwd=3)

# refit a model with residual squares
data_ = data.frame(X=data$X,y=residuals^2)
model_ = local_linear_estimator(data_, best_kernel)
fitted_ = sqrt(sapply(data$X, model_$predict))

par(new=TRUE)
plot(sorted$x, (fitted+1.96*fitted_)[sorted$ix],
     type="l",xlim=xlim,ylim=ylim,xlab="x",ylab="y",xaxt="n",yaxt="n",lwd=2, lty=2)

par(new=TRUE)
plot(sorted$x, (fitted-1.96*fitted_)[sorted$ix],
     type="l",xlim=xlim,ylim=ylim,xlab="x",ylab="y",xaxt="n",yaxt="n",lwd=2, lty=2)


par(new=TRUE)
polygon(c(sorted$x,rev(sorted$x)), c(
  (fitted+1.96*fitted_)[sorted$ix], rev((fitted-1.96*fitted_)[sorted$ix])),
  col=alpha("blue", 0.3))

legend("topright", legend=c("97.5% quantile","2.5% quantile", "fitted","data"),
       lwd=c(2,2,3,NA), pch=c(rep(NA,3),1), cex=2, lty=c(2,2,1,NA))

# dev.off()

get_hmatrix <- function(data, kernel){
  X = data$X
  n = length(X)
  ws = pairwise_inner_product(X, kernel)
  A = outer(X, X, function(x, y){y-x})
  H = matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
    wi = ws[i,]
    H[i,] = A%*%(wi*(X-X[i]))*wi
  }
  H = H/rowSums(H)
  return(H)
}
n = nrow(data)
H = get_hmatrix(data, best_kernel)
sigma2 = sum(residuals^2)/(n-2*sum(diag(H))+sum(diag(t(H)%*%H)))
sigma_x = sqrt(rowSums(H^2)*sigma2)

png("./figs/local_linear_regression_95_confint.png", width=1200, height=800)

xlim = c(min(data$X)-5, max(data$X)+5)
ylim = c(0, max(data$y)+1)
plot(data$X, data$y, type="p",xlim=xlim,ylim=ylim,xlab="x",ylab="y", pch=1)

par(new=TRUE)
sorted = sort(data$X, index.return=TRUE)
plot(sorted$x, fitted[sorted$ix],
     type="l",xlim=xlim,ylim=ylim,xlab="x",ylab="y",xaxt="n",yaxt="n",lwd=3)

par(new=TRUE)
plot(sorted$x, (fitted+1.96*sigma_x)[sorted$ix],
     type="l",xlim=xlim,ylim=ylim,xlab="x",ylab="y",xaxt="n",yaxt="n",lwd=2, lty=2)

par(new=TRUE)
plot(sorted$x, (fitted-1.96*sigma_x)[sorted$ix],
     type="l",xlim=xlim,ylim=ylim,xlab="x",ylab="y",xaxt="n",yaxt="n",lwd=2, lty=2)


par(new=TRUE)
polygon(c(sorted$x,rev(sorted$x)), c(
  (fitted+1.96*sigma_x)[sorted$ix], rev((fitted-1.96*sigma_x)[sorted$ix])),
  col=alpha("blue", 0.3))

legend("topright", legend=c("97.5% quantile","2.5% quantile", "fitted","data"),
       lwd=c(2,2,3,NA), pch=c(rep(NA,3),1), cex=2, lty=c(2,2,1,NA))

dev.off()