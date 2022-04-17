get_kernel <- function(type, bandwidth){
  if(type=="gaussian"){
    kernel <- function(x, y){
      if(!is.null(dim(x))){
        if(nrow(y)==1){
          dif = abs(sweep(x, 1, drop(y), "-"))/bandwidth
        }
      } else {
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

pairwise_inner_product <- function(X, kernel){
  n = nrow(X)
  grid = expand.grid(1:n,1:n)
  inds = grid$Var1<=grid$Var2
  grid = grid[inds,]
  func <- function(inds){
    ind_1 = inds[1]
    ind_2 = inds[2]
    return(kernel(X[ind_1,],X[ind_2,]))
  }
  upper = apply(grid, 1, func)
  m = matrix(0, nrow=n, ncol=n)
  m[inds] = upper
  m = m + t(m*upper.tri(m))
  return(m)
}

get_weights <- function(X, newx, kernel){
  weights = drop(kernel(X, newx))
  return(weights/sum(weights))
}

kernel_regression <- function(data, kernel_type="gaussian", bandwidth=1){
  X = data[,"X", drop=FALSE]
  y = data$y
  kernel = get_kernel(kernel_type, bandwidth)
  # feat = pairwise_inner_product(X, kernel)
  # params = solve(t(feat)%*%feat,feat%*%matrix(y))
  weights = pairwise_inner_product(X, kernel)
  weights = weights/rowSums(weights)
  fitted = drop(weights%*%y)
  predict <- function(newX){
    func <- function(newx){
      return(get_weights(X, matrix(newx, nrow=1), kernel))
    }
    new_weights = t(apply(newX, 1, func))
    return(drop(new_weights%*%matrix(y)))
  }
  return(list(fitted=fitted, predict=predict, kernel=kernel))
}

# noisy_data <- function(n, sigma=0.1){
#   thetas = runif(n, -pi, pi)
#   rads = runif(n, 0, 1)
#   X = rads * cbind(cos(thetas),sin(thetas))
#   y = rads + sigma*rnorm(n)
#   return(list(X=X,y=y))
# }

noisy_data <- function(n, sigma=0.1){
  X = matrix(seq(-1, 1, length.out=n))
  y = cos(3*pi*X) + sigma*rnorm(n)
  return(data.frame(X=X,y=y))
}

n = 100
mydata = noisy_data(n, sigma=0.1)
model = kernel_regression(mydata, kernel_type="gaussian", bandwith=1)

# plot(X,cex=y/max(y),type="p")
plot(X, y, type="p")

train_mse = mean((model$fitted-y)^2)

evaluate <- function(model, new_data){
  predicted = model$predict(new_data[,"X",drop=FALSE])
  return(mean((predicted-new_data$y)^2))
}

loocv <- function(model, data){
  X = data[,"X", drop=FALSE]
  y = data$y
  weights = pairwise_inner_product(X, model$kernel)
  weights = weights/rowSums(weights)
  cv = sum(((y-model$fitted)/(1-diag(weights)))^2)
  return(cv)
}

new_data = noisy_data(n)
test_mse = evaluate(model, new_data)

train_test_split <- function(data, train_size=0.75){
  n = nrow(data[,"X",drop=FALSE])
  permuted = sample(1:n, n, replace=FALSE)
  train_size = round(train_size*n)
  train_data = data[permuted[1:train_size],]
  test_data = data[permuted[(train_size+1):n],]
  return(list(train_data=train_data,test_data=test_data))
}

set.seed(1234)
sigma=0.1
data = noisy_data(n, sigma=sigma)
splits = train_test_split(data)
train_data = splits$train_data
test_data = splits$test_data
models = list()
mses = c()

bandwidths = c(0.01, 0.03, 0.1, 0.3, 1)
k = length(bandwidths)
for(i in 1:k){
  bwd = bandwidths[i]
  model = kernel_regression(train_data, bandwidth=bwd)
  models[[i]] = model
  mse = evaluate(model, test_data)
  mses = c(mses, mse)
}

cvs = c()
for(i in 1:k){
  bwd = bandwidths[i]
  model = kernel_regression(data, bandwidth=bwd)
  cv = loocv(model, data)
  cvs = c(cvs, cv)
}

grid = seq(-1,1,length.out=100)
ground_truth = cos(3*pi*grid)

mycolors = rainbow(1+k)

png(filename="./figs/wiggly_less_noisy.png", width=1200, height=800)

plot(grid, ground_truth, type="l",lty=1,
     xlim=c(-1,1), ylim=c(-2, 4),
     xlab="x", ylab="y",
     col=mycolors[1],
     lwd=3
)
par(new=TRUE)
plot(train_data$X, train_data$y, type="p", cex=3,
     xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,1), ylim=c(-2, 4))

for(i in 1:k){
  par(new=TRUE)
  plot(grid, models[[i]]$predict(matrix(grid)),
       type="l", lty=i+1,
       xlim=c(-1, 1), ylim=c(-2, 4),
       xaxt="n", yaxt="n", xlab="", ylab="",
       col=mycolors[i+1], lwd=5)
}

best_bwd = bandwidths[which.min(mses)]
legend("topright", legend=c("training data",
"ground truth", sapply(bandwidths, function(x){
  paste0("bandwidth=", x, ifelse(x==best_bwd," (best)",""), collapse="")
  })), col=c("black", mycolors), pch=c(1,rep(NA,k+1)),
  lty=c(NA,1:(1+k)), 
  lwd=c(NA,rep(5,1+k)), cex=2
)

dev.off()


set.seed(1234)
sigma=0.5
data = noisy_data(n, sigma=sigma)
splits = train_test_split(data)
train_data = splits$train_data
test_data = splits$test_data
models = list()
mses = c()

for(i in 1:length(bandwidths)){
  bwd = bandwidths[i]
  model = kernel_regression(train_data, bandwidth=bwd)
  models[[i]] = model
  mse = evaluate(model, test_data)
  mses = c(mses, mse)
}

cvs = c()
for(i in 1:k){
  bwd = bandwidths[i]
  model = kernel_regression(data, bandwidth=bwd)
  cv = loocv(model, data)
  cvs = c(cvs, cv)
}

grid = seq(-1,1,length.out=100)
ground_truth = cos(3*pi*grid)

mycolors = rainbow(1+k)

png(filename="./figs/wiggly_more_noisy.png", width=1200, height=800)

plot(grid, ground_truth, type="l",lty=1,
     xlim=c(-1,1), ylim=c(-2, 4),
     xlab="x", ylab="y",
     col=mycolors[1],
     lwd=3
)
par(new=TRUE)
plot(train_data$X, train_data$y, type="p", cex=3,
     xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,1), ylim=c(-2, 4))

for(i in 1:k){
  par(new=TRUE)
  plot(grid, models[[i]]$predict(matrix(grid)),
       type="l", lty=i+1,
       xlim=c(-1, 1), ylim=c(-2, 4),
       xaxt="n", yaxt="n", xlab="", ylab="",
       col=mycolors[i+1], lwd=5)
}

best_bwd = bandwidths[which.min(mses)]
legend("topright", legend=c("training data",
        "ground truth", sapply(bandwidths, function(x){
          paste0("bandwidth=", x, ifelse(x==best_bwd," (best)",""), collapse="")
        })), col=c("black", mycolors), pch=c(1,rep(NA,k+1)),
       lty=c(NA,1:(1+k)), 
       lwd=c(NA,rep(5,1+k)), cex=2
)

dev.off()


noisy_data <- function(n, sigma=0.1){
  X = matrix(seq(-1, 1, length.out=n))
  y = 4*X*(X^2-3/4) + sigma*rnorm(n)
  return(data.frame(X=X,y=y))
}

set.seed(1234)
sigma=0.1
data = noisy_data(n, sigma=sigma)
splits = train_test_split(data)
train_data = splits$train_data
test_data = splits$test_data
models = list()
mses = c()

for(i in 1:length(bandwidths)){
  bwd = bandwidths[i]
  model = kernel_regression(train_data, bandwidth=bwd)
  models[[i]] = model
  mse = evaluate(model, test_data)
  mses = c(mses, mse)
}

cvs = c()
for(i in 1:k){
  bwd = bandwidths[i]
  model = kernel_regression(data, bandwidth=bwd)
  cv = loocv(model, data)
  cvs = c(cvs, cv)
}

grid = seq(-1,1,length.out=100)
ground_truth = 4*grid*(grid^2-3/4)

mycolors = rainbow(1+k)

png(filename="./figs/smooth_less_noisy.png", width=1200, height=800)

plot(grid, ground_truth, type="l",lty=1,
     xlim=c(-1,1), ylim=c(-2, 4),
     xlab="x", ylab="y",
     col=mycolors[1],
     lwd=3
)
par(new=TRUE)
plot(train_data$X, train_data$y, type="p", cex=3,
     xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,1), ylim=c(-2, 4))

for(i in 1:k){
  par(new=TRUE)
  plot(grid, models[[i]]$predict(matrix(grid)),
       type="l", lty=i+1,
       xlim=c(-1, 1), ylim=c(-2, 4),
       xaxt="n", yaxt="n", xlab="", ylab="",
       col=mycolors[i+1], lwd=5)
}

best_bwd = bandwidths[which.min(mses)]
legend("topright", legend=c("training data",
      "ground truth", sapply(bandwidths, function(x){
        paste0("bandwidth=", x, ifelse(x==best_bwd," (best)",""), collapse="")
      })), col=c("black", mycolors), pch=c(1,rep(NA,k+1)),
     lty=c(NA,1:(1+k)), 
     lwd=c(NA,rep(5,1+k)), cex=2
)

dev.off()

set.seed(1234)
sigma=0.5
data = noisy_data(n, sigma=sigma)
splits = train_test_split(data)
train_data = splits$train_data
test_data = splits$test_data
models = list()
mses = c()

for(i in 1:length(bandwidths)){
  bwd = bandwidths[i]
  model = kernel_regression(train_data, bandwidth=bwd)
  models[[i]] = model
  mse = evaluate(model, test_data)
  mses = c(mses, mse)
}

cvs = c()
for(i in 1:k){
  bwd = bandwidths[i]
  model = kernel_regression(data, bandwidth=bwd)
  cv = loocv(model, data)
  cvs = c(cvs, cv)
}

grid = seq(-1,1,length.out=100)
ground_truth = 4*grid*(grid^2-3/4)

mycolors = rainbow(1+k)

png(filename="./figs/smooth_more_noisy.png", width=1200, height=800)

plot(grid, ground_truth, type="l",lty=1,
     xlim=c(-1,1), ylim=c(-2, 4),
     xlab="x", ylab="y",
     col=mycolors[1],
     lwd=3
)
par(new=TRUE)
plot(train_data$X, train_data$y, type="p", cex=3,
     xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,1), ylim=c(-2, 4))

for(i in 1:k){
  par(new=TRUE)
  plot(grid, models[[i]]$predict(matrix(grid)),
       type="l", lty=i+1,
       xlim=c(-1, 1), ylim=c(-2, 4),
       xaxt="n", yaxt="n", xlab="", ylab="",
       col=mycolors[i+1], lwd=5)
}

best_bwd = bandwidths[which.min(mses)]
legend("topright", legend=c("training data",
        "ground truth", sapply(bandwidths, function(x){
          paste0("bandwidth=", x, ifelse(x==best_bwd," (best)",""), collapse="")
        })), col=c("black", mycolors), pch=c(1,rep(NA,k+1)),
       lty=c(NA,1:(1+k)), 
       lwd=c(NA,rep(5,1+k)), cex=2
)

dev.off()
