GM <- function(X, y){
  n = dim(X)[1]; p = dim(X)[2]
  tau <- 1/sqrt(diag(solve(t(X) %*% X/n)))
  M <- rep(0, p)
  for(j in 1:p){
    Z <- rnorm(n)
    P <- X[, -j]%*%solve(t(X[, -j])%*%X[, -j])%*%t(X[, -j])
    c <- sqrt((sum(X[, j]^2) - sum((P%*%X[, j])^2))/(sum(Z^2) - sum((P%*%Z)^2)))
    
    fit <- glmnet(cbind(X[,j] + c*Z, X[,j] - c*Z, X[, -j]),y  ,
                   family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    # fit <- glm(y ~ cbind(X[,j] + tau[j]*Z, X[,j] - tau[j]*Z, X[, -j]) - 1, family = 'binomial')
    beta1 <- as.numeric( fit$beta)[1]
    beta2 <- as.numeric( fit$beta)[2]
    M[j] <- beta1*beta2*(tau[j]^2 + c^2)
  }
  select_index <- analys(M, abs(M), q)
  num_select <- length(select_index)
  return(list(select_index = select_index, num_select = num_select))
}
