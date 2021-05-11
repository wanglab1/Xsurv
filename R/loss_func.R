#'  define objective function and evaluation function
#'
#'
#'
#'
#'  @export Cox_lgb_obj

Cox_lgb_obj <- function(preds, dtrain) {
  y_true <- lightgbm::getinfo(dtrain, "label") #labels<-dtrain$label
  #print(labels)
  censor <- as.numeric(y_true > 0)
  ord <- order(y_true)
  ran <- rank(y_true)
  #print(ord)
  #foo<<-censor
  #compute the first gradient
  d <- censor[ord]  #status
  etas <- preds[ord] #linear predictor
  #print(etas)
  haz <- as.numeric(exp(etas)) #w[i]
  #print(haz)
  rsk <- rev(cumsum(rev(haz))) #W[i]
  P <- outer(haz, rsk, '/')

  P[upper.tri(P)] <- 0
  grad<- -(d - P %*% d)
  grad <- grad[ran]
  #print(grad)
  #derive hessian
  # H1=outer(haz,rsk^2,'/')
  # H1=t(t(H1)*rsk)
  H1 <- P
  H2 <- outer(haz ^ 2, rsk ^ 2, '/')
  H <- H1 - H2
  H[upper.tri(H)] <- 0
  hess <- H %*% d
  hess <- hess[ran]
  #hess=rep(0,length(grad))
  # Return the result as a list
  return(list(grad = grad, hess = hess))
}

### XGB C-index objective
cidx_xgb_obj <- function(preds, dtrain) {
  alpha <- 2
  y_true <- xgboost::getinfo(dtrain, 'label')
  censor <- as.numeric(y_true > 0)

  ord <- order(y_true)
  ran <- rank(y_true)
  #print(ord)
  #foo<<-censor
  #compute the first gradient
  d <- censor[ord]  #status
  etas <- preds[ord] #linear predictor
  n <- length(etas)
  temp1l <- matrix(etas[1:(n - 1)], nrow = n - 1, ncol = n - 1) -
    matrix(etas[2: n], nrow = n - 1, ncol = n - 1, byrow = 1)
  #print(temp1l)
  temp2l <- exp(alpha) ^ temp1l
  temp2l[lower.tri(temp2l)] <- 0
  temp3l <- (1 + temp2l)^2
  #print(temp3l)
  temp4l <- -alpha*temp2l / temp3l
  temp4l[lower.tri(temp4l)] <- 0
  #print(temp4l)
  temp4r <- -temp4l
  gl <- c(apply(temp4l, 1, sum), 0)*d
  gr <- c(0, -apply(temp4l, 2, sum))
  grad <- gl + gr
  grad <- grad[ran]
  #print(grad)
  #derive hessian
  # H1=outer(haz,rsk^2,'/')
  # H1=t(t(H1)*rsk)
  temp5l <- -(alpha*temp2l*temp3l-temp2l*2*temp3l*alpha*temp2l)/(1+temp2l)^4
  temp5l[lower.tri(temp5l)] <- 0
  hl <- c(apply(temp5l, 1, sum), 0)*d
  hr <- c(0,apply(temp5l, 2, sum))
  hess <- hl + hr
  hess <- hess[ran]
  #hess=rep(0,length(grad))
  # Return the result as a list
  return(list(grad = grad, hess = hess))
}


#'  define objective function and evaluation function
#'
#'
#'
#'
#'  @export cidx_lgb_obj
cidx_lgb_obj <- function(preds, dtrain) {
  alpha <- 2
  y_true <- lightgbm::getinfo(dtrain, 'label')
  censor <- as.numeric(y_true > 0)

  ord <- order(y_true)
  ran <- rank(y_true)
  #print(ord)
  #foo<<-censor
  #compute the first gradient
  d <- censor[ord]  #status
  etas <- preds[ord] #linear predictor
  n <- length(etas)
  temp1l <- matrix(etas[1:(n - 1)], nrow = n - 1, ncol = n - 1) -
    matrix(etas[2: n], nrow = n - 1, ncol = n - 1, byrow = 1)
  #print(temp1l)
  temp2l <- exp(alpha) ^ temp1l
  temp2l[lower.tri(temp2l)] <- 0
  temp3l <- (1 + temp2l)^2
  #print(temp3l)
  temp4l <- -alpha*temp2l / temp3l
  temp4l[lower.tri(temp4l)] <- 0
  #print(temp4l)
  temp4r <- -temp4l
  gl <- c(apply(temp4l, 1, sum), 0)*d
  gr <- c(0, -apply(temp4l, 2, sum))
  grad <- gl + gr
  grad <- grad[ran]
  #print(grad)
  #derive hessian
  # H1=outer(haz,rsk^2,'/')
  # H1=t(t(H1)*rsk)
  temp5l <- -(alpha*temp2l*temp3l-temp2l*2*temp3l*alpha*temp2l)/(1+temp2l)^4
  temp5l[lower.tri(temp5l)] <- 0
  hl <- c(apply(temp5l, 1, sum), 0)*d
  hr <- c(0,apply(temp5l, 2, sum))
  hess <- hl + hr
  hess <- hess[ran]
  #hess=rep(0,length(grad))
  # Return the result as a list
  return(list(grad = grad, hess = hess))
}

### XGB Cox deviance function
#'  define objective function and evaluation function
#'
#'
#'
#'
#'  @export Cox_xgb_deviance_func
#'
#'
Cox_xgb_deviance_func <- function(preds, dtrain) {
  y_true <- xgboost::getinfo(dtrain, "label") #labels<-dtrain$label
  censor <- as.numeric(y_true > 0) #not working!
  #foo<<-censor
  #compute the first gradient
  ord <- order(y_true)
  d <- censor[ord]  #status
  etas <- preds[ord] #linear predictor
  haz <- as.numeric(exp(etas)) #w[i]
  rsk <- rev(cumsum(rev(haz)))
  err <- -sum(d*(etas - log(rsk)))
  return(list(metric = "deviance", value = err))
}

### LGB Cox deviance function
#'  define objective function and evaluation function
#'
#'
#'
#'
#'  @export Cox_lgb_deviance_func
#'
Cox_lgb_deviance_func <- function(preds, dtrain) {
  y_true <- lightgbm::getinfo(dtrain, "label") #labels<-dtrain$label
  censor <- as.numeric(y_true > 0) #not working!
  #foo<<-censor
  #compute the first gradient
  ord <- order(y_true)
  d <- censor[ord]  #status
  etas <- preds[ord] #linear predictor
  haz <- as.numeric(exp(etas)) #w[i]
  rsk <- rev(cumsum(rev(haz)))
  err <- -sum(d*(etas - log(rsk)))
  return(list(metric = "deviance", value = err))
}

### XGB C-index deviance function
#'  define objective function and evaluation function
#'
#'
#'
#'
#'  @export cidx_xgb_deviance_func
#'
cidx_xgb_deviance_func <- function(preds, dtrain) {
  y_true <- xgboost::getinfo(dtrain, "label") #labels<-dtrain$label
  censor <- as.numeric(y_true > 0)
  #foo<<-censor
  #compute the first gradient
  ord<-order(y_true)
  d=censor[ord]  #status
  etas=preds[ord] #linear predictor
  haz<-as.numeric(exp(etas)) #w[i]
  rsk<-rev(cumsum(rev(haz)))
  err <- -sum(d*(etas-log(rsk)))
  return(list(metric = "deviance", value = err))
}

### LGB Cox deviance function
cidx_lgb_deviance_func <- function(preds, dtrain) {
  y_true <- lightgbm::getinfo(dtrain, "label") #labels<-dtrain$label
  censor <- as.numeric(y_true > 0)
  #foo<<-censor
  #compute the first gradient
  ord <- order(y_true)
  d <- censor[ord]  #status
  etas <- preds[ord] #linear predictor
  haz <- as.numeric(exp(etas)) #w[i]
  rsk <- rev(cumsum(rev(haz)))
  err <- -sum(d*(etas-log(rsk)))
  return(list(name = 'deviance', value = err, higher_better = F))
}

### XGB C-index
#'  define objective function and evaluation function
#'
#'
#'
#'
#'  @export cidx_xgb_func
#'
#'
cidx_xgb_func <- function(preds, dtrain) {
  y_true <- xgboost::getinfo(dtrain, "label")
  censor <- as.numeric(y_true > 0)
  surv_t <- survival::Surv(abs(y_true), censor)
  return(list(metric = 'cidx', value = concordance(surv_t ~ preds)$con))
}

### LGB C-index
#'  define objective function and evaluation function
#'
#'
#'
#'
#'  @export cidx_lgb_func
#'
cidx_lgb_func <- function(preds, dtrain) {
  y_true <- lightgbm::getinfo(dtrain, "label")
  censor <- as.numeric(y_true > 0)
  surv_t <- survival::Surv(abs(y_true), censor)
  return(list(name = 'cidx', value = concordance(surv_t ~ preds)$con, higher_better = T))
}
