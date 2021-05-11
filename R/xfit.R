#' Fitting xgb model
#'
#' This function allows you to fit survival data with a xgb model.
#' @param datax X data set
#' @param datay Y data set
#' @param method fitting method
#' @param nfolds number of folds
#' @param nround number of rounds
#' @param lambda l1 parameter
#' @param alpha l2 parameter
#' @param  eta learning rate
#' @param early_stopping_rounds force a stopping round
#' @import xgboost
#' @import survival
#' @import stats
#' @export xgb.sur
#'
#' @return a fitted xgb model
#'

xgb.sur<-function(datax,datay,method=c('defaut','pl','C'),nfolds=5,nround=NULL,
                  lambda=NULL,alpha=NULL,eta=NULL,early_stopping_rounds=NULL
)
{
  x_train=data.matrix(datax)
  y=as.data.frame(datay)
  method=match.arg(method)
  if(is.null(lambda))
    lambda=.01
  if(is.null(alpha))
    alpha=.01
  if(is.null(eta))
    eta=.01
  if(is.null(nround))
    nround=1000
  if(is.null(early_stopping_rounds))
    early_stopping_rounds=20

  tt<-length(x_train[,1])

  y_train_boost <-  2 * y$time * (y$status - .5) #make fisrt col status and second col time
  #y_train<-surv_time
  XDtrain <- xgboost::xgb.DMatrix(x_train, label = y_train_boost)
  if(method=='C')
    model<-xgboost::xgb.cv(list(objective = cidx_xgb_obj, eval_metric = cidx_xgb_func,
                                tree_method = 'hist', grow_policy = 'lossguide',
                                eta = eta, lambda = lambda, alpha = alpha, subsample = .5,
                                colsample_bytree = .5), XDtrain, nround = nround,
                           nfold = nfolds, verbose = F, early_stopping_rounds = early_stopping_rounds, maximize = T,
                           callbacks = list(cb.cv.predict(T)))
  else    model<-xgboost::xgb.cv(list(objective = 'survival:cox', eval_metric = cidx_xgb_func,
                                      tree_method = 'hist', grow_policy = 'lossguide',
                                      eta = eta, lambda = lambda, alpha = alpha, subsample = .5,
                                      colsample_bytree = .5), XDtrain, nround = nround,
                                 nfold = nfolds, verbose = F, early_stopping_rounds = early_stopping_rounds, maximize = T,
                                 callbacks = list(cb.cv.predict(T)))
  model
}



#' Fitting lgb model
#'
#' This function allows you to fit survival data with a lgb model.
#' @param datax X data set
#' @param datay Y data set
#' @param method fitting method
#' @param nfolds number of folds
#' @param nround number of rounds
#' @param lambda l1 parameter
#' @param alpha l2 parameter
#' @param  eta learning rate
#' @param early_stopping_rounds force a stopping round
#' @import xgboost
#' @import survival
#' @import stats
#' @export lgb.sur
#'
#' @return a fitted xgb model
#'
lgb.sur<-function(datax,datay,method=c('defaut','pl','C'),nfolds=5,nround=NULL,
                  lambda=NULL,alpha=NULL,eta=NULL,early_stopping_rounds=NULL
)
{
  x_train=data.matrix(datax)
  y=as.data.frame(datay)
  method=match.arg(method)
  if(is.null(lambda))
    lambda=.01
  if(is.null(alpha))
    alpha=.01
  if(is.null(eta))
    eta=.01
  if(is.null(nround))
    nround=1000
  if(is.null(early_stopping_rounds))
    early_stopping_rounds=20

  tt<-length(x_train[,1])

  y_train_boost <-  2 * y$time * (y$status - .5) #make fisrt col status and second col time
  #y_train<-surv_time
  LDtrain <- lgb.Dataset(x_train, label = y_train_boost)
  if(method=='C')
    model<-lightgbm::lgb.cv(list(objective = cidx_lgb_obj,
                                 eta = eta, lambda = lambda, alpha = alpha, subsample = .5,
                                 colsample_bytree = .5), LDtrain, nround = nround,eval=cidx_lgb_func,
                            nfold = nfolds, verbose = 0, early_stopping_rounds = early_stopping_rounds)
  else    model<-lightgbm::lgb.cv(list(objective = Cox_lgb_obj,
                                       eta = eta, lambda = lambda, alpha = alpha, subsample = .5,
                                       colsample_bytree = .5), LDtrain, nround = nround,eval=cidx_lgb_func,
                                  nfold = nfolds, verbose = 0, early_stopping_rounds = early_stopping_rounds)
  model
}

