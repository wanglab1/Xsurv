#CV for Xsurv
#'
#' This function can use cv to automatically tuning the paramter to fit the survival model for xgboost and lightgbm.
#' @param datax X data set
#' @param datay Y data set including time and event status
#' @param top_n number of top features for survival tree fitting
#' @param option model fitting option,defaut is xgb
#' @param method fitting metohd,defaut is 'pl' means using loss function:coxph likelihood
#' @param nfolds number of folds for crossvalidation
#' @param number of rounds
#' @param lambda l1 penalty parameter
#' @param alpha l2 penalty parameter
#' @param eta learning rate
#' @param early_stopping_rounds force a stopping round
#' @import xgboost
#' @import lightgbm
#' @import rpart
#' @import SHAPforxgboost
#' @import partykit
#' @import survival
#' @import dplyr
#' @export Xsurv.cv
#' @return a list object containing:model,cindex,tree,SHAP and risk
#'
#' @example R/example/example.R
#'


Xsurv.cv<-function(datax,datay,top_n=NULL,option=c('defaut','xgb','lgb','gbm','rf'),method=c('defaut','pl','C'),
                   search=c('rd','grid'),nfolds=5,cvfrac=0.8,rdtime=10,
                   nround=NULL,Lambda=NULL,Alpha=NULL,Eta=NULL,early_stopping_rounds=NULL,
                   cp=NULL,maxdpth=NULL)
{
  tt=nrow(datax)
  n_train=as.integer(tt*cvfrac)
  n_test=tt-n_train
  x_train=as.data.frame(datax[1:n_train,])
  x_test=data.matrix(datax[(n_train+1):tt,])
  y=as.data.frame(datay)
  yy_train<-y[1:n_train,]
  yy_test<-y[(n_train+1):tt,]
  yt=survival::Surv(y$time,y$status)
  d_train<-x_train
  y_train=yt[1:n_train]
  y_test=yt[(n_train+1):tt]
  d_train<-magrittr::`%<>%`(d_train,dplyr::mutate(yy = y_train))
  d_test=as.data.frame(x_test)
  yy<-y_train

  option=match.arg(option)
  search=match.arg(search)

  method=match.arg(method)
  sp_tree<-NULL
  sh<-NULL
  if(is.null(Lambda))
    Lambda=c(0.01,0.05,0.1,0.2,0.5)
  if(is.null(Alpha))
    Alpha=c(0.01,0.05,0.1,0.2,0.5)
  if(is.null(Eta))
    Eta=c(0.01,0.05,0.1,0.2,0.5)
  if(is.null(nround))
    nround=1000
  if(is.null(early_stopping_rounds))
    early_stopping_rounds=20

  lam=0
  alp=0
  et=0
  y_boost <-  2 * y$time * (y$status - .5)
  y_train_boost<-y_boost[1:n_train]
  y_test_boost <- y_boost[(n_train+1):tt]

  if(search=='grid')
  { #grid search
    cdx=0

    for(lambda in Lambda)
    { for (alpha in Alpha) {
      for(eta in Eta)
      {
        if(option=='lgb'){

          if(method=='C')
            model<-lgb.sur(x_train,yy_train,method = 'C',nfolds = nfolds,nround=nround,lambda=lambda,
                           alpha = alpha,eta=eta,
                           early_stopping_rounds = early_stopping_rounds)
          else   model<-lgb.sur(x_train,yy_train,nfolds = nfolds,nround=nround,lambda=lambda,
                                alpha = alpha,eta=eta,
                                early_stopping_rounds = early_stopping_rounds)


          y_lgcox_pred <- matrix(0, n_test, nfolds)
          for (i in 1:nfolds) {
            y_lgcox_pred[, i] <- -predict(model$boosters[[i]]$booster, x_test)

          }
          lgbcx<-rep(0,nfolds)
          for(i in 1:nfolds){
            aa=(as.matrix(y_lgcox_pred[,i]))
            lgbcx[i]<-survival::concordance(y_test ~ aa)$con


          }
          k=which.max(lgbcx)
          if(lgbcx[k]>cdx)
          {
            cdx<-lgbcx[k]
            mod<-model$boosters[[k]]$booster
            lam=lambda
            alp=alpha
            et=eta
          }




        }  else  {

          if(method=='C')
            model<-xgb.sur(x_train,yy_train,method = 'C',nfolds = nfolds,nround=nround,lambda=lambda,
                           alpha = alpha,eta=eta,
                           early_stopping_rounds = early_stopping_rounds)
          else      model<-xgb.sur(x_train,yy_train,nfolds = nfolds,nround=nround,lambda=lambda,
                                   alpha = alpha,eta=eta,
                                   early_stopping_rounds = early_stopping_rounds)

          y_xgbcox_pred <- matrix(0, n_test, nfolds)

          for (i in seq_len(5)) {

            y_xgbcox_pred[, i] <- -predict(model$models[[i]], x_test)

          }
          xgbcx<-rep(0,nfolds)

          for(i in 1:nfolds){
            aa=(as.matrix(y_xgbcox_pred[,i]))
            xgbcx[i]<-survival::concordance(y_test ~ aa)$con

          }
          k=which.max(xgbcx)
          if(xgbcx[k]>cdx)
          {
            cdx<-xgbcx[k]
            mod<-model$models[[k]]
            lam=lambda
            alp=alpha
            et=eta
          }
          print(cdx)




        }
      }
    }
    }
  } else {
    ##random search
    cdx=0
    for (rt in 1:rdtime) {

      lambda<-stats::runif(1,0.01,0.8)
      alpha<-stats::runif(1,0.01,0.8)
      eta<-stats::runif(1,0.01,0.8)
      if(option=='lgb'){

        if(method=='C')
          model<-lgb.sur(x_train,yy_train,method = 'C',nfolds = nfolds,nround=nround,lambda=lambda,
                         alpha = alpha,early_stopping_rounds = early_stopping_rounds)
        else   model<-lgb.sur(x_train,yy_train,nfolds = nfolds,nround=nround,lambda=lambda,
                              alpha = alpha,early_stopping_rounds = early_stopping_rounds)


        y_lgcox_pred <- matrix(0, n_test, nfolds)
        for (i in 1:nfolds) {
          y_lgcox_pred[, i] <- -predict(model$boosters[[i]]$booster, x_test)

        }
        lgbcx<-rep(0,nfolds)
        for(i in 1:nfolds){
          aa=(as.matrix(y_lgcox_pred[,i]))
          lgbcx[i]<-survival::concordance(y_test ~ aa)$con


        }
        k=which.max(lgbcx)
        if(lgbcx[k]>cdx)
        {
          cdx<-lgbcx[k]
          mod<-model$boosters[[k]]$booster
          lam=lambda
          alp=alpha
          et=eta
        }




      }  else  {

        if(method=='C')
          model<-xgb.sur(x_train,yy_train,method = 'C',nfolds = nfolds,nround=nround,lambda=lambda,
                         alpha = alpha,eta=eta,
                         early_stopping_rounds = early_stopping_rounds)
        else

          model<-xgb.sur(x_train,yy_train,nfolds = nfolds,nround=nround,lambda=lambda,
                         alpha = alpha,eta=eta,
                         early_stopping_rounds = early_stopping_rounds)


        y_xgbcox_pred <- matrix(0, n_test, nfolds)

        for (i in seq_len(5)) {

          y_xgbcox_pred[, i] <- -predict(model$models[[i]], x_test)

        }
        xgbcx<-rep(0,nfolds)

        for(i in 1:nfolds){
          aa=(as.matrix(y_xgbcox_pred[,i]))
          xgbcx[i]<-survival::concordance(y_test ~ aa)$con

        }
        k=which.max(xgbcx)
        if(xgbcx[k]>cdx)
        {
          cdx<-xgbcx[k]
          mod<-model$models[[k]]
          lam=lambda
          alp=alpha
          et=eta
        }

      }


    }






  }

  sp_tree<-sim_surv_tree(mod,x_train,yy_train,top_n,maxdpth,cp)
  xx_train=data.matrix(datax)
  sh=SHAPforxgboost::shap.plot.summary.wrap1(mod,xx_train,top_n = top_n)

  xrisk<-surv_risk_aut(mod,x_train,datax)
  y$risk<-factor(xrisk,levels=c('High Risk','Medium Risk','Low Risk'))
  kmrisk<-survival::survfit(Surv(time,status)~risk,data=y)

  risk<-list('fit'=kmrisk,'data'=y)
  ls<-list('model'=mod,'cindex'=cdx,'tree'=sp_tree,'SHAP'=sh,'risk'=risk,'lambda'=lam,
           'alpha'=alp,'eta'=et)
  ls
}
