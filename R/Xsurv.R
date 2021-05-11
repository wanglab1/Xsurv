#main function
#'
#' This function can fit the survival data with different model including xgboost,lightgbm,random forests, and gbm.
#' @param datax X data set
#' @param datay Y data set including time and event status
#' @param top_n number of top features for survival tree fitting(avaliable for xgb and lgb)
#' @param option model fitting option,defaut is xgb
#' @param method fitting metohd,defaut is 'pl' means using loss function:coxph likelihood
#' @param nfolds number of folds for crossvalidation
#' @param number of rounds
#' @param lambda l1 penalty parameter
#' @param alpha l2 penalty parameter
#' @param eta learning rate
#' @param early_stopping_rounds force a stopping round
#' @param gtree number of trees using in gbm model
#' @param ncores number of cores using in gbm
#' @param gfrac bag fraction in gbm with defaut 0.5
#' @param gsh shrinkage paramter in gbm with defaut 0.001
#' @param rfnsp number of splits in random forests
#' @import xgboost
#' @import lightgbm
#' @import rpart
#' @import SHAPforxgboost
#' @import partykit
#' @import survival
#' @import dplyr
#' @import magrittr
#' @export Xsurv
#' @return a list object containing:model,cindex,tree,SHAP and risk
#'
#' @example R/example/example.R


Xsurv<-function(datax,datay,top_n=NULL,option=c('defaut','xgb','lgb','gbm','rf'),method=c('defaut','pl','C'),nfolds=5,
                nround=NULL,lambda=NULL,alpha=NULL,eta=NULL,early_stopping_rounds=NULL,gtree=NULL,
                ncores=NULL,gfrac=NULL,gsh=NULL,rfnsp=NULL,cp=NULL,maxdpth=NULL)
{

  x_train=as.data.frame(datax)
  y=as.data.frame(datay)
  d_train<-x_train
  y_train=survival::Surv(y$time,y$status)
  d_train<-magrittr::`%<>%`(d_train,dplyr::mutate(yy = y_train))

  yy<-y_train

  option=match.arg(option)

  method=match.arg(method)
  sp_tree<-NULL
  sh<-NULL
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
  if(is.null(gtree))
    gtree=1000
  if(is.null(ncores))
    ncores=1
  if(is.null(gfrac))
    gfrac=0.5
  if(is.null(gsh))
    gsh=0.001
  if(is.null(rfnsp))
    rfnsp=10
  tt<-length(x_train[,1])
  y_train_boost <-  2 * y$time * (y$status - .5) #make fisrt col status and second col time
  #y_train<-surv_time
  if(option=='lgb'){

    if(method=='C')
      model<-lgb.sur(datax,datay,method = 'C',nfolds = nfolds,nround=nround,lambda=lambda,
                     alpha = alpha,eta=eta,
                     early_stopping_rounds = early_stopping_rounds)
    else   model<-lgb.sur(datax,datay,nfolds = nfolds,nround=nround,lambda=lambda,
                          alpha = alpha,eta=eta,
                          early_stopping_rounds = early_stopping_rounds)

    x_train=data.matrix(x_train)
    y_lgcox_pred <- matrix(0, tt, nfolds)
    for (i in 1:nfolds) {
      y_lgcox_pred[, i] <- -predict(model$boosters[[i]]$booster, x_train)

    }
    lgbcx<-rep(0,nfolds)
    for(i in 1:nfolds){
      aa=(as.matrix(y_lgcox_pred[,i]))
      lgbcx[i]<-survival::concordance(y_train ~ aa)$con


    }
    k=which.max(lgbcx)
    mod<-model$boosters[[k]]$booster
    cdx<-lgbcx[k]
    sp_tree<-sim_surv_tree(mod,x_train,datay,top_n,maxdpth,cp)
    sh=sh=SHAPforxgboost::shap.plot.summary.wrap1(mod,x_train,top_n = top_n)
    xrisk<-surv_risk_aut(mod,datax,datax)

  } else if(option=='gbm'){

    mod <- gbm::gbm(yy ~ .,       # formula
                    data = d_train,                 # dataset
                    distribution = "coxph",
                    n.trees =gtree ,              # number of trees
                    shrinkage = gsh,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                    bag.fraction = gfrac,        # subsampling fraction, 0.5 is probably best
                    train.fraction = 0.8,      # fraction of data for training, first train.fraction*N used for training
                    cv.folds = nfolds,              # do 5-fold cross-validation
                    keep.data = T,
                    verbose = T,
                    n.cores = 1
    )

    y_gbm_predict <- exp(gbm::predict.gbm(mod))

    cdx<-survival::concordance(yy~y_gbm_predict)$con
    xrisk<-surv_risk_aut_gbm(mod,datax,datax)


  } else if(option=='rf'){
    imp.d_train <- impute.rfsrc(yy ~ ., data = d_train, nsplit = rfnsp)
    ncl=ncol(d_train)
    imp.d_train[,ncl]=d_train[,ncl]


    mod <- randomForestSRC::rfsrc(yy ~ ., data = imp.d_train, nsplit = rfnsp,importance = TRUE)
    print(mod$n)
    y_rf_predict <- randomForestSRC::predict.rfsrc(object=mod)

    y_rf_predict1<-y_rf_predict$regrOutput$yy$predicted
    print(length(y_rf_predict1))
    cdx<-survival::concordance(yy ~ y_rf_predict1)$con

    xrisk<-surv_risk_aut_rf(mod,datax,datax)

    print('here')


  } else  {

    if(method=='C')
      model<-xgb.sur(datax,datay,method = 'C',nfolds = nfolds,nround=nround,lambda=lambda,
                     eta=eta,
                     alpha = alpha,early_stopping_rounds = early_stopping_rounds)
    else      model<-xgb.sur(datax,datay,nfolds = nfolds,nround=nround,lambda=lambda,
                             eta=eta,
                             alpha = alpha,early_stopping_rounds = early_stopping_rounds)

    y_xgbcox_pred <- matrix(0, tt, nfolds)
    x_train=data.matrix(x_train)
    for (i in seq_len(5)) {

      y_xgbcox_pred[, i] <- -predict(model$models[[i]], x_train)

    }
    xgbcx<-rep(0,nfolds)

    for(i in 1:nfolds){
      aa=(as.matrix(y_xgbcox_pred[,i]))
      xgbcx[i]<-survival::concordance(y_train ~ aa)$con

    }

    k=which.max(xgbcx)
    mod<-model$models[[k]]
    cdx<-xgbcx[k]
    sp_tree<-sim_surv_xgb_tree(mod,x_train,y,top_n,maxdpth,cp)
    sh=SHAPforxgboost::shap.plot.summary.wrap1(mod,x_train,top_n = top_n)
    xrisk<-surv_risk_aut(mod,datax,datax)



  }


  y$risk<-factor(xrisk,levels=c('High Risk','Medium Risk','Low Risk'))
  kmrisk<-survival::survfit(Surv(time,status)~risk,data=y)

  risk<-list('fit'=kmrisk,'data'=y)
  ls<-list('model'=mod,'cindex'=cdx,'tree'=sp_tree,'SHAP'=sh,'risk'=risk)
  ls
}



