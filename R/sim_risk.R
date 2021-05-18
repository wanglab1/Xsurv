#' Automatically define survival risk with xgb or lgb model
#'
#' Find the risk level for patient under model
#'
#' @param model fitted model (xgb or lgb)
#' @param  train_data train data (only covariates)
#' @param  test_data  test data (only covaraites)
#'
#' @import xgboost
#' @import class
#' @import stats
#' @import lightgbm
#' @import survival
#' @export surv_risk_aut
#'
#' @return a factor vector with three different risk levels
#'

surv_risk_aut<-function(model,train_data,test_data){
  train_data<-data.matrix(train_data)
  test_data<-data.matrix(test_data)



  pred_train<-as.data.frame(-predict(model,train_data))
  pred_test<-as.data.frame(-predict(model,test_data))

  cl1=stats::kmeans(pred_train,3)



  center=cl1$centers
  risklevel<-sort(center)



  r_test<-stats::kmeans(pred_train,centers = risklevel)
  ris_tran <- function(x) {
    k<-length(x)
    y<-rep(0,k)
    x<-as.numeric(x)
    for (i in 1:k) {
      if(x[i]==1) {y[i]<-'High Risk'}
      else if(x[i]==2) {y[i]<-'Medium Risk'}
      else y[i]<-'Low Risk'
    }



    return(y)
  }
  pred_train_lf<-leaf_tf(model,train_data)
  pred_test_lf<-leaf_tf(model,test_data)
  pred_train_lf<-as.data.frame(pred_train_lf)
  pred_test_lf<-as.data.frame(pred_test_lf)
  cl_pred<-class::knn(pred_train_lf,pred_test_lf,cl=r_test$cluster)


  pred_risk<-ris_tran(cl_pred)


  prisk<-factor(pred_risk,levels = c('High Risk','Medium Risk','Low Risk'))

  return(prisk)


}

#' Automatically define survival risk with gbm model
#'
#' Find the risk level for patient under a gbm model
#'
#' @param model fitted gbm model
#' @param  train_data train data (only covariates)
#' @param  test_data  test data (only covaraites)
#'
#' @import gbm
#' @import class
#' @import stats
#' @import survival
#' @export surv_risk_aut_gbm
#'
#' @return a factor vector with three different risk levels
#'
surv_risk_aut_gbm<-function(model,train_data,test_data){
  train_data<-as.data.frame(train_data)
  test_data<-as.data.frame(test_data)
  pred_train<-as.data.frame(-gbm::predict.gbm(model,train_data))
  pred_test<-as.data.frame(-gbm::predict.gbm(model,test_data))
  cl1=stats::kmeans(pred_train,3)



  center=cl1$centers
  risklevel<-sort(center)

  r_test<-stats::kmeans(pred_train,centers = risklevel)
  ris_tran <- function(x) {
    k<-length(x)
    y<-rep(0,k)
    x<-as.numeric(x)
    for (i in 1:k) {
      if(x[i]==1) {y[i]<-'High Risk'}
      else if(x[i]==2) {y[i]<-'Medium Risk'}
      else y[i]<-'Low Risk'
    }



    return(y)
  }



  cl_pred<-class::knn(pred_train,pred_test,cl=r_test$cluster)


  pred_risk<-ris_tran(cl_pred)


  prisk<-factor(pred_risk,levels = c('High Risk','Medium Risk','Low Risk'))

  return(prisk)


}

#' Automatically define survival risk with rf model
#'
#' Find the risk level for patient under a rf model
#'
#' @param model fitted rf model
#' @param  train_data train data (only covariates)
#' @param  test_data  test data (only covaraites)
#'
#' @import randomForestSRC
#' @import class
#' @import stats
#' @export surv_risk_aut_rf
#'
#' @return a factor vector with three different risk levels
#'
surv_risk_aut_rf<-function(model,train_data,test_data){
  train_data<-as.data.frame(train_data)
  test_data<-as.data.frame(test_data)
  pred_train<-as.data.frame(randomForestSRC::predict.rfsrc(mod, train_data)$regrOutput$yy$predicted)
  pred_test<-as.data.frame(randomForestSRC::predict.rfsrc(mod, test_data)$regrOutput$yy$predicted)
  cl1=stats::kmeans(pred_train,3)



  center=cl1$centers
  risklevel<-sort(center)

  r_test<-stats::kmeans(pred_train,centers = risklevel)
  ris_tran <- function(x) {
    k<-length(x)
    y<-rep(0,k)
    x<-as.numeric(x)
    for (i in 1:k) {
      if(x[i]==1) {y[i]<-'High Risk'}
      else if(x[i]==2) {y[i]<-'Medium Risk'}
      else y[i]<-'Low Risk'
    }



    return(y)
  }



  cl_pred<-class::knn(pred_train,pred_test,cl=r_test$cluster)


  pred_risk<-ris_tran(cl_pred)


  prisk<-factor(pred_risk,levels = c('High Risk','Medium Risk','Low Risk'))

  return(prisk)


}

#' Manually define survival risk with xgb or lgb model
#'
#' Find the risk level for patient under model
#'
#' @param model fitted model (xgb or lgb)
#' @param  train_data train data (only covariates)
#' @param  test_data  test data (only covaraites)
#' @param  upper upper percentil cutting point
#' @param  lower lower percentil cutting point
#' @import xgboost
#' @import class
#' @import stats
#' @import lightgbm
#' @import survival
#' @export surv_risk_m
#'
#' @return a factor vector with three different risk levels
#'

surv_risk_m<-function(model,train_data,test_data,upper,lower){
  train_data<-data.matrix(train_data)
  test_data<-data.matrix(test_data)
  pred_train<-as.data.frame(-predict(model,train_data))
  pred_test<-as.data.frame(-predict(model,test_data))

  ris_tran_tr <- function(x,a,b) {
    k<-length(x)
    y<-rep(0,k)
    x<-as.numeric(x)
    for (i in 1:k) {
      if(x[i]<=a) {y[i]<-'Low Risk'}
      else if(x[i]>=b) {y[i]<-'High Risk'}
      else y[i]<-'Medium Risk'
    }



    return(y)
  }
  h_mq=stats::quantile(exp(-pred_train[,1]),upper)

  l_mq=stats::quantile(exp(-pred_test[,1]),lower)

  cl1=ris_tran_tr(exp(-pred_train[,1]),l_mq,h_mq)



  center=cl1$centers
  risklevel<-sort(center)

  r_test<-stats::kmeans(pred_train,centers = risklevel)

  pred_train_lf<-leaf_tf(model,train_data)
  pred_test_lf<-leaf_tf(model,test_data)
  pred_train_lf<-as.data.frame(pred_train_lf)
  pred_test_lf<-as.data.frame(pred_test_lf)
  cl_pred<-class::knn(pred_train_lf,pred_test_lf,cl=r_test$cluster)




  prisk<-factor(pred_risk,levels = c('High Risk','Medium Risk','Low Risk'))

  return(prisk)


}

