#Survival time prediction
#'
#' This function allows you to transform survival prediction from risk score to survival time.
#' @param model xgboost or lightgbm model object
#' @param x_train X train data set
#' @param y_train Y trian data set
#' @param x_test X test data set
#' @param nc number of cutting points defaut is 25, and 20-30 is recommended
#' @param nq quantile of train data time in each cut group for prediction time
#' @import xgboost
#' @import lightgbm
#' @import stats
#' @import class
#' @export Xsurv_predict
#' @return predicted survival time
#'
#' @example R/example/example.R
#'
Xsurv_predict<-function(model,x_train,y_train,x_test,nc=25,nq=0.8)
{

  train_data<-data.matrix(x_train)
  n_train<-length(train_data[,1])
  test_data<-data.matrix(x_test)


  full_data<-rbind(train_data,test_data)
  n_fu<-length(full_data[,1])

  pred_lf<-as.data.frame(leaf_tf(model,full_data))

  pred_train_lf<-pred_lf[1:n_train,]
  pred_test_lf<-pred_lf[(n_train+1):n_fu,]



  l_pred=pred_train_lf
  l_pred_test=pred_test_lf

  cl25=stats::kmeans(l_pred,nc)
  clus<-cl25$cluster



  time_train<-as.data.frame(y_train$time)

  time_train$clus<-clus

  kpre_test=class::knn(l_pred,l_pred_test,cl = clus)

  pre_cl<-rep(0,nc)
  for(j in 1:nc)
  {
    ind<-which(time_train[,2]==j)

    #censored data included,hence nq should be something over 0.5,defaut is 0.8
    pre_cl[j]=quantile(time_train[ind,1],nq)

  }
  # pre_adj<-sort(pre_cl)








  surv_pred<-pre_cl[kpre_test]


  d_pred<-data.frame(time=surv_pred)

  d_pred

}

#Survival probabilty prediction
#'
#' This function allows you to predict a patient's survival curve
#' @param model xgboost or lightgbm model object
#' @param x_train X train data set
#' @param y_train Y trian data set
#' @param x_test X test data set
#' @param nc number of cutting points defaut is 25, and 20-30 is recommended
#' @param nq quantile of train data time in each cut group for prediction time
#' @import xgboost
#' @import lightgbm
#' @import stats
#' @import class
#' @export Xsurv_predict_sv
#' @return predicted survival probabilty
#'
Xsurv_predict_sv<-function(model,x_train,y_train,x_new,nc=25)
{

  train_data<-data.matrix(x_train)
  n_train<-length(train_data[,1])
  test_data<-data.matrix(x_new)


  full_data<-rbind(train_data,test_data)
  n_fu<-length(full_data[,1])

  pred_lf<-leaf_tf(model,full_data)

  pred_lf<-as.data.frame(pred_lf)
  pred_train_lf<-pred_lf[1:n_train,]
  pred_test_lf<-pred_lf[n_fu,]


  l_pred=pred_train_lf
  l_pred_test=pred_test_lf





  cl25=stats::kmeans(l_pred,nc)
  clus<-cl25$cluster


  time_train<-as.data.frame(y_train$time)

  time_train$clus<-clus

  kpre_test=class::knn(l_pred,l_pred_test,cl = clus)

  new.cl<-kpre_test



  ind<-which(time_train[,2]==new.cl)

  d_sv<-as.data.frame(y_train)

  d_sv=d_sv[ind,]

  fit_sv<-survival::survfit(Surv(time,status)~1,data=d_sv)

  fit_sv
}
