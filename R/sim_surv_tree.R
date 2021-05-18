#' Simple survival tree
#'
#' This function allows you to find risk level automatically with xgbor lgb object.
#' @param model lightgbm model object
#' @param x_data X data set
#' @param y_data Y data set
#' @param top_nc numbers of top features selected
#' @param  maxdp max depth for rpart tree
#' @param  cp complexity param for rpart tree
#' @import SHAPforxgboost
#' @import rpart
#' @import partykit
#' @import lightgbm
#' @import xgboost
#' @import stats
#' @export sim_surv_tree
#'
#' @return a list contains 1.rpart tree;2.ctree
#'

sim_surv_tree<-function(model,x_data,y_data,top_nc=NULL,maxdp=NULL,cp=NULL){
  y_data=as.data.frame(y_data)
  cnames<-colnames(x_data)
  xdata=data.matrix(x_data)
  if(is.null(top_nc)){top_nc=3}
  if(is.null(maxdp)){maxdp=3}


  imp<-SHAPforxgboost::shap.values(model,xdata)$mean_shap_score
  top_name<-names(imp)

  idx=rep(0,top_nc)
  for(i in 1:top_nc){
    idx[i]<-which(cnames==top_name[i])
  }

  x_tree<-as.data.frame(x_data[,idx])

  yt<-survival::Surv(y_data$time,y_data$status)
  fit<-rpart::rpart(yt~.,data=x_tree,control = rpart.control(maxdepth=maxdp,cp=cp))
  tfit<-partykit::as.party(fit)
  tfit1<-partykit::as.party(fit)
  plot(tfit1)
  tfit2<-partykit::ctree(yt~.,data=x_tree)
  ls<-list('tree1'=tfit1,'tree2'=tfit2)

  ls

}
