#' Transform feature space to high dimension leaf space for xgb model
#'
#' This function transforms the feature space to a high dimension leaf space
#' @param model xboost model object
#' @param datax X data set
#' @param datay Y data set
#' @import xgboost
#' @import stats
#'
#' @export leaf_xgb
#' @return a one hot-fashion feature space
leaf_xgb<-function(model,datax,datay){
  x_train=data.matrix(datax)
  y=as.data.frame(datay)
  y_train_boost <-  2 * y$time * (y$status - .5)
  XDtrain <- xgboost::xgb.DMatrix(x_train, label = y_train_boost)
  dt<-predict(model,XDtrain,predleaf=TRUE)
  ntree<-length(dt[1,])
  ndata<-length(x_train[,1])
  nl<-max(dt)
  leaf_ind<-matrix(0,ndata,nl)
  leaf_tree<-array(leaf_ind,dim = c(ndata,nl,ntree))
  for(i in 1:ntree)
  {
    leaf_ind<-matrix(0,ndata,nl)
    for (j in 1:ndata) {
      indx<-dt[j,i]
      leaf_ind[j,indx]<-1
    }
    leaf_tree[,,i]<-leaf_ind
  }
  lf_pred<-matrix(0,ndata,nl*ntree)
  for (k in 1:ndata) {
    for (j in 1:ntree) {
      a=1+nl*(j-1)
      b=nl*j

      lf_pred[k,a:b]<-leaf_tree[k,,j]
    }
  }
  lf_pred
}



#' Transform feature space to high dimension leaf space for xgb or lgb model
#'
#' This function transforms the feature space to a high dimension leaf space
#' @param model  model object
#' @param datax X data set
#' @param datay Y data set
#' @import xgboost
#' @import lightgbm
#' @import stats
#'
#' @export leaf_tf
#' @return a one hot-fashion feature space
leaf_tf<-function(model,datax){
  x_train=data.matrix(datax)

  dt<-predict(model,x_train,predleaf=TRUE)
  ntree<-length(dt[1,])
  ndata<-length(x_train[,1])
  nl<-max(dt)
  leaf_ind<-matrix(0,ndata,nl)
  leaf_tree<-array(leaf_ind,dim = c(ndata,nl,ntree))
  for(i in 1:ntree)
  {
    leaf_ind<-matrix(0,ndata,nl)
    for (j in 1:ndata) {
      indx<-dt[j,i]
      leaf_ind[j,indx]<-1
    }
    leaf_tree[,,i]<-leaf_ind
  }
  lf_pred<-matrix(0,ndata,nl*ntree)
  for (k in 1:ndata) {
    for (j in 1:ntree) {
      a=1+nl*(j-1)
      b=nl*j

      lf_pred[k,a:b]<-leaf_tree[k,,j]
    }
  }
  lf_pred
}
