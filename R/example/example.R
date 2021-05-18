#example
library(survival)
library(Xsurv)
data(lung)

mydata<-(lung[,-1])
mydata[,2]<-mydata[,2]-1
length(mydata[,1])
names(mydata)<-colnames(mydata)
datay_train<-mydata[1:180,c(1,2)]
datax_train<-mydata[1:180,-c(1,2)]
datay_test<-mydata[181:228,c(1,2)]
datax_test<-mydata[181:228,-c(1,2)]
xs<-Xsurv(datax_train,datay_train,top_n = 5,cp=0.01)

#xs<-Xsurv.cv(datax_train,datay_train,top_n=5)
xm<-xs$model
xtree<-xs$tree
x_ctree<-xtree$tree2
#plot(x_ctree)
shap=xs$SHAP

shap

risk=xs$risk
fit=risk$fit

#plot(fit)
#prediction

pre_time<-pre<-Xsurv_predict(xm,datax_train,datay_train,datax_test)

#predict survival probabilty

pre_x<-Xsurv_predict_sv(xm,datax_train,datay_train,datax_test[1,])

plot(pre_x)

