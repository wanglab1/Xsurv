# Xsurv
This package implements XGBoost and LightGBM in survival data with partiallikelihood function of cox ph model and smoothed C-index as loss functions.
gbm and random forests are also provided in Xsurv function to fit the model.
Additional survival information can be found with XGB and LGB models,including simple survival tree,risk level and SHAP importance. 
## Installation

Please install from Github:
``` r
devtools::install_github("topycyao/Xsurv")
```
## Examples
```{r}
library(Xsurv)
library(survival)
#fit the data 
data(lung)
View(lung)

mydata<-(lung[,-1])

datay_train<-mydata[,c(1,2)]
datax_train<-mydata[,-c(1,2)]

xs<-Xsurv(datax_train,datay_train,top_n = 5,cp=0.01) #cp is the complexity of the fitted survival tree

#other options can be used 
#xs_lgb<-Xsurv(datax_train,datay_train,option='lgb',method='C',top_n=5)
#xs_gbm<-Xsurv(datax_train,datay_train,option='gbm')
#xs_rf<-Xsurv(datax_train,datay_train,option='rf')

#model can be trained automatically by Xsurv.cv
#xs.cv<-Xsurv.cv(datax_train,datay_train,top_n = 5)

```
## model analysis


SHAP importance
```{r}
shap=xs$SHAP
shap
```
<p align="center">
  <img src = "https://github.com/topycyao/Xsurv/blob/master/docs%20/figures/shaplung.png" width="600" height="500">
</p>

Suvival tree

```{r}
xm<-xs$mod
xtree<-xs$tree
tr1<-xtree$tree1
plot(tr1)
```
<p align="center">
  <img src = "https://github.com/topycyao/Xsurv/blob/master/docs%20/figures/exampletree.png"  width="500" height="400">
</p>



Risk analysis
```{r}
library(survminer)
risk=xs$risk
#plot the kaplan-meier curve for different risk levels
risk_data=risk$data
fit=survival::survfit(Surv(time,status)~risk,data=risk_data)
ggsurvplot(fit,pval = TRUE,palette = c('coral','burlywood1','cadetblue1'),size=2,
                          legend=c(0.85,0.85),legend.title='',font.x=c(15,"plain","black"),
                           font.y=c(18,"plain","black"))
```
<p align="center">
  <img src = "https://github.com/topycyao/Xsurv/blob/master/docs%20/figures/kmrisk.png"  width="600" height="500" >
</p>

## Data generation
```{r}
sim_dat1<-Xsurv_sim_data(size=500,dim=20,lambda=2,vu=1,
                         c_rate=0.3)
sim_dat2<-Xsurv_sim_data(size=500,dim=20,lambda=2,vu=2,
                         c_rate=0.3)    
beta=c(rep(1,5),rep(0,15))
sim_dat3<-Xsurv_sim_data(size=500,dim=20,lambda=2,vu=2,beta=beta,
                         c_rate=0.3) 
                         
```
## Reference

Li,K. et al. Efficient gradient boosting for prognostic biomaker discovery
