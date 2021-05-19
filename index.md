# Xsurv <img src="docs/figures/Xlogo.png" align="right" alt="" width="120" />

# `Xsurv` --- Efficient grdient boosting for survival data  in R

## Why `Xsurv`?
- **Work with any type of survival data.** 
- **One of the simplest syntax.** 
- **Popular machine learning methods included** Xgboost, Lightgbm, random forests, and gbm.
- **No parameters tuning bothering with cross-validation**: everything can be done automatically! - **Robust model determined risk levels.**
- **Know the important  features associated with survival outcome**:clearly presented by SHAP importance plot. 

### Let's Start
```r
devtools::install_github("topycyao/Xsurv")          # Install the package
```



---

## Table of Contents

- [Introduction](#introduction)
- [Contact](#contact)
- [Installation](#installation)
- [Features](#features)
- [Quick Start](#quick-start)
- [Bugs and Issues](#bugs-and-issues)
- [Citation](#citation)

## Introduction

`Xsurv` is a useful and timely computational tool to identify candidate biomarkers that can predict or modulate patient prognosis, which may facilitate cancer translational and clinical research. 


## Contact
Create a ticket with a bug or question on [GitHub Issues](https://github.com/topycyao/Xsurv/issues) to help you and enrich it with your experience. 


## Installation
### Latest release on GitHub
```{r}
install.packages('devtools') # Ignore this if devtools is already installed.
devtools::install_github('topycyao/Xsurv')
```
## Features

1.  Easy manipulation of survival data:

    + The only thing you need to do is to divide data into two parts:covariates (x) and survival outcomes ;
  
    + Support different algorithms and keep updating;

    + No effort needed to tune your model even you have no experience:simply run Xsurv.cv and everything is done;
    

2. Prognostic biomarker discovery analysis made simple:

    + Directly know the top n features in fitted model;

    + A re-construct survival tree with important features helps to understand;
  

3. Model determined risk levels and robust predictions for survival probabilty: 

    
    
## Quick start

### Load the package

```{r}
library(Xsurv) #Load Xsurv into R

```
### Xsurv can help to generate survival data
```{r}
sim_dat<-Xsurv_sim_data(size=500,dim=20,lambda=2,vu=1, 
                                    c_rate=0.3)   # A data set is generated with sample size =500

#Covariates and survival outcome should be separted before fitting to Xsurv models
sim_x<-sim_dat[,1:20] # The first 20 (equal to dimension of covariates) columns
sim_y<-sim_dat[,c(21,22)] # The last 2 columns
```
### Quick fit from data
```{r}

fit<-Xsurv.cv(sim_x,sim_y,top_n=5)

```

## Bugs and Issues
All bug reports, documentation improvements, enhancements and ideas are appreciated. Just let us know via [GitHub](https://github.com/topycyao/Xsurv/issues).


## Citation

Li,K.et al. (2021). Efficient gradient boosting for prognostic biomarker discovery.


