#simulation
#'
#' This function allows you to find risk level automatically with lgbcv object.
#' @param size sample size of the simulation data
#' @param dim dimension of parameters
#' @param lambda scale parameter
#' @param vu shape parameter,defaut is 1 for exponential distribution
#' @param beta coefficients
#' @param c_rate approximate censor rate
#'
#' @import stats
#' @export Xsurv_sim_data
#' @return a simulated dataframe
#'
#' @example R/example/example_sim.R
#'
Xsurv_sim_data<-function(size,dim,lambda=NULL,vu=1,
                         beta=NULL,c_rate){
  n=size
  p=dim
  V<-diag(p)

  x <- matrix(stats::rnorm(n*p), n, p)

  if(is.null(beta)){

    beta=rep(1,p)
  }
  if(is.null(lambda)){

    lambda=2
  }
  if(vu<0){

    stop('vu must be larger than 0')
  }
  mu<-exp(drop(x %*% beta))*lambda
  #vu=1 is exponential distribution,vu>1 is weibull distribution
  real_time<--(log(runif(n)))/((mu)^(1/vu))
  if(c_rate<1)
  {lam_c<-(lambda*c_rate/(1-c_rate))^vu}

  mu_c<-exp(drop(x %*% beta))*lam_c
  fl_time<--(log(runif(n)))/((mu_c)^(1/vu))
  status<-as.numeric(real_time<=fl_time)
  if(c_rate>=1) status=rep(0,size)
  time<-pmin(real_time,fl_time)
  sim_dat<-data.frame(time,status)
  x<-as.data.frame(x)
  sim_dat<-cbind(x,sim_dat)
  sim_dat
  }
