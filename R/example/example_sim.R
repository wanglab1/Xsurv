#
sim_dat<-Xsurv_sim_data(size=500,dim=20,lambda=2,vu=1,
                                    c_rate=0.3)


sim_x<-sim_dat[,1:20]
sim_y<-sim_dat[,c(21,22)]

fit<-Xsurv.cv(sim_x,sim_y,top_n=5)
