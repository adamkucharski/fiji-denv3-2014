thetatabA=NULL
theta_inittabA=NULL
s_trace_tab0_C=NULL
s_trace_tab0_A=NULL
c_trace_tab0=NULL
r_trace_tab0_C=NULL
r_trace_tab0_A=NULL
x_trace_tab0_C=NULL
x_trace_tab0_A=NULL

load(paste("outputs/outputR",country.name,"_",epi.name,pick_posterior,"_ELISA_",use.ELISA.data,exclude.p,".RData",sep=""))

thetatab=cbind(data.frame(thetatab),data.frame(thetaAlltab[,iiH,]))
theta_inittab=data.frame(theta_initAlltab[,iiH,])

mcmc_samples=length(sim_liktab)
maxB=sum(sim_liktab!=-Inf)/mcmc_samples
minB=mcmc.burn*maxB
picks0=c(round(minB*mcmc_samples):round(maxB*mcmc_samples))

sim_likOut = sim_liktab[picks0]
thetatabA=rbind(thetatabA,thetatab[picks0,])
theta_inittabA=rbind(theta_inittabA,theta_inittab[picks0,])

s_trace_tab0_C=rbind(s_trace_tab0_C,s_trace_tabC[picks0,iiH,])
s_trace_tab0_A=rbind(s_trace_tab0_A,s_trace_tabC[picks0,iiH,])
r_trace_tab0_C=rbind(r_trace_tab0_C,r_trace_tabC[picks0,iiH,])
r_trace_tab0_A=rbind(r_trace_tab0_A,r_trace_tabA[picks0,iiH,])
c_trace_tab0 = rbind(c_trace_tab0,c_trace_tab[picks0,iiH,])
x_trace_tab0_C=rbind(x_trace_tab0_C,x_trace_tabC[picks0,iiH,])
x_trace_tab0_A=rbind(x_trace_tab0_A,x_trace_tabA[picks0,iiH,])

picks=c(1:length(thetatabA[,1]))

pick.max = picks[sim_likOut[picks]==max(sim_likOut[picks])][1]

thetatab=thetatabA
theta_inittab=theta_inittabA
c_trace_tab=c_trace_tab0
s_trace_tabC=s_trace_tab0_C
s_trace_tabA=s_trace_tab0_A
r_trace_tabC=r_trace_tab0_C
r_trace_tabA=r_trace_tab0_A
x_trace_tabC=x_trace_tab0_C
x_trace_tabA=x_trace_tab0_A

