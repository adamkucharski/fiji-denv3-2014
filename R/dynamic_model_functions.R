# - - - - - - - - - - - - - - - - - - - - - - - - 
# Main model functions
# Github version
#
# Kucharski AJ et al. Transmission dynamics of Zika virus in island populations: a modelling analysis of the 2013-14 French Polynesia outbreak.
# PLoS Negl Trop Dis. 2016:10(5)
# - - - - - - - - - - - - - - - - - - - - - - - - 


calculate_r0 <- function(th_in,sus_c=1,sm_c=1,b_vary=1){
  
  # Rate humans get infected -- FORMULATION WITH ONE MOSQUITO POP AND DIFFERENT BITING RATES
  b_hv = b_vary * th_in$beta ; b_hh = 0
  
  # Rate vectors get infected
  b_vh = b_hv * th_in$beta_v 
  b_vv = 0

  rr_hh=rep(0,length(b_vary)); rr_vv = rr_hh;
  
  rr_hv = (sus_c*b_hv/th_in$mu_v)*(th_in$v_exp/(th_in$v_exp+th_in$mu_v)) 
  rr_vc = (sm_c*b_vh/th_in$r_inf)

  rr_post = NULL
  for(ii in 1:length(b_vary)){
    rr_post=c(rr_post,( max(Re(eigen(matrix(c(rr_hh[ii],rr_hv[ii],rr_vc[ii],rr_vv[ii]),nrow=2))$values)))  )
  }
  
  return( list(rr_out=rr_post,rr_mat = matrix(c(rr_hh[1],rr_hv[1],rr_vc[1],rr_vv[1]),nrow=2,byrow=T)) )
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

ReportC<-function(ct, theta,repSS){
  mu00=ct
  mu01=sapply(mu00,function(x){max(x,0)})
  

  if(cutt.date<=swap.date){
  sapply(mu01,function(x){rnbinom(1, mu=theta[["repR"]]*x,size=1/theta[["repvol"]])})
  }else{
    c(sapply(mu01[1:(theta[["rep_drop"]])] ,function(x){rnbinom(1, mu=theta[["repR"]]*x,size=1/theta[["repvol"]])}),rep(-1,exclude.p-1),
      sapply(mu01[(theta[["rep_drop"]]+exclude.p):length(mu01)] ,function(x){rnbinom(1, mu=theta[["repRA"]]*x,size=1/theta[["repvolA"]])})
    )
  }

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

LikelihoodFunction<-function(y, c1, c2, theta,iiN,repSS=NULL){
  muAll=sapply(c1,function(x){max(0,x)}) + sapply(c2,function(x){max(0,x)})

  # Add together early and late fitting:
  if(cutt.date<=swap.date){
    sum( dnbinom(y, mu=theta[["repR"]]*muAll,size=1/theta[["repvol"]],log=T) )
  }else{
    sum( dnbinom(y[1:(theta[["rep_drop"]])], mu=theta[["repR"]]*muAll[1:(theta[["rep_drop"]])],size=1/theta[["repvol"]],log=T) ) + 
    sum( dnbinom(y[(theta[["rep_drop"]]+exclude.p):length(y)], mu=theta[["repRA"]]*muAll[(theta[["rep_drop"]]+exclude.p):length(y)],size=1/theta[["repvolA"]],log=T) )
    
  }

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


SampleTheta<-function(theta_in, theta_init_in,m,covartheta,covartheta_init,singleI=NULL){
  
  # sample new parameters from nearby: 
  mean_vector_theta = theta_in
  mean_vector_theta0=mean_vector_theta
  theta_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta0), covartheta)))
  names(theta_star)=names(theta_in)
  
  if(sum(names(theta_star)=="repR")>0){ # check theta contains this vector
    theta_star[["repR"]]=min(theta_star[["repR"]],2-theta_star[["repR"]]) # Ensure reporting between zero and 1
  }
  
  if(sum(names(theta_star)=="beta_v_amp")>0){
    theta_star[["beta_v_amp"]]=min(theta_star[["beta_v_amp"]],2-theta_star[["beta_v_amp"]]) # Ensure amplitude between zero and 1
  }
  
  if(sum(names(theta_star)=="beta_c_base")>0){
    theta_star[["beta_c_base"]]=min(theta_star[["beta_c_base"]],2-theta_star[["beta_c_base"]]) # Ensure amplitude between zero and 1
    theta_star[["beta_c_base2"]]=min(theta_star[["beta_c_base2"]],2-theta_star[["beta_c_base2"]]) # Ensure amplitude between zero and 1
  }
  
  if(sum(names(theta_star)=="beta2")>0){
    theta_star[["beta2"]]=min(theta_star[["beta2"]],2-theta_star[["beta2"]]) # Ensure beta is between zero and 1
    theta_star[["beta3"]]=min(theta_star[["beta3"]],2-theta_star[["beta3"]]) # Ensure beta is between zero and 1
  }

  if(!is.null(singleI)){
    mean_vector_theta_init = theta_init_in
    
    theta_init_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta_init), covartheta_init)))
    names(theta_init_star)=names(theta_init_in)
    
    theta_init_star[["im_initC"]] = min(theta_init_star[["im_initC"]],1-theta_init_star[["im_initC"]])
    theta_init_star[["em_initC"]] = theta_init_star[["im_initC"]] 
    theta_init_star[["sm_initC"]] = 1-theta_init_star[["im_initC"]]-theta_init_star[["em_initC"]]
    #theta_init_star[["im_initA"]] = theta_init_star[["im_initC"]] ; theta_init_star[["em_initA"]] = theta_init_star[["im_initA"]] # FIX SAME INITIAL CONDITIONS IN BOTH  C and A
    theta_init_star[["im_initA"]] = min(theta_init_star[["im_initA"]],1-theta_init_star[["im_initA"]])
    theta_init_star[["em_initA"]] = theta_init_star[["im_initA"]]
    theta_init_star[["sm_initA"]] = 1-theta_init_star[["im_initA"]]-theta_init_star[["em_initA"]]

    if(singleI==1){ # Resample recovereds
      theta_init_star[["r_initC"]]= min(theta_init_star[["r_initC"]], 2*theta_in[["npopC"]] - theta_init_star[["r_initC"]]) # Bound at population size
      theta_init_star[["r_initA"]]= min(theta_init_star[["r_initA"]], 2*theta_in[["npopA"]] - theta_init_star[["r_initA"]]) # Bound at population size
    }
    theta_init_star[["e_initC"]] = theta_init_star[["i1_initC"]]

    theta_init_star[["e_initA"]] = theta_init_star[["i1_initA"]]
    
    theta_init_star[["s_initC"]] = (theta_in[["npopC"]]-theta_init_star[["i1_initC"]]-theta_init_star[["e_initC"]]-theta_init_star[["r_initC"]])
    theta_init_star[["s_initA"]] = (theta_in[["npopA"]]-theta_init_star[["i1_initA"]]-theta_init_star[["e_initA"]]-theta_init_star[["r_initA"]])

    
    theta_init_star=theta_init_star
  }else{
    theta_init_star=theta_init_in
  }
  return(list(thetaS=theta_star,theta_initS=theta_init_star))
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

ComputeProbability<-function(sim_likelihood,sim_likelihood_star,thetatab,theta_star,thetaItab,theta_Istar,thetaAllstar,itertab,c_case_ratioStar,c_case_ratioTab){
  
  # sim_liktab[m],sim_marg_lik_star,thetatab=thetatab[m,]; 
  
  # Include priors - Note have prior on Amplitude now as well
  p_theta_star = priorInf(1/theta_star[["r_inf"]])*priorExp(1/theta_star[["r_exp"]])*priorVEx(1/theta_star[["v_exp"]])*priorMuV(1/theta_star[["mu_v"]])*priorAmplitude(theta_star[["beta_v_amp"]])
  p_theta = priorInf(1/thetatab[["r_inf"]])*priorExp(1/thetatab[["r_exp"]])*priorVEx(1/thetatab[["v_exp"]])*priorMuV(1/thetatab[["mu_v"]])*priorAmplitude(thetatab[["beta_v_amp"]])

  # Calculate probability of correct children/adult reporting
  p_ratio_star = 1 
  p_ratio_base = 1 

  # P(theta | theta_star)
  q_theta_given_theta_star = 1
  q_theta_star_given_theta = 1
  
  val = exp((sim_likelihood_star-sim_likelihood))*(p_ratio_star/p_ratio_base)*(p_theta_star/p_theta)*(q_theta_given_theta_star/q_theta_star_given_theta) 
  min(val, 1)
  
  
}


sigmf <- function(x,x0,k){exp(-(x-x0)/k)/(1+exp(-(x-x0)/k))^2}

sigmd1 <- function(x,x0){as.numeric(x>x0)} #as.numeric(x>x0)

# Seasonality function - note that day 0 is assumed to be as.Date("2013-10-28")

seasonal_f <- function(x,date0,theta){
  yy = (1 + theta[["beta_v_mask"]]*theta[["beta_v_amp"]]*sin(((x+date0)/365- season.shift )*2*pi)) # Pin to Feb max -- CALC as.Date("2013-11-04")+0.048*2*pi*365
  yy
}


decline_f <- function(x,date0,theta){

  pp = control.shift + control.range /( 1 + exp(-10*theta[["beta_c_mid"]])) 
  yy = 1 - theta[["beta_c_mask"]]*theta[["beta_c_base"]]/(1+exp(-10*theta[["beta_c_grad"]]*((x+date0)/365-pp)))  
  yy
}



Simulate_model2<-function(NN,dt=0.1, theta, theta_init, y.vals,time.vals,repTab,date_list,plotA=FALSE,simzetaA=NULL,fitplot=FALSE,locationI){

  output <- Deterministic_modelR(1,dt,theta, theta_init, y.vals,time.vals,repTab,locationI)

  if(plotA==T){
    case_count=output$C_trace
    
    # Reporting model
    xmax = max(min(date_list)+time.vals) 
    #par(mfrow=c(2,1))
    par(mar = c(5,5,1,3),mfrow=c(2,1))
    mu1=case_count
    plot(date_list,y.vals,ylim=c(0,1.2*max(y.vals,mu1*theta[["repR"]])),xlim=c(min(date_list),xmax))
    for(kk in 1:NN){
      lines(min(date_list)+time.vals-min(time.vals),ReportC(mu1,theta,repTab),col=rgb(0,0,1,0.5))

    }

    par(new=TRUE)

    xxmax = max(time.vals)  # Adjust because seasonality starts at 0, observations start at 1
    xxmin = min(time.vals)
    xxlength = xxmax - xxmin +1

    
    plot(min(date_list)+c(xxmin:xxmax) -xxmin   ,seasonal_f(c(xxmin:xxmax),date0 = theta[["shift_date"]],theta),type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1.5),xlim=c(min(date_list),xmax))
    lines(min(date_list)+c(xxmin:xxmax) -xxmin  ,sapply(c(xxmin:xxmax),function(x){decline_f(x,date0 = theta[["shift_date"]],theta)}) ,type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="")

    axis(4,col="red",col.axis="red")
    mtext("Transmission rate", side=4, line=2,col="red",cex=1) # Label for 2nd axis
    
    # Plot susceptibles etc.
    plot(min(date_list)+time.vals-7,output$X_traceC,ylim=c(0,1),xlim=c(min(date_list),xmax),type="l")
    lines(min(date_list)+time.vals-7,output$S_traceC/theta[["npopC"]],ylim=c(0,1),xlim=c(min(date_list),xmax),type="l",col="red")
    lines(min(date_list)+time.vals-7,output$S_traceC/theta[["npopA"]],ylim=c(0,1),xlim=c(min(date_list),xmax),type="l",col="blue")
    
    
    dev.copy(pdf,paste("plotsD/Simulate0.pdf",sep=""),width=10,height=6)
    dev.off()
    
    
  }
  
  if(fitplot==T){
    mu1=output$C_trace
    list(report_traj=ReportC(mu1,theta,repTab))
  }

  
}
  


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Output vector-borne simulated deterministic likelihood - THIS IS USED IN MODEL

# Simulation model ODE
simulate_deterministic <- function(theta, init.state, times) {
  SIR_ode <- function(time, state, theta) {
    
    ## extract parameters
    beta_h1 <- theta[["beta"]] *  seasonal_f(time,date0=theta[["shift_date"]],theta) * decline_f(time,date0=theta[["shift_date"]],theta) # Depends on load_timeseries_data_values
    beta_h3 <-  beta_h1 # *  seasonal_f(time,date0=0,theta) #theta[["beta3"]]
    beta_h2 <- beta_h1 # Reduced transmission between groups
    beta_v1 <-  theta[["beta_v"]] * beta_h1 # Scale with human contact
    beta_v2 <-  beta_v1 # Scale with human contact
    beta_v3 <-  beta_v1 # Scale with human contact
    NsizeC <- theta[["npopC"]]
    NsizeA <- theta[["npopA"]]
    delta_v  <- theta[["mu_v"]]
    alpha_v <- theta[["v_exp"]]
    alpha_h <- theta[["r_exp"]]
    gamma <- theta[["r_inf"]]
    
    ## states
    SC <- state[["s_initC"]]
    EC <- state[["e_initC"]]
    IC <- state[["i_initC"]]
    RC <- state[["r_initC"]]
    CC <- state[["c_initC"]] 
    SMC <- state[["sm_initC"]]
    EMC <- state[["em_initC"]]
    IMC <- state[["im_initC"]]
    
    SA <- state[["s_initA"]]
    EA <- state[["e_initA"]]
    IA <- state[["i_initA"]]
    RA <- state[["r_initA"]]
    CA <- state[["c_initA"]] 
    #SMA <- state[["sm_initA"]]
    #EMA <- state[["em_initA"]]
    #IMA <- state[["im_initA"]]
    
    # Allow potential for extinction if not enough infected humans

    ICpos = sigmd1(IC,1) # Need at least one infective
    IApos = sigmd1(IA,1) # Need at least one infective
    Ipos = sigmd1(IC+IA,1) # Need at least one infective
    
    # Human populations-- FORMULATION WITH ONE MOSQUITO POP AND DIFFERENT BITING RATES

    dSC  =  - ICpos*SC*(beta_h1*IMC) # EDITED ** HOMOGENEOUS MOSQUITOES
    dEC  =  ICpos*SC*(beta_h1*IMC) - alpha_h*EC  # EDITED ** HOMOGENEOUS MOSQUITOES
    dIC  = alpha_h*EC  - gamma*IC
    dRC  = gamma*IC
    dCC  = alpha_h*EC  # - assume reported when infectious
    
    dSA  = - IApos*SA*(beta_h3*IMC)      # EDITED ** HOMOGENEOUS MOSQUITOES
    dEA  =  IApos*SA*(beta_h3*IMC) - alpha_h*EA       # EDITED ** HOMOGENEOUS MOSQUITOES
    dIA  = alpha_h*EA  - gamma*IA
    dRA  = gamma*IA
    dCA  = alpha_h*EA
    
    # Mosquito populations -- ONLY USE ONE MOSQUITO POPULATION

    dSMC = delta_v - Ipos*SMC*(beta_v1*IC + beta_v3*IA)/(NsizeC + NsizeA) - delta_v*SMC
    dEMC = Ipos*SMC*(beta_v1*IC + beta_v3*IA)/(NsizeC + NsizeA) - (delta_v+alpha_v)*EMC 
    dIMC = alpha_v*EMC - delta_v*IMC

    return(list(c(dSC,dEC,dIC,dRC,dCC,dSMC,dEMC,dIMC,
                  dSA,dEA,dIA,dRA,dCA))) #,dSMA,dEMA,dIMA
  }
  
  # put incidence at 0 in init.state
  traj <- as.data.frame(ode(init.state, times, SIR_ode, theta, method = "ode45"))

  return(traj)
}

Deterministic_modelR<-function(iiN,dt,theta, theta_init, y.vals,time.vals,repTN,locationI){
  
  sim.vals = seq(0,max(time.vals)-min(time.vals),7) + 7 # These values tell how to match with data points

  init1=c(
    s_initC=theta_init[["s_initC"]],e_initC=theta_init[["i1_initC"]],i_initC=theta_init[["i1_initC"]],r_initC=theta_init[["r_initC"]],c_initC=0,
    sm_initC=theta_init[["sm_initC"]],em_initC=theta_init[["em_initC"]],im_initC=theta_init[["im_initC"]],
    s_initA=theta_init[["s_initA"]],e_initA=theta_init[["i1_initA"]],i_initA=theta_init[["i1_initA"]],r_initA=theta_init[["r_initA"]],c_initA=0)

  # Output simulation data
  output <- simulate_deterministic(theta,init1,seq(0,max(sim.vals),dt) )
  
  S_trajC <- output[match(sim.vals,output$time),"s_initC"]
  S_trajA <- output[match(sim.vals,output$time),"s_initA"]
  X_trajC <- output[match(sim.vals,output$time),"sm_initC"]
  X_trajA <- X_trajC #output[match(sim.vals,output$time),"sm_initA"] # NOTE THIS IS DEPRECATED
  R_trajC <- output[match(sim.vals,output$time),"r_initC"]
  R_trajA <- output[match(sim.vals,output$time),"r_initA"]
  cases1 <- output[match(sim.vals,output$time),"c_initC"]
  cases2 <- output[match(sim.vals,output$time),"c_initA"]
  casecountC <- cases1-c(0,cases1[1:(length(sim.vals)-1)])
  casecountA <- cases2-c(0,cases2[1:(length(sim.vals)-1)])
  
  casecount <- casecountC + casecountA
  seroPC_1 <- theta_init[["r_initC"]]/theta[["npopC"]]
  seroPA_1 <- theta_init[["r_initA"]]/theta[["npopA"]]
  seroPC_2 <- tail(R_trajC,1)/theta[["npopC"]]
  seroPA_2 <- tail(R_trajA,1)/theta[["npopA"]]
  
  # Calculate effective R and constrain to be <1
  likesero1 = (ifelse(locationI=="Central2014",dbinom(n_Luminex_C_D3[1], size=n_Luminex_C_D3[3], prob=seroPC_1, log = T),0) +
                 ifelse(locationI=="Central2014",dbinom(n_Luminex_A_D3[1], size=n_Luminex_A_D3[3], prob=seroPA_1, log = T),0))*theta[["sero_lik1"]]
  liksero2 = (ifelse(locationI=="Central2014",dbinom(n_Luminex_C_D3[2], size=n_Luminex_C_D3[3], prob=seroPC_2, log = T),0) +
               ifelse(locationI=="Central2014",dbinom(n_Luminex_A_D3[2], size=n_Luminex_A_D3[3], prob=seroPA_2, log = T),0) )*theta[["sero_lik2"]]

  likcases = LikelihoodFunction(y.vals,casecountC[1:length(y.vals)],casecountA[1:length(y.vals)], theta,1,repTN)
  
  # Check if there is a second epidemic (i.e. at least one reported case). 60 is cutoff for 2015
  likcases2015 = -1e10*((theta[["repR"]] * max(casecountC[60:length(casecountC)] + casecountA[60:length(casecountC)]) > 1) ) #& min(casecountC+casecountC)>1 ) # Condition on persistence
  
  likelihood = likcases + likesero1 + liksero2 + theta[["no2015"]]*likcases2015
  
  # Avoid infinity in MCMC sum over regions
  if(likelihood == -Inf | is.na(likelihood)){likelihood=-1e10}
  
  return(list(C_trace=casecount,C_traceC=casecountC,C_traceA=casecountA,S_traceC=S_trajC,S_traceA=S_trajA,R_traceC=R_trajC,R_traceA=R_trajA,X_traceC=X_trajC,X_traceA=X_trajA,lik=likelihood)) #,
  
}
