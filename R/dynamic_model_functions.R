# - - - - - - - - - - - - - - - - - - - - - - - - 
# Main model functions
# Github version
#
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - - - 


# - - - - - -
# Start functions adapted from Mordecai et al (2017) PLOS NTD. DOI: journal.pntd.0005568

briere <- function(t, c, Tm, T0){
  b = (c*t*(t-T0)*sqrt(Tm-t)) 
  b[!(t>T0 & t<Tm)] = 0 
  
  b
}

briere.trunc <- function(t, c, Tm, T0){
  b = (c*t*(t-T0)*sqrt(Tm-t)) 
  b[!(t>T0 & t<Tm)] = 0 
  #b[b>1]=1
  b
}

quad.2 <- function(t, T0, Tm, qd){
  b = qd*(t-T0)*(t-Tm)
  b[!(t>T0 & t<Tm)] = 0 
  b
}

quad.2.trunc <- function(t, T0, Tm, qd, lim=0.0001){
  b = qd*(t-T0)*(t-Tm)
  b[!(t>T0 & t<Tm)] = 0 
  #b[b>1]=1
  b
}

# End functions adapted from Mordecai et al (2017) PLOS NTD. DOI: journal.pntd.0005568
# - - - - - -


# - - - - - -
# R0 matrix calculation
calculate_r0 <- function(th_in,sus_c=1,sm_c=1,t_vary=1,controlT=1){
  
  # DEBUG   sus_c=s_pickH[cosPick]; sm_c=x_pickC[cosPick]; th_in = thetatab[b_ii,]; t_vary=time.V

  rain_scale = th_in$beta_v_mask*carrying_f(t_vary,0,th_in) + (1-th_in$beta_v_mask)
  
  
  temp_t =  th_in$beta_v_mask*seasonal_f(t_vary,th_in) + (1-th_in$beta_v_mask)*26
  
  contact_rate = th_in$beta_v * bite_temp(temp_t) * (controlT*decline_f(t_vary,date0=th_in$shift_date,th_in) +(1-controlT))
  
  delta_v = mortality_rate_temp(temp_t) *th_in$mu_v
  exp_v = EI_rate_temp(temp_t)*th_in$v_exp
    
  # par(mfrow=c(2,1)); plot(temp_t); plot(contact_rate)
  
  # Rate humans get infected -- FORMULATION WITH ONE MOSQUITO POP AND DIFFERENT BITING RATES
  b_hv =  prob_to_h_temp(temp_t) * th_in$m_density * contact_rate * density_vary(temp_t) * rain_scale
  b_hh = 0
  
  # Rate vectors get infected
  b_vh = prob_to_v_temp(temp_t) * contact_rate  
  b_vv = 0

  rr_hh=rep(0,length(t_vary)); rr_vv = rr_hh;
  
  rr_hv = (sus_c*b_hv/delta_v)*(exp_v/(exp_v+delta_v)) 
  rr_vc = (sm_c*b_vh/th_in$r_inf)

  rr_post = NULL
  for(ii in 1:length(t_vary)){
    rr_post=c(rr_post,( max(Re(eigen(matrix(c(rr_hh[ii],rr_hv[ii],rr_vc[ii],rr_vv[ii]),nrow=2))$values)))  )
  }
  
  return( list(rr_out=rr_post,rr_mat = matrix(c(rr_hh[1],rr_hv[1],rr_vc[1],rr_vv[1]),nrow=2,byrow=T)) )
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Case reporting

ReportC<-function(ct, theta,repSS,yprop){
  
  #DEBUG:  ct = mu1; repSS = "susp"; yprop = y.vals.prop
  
  mu01=sapply(ct,function(x){max(x,0)})
  
  if(repSS=="lab"){ mu_output = sapply(yprop*mu01,function(x){rnbinom(1, mu=theta[["repR"]]*x,size=1/theta[["repvol"]] )}) }
  if(repSS=="sus"){ mu_output = sapply((1-yprop)*mu01,function(x){rnbinom(1, mu=theta[["repRA"]]*x,size=1/theta[["repvol"]] )}) }
  
  # if(cutt.date<=swap.date){
  # sapply(mu01,function(x){rnbinom(1, mu=theta[["repR"]]*x,size=1/theta[["repvol"]])})
  # }else{
  #   c(sapply(mu01[1:(theta[["rep_drop"]])] ,function(x){rnbinom(1, mu=theta[["repR"]]*x,size=1/theta[["repvol"]])}),rep(-1,exclude.p-1),
  #     sapply(mu01[(theta[["rep_drop"]]+exclude.p):length(mu01)] ,function(x){rnbinom(1, mu=theta[["repRA"]]*x,size=1/theta[["repvolA"]])})
  #   )
  # }
  mu_output

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Surveillance data likelihood

LikelihoodFunction<-function(y, c1, c2, theta,iiN,repSS=NULL,yprop){
  
  muAll=sapply(c1,function(x){max(0,x)}) + sapply(c2,function(x){max(0,x)})

  # # Add together early and late fitting:
  # if(cutt.date<=swap.date){
  #   sum( dnbinom(y, mu=theta[["repR"]]*muAll,size=1/theta[["repvol"]],log=T) )
  # }else{
  #   sum( dnbinom(y[1:(theta[["rep_drop"]])], mu=theta[["repR"]]*muAll[1:(theta[["rep_drop"]])],size=1/theta[["repvol"]],log=T) ) + 
  #   sum( dnbinom(y[(theta[["rep_drop"]]+exclude.p):length(y)], mu=theta[["repRA"]]*muAll[(theta[["rep_drop"]]+exclude.p):length(y)],size=1/theta[["repvolA"]],log=T) )
  #   
  # }

  if(repSS=="lab"){ lik_output = ( dnbinom(y, mu=theta[["repR"]]*muAll*yprop,size=1/theta[["repvol"]],log=T) )}
  if(repSS=="sus"){ lik_output = ( dnbinom(y, mu=theta[["repRA"]]*muAll*(1-yprop),size=1/theta[["repvol"]],log=T) )}

  lik_output

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Resample theta

SampleTheta<-function(theta_in, theta_init_in,m,covartheta,covartheta_init,singleI=NULL,pmask){
  
  #DEBUG: theta_in = thetaAlltab[m,iiH,]; theta_init_in = theta_initAlltab[m,iiH,]; covartheta=0*cov_matrix_thetaA; covartheta_init = 0*cov_matrix_theta_init; singleI=kk
  #DEBUG: theta_in = thetatab[m,]; theta_init_in = theta_initAlltab[m,iiH,]; covartheta=0*cov_matrix_thetaA; covartheta_init = 0*cov_matrix_theta_init; singleI=kk
  
  # sample new parameters from nearby: 
  mean_vector_theta = theta_in

  if(sum(names(mean_vector_theta)=="shift_date")>0 & !is.null(singleI)){
    mean_vector_theta[["shift_date"]]=10 # To avoid errors
  }
  
  mean_vector_theta0 = mean_vector_theta 
  theta_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta0), covartheta)))
  names(theta_star)=names(theta_in)
  
  if(sum(names(theta_star)=="repR")>0){ # check theta contains this vector
    theta_star[["repR"]]=min(theta_star[["repR"]],2-theta_star[["repR"]]) # Ensure reporting between zero and 1
  }
  
  
  if(sum(names(theta_star)=="prop_at_risk")>0){ # check theta contains this vector
    theta_star[["prop_at_risk"]]=min(theta_star[["prop_at_risk"]],2-theta_star[["prop_at_risk"]]) # Ensure at risk group between zero and 1
  }
  
  if(sum(names(theta_star)=="beta_v_amp")>0){
    theta_star[["beta_v_amp"]]=min(theta_star[["beta_v_amp"]],200-theta_star[["beta_v_amp"]]) # Ensure rainfall effect between zero and 100
  }
  
  if(sum(names(theta_star)=="beta_c_grad")>0){
    theta_star[["beta_c_grad"]]=min(theta_star[["beta_c_grad"]],2000-theta_star[["beta_c_grad"]]) # Ensure control gradient between zero and 1000
  }
  
  if(sum(names(theta_star)=="beta_c_base")>0){
    theta_star[["beta_c_base"]]=min(theta_star[["beta_c_base"]],2-theta_star[["beta_c_base"]]) # Ensure amplitude between zero and 1
    theta_star[["beta_c_base2"]]=min(theta_star[["beta_c_base2"]],2-theta_star[["beta_c_base2"]]) # Ensure amplitude between zero and 1
  }

  if(!is.null(singleI)){
    theta_star[pmask] = theta_in[pmask] # Replace masked values (some may be negative)
  }
  
  if(!is.null(singleI)){
    mean_vector_theta_init = theta_init_in
    
    theta_init_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta_init), covartheta_init)))
    names(theta_init_star)=names(theta_init_in)
    
    theta_init_star[["im_initC"]] = min(theta_init_star[["im_initC"]],theta_star[["npopM"]] -theta_init_star[["im_initC"]])
    theta_init_star[["em_initC"]] = theta_init_star[["im_initC"]] 
    theta_init_star[["sm_initC"]] = theta_star[["npopM"]] - theta_init_star[["im_initC"]]-theta_init_star[["em_initC"]]
    #theta_init_star[["im_initA"]] = theta_init_star[["im_initC"]] ; theta_init_star[["em_initA"]] = theta_init_star[["im_initA"]] # FIX SAME INITIAL CONDITIONS IN BOTH  C and A
    theta_init_star[["im_initA"]] = min(theta_init_star[["im_initA"]],theta_star[["npopM"]] -theta_init_star[["im_initA"]])
    theta_init_star[["em_initA"]] = theta_init_star[["im_initA"]]
    #theta_init_star[["sm_initA"]] = theta_star[["npopM"]] -theta_init_star[["im_initA"]]-theta_init_star[["em_initA"]]

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

ComputeProbability<-function(sim_likelihood,sim_likelihood_star,thetatab,theta_star,thetaItab,theta_Istar,thetaAllstar,thetaAlltab,itertab,c_case_ratioStar,c_case_ratioTab){
  
  # sim_likelihood=sim_liktab[m]; sim_likelihood_star=sim_marg_lik_star; thetatab=thetatab[m,]; 
  
  # Include priors - Note have prior on Amplitude now as well
  p_theta_star = priorInf(1/theta_star[["r_inf"]])*priorExp(1/theta_star[["r_exp"]])*priorVEx(1/theta_star[["v_exp"]])*priorMuV(1/theta_star[["mu_v"]])*priorBeta_v(theta_star[["beta_v"]]) *priorDensity(theta_star[["m_density"]]) 
  p_theta = priorInf(1/thetatab[["r_inf"]])*priorExp(1/thetatab[["r_exp"]])*priorVEx(1/thetatab[["v_exp"]])*priorMuV(1/thetatab[["mu_v"]])*priorBeta_v(thetatab[["beta_v"]]) * priorDensity(thetatab[["m_density"]]) 

  # Calculate acceptance probability
  val = exp((sim_likelihood_star-sim_likelihood))*(p_theta_star/p_theta)
  min(val, 1)
  
  
}


sigmf <- function(x,x0,k){exp(-(x-x0)/k)/(1+exp(-(x-x0)/k))^2}

sigmd1 <- function(x,x0){as.numeric(x>x0)} #as.numeric(x>x0)

# Seasonality function
seasonal_f <- function(x,theta){
  #yy = (1 + theta[["beta_v_mask"]]*theta[["beta_v_amp"]]*sin(((x+date0)/365- season.shift )*2*pi)) 
  #yy = (1-theta[["temp_balance"]])*seasonaltemp(x,theta_fit1) + theta[["temp_balance"]]*seasonaltemp(x,theta_fit2)
  yy = seasonaltemp(x,theta_fit1) 
  yy
}

# Carrying capacity:
carrying_f <- function(x,date0,theta){
  #yy = (1 + theta[["beta_v_mask"]]*theta[["beta_v_amp"]]*sin(((x+date0)/365- season.shift )*2*pi)) 
  #yy = theta[["beta_v_mid"]]*(1 + theta[["beta_v_amp"]]*(seasonalrain(x,theta_fitRain)-84)/(220-84))  # Scale between 0 and 1
  #yy = theta[["beta_v_mid"]]*((seasonalrain(x,theta_fitRain)-84)/(220-84))  # Scale between 0 and 1
  #yy = yy[yy>0]
  #yy = 1+theta[["beta_v_amp"]]*(seasonalrain(x,theta_fitRain)-221.2756) # Scale to max value
  #yy = (1 + theta[["beta_v_amp"]]*(seasonalrain(x,theta_fitRain)-221.2756)/136.7567) # Scale to oscillate around 1
  
  #theta[["beta2"]]=1; theta[["beta_v_amp"]]=1;
  
  yy = 1/(1+1/(0 + theta[["beta_v_amp"]]*seasonalrain(x,theta_fitRain)/222.4374)) # Scale to oscillate around 1 # theta[["beta2"]]
  
  yy = yy[yy>0]
  yy
}


decline_f <- function(x,date0,theta){
  pp = control.shift + control.range /( 1 + exp(-10*(theta[["beta_c_mid"]] -1)))  # made sure drop after period control starts
  
  yy = 1 - theta[["beta_c_mask"]]*theta[["beta_c_base"]]/(1+exp(-theta[["beta_c_grad"]]*((x+date0)/365-pp )))  
  yy
}

#xx1=0
#start.date + 365* (control.shift + control.range /( 1 + exp(-10*(xx1-1))) )  # made sure drop after period control starts

Simulate_model2<-function(NN,dt=0.1, theta, theta_init, y.vals,y.vals2, y.vals.prop,time.vals,repTab,date_list,plotA=FALSE,simzetaA=NULL,fitplot=FALSE,locationI){

  # DEBUG
  # theta=c(theta,thetaAll[iiH,]);  theta_init=theta_initAll[iiH,]; repTab=repTN; locationI=locationtab[iiH] ; sero_lik=1; NN=5
  
  
  output <- Deterministic_modelR(1,dt,theta, theta_init, y.vals,y.vals2, y.vals.prop, time.vals,repTab,locationI)

  print( output$lik)
  
  if(plotA==T){
    
    case_count=output$C_trace
    
    # Reporting model
    xmax = max(min(date_list)+time.vals) 
    #par(mfrow=c(2,1))
    par(mar = c(4,4,1,3),mfrow=c(2,1))
    mu1=case_count
    case_actual = (y.vals+ y.vals2)
    plot(date_list,case_actual,xlim=c(min(date_list),xmax),pch=19,ylab="cases",ylim=c(0,2000)) #,ylim=c(0,2000) #
    #points(date_list,y.vals2,ylim=c(0,1000),xlim=c(min(date_list),xmax))
    
    for(kk in 1:NN){
      #lines(min(date_list)+time.vals-min(time.vals),ReportC(mu1,theta,repSS = "lab",yprop = y.vals.prop),col=rgb(0,0,1,0.5))
      case_out =( ReportC(mu1,theta,repSS = "lab",yprop = y.vals.prop) + ReportC(mu1,theta,repSS = "sus",yprop = y.vals.prop) ) #%>% log()
      case_all =( theta[["repR"]]*mu1 )#%>% log()
      lines(min(date_list)+time.vals-min(time.vals),case_out,col=rgb(0,0.5,1,0.5))
      lines(min(date_list)+time.vals-min(time.vals),case_all,col=rgb(0,0,1,0.5),lty=2)
      
    }
    polygon(c(as.Date("2014-03-08"),as.Date("2014-03-08"),as.Date("2014-03-22"),as.Date("2014-03-22")),c(-1,1e4,1e4,-1),col=rgb(1,0,0,0.3),lty=0)
    
    #lines(c(as.Date("2014-02-15"),as.Date("2014-02-15")),c(-100,100),col=rgb(1,0,0),lty=1)
    

    par(new=TRUE)

    xxmax = max(time.vals)  # Adjust because seasonality starts at 0, observations start at 1
    xxmin = min(time.vals)
    xxlength = xxmax - xxmin +1

    plot(min(date_list)+c(xxmin:xxmax) -xxmin   ,seasonal_f(c(xxmin:xxmax),theta),type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(20,30),xlim=c(min(date_list),xmax))
    lines(min(date_list)+c(xxmin:xxmax) -xxmin  ,sapply(c(xxmin:xxmax),function(x){20+5*decline_f(x,date0 = theta[["shift_date"]],theta)}) ,type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="")

    axis(4,col="red",col.axis="red")
    mtext("Transmission rate", side=4, line=2,col="red",cex=1) # Label for 2nd axis
    
    # Plot susceptibles etc.
    plot(min(date_list)+time.vals-7,output$X_traceA/theta[["npopM"]],xlim=c(min(date_list),xmax),ylim=c(0,2),type="l") #ylim=c(0,3)
    lines(min(date_list)+time.vals-7,output$X_traceC/theta[["npopM"]],ylim=c(0,1),xlim=c(min(date_list),xmax),type="l",col="black",lty=2)
    lines(min(date_list)+time.vals-7,theta[["prop_at_risk"]]*(1-output$S_traceC/theta[["npopC"]]),ylim=c(0,1),xlim=c(min(date_list),xmax),type="l",col="red")
    lines(min(date_list)+time.vals-7,theta[["prop_at_risk"]]*(1-output$S_traceA/theta[["npopA"]]),ylim=c(0,1),xlim=c(min(date_list),xmax),type="l",col="red",lty=2)
    
    
    
    rr_vals = calculate_r0(theta,sus_c=(output$S_traceC+output$S_traceA)/(theta[["npopC"]]+theta[["npopA"]]),sm_c=output$X_traceC,t_vary=time.vals,controlT=1)$rr_out
    
    lines(min(date_list)+time.vals-7,rr_vals,ylim=c(0,1),xlim=c(min(date_list),xmax),type="l",col="blue")
    
    
    dev.copy(pdf,paste("plots/Simulate0.pdf",sep=""),width=10,height=6)
    dev.off()
    
    
  }
  
  if(fitplot==T){
    #mu1=output$C_trace
    #list(report_traj=ReportC(mu1,theta,repTab))
  }

  
}
  


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Output vector-borne simulated deterministic likelihood - THIS IS USED IN MODEL

# Simulation model ODE
simulate_deterministic <- function(theta, init.state, times) {
  SIR_ode <- function(time, state, theta) {

    # Define temperature and variable carrying capacity
    temp_t <-  theta[["beta_v_mask"]]*seasonal_f(time,theta) + (1-theta[["beta_v_mask"]])*25
    density_scale <-  theta[["beta_v_mask"]]*carrying_f(time,0,theta) + (1-theta[["beta_v_mask"]])

    # Define vector-specific parameters
    contact_rate <- theta[["beta_v"]] * bite_temp(temp_t)  #
    mosquito_to_human <- prob_to_h_temp(temp_t) * contact_rate * theta[["m_density"]] * density_vary(temp_t) * density_scale * decline_f(time,date0=theta[["shift_date"]],theta)     # scale by density * rain_scale
    human_to_mosquito <- prob_to_v_temp(temp_t) * contact_rate 
    delta_v  <- mortality_rate_temp(temp_t) *theta[["mu_v"]]
    recruit_v <-  delta_v 
    alpha_v <- EI_rate_temp(temp_t) *theta[["v_exp"]] # EIP

    # Define human parameters
    beta_h1 <- mosquito_to_human
    beta_v1 <-  human_to_mosquito
    NsizeC <- theta[["npopC"]]
    NsizeA <- theta[["npopA"]]
    
    alpha_h <- theta[["r_exp"]] #IIP
    gamma <- theta[["r_inf"]]
    
    # Load initial conditions
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
    
    # Allow potential for extinction if not enough infected humans (deprecated)

    ICpos = 1 #sigmd1(IC,1) # Need at least one infective
    IApos = 1 #sigmd1(IA,1) # Need at least one infective
    Ipos = 1 #sigmd1(IC+IA,1) # Need at least one infective
    
    Nmsize = SMC + EMC + IMC
    
    # ODEs for human populations-- FORMULATION WITH ONE MOSQUITO POP

    dSC  =  - ICpos*SC*(beta_h1*IMC) 
    dEC  =  ICpos*SC*(beta_h1*IMC) - alpha_h*EC  
    dIC  = alpha_h*EC  - gamma*IC
    dRC  = gamma*IC
    dCC  = alpha_h*EC  # - assume reported when infectious
    
    dSA  = - IApos*SA*(beta_h1*IMC)     
    dEA  =  IApos*SA*(beta_h1*IMC) - alpha_h*EA     
    dIA  = alpha_h*EA  - gamma*IA
    dRA  = gamma*IA
    dCA  = alpha_h*EA  # - assume reported when infectious
    
    # ODEs for  mosquito populations -- ONLY USE ONE MOSQUITO POPULATION

    dSMC = recruit_v - Ipos*SMC*(beta_v1*IC + beta_v1*IA)/(NsizeC + NsizeA) - delta_v*SMC  #larvae_dev*SMA 
    dEMC = Ipos*SMC*(beta_v1*IC + beta_v1*IA)/(NsizeC + NsizeA) - (delta_v+alpha_v)*EMC 
    dIMC = alpha_v*EMC - delta_v*IMC

    # Output values
    return(list(c(dSC,dEC,dIC,dRC,dCC,dSMC,dEMC,dIMC,
                  dSA,dEA,dIA,dRA,dCA)))
  }
  
  # Put incidence at 0 in init.state
  traj <- as.data.frame(ode(init.state, times, SIR_ode, theta, method = "ode45"))

  return(traj)
}

Deterministic_modelR<-function(iiN,dt,theta, theta_init, y.vals, y.vals2, y.vals.prop, time.vals,repTN,locationI){
  
  sim.vals = seq(0,max(time.vals)-min(time.vals),7) + 7 # These values tell how to match with data points

  # Set initial conditions
  init1=c(
    s_initC=theta_init[["s_initC"]],e_initC=theta_init[["i1_initC"]],i_initC=theta_init[["i1_initC"]],r_initC=theta_init[["r_initC"]],c_initC=0,
    sm_initC=theta_init[["sm_initC"]],em_initC=theta_init[["em_initC"]],im_initC=theta_init[["im_initC"]],
    s_initA=theta_init[["s_initA"]],e_initA=theta_init[["i1_initA"]],i_initA=theta_init[["i1_initA"]],r_initA=theta_init[["r_initA"]],c_initA=0)

  # Output simulation data
  output <- simulate_deterministic(theta,init1,seq(0,max(sim.vals),dt) )
  
  S_trajC <- output[match(sim.vals,output$time),"s_initC"]
  S_trajA <- output[match(sim.vals,output$time),"s_initA"]
  X_trajC <- output[match(sim.vals,output$time),"sm_initC"]
  X_trajA <- output[match(sim.vals,output$time),"sm_initC"] + output[match(sim.vals,output$time),"em_initC"] + output[match(sim.vals,output$time),"im_initC"]
  R_trajC <- output[match(sim.vals,output$time),"r_initC"]
  R_trajA <- output[match(sim.vals,output$time),"r_initA"]
  cases1 <- output[match(sim.vals,output$time),"c_initC"]
  cases2 <- output[match(sim.vals,output$time),"c_initA"]
  casecountC <- cases1-c(0,cases1[1:(length(sim.vals)-1)])
  casecountA <- cases2-c(0,cases2[1:(length(sim.vals)-1)])
  
  # Compile case counts and serological data
  casecount <- casecountC + casecountA
  seroPC_1 <- (theta_init[["r_initC"]]/theta[["npopC"]]) *theta[["prop_at_risk"]]
  seroPA_1 <- (theta_init[["r_initA"]]/theta[["npopA"]]) *theta[["prop_at_risk"]]
  seroPC_2 <- (tail(R_trajC,1)/theta[["npopC"]]) *theta[["prop_at_risk"]]
  seroPA_2 <- (tail(R_trajA,1)/theta[["npopA"]]) *theta[["prop_at_risk"]]
  
  # Calculate serology and surveillance likelihoods
  liksero1 = (ifelse(locationI=="Central",dbinom(n_Luminex_C_D3[1], size=n_Luminex_C_D3[3], prob=seroPC_1, log = T),0) +
                 ifelse(locationI=="Central",dbinom(n_Luminex_A_D3[1], size=n_Luminex_A_D3[3], prob=seroPA_1, log = T),0))*theta[["sero_lik1"]]
  liksero2 = (ifelse(locationI=="Central",dbinom(n_Luminex_C_D3[2], size=n_Luminex_C_D3[3], prob=seroPC_2, log = T),0) +
               ifelse(locationI=="Central",dbinom(n_Luminex_A_D3[2], size=n_Luminex_A_D3[3], prob=seroPA_2, log = T),0) )*theta[["sero_lik2"]]

  likcasesLab = LikelihoodFunction(y.vals,casecountC[1:length(y.vals)],casecountA[1:length(y.vals)], theta,1,repSS="lab",y.vals.prop[1:length(y.vals)]) 
  likcasesDLI = LikelihoodFunction(y.vals2,casecountC[1:length(y.vals)],casecountA[1:length(y.vals)], theta,1,repSS="sus",y.vals.prop[1:length(y.vals)])
  
  likelihood = sum(likcasesLab) + sum(likcasesDLI) + liksero1 + liksero2 
  
  write.csv(rbind(n_Luminex_C_D3,n_Luminex_A_D3),paste("plots/check_serology",n_Luminex_A_D3[1],".csv",sep=""))
  
  #print(likcasesLab)
  
  # Avoid infinity in MCMC loop
  if(likelihood == -Inf | is.na(likelihood)){likelihood=-1e10}
  
  return(list(C_trace=casecount,C_traceC=casecountC,C_traceA=casecountA,S_traceC=S_trajC,S_traceA=S_trajA,R_traceC=R_trajC,R_traceA=R_trajA,X_traceC=X_trajC,X_traceA=X_trajA,lik=likelihood)) #,
  
}
