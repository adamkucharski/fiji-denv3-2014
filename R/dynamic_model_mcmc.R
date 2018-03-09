# - - - - - - - - - - - - - - - - - - - - - - - 
# Transmission modelling: run MCMC
# Author: Adam Kucharski
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - - 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN MCMC Model

run_transmission_mcmc <- function(MCMC.runs = 10){

  source("R/dynamic_model_characteristics.R",local=F)
  
  #multichain <- c(1:4) # run in parallel
  
  multichain=4
  
  #foreach(iiM=multichain) %dopar% {  # Loop over scenarios with parallel MCMC chains
  for(iiM in multichain){
    # - - - - - - - - - - - 
    # Initialise ICs 
    
    thetaR_IC <- read.csv(paste("data/thetaR_IC_",country.name,"_only_ELISA",use.ELISA.data,".csv",sep=""),stringsAsFactors=FALSE)
    
    # Outbreak-specific parameters
    thetaAll=data.frame(beta=rep(NA,locnn),
                        beta_c_mid=rep(NA,locnn),
                        beta_c_base=rep(NA,locnn),
                        beta_c_mid2=rep(NA,locnn),
                        beta_c_base2=rep(NA,locnn),
                        beta_c_constrain=rep(NA,locnn),
                        sero_lik1=rep(NA,locnn),
                        sero_lik2=rep(NA,locnn),
                        no2015=rep(0,locnn), # default is to allow 2nd wave
                        repR=rep(NA,locnn),
                        repRA=rep(NA,locnn),
                        rep_drop=rep(NA,locnn),
                        repvol=rep(0.05,locnn),
                        repvolA=rep(0.02,locnn),
                        t_start=0,
                        npop=rep(NA,locnn),npopC=rep(NA,locnn),npopA=rep(NA,locnn),
                        shift_date=rep(NA,locnn))
    
    # Global parameters
    theta = c(r_exp=1/prior_p_Exp[1], # latent (h) # from Chan et al at 30C
              r_inf=1/prior_p_Inf[1], # inf (h)
              v_exp=1/prior_p_VEx[1], # latent (v) # from Chan et al at 30C
              mu_v=1/prior_p_MuV[1], # lifespan
              beta2=as.numeric(thetaR_IC[thetaR_IC$param=="beta_h_2",2]), # baseline after control - jointly fitted
              beta3=as.numeric(thetaR_IC[thetaR_IC$param=="beta_h_3",2]), # DEPRECATED
              beta_v=as.numeric(thetaR_IC[thetaR_IC$param=="beta_v",2]), # Relative mixing M-to-H
              b_rate=15.47/(365*1000),
              d_rate=4.93/(365*1000),
              beta_v_amp=as.numeric(thetaR_IC[thetaR_IC$param=="beta_v_amp",2]),
              beta_v_mid=as.numeric(thetaR_IC[thetaR_IC$param=="beta_v_mid",2])
    )
    
    # Initial conditions
    theta_initAll=data.frame(s_initC=rep(NA,locnn),e_initC=rep(NA,locnn),
                             i1_initC=rep(NA,locnn),r_initC=rep(NA,locnn),
                             sm_initC=rep(NA,locnn),em_initC=rep(NA,locnn),
                             im_initC=rep(NA,locnn),
                             s_initA=rep(NA,locnn),e_initA=rep(NA,locnn),
                             i1_initA=rep(NA,locnn),r_initA=rep(NA,locnn),
                             sm_initA=rep(NA,locnn),em_initA=rep(NA,locnn),
                             im_initA=rep(NA,locnn))
    
    # Loop to set up initial conditions (as may be multiple timeseries)
    
    time.store = list()
    
    for(iiH in itertab){
      
      # Load data for specific place + set initial conditions
      source("R/load_timeseries_data.R",local=TRUE)
      repTN=rR.vals/totalSS[iiH]
      
      time.store[[iiH]]=list(t_store = time.vals,y_store = y.vals,y_store2 = y.vals2, y.vals.prop)
      
      c1=(names(thetaR_IC)==locationtab[iiH])
      
      # Scale age groups by total size
      popsize=thetaR_IC[thetaR_IC$param=="npop",c1]
      popsizeC=thetaR_IC[thetaR_IC$param=="npopC",c1]
      popsizeA=thetaR_IC[thetaR_IC$param=="npopA",c1]
      popsizeTot=(popsizeC+popsizeA)
      popsizeC=round(popsize* popsizeC/popsizeTot); popsizeA=round(popsize* popsizeA/popsizeTot)
      
      thetaAll[iiH,"beta"]=as.numeric(thetaR_IC[thetaR_IC$param=="beta_h",c1])
      thetaAll[iiH,"beta_c_mid"]=as.numeric(thetaR_IC[thetaR_IC$param=="beta_control_mid",c1])
      thetaAll[iiH,"beta_c_base"]=as.numeric(thetaR_IC[thetaR_IC$param=="beta_control_base",c1])
      thetaAll[iiH,"beta_c_grad"]=as.numeric(thetaR_IC[thetaR_IC$param=="beta_control_grad",c1]) # Gradient of control function
      thetaAll[iiH,"beta_c_mask"]=as.numeric(thetaR_IC[thetaR_IC$param=="beta_control_mask",c1]) # Mask control measure (i.e. turn off)
      thetaAll[iiH,"beta_c_constrain"]=as.numeric(thetaR_IC[thetaR_IC$param=="beta_control_mask2",c1]) # Use to constrain fitting of midpoint in 2014 or not
      thetaAll[iiH,"beta_v_mask"]=as.numeric(thetaR_IC[thetaR_IC$param=="beta_v_mask",c1]) # Mask sesaonaltiy (i.e. turn off)
      thetaAll[iiH,"sero_lik1"]=1 # fit PRE serology
      thetaAll[iiH,"sero_lik2"]=1 # fit POST serology
      thetaAll[iiH,"npop"]=popsize; thetaAll[iiH,"npopC"]=popsizeC; thetaAll[iiH,"npopA"]=popsizeA
      thetaAll[iiH,"repR"]=as.numeric(thetaR_IC[thetaR_IC$param=="rep",c1])
      thetaAll[iiH,"repRA"]=as.numeric(thetaR_IC[thetaR_IC$param=="repA",c1])
      thetaAll[iiH,"rep_drop"]=rep.drop.date
      thetaAll[iiH,"repvol"]=as.numeric(thetaR_IC[thetaR_IC$param=="repvol",c1])
      thetaAll[iiH,"repvolA"]=as.numeric(thetaR_IC[thetaR_IC$param=="repvolA",c1])
      thetaAll[iiH,"shift_date"]=shift_to_nov
      
      # Set up initial parameters
      
      initial_inf=as.numeric(thetaR_IC[thetaR_IC$param=="inf0",c1]) #
      init_vec=as.numeric(thetaR_IC[thetaR_IC$param=="vec0",c1])/2
      
      theta_initAll[iiH,"r_initC"]=thetaAll[iiH,"npopC"] * (n_Luminex_C_D3[1]/n_Luminex_C_D3[3])
      theta_initAll[iiH,"r_initA"]=thetaAll[iiH,"npopA"] * (n_Luminex_A_D3[1]/n_Luminex_A_D3[3])
      
      theta_initAll[iiH,"e_initC"]=initial_inf; theta_initAll[iiH,"i1_initC"]=initial_inf
      theta_initAll[iiH,"em_initC"]=init_vec; theta_initAll[iiH,"im_initC"]=init_vec
      theta_initAll[iiH,"e_initA"]=initial_inf; theta_initAll[iiH,"i1_initA"]=initial_inf
      theta_initAll[iiH,"em_initA"]=init_vec; theta_initAll[iiH,"im_initA"]=init_vec
      
      theta_initAll[iiH,"s_initC"]=thetaAll[iiH,"npopC"]-theta_initAll[iiH,"i1_initC"]-theta_initAll[iiH,"e_initC"]-theta_initAll[iiH,"r_initC"]
      theta_initAll[iiH,"sm_initC"]=1-theta_initAll[iiH,"em_initC"]-theta_initAll[iiH,"im_initC"]
      
      theta_initAll[iiH,"s_initA"]=thetaAll[iiH,"npopA"]-theta_initAll[iiH,"i1_initA"]-theta_initAll[iiH,"e_initA"]-theta_initAll[iiH,"r_initA"]
      theta_initAll[iiH,"sm_initA"]=1-theta_initAll[iiH,"em_initA"]-theta_initAll[iiH,"im_initA"]
      
    }
    
    
    # - - - - - - - - - - - 
    # Specific different model types
    # - - - - - - - - - - -
    
    # Turn on/off 2014 control and seasonality
    # 1: SIR model cases  2: SIR model serology and cases  3: SIR + climate  4: SIR + climate + control
    thetaAll[1,"beta_c_constrain"]=1
    #if(use.ELISA.data==T){thetaAll[1,"beta"]=1.3*thetaAll[1,"beta"]} # Adjust for higher immunity in ELISA data
    if(iiM==1){ thetaAll[1,"beta_c_mask"]=0 ; thetaAll[1,"beta_v_mask"]=0 ; thetaAll[1,"sero_lik1"] = 0; thetaAll[1,"sero_lik2"] = 0; thetaAll[1,"beta"]=0.2; theta[["beta_v"]]=10 }
    if(iiM==2){ thetaAll[1,"beta_c_mask"]=0 ; thetaAll[1,"beta_v_mask"]=0; thetaAll[1,"beta"]=0.2;  theta[["beta_v"]]=10 }
    if(iiM==3){ thetaAll[1,"beta_c_mask"]=0; theta[["beta_v_amp"]]=0.3 } #; theta[["beta_v_amp"]]=0.9  } #; thetaAll[1,"beta"]=0.05 } #; ; thetaAll[1,"beta_v"]=10 }
    if(iiM==4){theta[["beta_v_amp"]]=0.6 }
    
    # Covariance matrices - Add theta and thetaAll together in MCMC runs
    nparam=length(theta) 
    npc=rep(1,nparam)
    pmask=c("beta3") # ** THIS FIXES UNIVERSAL PARAMETERS ** CAN TURN OFF AGE MIXING HERE
    npc[match(pmask,names(theta))]=0
    cov_matrix_theta0 = diag(npc)
    
    nparamA=length(thetaAll[1,])
    npcA=rep(1,nparamA)
    pmask=c("t_start","npop","npopC","npopA","shift_date","beta_c_mask","beta_c_constrain","beta_v_mask","rep_drop","sero_lik1","sero_lik2") # ** THIS FIXES OUBTREAK-SPECIFIC PARAMETERS **
    npcA[match(pmask,names(thetaAll[1,]))]=0
    cov_matrix_thetaA0 = diag(npcA)
    
    lmv=length(theta_initAll[1,])
    npcov_init=rep(1,lmv)
    pmaskInit=c("s_initC","e_initC","sm_initC","em_initC","s_initA","e_initA","sm_initA","em_initA")
    npcov_init[match(pmaskInit,names(theta_initAll[1,]) )]=0
    cov_matrix_theta_init0 = diag(npcov_init)
    
    
    # Quick simulation to check looks OK
    if(length(multichain)==1){
      par(mfrow=c(1,1),mar=c(4,4,1,1),mgp=c(2,0.7,0))
      repTN=1
      iiH=1; source("R/load_timeseries_data.R",local=TRUE)
      time.valsSim=time.vals
      Simulate_model2(5, dt, c(theta,thetaAll[iiH,]), theta_initAll[iiH,], y.vals,y.vals2, y.vals.prop,time.valsSim,repTN,date_list,plotA=TRUE,locationI=locationtab[iiH])
      
    }
    
    
    # - - - - - - - - - - - 
    # Set up matrices for MCMC run
    # - - - - - - - - - - -
    
    m = 1 #set initial MCMC
    thetatab=matrix(NA,nrow=(MCMC.runs+1),ncol=length(theta))
    colnames(thetatab)=names(theta)
    thetatab[1,]=theta
    
    thetaAlltab=array(NA, dim=c(MCMC.runs+1,locnn,length(thetaAll[1,])),dimnames=list(NULL,NULL,names(thetaAll)))
    thetaAlltab[1,,]=as.matrix(thetaAll)
    
    theta_initAlltab=array(NA, dim=c(MCMC.runs+1,locnn,length(theta_initAll[1,])),dimnames=list(NULL,NULL,names(theta_initAll)))
    theta_initAlltab[1,,]=as.matrix(theta_initAll)
    
    sim_liktab=rep(-Inf,(MCMC.runs+1))
    accepttab=rep(NA,(MCMC.runs))
    max.length = max(length(time.vals)) + 1 # Need to store enough values
    c_trace_tab=array(NA, dim=c(MCMC.runs+1,locnn,max.length)) # Note max length
    s_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    c_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    s_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    c_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    r_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    r_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    x_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    x_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
    
    # - - - - - - - - - - - 
    #RUN MCMC:
    # - - - - - - - - - - - 
    
    for (m in 1:MCMC.runs){
      
      if(m==1){
        epsilon0=0.001
        cov_matrix_theta=epsilon0*cov_matrix_theta0
        cov_matrix_thetaA=epsilon0*cov_matrix_thetaA0
        cov_matrix_theta_init=epsilon0*cov_matrix_theta_init0
      }else{
        epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^m)),1e-6) # Stop epsilon getting too big or small
        cov_matrix_theta=epsilon0*cov_matrix_theta0
        cov_matrix_thetaA=epsilon0*cov_matrix_thetaA0
        cov_matrix_theta_init=epsilon0*cov_matrix_theta_init0
      }
      
      # Resample global theta
      if(m>1){
        output_theta = SampleTheta(thetatab[m,],theta_initAlltab[1,1,],m,cov_matrix_theta,cov_matrix_theta_init,singleI=NULL,pmask) #sample nearby global parameter space
        theta_star=output_theta$thetaS
      }else{
        theta_star=thetatab[m,]
      }
      
      # Resample outbreak-specific parameters
      sim_marg_lik_star=0
      thetaAllstar=0*thetaAlltab[m,,]
      theta_initAllstar=0*theta_initAlltab[m,,]
      sTraceCStar=0*s_trace_tabC[m,,]; sTraceAStar=0*s_trace_tabA[m,,]
      cTraceStar=0*c_trace_tab[m,,]
      rTraceCStar=0*r_trace_tabC[m,,]; rTraceAStar=0*r_trace_tabA[m,,]
      cTraceCStar=0*c_trace_tabC[m,,]; cTraceAStar=0*c_trace_tabA[m,,]
      xTraceCStar=0*x_trace_tabC[m,,]; xTraceAStar=0*x_trace_tabA[m,,]
      
      for(kk in itertabM){
        iiH=kk
        time.vals = time.store[[iiH]]$t_store
        y.vals = time.store[[iiH]]$y_store
        repTN=1
        
        # Resample specific theta - note this includes masked variables
        if(m==1){ # Don't resample on 1st step
          output_H = SampleTheta(thetaAlltab[m,iiH,],theta_initAlltab[m,iiH,],m,0*cov_matrix_thetaA,0*cov_matrix_theta_init,singleI=kk,pmask) #sample parameter space
        }else{
          output_H = SampleTheta(thetaAlltab[m,iiH,],theta_initAlltab[m,iiH,],m,cov_matrix_thetaA,cov_matrix_theta_init,singleI=kk,pmask) #sample parameter space
        } 
        
        thetaA_star=output_H$thetaS
        theta_init_star=output_H$theta_initS
        
        # Run model simulation
        output1= Deterministic_modelR(1, dt, c(theta_star,thetaA_star), theta_init_star, y.vals,y.vals2, y.vals.prop,time.vals,repTN,locationI=locationtab[iiH])
        sim_marg_lik_star=sim_marg_lik_star+output1$lik
        
        #Store vales
        thetaAllstar[iiH,]=thetaA_star
        theta_initAllstar[iiH,]=theta_init_star
        sTraceCStar[iiH,(1:length(output1$S_traceC))]=output1$S_traceC
        sTraceAStar[iiH,(1:length(output1$S_traceA))]=output1$S_traceA
        cTraceStar[iiH,(1:length(output1$C_trace))]=output1$C_trace
        rTraceCStar[iiH,(1:length(output1$R_traceC))]=output1$R_traceC
        rTraceAStar[iiH,(1:length(output1$R_traceA))]=output1$R_traceA
        cTraceCStar[iiH,(1:length(output1$C_traceC))]=output1$C_traceC
        cTraceAStar[iiH,(1:length(output1$C_traceA))]=output1$C_traceA
        xTraceCStar[iiH,(1:length(output1$X_traceC))]=output1$X_traceC
        xTraceAStar[iiH,(1:length(output1$X_traceA))]=output1$X_traceA
      } # end loop over regions
      
      # child/adult case ratio - NOT USED IN LIKELIHOOD
      #c_case_ratioStar = sum(cTraceCStar[2,],na.rm = T) /sum(cTraceStar[2,],na.rm = T)
      #c_case_ratioTab = if(m>1){sum(c_trace_tabC[m,2,],na.rm = T) /sum(c_trace_tab[m,2,],na.rm = T)}else{c_case_ratioStar}
      
      # Calculate probability function
      output_prob = ComputeProbability(sim_liktab[m],sim_marg_lik_star,thetatab[m,],theta_star,theta_initAlltab[m,,],theta_initAllstar,thetaAllstar,itertab,c_case_ratioStar,c_case_ratioTab) 
      
      # Update parameter values
      
      if(runif(1,0,1) < output_prob){
        thetatab[m+1,] = theta_star
        thetaAlltab[m+1,,] = thetaAllstar
        theta_initAlltab[m+1,,] = theta_initAllstar
        s_trace_tabC[m+1,,]=sTraceCStar; s_trace_tabA[m+1,,]=sTraceAStar
        r_trace_tabC[m+1,,]=rTraceCStar; r_trace_tabA[m+1,,]=rTraceAStar
        c_trace_tab[m+1,,]=cTraceStar
        c_trace_tabC[m+1,,]=cTraceCStar; c_trace_tabA[m+1,,]=cTraceAStar
        x_trace_tabC[m+1,,]=xTraceCStar; x_trace_tabA[m+1,,]=xTraceAStar
        sim_liktab[m+1] = sim_marg_lik_star
        accepttab[m]=1
        
      }else{
        thetatab[m+1,] = thetatab[m,]
        thetaAlltab[m+1,,] = thetaAlltab[m,,]
        theta_initAlltab[m+1,,] = theta_initAlltab[m,,]
        s_trace_tabC[m+1,,]=s_trace_tabC[m,,]; s_trace_tabA[m+1,,]=s_trace_tabA[m,,]
        r_trace_tabC[m+1,,]=r_trace_tabC[m,,]; r_trace_tabA[m+1,,]=r_trace_tabA[m,,]
        c_trace_tab[m+1,,]=c_trace_tab[m,,]
        c_trace_tabC[m+1,,]=c_trace_tabC[m,,]; c_trace_tabA[m+1,,]=c_trace_tabA[m,,]
        x_trace_tabC[m+1,,]=x_trace_tabC[m,,]; x_trace_tabA[m+1,,]=x_trace_tabA[m,,]
        sim_liktab[m+1] = sim_liktab[m]
        accepttab[m]=0
      }
      
      if(m<20){
        accept_rate=0.234
      }else{
        accept_rate=sum(accepttab[1:m])/m
      }

      # Save outputs every 1000 iterations
      
      if(m %% min(MCMC.runs,1000) == 0){
        print(c(m,accept_rate,sim_liktab[m],epsilon0))
        save(sim_liktab,accepttab,s_trace_tabC,s_trace_tabA,c_trace_tab,r_trace_tabC,r_trace_tabA,x_trace_tabC,x_trace_tabA,thetatab,thetaAlltab,theta_initAlltab,file=paste("outputs/outputR",country.name,"_",epi.name,iiM,"_ELISA_",use.ELISA.data,exclude.p,".RData",sep=""))
      }
      
    } # End MCMC run
    
  } # End multichains
  
}
