# - - - - - - - - - - - - - - - - - - - - - - - 
# Transmission modelling plot posteriors
# Author: Adam Kucharski
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - - 

# Need to debug all the reporting functions


mcmc.burn=0.2


c.text <- function(x,sigF=3){

  
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

# - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot posteriors

plot_posteriors<-function( p_pick=4 ){
  
  paramA=NULL
  
  locnnLoop=1
  iiH = 1

  # Import posteriors
  pick_posterior = p_pick
  source("R/load_posterior_single.R",local=TRUE)
  
  # PLOT HISTOGRAMS
  
  par(mfrow=c(3,6),mar = c(3,3,1,1),mgp=c(2,0.7,0),las=0)
  
  if(iiH==1){
    colW='grey' 
    colM=rgb(0.2,0.4,1)
    colB='white'
    baW=0.3
    breaks0=seq(0,50,1)
    hist(1/thetatab$v_exp[picks],breaks=breaks0,xlab=expression("incubation period (v)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,30))
    curve(priorVEx(x), col="red", lwd=2, add=TRUE, yaxt="n")

    hist(1/thetatab$r_exp[picks],breaks=breaks0,xlab=expression("incubation period (h)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,20))
    curve(priorExp(x), col="red", lwd=2, add=TRUE, yaxt="n")
    #lines(density(1/thetatab$r_exp[picks],bw=baW), col=colM,lwd=2)
    
    hist(1/thetatab$r_inf[picks],breaks=breaks0,xlab=expression("infectious period (h)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,20))
    curve(priorInf(x), col="red", lwd=2, add=TRUE, yaxt="n")
    #lines(density(1/thetatab$r_inf[picks],bw=baW), col=colM,lwd=2)
    
    hist(1/thetatab$mu_v[picks],breaks=breaks0,xlab=expression("lifespan (v)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,20))
    curve(priorMuV(x), col="red", lwd=2, add=TRUE, yaxt="n")
  
  }

  brekN=15 

  hist(thetatab$m_density[picks],xlab=expression('density'),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
  curve(priorDensity(x), col="red", lwd=2, add=TRUE, yaxt="n")

  hist(thetatab$beta_v[picks],xlab=expression(beta[v]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
  curve(priorBeta_v(x), col="red", lwd=2, add=TRUE, yaxt="n")
  
  hist(thetatab$beta_v_amp[picks],xlab="rainfall amplitude",main=NULL,border=colB,col=colW,prob=TRUE)

  hist(thetatab$beta_c_grad[picks],xlab=expression('a'[1]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
  hist(thetatab$beta_c_base[picks],xlab=expression('a'[2]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
  
  control_mid = ( control.range /( 1 + exp(-10*(thetatab$beta_c_mid[picks] -1))))*365
  
  hist(control_mid,xlab=expression('a'[tau]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
  
  hist(thetatab$repR[picks],xlab="propn cases reported (lab)",main=NULL,border=colB,col=colW,prob=TRUE)
  hist(thetatab$repRA[picks],xlab="propn cases reported (DLI)",main=NULL,border=colB,col=colW,prob=TRUE)
  hist(thetatab$repvol[picks],xlab="reporting dispersion",main=NULL,border=colB,col=colW,prob=TRUE)


  hist(theta_inittab$i1_initC[picks],xlab=expression('I'[hc]^0),main=NULL,border=colB,col=colW,prob=TRUE)
  hist(theta_inittab$i1_initA[picks],xlab=expression('I'[ha]^0),main=NULL,border=colB,col=colW,prob=TRUE)
  hist(theta_inittab$r_initC[picks]/thetatab$npopC[picks],xlab=expression('R'[hc]^0),main=NULL,border=colB,col=colW,prob=TRUE)
  hist(theta_inittab$r_initA[picks]/thetatab$npopA[picks],xlab=expression('R'[ha]^0),main=NULL,border=colB,col=colW,prob=TRUE)
  hist(theta_inittab$im_initC[picks],xlab=expression('I'[v]^0),main=NULL,border=colB,col=colW,prob=TRUE)
  
  dev.copy(pdf,paste("plots/Figure_S4_posteriors_",use.ELISA.data,"_",p_pick,".pdf",sep=""),width=10,height=6) #,locationtab[iiH],
  dev.off()
  
  # Output max likelihood theta for future runs
  max.init.names = c("i1_initC","im_initC")
  max.out.theta.init = theta_inittab[pick.max,max.init.names]

  max.names = c("beta","beta_v","beta_v_amp","beta_v_mask","beta_c_mid","beta_c_grad","beta_c_base","beta_c_mask","m_density","repR","repRA","repvol","npop","npopC","npopA")
  max.out.theta = thetatab[pick.max,max.names]
  max.param = cbind(c(max.names,max.init.names),c(max.out.theta,max.out.theta.init))
  
  write.csv(max.param,"plots/Table_5_params_best.csv")


  # PLOT Parameter correlations
  param.names = c("m_density","beta_v","beta_v_amp","beta_c_grad","beta_c_base","beta_c_mid"); # beta_c_base is base level
  param.labels = c("density","contact rate",expression('rain amp'),expression('a'[1]),expression('a'[2]),expression('a'[tau]))
  
  par(mfcol=c(length(param.names),length(param.names)))
  par(mar = c(3,3,1,1),mgp=c(1.8,0.5,0))
  
  thetatab0 = thetatab %>% data.frame()
  sample.p = sample(length(thetatab0$beta),1000,replace=T)
  thinner.theta=thetatab0[sample.p,]

    for(ii in 1:length(param.names)){
      for(jj in 1:length(param.names)){
        if(ii<=jj){
          if(ii == jj){
            hist(thetatab0[[param.names[ii]]],xlab=param.labels[ii],main=NULL,freq=F) #paste("ESS=",round(effectiveSize(thetatab0[[param.names[ii]]])))
          }else{
            plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
          }
        }else{
          plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.3, xlab=param.labels[ii], ylab=param.labels[jj])
          points(median(thinner.theta[[param.names[ii]]]),median(thinner.theta[[param.names[jj]]]),col="orange",cex=1.5,lwd=2)

        }

      }
    }
    
  dev.copy(png,paste("plots/Figure_S5_correlations",use.ELISA.data,"_",p_pick,".png",sep=""),units="cm",width=20,height=20,res=200)
  dev.off()

  # Compile beta matrix. Structure:  beta; 0.5 * beta * beta2  ; 0.5 * beta

  source("R/load_timeseries_data.R",local=TRUE)
  wks=length(time.vals)
  
  # Include bootstrap reporting uncertainty - denominator is total people infected

  repBTS=sapply(sample(picks,1000,replace=T),function(x){rr.bt=(thetatab$npop[x]-s_trace_tabC[x,wks]-s_trace_tabA[x,wks]); rnbinom(1,mu=(rr.bt*thetatab[x,"repR"]),size=1/thetatab[x,"repvol"])/rr.bt })
  repBTS2=sapply(sample(picks,1000,replace=T),function(x){rr.bt=(thetatab$npop[x]-s_trace_tabC[x,wks]-s_trace_tabA[x,wks]); rnbinom(1,mu=(rr.bt*thetatab[x,"repRA"]),size=1/thetatab[x,"repvol"])/rr.bt })
  
  param1=cbind(
    c.text(1/thetatab$r_exp[picks],2),
    c.text(1/thetatab$r_inf[picks],2),
    c.text(1/thetatab$v_exp[picks],2),
    c.text(1/thetatab$mu_v[picks],2),
    c.text(thetatab$beta_v[picks],2),
    c.text(thetatab$m_density[picks],2),
    c.text(thetatab$beta_v_amp[picks],2),
    
    c.text(thetatab$beta_c_grad[picks],2),
    c.text(thetatab$beta_c_base[picks],2),
    c.text(control_mid,2),
    c.text(thetatab$repR[picks],2),
    c.text(thetatab$repRA[picks],2),
    c.text(thetatab$repvol[picks],2),
    
    c.text(theta_inittab$i1_initC[picks],2),
    c.text(theta_inittab$i1_initA[picks],2),
    c.text(r_trace_tabC[picks,1]/thetatab$npopC[picks],2),
    c.text(r_trace_tabA[picks,1]/thetatab$npopA[picks],2),
    c.text(theta_inittab$im_initC[picks],2),

    c.text(100*repBTS,2),
    c.text(100*repBTS2,2)
    )
  
  rownames(param1)=c(locationtab[iiH])
  colnames(param1)=c("1/nu_h","1/gamma","1/nu_v","1/delta","alpha","m","K","a_1","a_2","a_tau","r_lab","r_DLI","rho","I_HC(0)","I_HA(0)","R_HC(0)","R_HA(0)","I_v(0)","propn reported lab (%)","propn reported DLI (%)")
  
  paramA=param1
  
  write.csv(t(paramA),paste("plots/Table_5_params_part1_",use.ELISA.data,"_",p_pick,".csv",sep=""))

}


# - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot Mordecai et al. data and weather

plot_weather_and_control <- function(){
  

  # Set up vectors
  dataP=NULL
  simM=F
  p_pick = 4
  labelN=1

  title.lab = "" 
  
  iiH = 1

  pick_posterior = p_pick
  source("R/load_posterior_single.R",local=TRUE)
  source("R/load_timeseries_data.R",local=TRUE)
  
  xRange = c(as.Date("2013-11-04"),as.Date("2015-10-16"))
  date_listSeason = seq(min(xRange),max(xRange),7)
  
  locnnLoop=1
  ll=1
  
  par(mfrow=c(2,5),mgp=c(2,0.7,0),mar = c(3,3,1,1))
  xRange = c(as.Date("2013-10-28"),as.Date("2014-11-01")) 
  weather.date = as.Date(weather.data$Date) + 15 # Select weather data
  
  tt_actualD = (seq(start.date-7,start.date+400,1) )
  tt_actualN = (tt_actualD - (start.date-7)) %>% as.numeric() # Convert to numeric

  # - - 
  # Plot temperature
  plot(weather.data.daily$date,xaxs="i",weather.data.daily$min_air_temp/2+weather.data.daily$max_air_temp/2,type="l",lty=1,col="orange",xlab="2013/14",ylab="temperature (°C)",lwd=1,ylim=c(21,31),xlim=xRange)
  #lines(weather.data.daily$date,weather.data.daily$max_air_temp,type="l",lty=1,col="red",lwd=1)
  
  lines(tt_actualD,seasonaltemp(tt_actualN,theta_fit1))
  
  title(main=LETTERS[ll],adj=0); ll=ll+1
  
  # - - 
  # Plot rainfall
  plot(weather.date,weather.data$Rain_av,type="l",lty=1,xaxs="i",col="blue",xlab="2013/14",ylim=c(0,500),ylab="rainfall (mm) ",lwd=0,xlim=xRange)
  for(ii in 1:(length(weather.data$Rain_av)/12)-1){ 
    lines(weather.date+365*ii,weather.data$Rain_av,type="l",lty=1,col=rgb(0,0,1,0.2),xaxs="i",lwd=1)
  }
  lines(weather.date,weather.data$Rain_av,type="l",lty=1,col=rgb(0,0,1),xaxs="i",lwd=2)

  title(main=LETTERS[ll],adj=0); ll=ll+1
  

  # Define weather range
  mean_temp = weather.data.daily$min_air_temp/2+weather.data.daily$max_air_temp/2

  # Plot vector parameters
  temp_plot = seq(20,30,0.1)
  weather_range = c(min(mean_temp[weather.data.daily$date<=as.Date("2014-10-01")]),max(mean_temp[weather.data.daily$date<=as.Date("2014-10-01")]))
  
  bite_temp(temp_plot)
  EI_rate_temp(temp_plot)
  MD_rate_temp(temp_plot)
  eggs_per_female_temp(temp_plot)
  eggs_to_adult_temp(temp_plot)
  prob_to_h_temp(temp_plot)
  prob_to_v_temp(temp_plot)
  
  # Lifespan
  plot(temp_plot,prior_p_MuV[1]/mortality_rate_temp(temp_plot),type="l",xlab="temperature",ylab="lifespan",ylim=c(6,8.5))
  lines(c(weather_range[2],weather_range[2]),c(-1,200),col="blue",lty=2); lines(c(weather_range[1],weather_range[1]),c(-1,200),col="blue",lty=2)
  title(main=LETTERS[ll],adj=0); ll=ll+1

  # EIP
  plot(temp_plot,prior_p_VEx[1]/EI_rate_temp(temp_plot),type="l",xlab="temperature",ylab="EIP",ylim=c(6,20))
  lines(c(weather_range[2],weather_range[2]),c(-1,200),col="blue",lty=2); lines(c(weather_range[1],weather_range[1]),c(-1,200),col="blue",lty=2)
  title(main=LETTERS[ll],adj=0); ll=ll+1

  # Biting rate
  plot(temp_plot,prior_betaV[1]*bite_temp(temp_plot),type="l",xlab="temperature",ylab="biting rate (per day)",ylim=c(0,0.5))
  lines(c(weather_range[2],weather_range[2]),c(-1,200),col="blue",lty=2); lines(c(weather_range[1],weather_range[1]),c(-1,200),col="blue",lty=2)
  title(main=LETTERS[ll],adj=0); ll=ll+1
  
  # Probability to human
  plot(temp_plot,prob_to_h_temp(temp_plot),type="l",xlab="temperature",ylab=expression('p'[hv]),ylim=c(0,1))
  lines(c(weather_range[2],weather_range[2]),c(-1,200),col="blue",lty=2); lines(c(weather_range[1],weather_range[1]),c(-1,200),col="blue",lty=2)
  title(main=LETTERS[ll],adj=0); ll=ll+1
  
  # Probability to vector
  plot(temp_plot,prob_to_v_temp(temp_plot),type="l",xlab="temperature",ylab=expression('p'[vh]),ylim=c(0,1))
  lines(c(weather_range[2],weather_range[2]),c(-1,200),col="blue",lty=2); lines(c(weather_range[1],weather_range[1]),c(-1,200),col="blue",lty=2)
  title(main=LETTERS[ll],adj=0); ll=ll+1


  # Density
  plot(temp_plot,density_vary(temp_plot),type="l",xlab="temperature",ylab="normalised density",ylim=c(0,1.5))
  lines(c(weather_range[2],weather_range[2]),c(-1,200),col="blue",lty=2); lines(c(weather_range[1],weather_range[1]),c(-1,200),col="blue",lty=2)
  title(main=LETTERS[ll],adj=0); ll=ll+1
  
  # - - -
  # Plot carrying capacity function
  rainfall_range = seq(0,400,10)
  theta_b_amp = 0.01; yy1 = 1/(1+1/(0 + theta_b_amp*rainfall_range/222.4374))
  theta_b_amp = 1; yy2 = 1/(1+1/(0 + theta_b_amp*rainfall_range/222.4374))
  theta_b_amp = 100; yy3 = 1/(1+1/(0 + theta_b_amp*rainfall_range/222.4374))
  
  plot(rainfall_range,yy1/max(yy1),type="l",lty=1,col="blue",xlab="rainfall (mm)",ylim=c(0,1),ylab="relative density")
  lines(rainfall_range,yy2/max(yy2),lty=2,col="blue")
  lines(rainfall_range,yy3/max(yy3),lty=3,col="blue")
  title(main=LETTERS[ll],adj=0); ll=ll+1
  
  # - - 
  # Plot control
  # - ADD, date_listSeason
  thetaBASE = c(beta_v_mask=1,beta_c_mask=1,beta_c_base=0.5,beta_c_grad=100,beta_c_mid=0.9)
  plot(date_listSeason,decline_f(as.numeric(date_listSeason-min(date_listSeason)+7),date0=0,thetaBASE),type="l",xaxs="i",xlim=xRange,xlab="2013/14",ylab="relative transmission",col=rgb(0,0.6,0.3),lwd=2,ylim=c(0,2))
  lines(c(min(xRange),max(xRange)),c(1,1),lty=2)
  polygon(c(as.Date("2014-03-08"),as.Date("2014-03-08"),as.Date("2014-03-22"),as.Date("2014-03-22")),c(-1,1e4,1e4,-1),col=rgb(1,0,0,0.2),lty=0)
  title(main=LETTERS[ll],adj=0); ll=ll+1

  
  dev.copy(pdf,paste("plots/Figure_S3_illustrate_control.pdf",sep=""),width=10,height=4)#,height=8)
  dev.off()
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - 
# Model comparison

model_comparison <- function( compareN = c(2,3,4) ){

  # Parameter numbers for AIC
  param.init = 5
  param.all = 7
  paramlist = c(0,0,1,1+3)
  
  aic.tab = NULL
  dic.tab = NULL
  lik.tab = NULL
  dev.tab = NULL
  
  # Load time series dataset - need to initial timeser as global (from main model.R)
  source("R/load_timeseries_data.R",local=TRUE)
  
  for(p_pick in compareN){ # NOTE CURRENTLY JUST 2013/14
    
    # Import multiple chains if available
    pick_posterior = p_pick
    source("R/load_posterior_single.R",local=TRUE)
    
    # Calculate BIC for each model
    deviance.at.max = - 2*max(sim_likOut)
    aic.model = 2*(param.init + param.all + paramlist[p_pick]) + deviance.at.max
    aic.tab = c(aic.tab,aic.model)
    
    lik.tab = c(lik.tab,max(sim_likOut))

    # Evaluate at maximum likelihood as a quick check
    output1 = Deterministic_modelR(1,dt, thetatab[pick.max,], theta_inittab[pick.max,], y.vals,y.vals2,y.vals.prop,time.vals,repTN,locationI=locationtab[1])
    print( output1$lik)
    
    # Calculate median parameters and incorporate initial conditions -  use median given correlation structure in posterior
    theta_init_mean = apply(theta_initAlltab[picks,iiH,],2,median)
    theta_mean = apply(thetatab,2,median)
    
    theta_init_mean[["s_initC"]] = theta_mean[["npopC"]] - theta_init_mean[["e_initC"]] - theta_init_mean[["i1_initC"]] - theta_init_mean[["r_initC"]]
    theta_init_mean[["s_initA"]] = theta_mean[["npopA"]] - theta_init_mean[["e_initA"]] - theta_init_mean[["i1_initA"]] - theta_init_mean[["r_initA"]]
    theta_init_mean[["sm_initC"]] = 1 - theta_init_mean[["em_initC"]] - theta_init_mean[["im_initC"]]
    
    output1 = Deterministic_modelR(1,dt, theta_mean, theta_init_mean, y.vals,y.vals2,y.vals.prop,time.vals,repTN,locationI=locationtab[1])
    loglik_theta_bar = output1$lik
    
    # Calculate DIC for each model
    deviance.at.post.mean = -2*loglik_theta_bar 
    effective.param = var(-2*sim_likOut)/2

    dic.calc = deviance.at.post.mean + 2*effective.param
    
    dic.tab = c(dic.tab,dic.calc)
    dev.tab = c(dev.tab,deviance.at.post.mean)
  }
  
  # Compile outputs
  aic.comp = aic.tab - min(aic.tab)
  dic.comp = dic.tab - min(dic.tab)
  name.models = c("SEIR","SEIR_climate","SEIR_climate_control") #,"SEIR_post","SEIR_climate_post","SEIR_climate_control_post")
  write.csv( cbind(name.models,signif(cbind(lik.tab,aic.tab,aic.comp,dic.tab,dic.comp),4)) ,paste("plots/Table_S3_DIC_table",use.ELISA.data,".csv",sep=""))
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot timeseries for DENV 3 

plot_figure_2014_dengue3 <- function(p_pick=4,simM=FALSE,Fmask=FALSE,long_time=F,DoubleFit=F){
  
  # 
  # Set up vectors
  # p_pick = 4; simM=FALSE; Fmask=FALSE; long_time = F; DoubleFit =F
  dataP=NULL
  simM=F
  
  locnnLoop=1

  layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,5,5), 4,byrow=T) )
  par(mgp=c(1.7,0.5,0),mar = c(3,3,1,3),mai=c(0.2,0.5,0.2,0.4),las=0)
  
  labelN=1

  title.lab = "" 
  
  iiH = 1 # PLOT 2013/14
  
  # Import multiple chains if available
  pick_posterior = p_pick
  source("R/load_posterior_single.R",local=TRUE)
  
  # Load time series dataset - need to initial timeser as global (from main model.R)
  source("R/load_timeseries_data.R",local=TRUE)

  repTab=rR.vals/totalSS[iiH]
  tMax <- length(time.vals) #length(y.vals)
  
  btsp=1000
  PickList = sample(picks,btsp,replace=T)
  
  cvector=matrix(NA,nrow=btsp,ncol=tMax)
  cvector_lab=matrix(NA,nrow=btsp,ncol=tMax)
  cvectorALL=matrix(NA,nrow=btsp,ncol=tMax)
  svectorC=matrix(NA,nrow=btsp,ncol=tMax)
  svectorA=matrix(NA,nrow=btsp,ncol=tMax)
  
  for(ii in 1:btsp){
    pick= PickList[ii]
    cvector[ii,]= ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,],repSS="sus",y.vals.prop) #+ ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,],repSS="lab",y.vals.prop) # Edit for both datasets
    cvector_lab[ii,]= ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,],repSS="lab",y.vals.prop) #+ ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,],repSS="lab",y.vals.prop) # Edit for both datasets
    cvectorALL[ii,]= ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,],repSS="sus",y.vals.prop) + ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,],repSS="lab",y.vals.prop)
    svectorC[ii,]= thetatab[pick,"prop_at_risk"]*r_trace_tabC[pick,1:tMax]/thetatab[pick,]$npopC # Proportion immune C 
    svectorA[ii,]= thetatab[pick,"prop_at_risk"]*r_trace_tabA[pick,1:tMax]/thetatab[pick,]$npopA # Proportion immune A
  }
  
  medP=apply(cvector,2,function(x){median(x)})
  ciP1=apply(cvector,2,function(x){quantile(x,0.025)})
  ciP2=apply(cvector,2,function(x){quantile(x,0.975)})
  ciP150=apply(cvector,2,function(x){quantile(x,0.25)})
  ciP250=apply(cvector,2,function(x){quantile(x,0.75)})

  medP_lab=apply(cvector_lab,2,function(x){median(x)})
  ciP1_lab=apply(cvector_lab,2,function(x){quantile(x,0.025)})
  ciP2_lab=apply(cvector_lab,2,function(x){quantile(x,0.975)})
  
  
  medP_All=apply(cvectorALL,2,function(x){median(x)})
  ciP1_All=apply(cvectorALL,2,function(x){quantile(x,0.025)})
  ciP2_All=apply(cvectorALL,2,function(x){quantile(x,0.975)})
  ciP150_All=apply(cvectorALL,2,function(x){quantile(x,0.25)})
  ciP250_All=apply(cvectorALL,2,function(x){quantile(x,0.75)})
  
  medP_RC=apply(svectorC,2,function(x){median(x)}); medP_RC=c(medP_RC,rep(tail(medP_RC,1),150)) # EXTEND IF NEEDED
  ciP1_RC=apply(svectorC,2,function(x){quantile(x,0.025)}); ciP1_RC=c(ciP1_RC,rep(tail(ciP1_RC,1),150)) # EXTEND IF NEEDED
  ciP2_RC=apply(svectorC,2,function(x){quantile(x,0.975)}); ciP2_RC=c(ciP2_RC,rep(tail(ciP2_RC,1),150)) # EXTEND IF NEEDED
  
  medP_RA=apply(svectorA,2,function(x){median(x)}); medP_RA=c(medP_RA,rep(tail(medP_RA,1),150)) # EXTEND IF NEEDED
  ciP1_RA=apply(svectorA,2,function(x){quantile(x,0.025)}); ciP1_RA=c(ciP1_RA,rep(tail(ciP1_RA,1),150)) # EXTEND IF NEEDED
  ciP2_RA=apply(svectorA,2,function(x){quantile(x,0.975)}); ciP2_RA=c(ciP2_RA,rep(tail(ciP2_RA,1),150)) # EXTEND IF NEEDED
  
  # - - - - - - - 
  # PLOT RESULTS

  date_list_all = seq(min(date_list),max(date_list)+1000,7) # This expands time.vals to full list
  date_plot=min(date_list)+ 0 - 7 

  xRange = c(as.Date("2013-11-04"),as.Date("2015-10-16")); xRange2 = xRange # DATES OF TIMESERIES
  if(long_time==T){
    xRangeTimeS = c(as.Date("2013-11-04"),as.Date("2015-10-16")) #as.Date("2015-12-16")) # DATE OF TIME SERIES as.Date("2014-09-01"))
  }else{
    xRangeTimeS = c(as.Date("2013-11-04"),as.Date("2014-08-30")) #as.Date("2014-09-01")# DATE OF TIME SERIES as.Date("2014-09-01"))
  }

  xSelect=(1:length(date_list))
  xSelect2= 1:tMax 
  xSelect2DROP1 = xSelect2[1:length(xSelect2)] #:(thetaAlltab[1,1,"rep_drop"])] # Adjust for pre- reporting drop
  xSelect2DROP2 = date_list>as.Date("2014-01-13")
  datacol=rgb(0.4,0.4,0.4)
  
  # PLOT TIMESERIES
  
  #par(mfrow=c(2,1),mgp=c(2,0.7,0),mar = c(3,4,1,4))
  
  ylimmax=1.1*max(c(y.vals,ciP2[(1:(length(date_list)+5))])) # DEBUG ylimmax XX 
  plot(date_list_all[xSelect2],medP[xSelect2],ylim=c(0,2000),xlab="",ylab=ifelse(iiH==1,"cases","confirmed cases"),pch=19,cex=1,col='white',yaxs="i",xlim=xRangeTimeS)
  

  if(DoubleFit==TRUE){ # Show model outputs?

    #polygon(c(date_list_all[xSelect2DROP1],rev(date_list_all[xSelect2DROP1])),c(ciP150[xSelect2DROP1],rev(ciP250[xSelect2DROP1])),lty=0,col=rgb(0,0.3,1,0.3))
    polygon(c(date_list_all[xSelect2DROP1],rev(date_list_all[xSelect2DROP1])),c(ciP1[xSelect2DROP1],rev(ciP2[xSelect2DROP1])),lty=0,col=rgb(0,0.3,1,0.2))
    lines(date_list_all[xSelect2DROP1],medP[xSelect2DROP1],type="l",col=rgb(0,0.3,1),lty=2,xaxt="n",yaxt="n",xlab="",ylab="")
    
    #polygon(c(date_list_all[xSelect2DROP2],rev(date_list_all[xSelect2DROP2])),c(ciP150[xSelect2DROP2],rev(ciP250[xSelect2DROP2])),lty=0,col=rgb(0,0.3,1,0.3))
    polygon(c(date_list_all[xSelect2DROP1],rev(date_list_all[xSelect2DROP1])),c(ciP1_lab[xSelect2DROP1],rev(ciP2_lab[xSelect2DROP1])),lty=0,col=rgb(0,0.3,1,0.2))
    lines(date_list_all[xSelect2DROP1],medP_lab[xSelect2DROP1],type="l",col=rgb(0,0.3,1),xaxt="n",yaxt="n",xlab="",ylab="")
    
    points(date_list,y.vals[date_list>"2014-01-13"],ylim=c(0,ylimmax),pch=19,cex=1,col='black') # Lab fitted points
    points(date_list,y.vals2,ylim=c(0,ylimmax),pch=1,cex=1,col='black') # DLI fitted points
  }


  if(DoubleFit==FALSE){ # Show model outputs?
    
    polygon(c(date_list_all[xSelect2DROP1],rev(date_list_all[xSelect2DROP1])),c(ciP150_All[xSelect2DROP1],rev(ciP250_All[xSelect2DROP1])),lty=0,col=rgb(0,0.3,1,0.3))
    
    polygon(c(date_list_all[xSelect2DROP1],rev(date_list_all[xSelect2DROP1])),c(ciP1_All[xSelect2DROP1],rev(ciP2_All[xSelect2DROP1])),lty=0,col=rgb(0,0.3,1,0.2))
    lines(date_list_all[xSelect2DROP1],medP_All[xSelect2DROP1],type="l",col=rgb(0,0.3,1),xaxt="n",yaxt="n",xlab="",ylab="")
    
    lines(date_list,y.vals+y.vals2,ylim=c(0,ylimmax),lwd=1.5,col='black') # Lab fitted points
    #points(date_list,y.vals+y.vals2,ylim=c(0,ylimmax),pch=19,cex=1,col='black') # Lab fitted points
    points(date_list,y.vals,ylim=c(0,ylimmax),pch=19,cex=1,col='black') # Lab fitted points
    points(date_list[xSelect2DROP2],y.vals2[xSelect2DROP2],ylim=c(0,ylimmax),pch=1,cex=1,col='black') # DLI fitted points in later stages
    
  }

  # Plot clean up campaign dates
  polygon(c(as.Date("2014-03-08"),as.Date("2014-03-08"),as.Date("2014-03-22"),as.Date("2014-03-22")),c(-1,1e4,1e4,-1),col=rgb(1,0,0,0.3),lty=0)

  title(main=LETTERS[1],adj=0) #paste(locationtabF[iiH],sep=""))

  
  # - - - 
  # PLOT SEROCONVERSION -- DEFINE BASED ON Luminex?
  
  test1C=binom.test(x=n_Luminex_C_D3[1],n=n_Luminex_C_D3[3]); test2C=binom.test(x=n_Luminex_C_D3[2],n=n_Luminex_C_D3[3])
  test1A=binom.test(x=n_Luminex_A_D3[1],n=n_Luminex_A_D3[3]); test2A=binom.test(x=n_Luminex_A_D3[2],n=n_Luminex_A_D3[3])


  col2=rgb(1,0.5,0,1)
  col2a=rgb(1,0.5,0,0.2)
  plot(date_list_all[xSelect2],medP_RC[xSelect2],type="l",col=NULL,xlab="",ylab="proportion seropositive",ylim=c(0,1),xlim=xRange2)
  
  if(Fmask==FALSE){
    polygon(c(date_list_all[xSelect2],rev(date_list_all[xSelect2])),c(ciP1_RC[xSelect2],rev(ciP2_RC[xSelect2])),lty=0,col=col2a)
    lines(date_list_all[xSelect2],medP_RC[xSelect2],lty=2,col=col2,xaxt="n",yaxt="n",xlab="",ylab="")
    
    polygon(c(date_list_all[xSelect2],rev(date_list_all[xSelect2])),c(ciP1_RA[xSelect2],rev(ciP2_RA[xSelect2])),lty=0,col=col2a)
    lines(date_list_all[xSelect2],medP_RA[xSelect2],type="l",col=col2,xaxt="n",yaxt="n",xlab="",ylab="")
  }
  
  shiftA=-1
  if(iiH==1){
    max.date.sero = max(date_list_all[xSelect2]) # max(xRange2)
    # Children sero in 2013/15
    points(max.date.sero-shiftA,n_Luminex_C_D3[2]/n_Luminex_C_D3[3],pch=1,col="red")
    lines(c(max.date.sero,max.date.sero)-shiftA,c(test2C$conf.int[1],test2C$conf.int[2]),pch=19,col="red")
    points(min(date_list_all)-shiftA,n_Luminex_C_D3[1]/n_Luminex_C_D3[3],pch=1,col="red")
    lines(c(min(date_list_all),min(date_list_all))-shiftA,c(test1C$conf.int[1],test1C$conf.int[2]),pch=19,col="red")
    # Adults sero in 2013/15
    points(max.date.sero,n_Luminex_A_D3[2]/n_Luminex_A_D3[3],pch=19,col="red")
    lines(c(max.date.sero,max.date.sero),c(test2A$conf.int[1],test2A$conf.int[2]),pch=19,col="red")
    points(min(date_list_all),n_Luminex_A_D3[1]/n_Luminex_A_D3[3],pch=19,col="red")
    lines(c(min(date_list_all),min(date_list_all)),c(test1A$conf.int[1],test1A$conf.int[2]),pch=19,col="red")
  }
  
  
  
  title(main=LETTERS[2],adj=0)
  
  
  # Plot seasonality function -- only for 2013/14
  date_listSeason = seq(min(xRange),max(xRange),7)
  t.start = 0
  time.V = (1:length(date_listSeason))*7 #time.vals + t.start
  if(length(time.V)>tMax){extra_date = (tMax+1):length(time.V)}else{extra_date = NULL}
  cosPick = 1:length(time.V)
  
  xRange2 = xRange #; xRange2[2] = xRange[2] + length(extra_date)*7 - 7*3
  
  plotCosR0 = NULL; plotReduce = NULL;plotCosR0vary = NULL # plotBoth = NULL; 
  plotR0store = NULL; plotRRstore = NULL; plotRCVstore = NULL; plotRAVstore = NULL; plotRVCstore = NULL
  btstrap = sample(picks,1000,replace=T)
  plotSuS_C = NULL;     plotSuS_A = NULL;     plotSuS_V = NULL
  
  btstrap = 100
  
  for(ii in 1:btstrap){
    b_ii = PickList[ii]
    #beta_ii = seasonal_f(time.V,date0= 0,thetatab[b_ii, ]) # INCLUDE SHIFT START WEEK LATER
    
    beta_control = sapply(time.V,function(x){decline_f(x,date0= thetatab[b_ii,"shift_date" ],thetatab[b_ii, ])})
    
    # Compile individual beta
    plotReduce=rbind(plotReduce,  beta_control)
    #plotCosR0=rbind(plotCosR0,  beta_ii)
    #plotBoth=rbind(plotBoth,  beta_ii*beta_control)
    
    #plotCosR0vary=rbind(plotCosR0vary,  1-min(beta_ii)/max(beta_ii))
    
    #Compile R/R0
    t.rpick = length(time.V)
    s_pickH = c(s_trace_tabC[b_ii,1:tMax]+s_trace_tabA[b_ii,1:tMax],rep(s_trace_tabC[b_ii,1:tMax]+s_trace_tabA[b_ii,tMax],length(extra_date)) )/(thetatab$npopC[b_ii] + thetatab$npopA[b_ii]) # Note this is reporting week
    x_pickC = c(x_trace_tabC[b_ii,1:tMax],rep(x_trace_tabC[b_ii,tMax],length(extra_date)) ) # Mosquito susceptibility =1
    x_pickA = c(x_trace_tabA[b_ii,1:tMax],rep(x_trace_tabA[b_ii,tMax],length(extra_date)) ) # Mosquito susceptibility =1
    
    # Susceptible children and adults
    s_pickC = c(s_trace_tabC[b_ii,1:tMax],rep(s_trace_tabC[b_ii,tMax],length(extra_date)) )/thetatab$npopC[b_ii] # Note this is reporting week
    s_pickA = c(s_trace_tabA[b_ii,1:tMax],rep(s_trace_tabA[b_ii,tMax],length(extra_date)) )/thetatab$npopA[b_ii] # Note this is reporting week
    
    # NOTE DEBUG -- HAVE REMOVED MOSQUITOES
    output_rr = calculate_r0(th_in=thetatab[b_ii,],sus_c=s_pickH[cosPick],sm_c=x_pickC[cosPick]/x_pickA[cosPick],t_vary=time.V,controlT=1)
    output_rr_nocontrol = calculate_r0(th_in=thetatab[b_ii,],sus_c=1+0*s_pickH[cosPick],sm_c=1+0*x_pickC[cosPick],t_vary=time.V,controlT=0)
    
    r0_post = output_rr_nocontrol$rr_out
    rr_post = output_rr$rr_out
    
    plotR0store=rbind(plotR0store,  r0_post)
    plotRRstore=rbind(plotRRstore,  rr_post)
    plotSuS_C = rbind(plotSuS_C,s_pickC[cosPick])
    plotSuS_A = rbind(plotSuS_A,s_pickA[cosPick])
    plotSuS_V = rbind(plotSuS_V,x_pickC[cosPick])
    
    # Store H-V and V-H reproduction numbers
    plotRCVstore=rbind(plotRCVstore,  output_rr$rr_mat[1,2])
    plotRVCstore=rbind(plotRVCstore,  output_rr$rr_mat[2,1])
    
  }
  
  plotReduceM = apply(plotReduce,2,c.nume)  #plotCosMinMax = apply(plotCosR0vary,2,c.nume) ; plotCosMR0 = apply(plotCosR0,2,c.nume) ; 
  plotRR_both = apply(plotRRstore,2,c.nume) 
  plotR0_both = apply(plotR0store,2,c.nume) 
  
  # - - - - - - - - - - - - - - - - 
  # Output R0 calcs
  # Proportion of year above 1
  
  # Calculate variation in seasonality too...
  
  meanR0 = apply(plotR0_both,1,mean) %>% signif(3)
  maxR0 = apply(plotR0_both,1,max) %>% signif(3)
  weeksR0 = (apply(plotR0_both>1,1,sum)/(1+as.numeric(max(date_listSeason)-min(date_listSeason))/7)) %>% signif(3)
  reduceProp = 1-c(tail(plotReduceM[1,],1),tail(plotReduceM[2,],1),tail(plotReduceM[3,],1)) %>% signif(3)
  #seasonProp = plotCosMinMax %>% signif(3)
  
  write.csv(
    cbind(c("meanR0","maxR0","controlreduce","R_cv","R_vh"),
    rbind(
    paste(meanR0[1]," (",meanR0[2],"-",meanR0[3],")",sep=""), # mean R0
    paste(maxR0[1]," (",maxR0[2],"-",maxR0[3],")",sep=""), # max R0
    #paste(weeksR0[1]," (",weeksR0[2],"-",weeksR0[3],")",sep=""), # weeks R0>1
    #paste(seasonProp[1]," (",seasonProp[2],"-",seasonProp[3],")",sep=""), # seasonal reduction
    paste(reduceProp[1]," (",reduceProp[3],"-",reduceProp[2],")",sep=""), # control reduction
    c.text(plotRCVstore), # R to H
    c.text(plotRVCstore) # R to V
    )),paste("plots/Table_5_params_part2_",use.ELISA.data,"_",p_pick,".csv",sep=""))


  # - - 
  # Plot interventions 
  
  #par(new=TRUE)
  #plot(date_listSeason,plotCosMR0[3,],type="l",col=rgb(0,0,0,0),ylim=c(0,2),xaxt="n",yaxt="n",xlim=xRange2,xlab="",ylab="")
  plot(date_listSeason,plotReduceM[3,],type="l",col=rgb(0,0,0,0),ylim=c(0,2.5),xlim=xRangeTimeS,xlab="",ylab=expression(paste("reproduction number, relative transmission",sep="")))
  
  # Plot clean up campaign dates
  polygon(c(as.Date("2014-03-08"),as.Date("2014-03-08"),as.Date("2014-03-22"),as.Date("2014-03-22")),c(-1,1e4,1e4,-1),col=rgb(1,0,0,0.3),lty=0)

  season_cut = date_listSeason
  plotReduceM = plotReduceM[,1:length(date_listSeason)]
  
  # Plot control measures
  if(p_pick>=4){
    polygon(c(season_cut,rev(season_cut)),c(plotReduceM[2,],rev(plotReduceM[3,])),lty=0,col=rgb(0,0.6,0.3,0.2))
    lines(season_cut,plotReduceM[1,],type="l",col=rgb(0,0.6,0.3),xlab="",ylab="")
    
  }
  
  # Plot R=1 line
  lines(c(min(date_listSeason),max(date_listSeason)),c(1,1),col="black",lty=3)
  
  # Plot basic reproduction number
  polygon(c(date_listSeason,rev(date_listSeason)),c(plotRR_both[2,],rev(plotRR_both[3,])),lty=0,col=rgb(0,0,1,0.2))
  lines(date_listSeason,plotRR_both[1,],type="l",col=rgb(0,0,1),xlab="",ylab="")
  
  # Plot effective reproduction number
  polygon(c(date_listSeason,rev(date_listSeason)),c(plotR0_both[2,],rev(plotR0_both[3,])),lty=0,col=rgb(0,0,0,0.2))
  lines(date_listSeason,plotR0_both[1,],type="l",col=rgb(0,0,0),xlab="",ylab="")

  title(main=LETTERS[3],adj=0)

  # - - 
  # Plot weather patterns - NEED TO UPDATE
  
  # Define data
  weather.date = as.Date(weather.data$Date) + 15
  
  t_numeric = as.numeric(seq(start.date,start.date+800,1) )
  temp_plot = seasonaltemp(t_numeric,theta_fit1)
  temp_actl = weather.data$Av_temp
  rain_actl = weather.data$Rain_av
  rain_plot = seasonalrain(t_numeric,theta_fitRain) 
  
  # Normalise values?
  tt_actualD = (seq(start.date-7,start.date+700,1) )
  tt_actualN = (tt_actualD - (start.date-7)) %>% as.numeric() # Convert to numeric
  
  plot(weather.data.daily$date,xaxs="i",weather.data.daily$min_air_temp/2+weather.data.daily$max_air_temp/2,type="l",lty=1,col="orange",xlab="2013/14",ylab="temperature (°C)",lwd=1,ylim=c(21,31),xlim=xRange)
  #lines(weather.data.daily$date,weather.data.daily$max_air_temp,type="l",lty=1,col="red",lwd=1)
  lines(tt_actualD,seasonaltemp(tt_actualN,theta_fit1),col="black")

  title(main=LETTERS[4],adj=0)

  par(new=TRUE)
  plot(weather.date,rain_actl,col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,500),xlim=xRange,type="l")
  #points(weather.date,rain_actl,col=rgb(0,0,1),xaxs="i",cex=1,pch=19)
  axis(4,col="blue",col.axis="blue")
  mtext("rainfall (mm)", side=4, line=1.7,col="blue",cex=0.7) # Label for 2nd axis

  # - - 
  # Plot basic reproduction number
  plot(date_listSeason,plotR0_both[1,],type="l",col="white",xlab="",ylab=expression(paste(R[0],sep="")),xaxs="i",ylim=c(0.3,2.5),xlim=xRange2)
  polygon(c(date_listSeason,rev(date_listSeason)),c(plotR0_both[2,],rev(plotR0_both[3,])),lty=0,col=rgb(0,0,0,0.2))
  lines(date_listSeason,plotR0_both[1,],type="l",col=rgb(0,0,0),xlab="",ylab="")
  
  # Plot R=1 line
  lines(c(min(date_listSeason),max(date_listSeason)),c(1,1),col="black",lty=3)
  

  title(main=LETTERS[5],adj=0)

  
  dev.copy(pdf,paste("plots/Figure_5_",use.ELISA.data,"_",p_pick,".pdf",sep=""),width=8,height=6)#,height=8)
  dev.off()

  
}