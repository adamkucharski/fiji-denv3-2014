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

plot_posteriors<-function(p_pick=4){
  
  paramA=NULL
  
  locnnLoop=1
    
    for(iiH in 1:locnnLoop){

      # Import multiple chains if available
      #source("load_posteriors.R",local=TRUE)
      pick_posterior = p_pick
      source("R/load_posterior_single.R",local=TRUE)
      
      max_season = 2 # peak seasonality
      
      r0_vtab = (thetatab$beta[picks]/thetatab$r_inf[picks])*(thetatab$beta_v[picks]/thetatab$mu_v[picks])*(thetatab$v_exp[picks]/(thetatab$v_exp[picks]+thetatab$mu_v[picks]))*max_season
      
      # PLOT HISTOGRAMS
      
      par(mfrow=c(3,6),mar = c(3,3,1,1),mgp=c(2,0.7,0))
      
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

      brekN=15 #,breaks=brekN,

      hist(thetatab$beta[picks],xlab=expression(beta[h]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      curve(priorBetaH2M(x), col="red", lwd=2, add=TRUE, yaxt="n")
      
      #hist(thetatab$beta_v[picks],xlab=expression('a'[v]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      hist(thetatab$beta_v[picks],xlab=expression(beta[v]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      curve(priorBetaM2H(x), col="red", lwd=2, add=TRUE, yaxt="n")
      
      #hist(thetatab$beta_v_amp[picks],xlab="rainfall effect",prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      
      
      hist(thetatab$beta_c_grad[picks],xlab=expression('a'[1]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      
      hist(thetatab$beta_c_base[picks],xlab=expression('a'[2]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      hist(thetatab$beta_c_mid[picks],xlab=expression('a'[tau]),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      
      hist(thetatab$repR[picks],xlab="propn cases reported (lab)",main=NULL,border=colB,col=colW,prob=TRUE)
      hist(thetatab$repRA[picks],xlab="propn cases reported (DLI)",main=NULL,border=colB,col=colW,prob=TRUE)
      
      hist(thetatab$repvol[picks],xlab="reporting dispersion",main=NULL,border=colB,col=colW,prob=TRUE)
      #hist(thetatab$repvolA[picks],xlab="reporting dispersion (2)",main=NULL,border=colB,col=colW,prob=TRUE)

      hist(thetatab$recruit_m[picks],xlab="recruitment rate",main=NULL,border=colB,col=colW,prob=TRUE)
      
      
      #hist(thetatab$repvol[picks],xlab=expression(phi),main=NULL,col=rgb(0.5,0.8,1),prob=TRUE)
      hist(theta_inittab$i1_initC[picks],xlab=expression('I'[hc]^0),main=NULL,border=colB,col=colW,prob=TRUE)
      hist(theta_inittab$i1_initA[picks],xlab=expression('I'[ha]^0),main=NULL,border=colB,col=colW,prob=TRUE)
      #hist(theta_inittab$r_initC[picks],xlab=expression('R'[hc]^0),main=NULL,border=colB,col=colW,prob=TRUE)
      #hist(theta_inittab$r_initA[picks],xlab=expression('R'[ha]^0),main=NULL,border=colB,col=colW,prob=TRUE)
      hist(theta_inittab$im_initC[picks],xlab=expression('I'[v]^0),main=NULL,border=colB,col=colW,prob=TRUE)
      
      dev.copy(pdf,paste("plots/Figure_S4_posteriors_",use.ELISA.data,exclude.p,"_",p_pick,".pdf",sep=""),width=10,height=6) #,locationtab[iiH],
      dev.off()
      
      # Output theta for future runs
      max.init.names = c("i1_initC","im_initC")
      max.out.theta.init = theta_inittab[pick.max,max.init.names]

      max.names = c("beta","beta3","beta3","beta_v","beta_v_amp","beta3","beta3","beta_c_mid","beta_c_grad","beta2","beta3","beta3","repR","repRA","rep_drop","repvol","repvolA","r_inf")
      max.out.theta = thetatab[pick.max,max.names]
      max.param = cbind(c(max.names,max.init.names),c(max.out.theta,max.out.theta.init))
      #write.csv(max.param,"plots/Table_5_params_part1.csv")
      

      # PLOT Parameter correlations
      param.names = c("beta","beta_v","recruit_m","beta_c_grad","beta_c_base","beta_c_mid"); # beta_c_base is base level
      param.labels = c(expression(beta[h]),expression(beta[v]),expression('a'[0]),expression('a'[1]),expression('a'[2]),expression('a'[tau]))

      par(mfcol=c(length(param.names),length(param.names)))
      par(mar = c(3,3,1,1),mgp=c(1.8,0.5,0))
      
      thetatab0 = thetatab %>% data.frame()
      sample.p = sample(length(thetatab0$beta),1000,replace=T)
      thinner.theta=thetatab0[sample.p,]

        for(ii in 1:length(param.names)){
          for(jj in 1:length(param.names)){
            if(ii<=jj){
              if(ii == jj){
                hist(thetatab0[[param.names[ii]]],xlab=param.labels[ii],main=NULL) #paste("ESS=",round(effectiveSize(thetatab0[[param.names[ii]]])))
              }else{
                plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
              }
            }else{
              plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.3, xlab=param.labels[ii], ylab=param.labels[jj])
              points(median(thinner.theta[[param.names[ii]]]),median(thinner.theta[[param.names[jj]]]),col="orange",cex=1.5,lwd=2)

            }

          }
        }
        
      dev.copy(png,paste("plots/Figure_S5_correlations",use.ELISA.data,exclude.p,"_",p_pick,".png",sep=""),units="cm",width=20,height=20,res=200)
      dev.off()

      # Compile beta matrix. Structure:  beta; 0.5 * beta * beta2  ; 0.5 * beta


      source("R/load_timeseries_data.R",local=TRUE)
      wks=length(time.vals)
      
      # Include bootstrap reporting uncertainty - denominator is total people infected

      repBTS=sapply(sample(picks,1000,replace=T),function(x){rr.bt=(thetatab$npop[x]-s_trace_tabC[x,wks]-s_trace_tabA[x,wks]); rnbinom(1,mu=(rr.bt*thetatab[x,"repR"]),size=1/thetatab[x,"repvol"])/rr.bt })
      repBTS2=sapply(sample(picks,1000,replace=T),function(x){rr.bt=(thetatab$npop[x]-s_trace_tabC[x,wks]-s_trace_tabA[x,wks]); rnbinom(1,mu=(rr.bt*thetatab[x,"repRA"]),size=1/thetatab[x,"repvolA"])/rr.bt })
      
      param1=cbind(
        c.text(r0_vtab,2),
        c.text(100*repBTS,2),
        c.text(100*repBTS2,2),
        c.text(theta_inittab$i1_initC[picks],2),
        c.text(theta_inittab$i1_initA[picks],2),
        c.text(100*(1-s_trace_tabC[picks,wks]/thetatab$npopC[picks]),2),
        c.text(100*(r_trace_tabC[picks,1]/thetatab$npopC[picks]),2),
        c.text(thetatab$beta[picks],2),
        c.text(thetatab$beta_v[picks],2),
        c.text(thetatab$beta3[picks]*thetatab$beta2[picks],2),
        c.text(thetatab$beta3[picks],2),
        c.text(1/thetatab$r_exp[picks],2),
        c.text(1/thetatab$r_inf[picks],2),
        c.text(1/thetatab$v_exp[picks],2),
        c.text(1/thetatab$mu_v[picks],2),
        c.text(thetatab$repR[picks],2),
        c.text(thetatab$repvol[picks],2),
        max(sim_liktab)
        )
      
      rownames(param1)=c(locationtab[iiH])
      colnames(param1)=c("R0","propn reported (%)","propn reported 2 (%)","I_HC(0)","I_HA(0)","final size","prop_imm","beta_h","beta_v","beta2","beta3","alpha_h","gamma","alpha_v","delta","r","phi","max_lik")
      
      paramA=rbind(paramA,param1)

    }
  
  write.csv(t(paramA),paste("plots/Table_5_params_part1_",use.ELISA.data,exclude.p,"_",p_pick,".csv",sep=""))

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
  
  locnnLoop=1
  par(mfrow=c(2,2),mgp=c(2,0.7,0),mar = c(3,4,1,3))
  


  # PLOT WEATHER
  xRange = c(as.Date("2013-11-01"),as.Date("2014-11-01")) # DATE OF SAMPLE COLLECTION
  weather.date = as.Date(weather.data$Date) + 15
  
  plot(weather.date,weather.data$Rain_av,type="l",lty=1,col="blue",xlab="",ylim=c(0,500),ylab="rainfall (mm) ",lwd=0,xlim=xRange)
  for(ii in 1:(length(weather.data$Rain_av)/12)-1){ 
    lines(weather.date+365*ii,weather.data$Rain_av,type="l",lty=1,col=rgb(0,0,1,0.2),xaxs="i",lwd=1)
  }
  lines(weather.date,weather.data$Rain_av,type="l",lty=1,col=rgb(0,0,1),xaxs="i",lwd=2)
  
  # Overlap MAP estimate of temperature function?
  #lines(date_listSeason,500*(plotCosMR0[1,]-0.55),type="l",col=rgb(0,0,0),ylim=c(0.5,2),xaxs="i",lwd=3)
  
  
  par(new=TRUE)
  plot(weather.date,weather.data$Av_temp,type="l",lty=1,col="red",xaxt="n",yaxt="n",xlab="",ylab="",lwd=0,ylim=c(21,28),xlim=xRange)
  for(ii in 1:(length(weather.data$Av_temp)/12)-1){
    lines(weather.date+365*ii,weather.data$Av_temp,type="l",lty=1,col=rgb(1,0,0,0.2),xaxs="i",lwd=1)
  }
  lines(weather.date,weather.data$Av_temp,type="l",lty=1,col=rgb(1,0,0),xaxs="i",lwd=2)
  
  axis(4,col="red",col.axis="red")
  mtext("temperature (°C)", side=4, line=2,col="red",cex=0.8) # Label for 2nd axis
  
  title(main=LETTERS[1],adj=0)
  
  # PLOT MORDECAI ET AL FUNCTION
  plot(temp.R0.data$aegy.temps.DTR8,temp.R0.data$R0.rel,type="l",xlab="temperature (°C)",ylab=expression("relative R"[0]))
  min(weather.data$Av_temp)
  weather2014.relative.R0[1]
  lines(c(max(weather2014$Av_temp),max(weather2014$Av_temp)),c(-1,2),col="blue")
  lines(c(min(weather2014$Av_temp),min(weather2014$Av_temp)),c(-1,2),col="blue")
  title(main=LETTERS[2],adj=0)
  
  # PLOT SEASONAL FUNCTION
  date0 = (as.Date("2013-10-28")-date_list[1]) %>% as.numeric() # Shift back as simulation starts from 2013-10-28
  date_listSeason = seq(min(xRange),max(xRange),7)
  thetaBASE = c(beta_v_mask=1,beta_v_amp=((1-rel.temp)/(1+rel.temp)),beta_c_mask=1,beta_c_base=0.5,beta_c_grad=10,beta_c_mid=0.38)
 
  plot(date_listSeason,seasonal_f(as.numeric(date_listSeason-min(date_listSeason)+7),date0,thetaBASE),type="l",xlab="date",ylab="relative transmission",col="red",lwd=2,ylim=c(0,2))
  lines(c(min(xRange),max(xRange)),c(1,1),lty=2)
  title(main=LETTERS[3],adj=0)
  
  # PLOT CONTROL
  plot(date_listSeason,decline_f(as.numeric(date_listSeason-min(date_listSeason)+7),date0,thetaBASE),type="l",xlab="date",ylab="relative transmission",col=rgb(0,0.6,0.3),lwd=2,ylim=c(0,2))
  lines(c(min(xRange),max(xRange)),c(1,1),lty=2)
  polygon(c(as.Date("2014-03-08"),as.Date("2014-03-08"),as.Date("2014-03-22"),as.Date("2014-03-22")),c(-1,1e4,1e4,-1),col=rgb(1,0,0,0.2),lty=0)
  title(main=LETTERS[4],adj=0)
  
  dev.copy(pdf,paste("plots/Figure_S3_illustrate_control.pdf",sep=""),width=8,height=6)#,height=8)
  dev.off()
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - 
# Model comparison

model_comparison <- function( compareN = c(2:4) ){

  # Need to update these
  param.init = 5
  param.all = 7
  paramlist = c(0,0,2,2+3,0,0,2,2+3)
  
  aic.tab = NULL
  dic.tab = NULL
  lik.tab = NULL
  
  # Load time series dataset - need to initial timeser as global (from main model.R)
  source("R/load_timeseries_data.R",local=TRUE)
  
  for(p_pick in compareN){ # NOTE CURRENTLY JUST 2013/14
    
    # Import multiple chains if available
    pick_posterior = p_pick
    source("R/load_posterior_single.R",local=TRUE)
    aic.model = 2*(param.init + param.all + paramlist[p_pick]) - 2*max(sim_likOut)
    aic.tab = c(aic.tab,aic.model)
    
    lik.tab = c(lik.tab,max(sim_likOut))

    # Validate likelihood
    output1 = Deterministic_modelR(1,dt, thetatab[pick.max,], theta_inittab[pick.max,], y.vals,y.vals2,y.vals.prop,time.vals,repTN,locationI=locationtab[1])
    print( output1$lik)
    
    output1 = Deterministic_modelR(1,dt, apply(thetatab,2,mean), apply(theta_initAlltab[picks,iiH,],2,mean), y.vals,y.vals2,y.vals.prop,time.vals,repTN,locationI=locationtab[1])
    loglik_theta_bar = output1$lik
    
    deviance.at.post.mean = -2*loglik_theta_bar #-2*max(sim_likOut)
    effective.param = var(-2*sim_likOut)/2
    dic.calc = deviance.at.post.mean + effective.param
    dic.tab = c(dic.tab,deviance.at.post.mean)
  }
  
  aic.comp = aic.tab - min(aic.tab)
  dic.comp = dic.tab - min(dic.tab)
  name.models = c("SEIR","SEIR_climate","SEIR_climate_control") #,"SEIR_post","SEIR_climate_post","SEIR_climate_control_post")
  write.csv( cbind(name.models,signif(cbind(lik.tab,aic.tab,aic.comp,dic.tab,dic.comp),4)) ,paste("plots/Table_S3_DIC_table",use.ELISA.data,exclude.p,".csv",sep=""))
  
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
  par(mgp=c(1.7,0.5,0),mar = c(3,3,1,3),mai=c(0.2,0.5,0.2,0.4))
  
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
    svectorC[ii,]= r_trace_tabC[pick,1:tMax]/thetatab[pick,]$npopC # Proportion immune C thetatab[pick,"prop_at_risk"]*
    svectorA[ii,]= r_trace_tabA[pick,1:tMax]/thetatab[pick,]$npopA # Proportion immune A
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
    )),paste("plots/Table_5_params_part2_",use.ELISA.data,exclude.p,"_",p_pick,".csv",sep=""))


  # - - 
  # Plot interventions 
  
  #par(new=TRUE)
  #plot(date_listSeason,plotCosMR0[3,],type="l",col=rgb(0,0,0,0),ylim=c(0,2),xaxt="n",yaxt="n",xlim=xRange2,xlab="",ylab="")
  plot(date_listSeason,plotReduceM[3,],type="l",col=rgb(0,0,0,0),ylim=c(0,2),xlim=xRangeTimeS,xlab="",ylab=paste("R and relative transmission",sep=""))
  
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
  # Plot weather patterns
  
  # Define data
  weather.date = as.Date(weather.data$Date) + 15
  
  t_numeric = as.numeric(seq(start.date-100,start.date+800,1) - start.date)
  temp_plot = seasonaltemp(t_numeric,theta_fit)
  temp_actl = weather.data$Av_temp
  rain_actl = weather.data$Rain_av
  rain_plot = seasonalrain(t_numeric,theta_fitRain) 
  
  # Normalise values?
  temp_plot = temp_plot
  temp_actl = temp_actl
  #rain_actl = rain_actl/mean(rain_actl)
  
  plot(t_numeric+start.date,t_numeric,type="l",lty=1,col="white",xlab="",ylim=c(21,27),ylab="temperature (°C) ",lwd=0,xaxs="i",xlim=xRange2)
  #for(ii in 1:(length(weather.data$Rain_av)/12)-1){ 
  #points(weather.date+365*ii,weather.data$Av_temp,lty=1,col=rgb(1,0,0,0.3),xaxs="i",pch=1,cex=1)
  #}
  points(weather.date,temp_actl,col=rgb(1,0,0),xaxs="i",pch=19,cex=1)
  lines(t_numeric+start.date,temp_plot,col="red",lty=1 )
  title(main=LETTERS[4],adj=0)
  
  points(weather.date,rain_actl,col=rgb(0,0,1),xaxs="i",cex=1)
  
  par(new=TRUE)
  plot(weather.date,rain_actl,col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,500),xlim=xRange)
  #points(weather.date,rain_actl,col=rgb(0,0,1),xaxs="i")
  lines(t_numeric+start.date,rain_plot,col="blue",lty=1 )
  axis(4,col="blue",col.axis="blue")
  mtext("rainfall (mm)", side=4, line=1.7,col="blue",cex=0.7) # Label for 2nd axis
  
  
  
  # - - 
  # Plot basic reproduction number
  plot(date_listSeason,plotR0_both[1,],type="l",col="white",xlab="",ylab=expression(paste(R[0]," and normalised rainfall",sep="")),xaxs="i",ylim=c(0.3,2),xlim=xRange2)
  polygon(c(date_listSeason,rev(date_listSeason)),c(plotR0_both[2,],rev(plotR0_both[3,])),lty=0,col=rgb(0,0,0,0.2))
  lines(date_listSeason,plotR0_both[1,],type="l",col=rgb(0,0,0),xlab="",ylab="")
  
  # Plot R=1 line
  lines(c(min(date_listSeason),max(date_listSeason)),c(1,1),col="black",lty=3)
  

  title(main=LETTERS[5],adj=0)

  
  dev.copy(pdf,paste("plots/Figure_5_",use.ELISA.data,exclude.p,"_",p_pick,".pdf",sep=""),width=8,height=6)#,height=8)
  dev.off()

  
}