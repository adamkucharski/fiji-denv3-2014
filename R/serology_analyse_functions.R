# - - - - - - - - - - - - - - - - - - - - - - - 
# Functions for Fiji serological analysis
# Author: Adam Kucharski
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - -

# Note: environmental variables are described in the 'fit_ELISA_univariable()' function
# AGE_2015 is used to identify age band (AGE_GROUP) in code

c_quantile1 <- function(x){as.numeric(quantile(x,0.025))}
c_quantile2 <- function(x){as.numeric(quantile(x,0.975))}

c_50_1 <- function(x){as.numeric(quantile(x,0.25))}
c_50_2 <- function(x){as.numeric(quantile(x,0.75))}

binom_output <- function(xx,nn){
  htest = binom.test(xx,nn, p = 1,conf.level=0.95)
  meanA=100*htest$estimate %>% as.numeric()  %>% signif(digits=3)
  conf1=100*htest$conf.int[1] %>% signif(digits=3)
  conf2=100*htest$conf.int[2] %>% signif(digits=3)
  c(paste(meanA," (",conf1,"-",conf2,")",sep="") )
}


# - - - - - - - - - - - - - - - - - - - - 
# Regression model for clusters with Luminex data - MAIN SET UP CODE
# - - - - - - - - - - - - - - - - - - - - 

fit_ELISA <- function(ELISA.fit = T, plot_venn=F){
  inputs =   data.frame(read.csv("data/serology_inputs.csv",stringsAsFactors = F))

  # - - - - - - - - -
  # Plot change in serology with mixture model
  
  if(ELISA.fit ==T){
    change.titre = inputs$ELISA2-inputs$ELISA1
    change.titre = change.titre[change.titre > -9] # Remove outliers
  }else{
    denv3.mia.1 = inputs$DENV3/1000 
    denv3.mia.2 = inputs$DENV3B/1000
    change.titre = denv3.mia.2 - denv3.mia.1
    change.titre = change.titre[change.titre>-5] # Remove outliers
    denv3.neut.1 = inputs$D3s13; denv3.neut.1 = denv3.neut.1[!is.na(denv3.neut.1)] #inputs[,c("D1s13","D2s13","D3s13","D4s13")]
    denv3.neut.2 = inputs$D3s15; denv3.neut.2 = denv3.neut.2[!is.na(denv3.neut.2)] #inputs[,c("D1s15","D2s15","D3s15","D4s15")]
  }

  param = c(a0 = 0, b0=1, a1=2, b1=3, lambda1 = 1)
  result2 = optim(param, fitmixture2, method="L-BFGS-B", val=change.titre,lower=c(rep(0,4),-10), hessian=FALSE) #, control=list(trace=1))
  
  parfit = result2$par #param #
  
  # Fit 2 mixture model
  
  x.titre = seq(-10,30,1)
  y.error = (1-transform.lambda(parfit[["lambda1"]]))*dnorm(x.titre, mean = 0, sd = parfit[["b0"]])
  y.signl1 = transform.lambda(parfit[["lambda1"]])*dgamma(x.titre, shape = parfit[["a1"]], scale = parfit[["b1"]]) 
  y.signl2 = 0
  
  p.signal = (y.signl1+y.signl2)/(y.error+y.signl1+y.signl2)
  inputsFit = inputs
  titre.weights = sapply(inputsFit$RiseT,function(x){ p.signal[match(x,x.titre)]   })

  # << Need to here to set up for remaining files << <<

  # Calculate R.squared of fits
  if(ELISA.fit ==T){
    xx <- seq(-10,20,0.1)
    all.model = transform.lambda(parfit[["lambda1"]])*dgamma(xx, shape = parfit[["a1"]], scale = parfit[["b1"]]) + (1-transform.lambda(parfit[["lambda1"]]))*dnorm(xx, mean = 0, sd = parfit[["b0"]])
    x.vals = names(table(change.titre)) %>% as.numeric()
    
    y.modl = all.model[match(x.vals,xx)]
    y.vals = as.numeric(table(change.titre)); y.tot=sum(y.vals); y.valsN=y.vals/y.tot
    r.square = 1-(sum((y.valsN-y.modl )^2)/sum((y.valsN-mean(y.valsN))^2))
  }else{
    xx <- seq(-10,20,0.1)
    all.model = transform.lambda(parfit[["lambda1"]])*dgamma(xx, shape = parfit[["a1"]], scale = parfit[["b1"]]) + (1-transform.lambda(parfit[["lambda1"]]))*dnorm(xx, mean = 0, sd = parfit[["b0"]])
    x.vals = names(table(round(change.titre,0))) %>% as.numeric()
    
    y.modl = all.model[match(x.vals,xx)]
    y.vals = as.numeric(table(round(change.titre,0))); y.tot=sum(y.vals); y.valsN=y.vals/y.tot
    r.square = 1-(sum((y.valsN-y.modl )^2)/sum((y.valsN-mean(y.valsN))^2))
    
    # - - - - - 
    # Plot MIA and neut data
    par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.7,0.5,0))
    
    # Plot Neut distribution
    
    # Put data on log scale
    denv3.neut.1.plot = denv3.neut.1; denv3.neut.1.plot[denv3.neut.1.plot==0]=5; denv3.neut.1.plot = log(as.numeric(denv3.neut.1.plot)/5,2)
    denv3.neut.2.plot = denv3.neut.2; denv3.neut.2.plot[denv3.neut.2.plot==0]=5; denv3.neut.2.plot = log(as.numeric(denv3.neut.2.plot)/5,2)

    hist(denv3.neut.1.plot,breaks=seq(-0.5,6.5,1),border="white",ylim=c(0,26),main="",xlab="DENV-3 neutralisation titre",xaxt="n",col=rgb(1,0.5,0,1),yaxs="i",ylab="number of samples")
    hist(denv3.neut.2.plot,freq=T,breaks=seq(-0.5,6.5,1),add=T,border="white",col=rgb(0,0,1,0.5))
    lines(c(1.5,1.5),c(0,100),lty=2)
    axis(1, at=c(0:8), labels=c(0,5*2^c(1:8)))
    title(adj=0,main=LETTERS[1])

    # Plot MIA distribution
    hist(denv3.mia.1*1000,freq=T,breaks=seq(-2000,26000,1000),col=rgb(1,0.5,0,1),border="white",main="",xlab="DENV-3 MIA value",ylab="number of samples",xlim=c(-2000,26000),xaxs="i",yaxs="i")
    hist(denv3.mia.2*1000,freq=T,breaks=seq(-2000,26000,1000),add=T,border="white",col=rgb(0,0,1,0.5))
    title(adj=0,main=LETTERS[2])
    
    txtSize=1
    text(x=1e4,y=40,adj=0,labels="2013",col=rgb(1,0.5,0),cex=txtSize)
    text(x=1e4,y=32,adj=0,labels="2015",col=rgb(0,0,1),cex=txtSize)
    
    dev.copy(pdf,paste("plots/Figure_S2_MIA_analysis.pdf",sep=""),width=8,height=4,useDingbats=FALSE)
    dev.off()
    
  }

  print(c("R2",r.square))
  
  if(ELISA.fit==F){return()}
  
  # - - - - - - - - - - - - - 
  # Calculate Kappa between seropositive Luminex and ELISA
  countEL = inputs; countEL$DconvertM = titre.weights
  countEL$sero1 = (countEL$ELISA1>=11) %>% as.numeric() # compare cutoff
  countEL$sero2 = (countEL$ELISA2>=11) %>% as.numeric() # compare cutoff
  
  # At least one MIA positive
  L1_positive = (rowSums(countEL[,c("DENV1P","DENV2P","DENV3P","DENV4P")])>0) %>% as.numeric()
  L2_positive = (rowSums(countEL[,c("DENV1BP","DENV2BP","DENV3BP","DENV4BP")])>0) %>% as.numeric()

  comp.1 = c(countEL$sero1,countEL$sero2); comp.2 = c(as.numeric(L1_positive),as.numeric(L1_positive))
  
  k_outEL = kappa2(cbind(comp.1,comp.2), "unweighted")  # Cohen's kappa
  k_accEL = ( sum(comp.1==1 & comp.2==1) +  sum(comp.1==0 & comp.2==0) )/length(comp.1)  # Accuracy
  
  # - - - - - - - - - - - - - 
  # Calculate Kappa between neutralisation and MIA for DENV-3
  countNL = inputs[!is.na(inputs$D3s15),];
  
  # At least one Luminex
  L1_positive_3 = ((countNL[,c("DENV3P")])>0) %>% as.numeric(); L2_positive_3 = ((countNL[,c("DENV3BP")])>0) %>% as.numeric()
  N1_positive_3 = (countNL$D3s13>=20) %>% as.numeric(); N2_positive_3 = (countNL$D3s15>=20) %>% as.numeric()

  comp.1N = c(N1_positive_3,N2_positive_3); comp.2N = c(as.numeric(L1_positive_3),as.numeric(L2_positive_3))

  k_outNL3 = kappa2(cbind(comp.1N,comp.2N), "unweighted")  # Cohen's kappa
  k_accNL3 = ( sum(comp.1N==1 & comp.2N==1) +  sum(comp.1N==0 & comp.2N==0) )/length(comp.1N)  # Accuracy
  
  # - - - - - - - - - - - - - 
  # Calculate Kappa between neutralisation, MIA and ELISA for Any dengue
  
  # At least one Luminex
  L1_positive_NT = (rowSums(countNL[,c("DENV1P","DENV2P","DENV3P","DENV4P")])>0) %>% as.numeric()
  L2_positive_NT = (rowSums(countNL[,c("DENV1BP","DENV2BP","DENV3BP","DENV4BP")])>0) %>% as.numeric()

  N1_positive_Any = (rowSums(countNL[,c("D1s13","D2s13","D3s13","D4s13")]>=20)>0) %>% as.numeric()
  N2_positive_Any = (rowSums(countNL[,c("D1s15","D2s15","D3s15","D4s15")]>=20)>0) %>% as.numeric()
  
  E1_positive_Any = (countNL$sero1==1) %>% as.numeric()
  E2_positive_Any = (countNL$sero2==1) %>% as.numeric()
  
  comp.1Any = c(N1_positive_Any,N2_positive_Any); comp.2Any = c(as.numeric(L1_positive_NT),as.numeric(L2_positive_NT)); comp.3Any = c(as.numeric(E1_positive_Any),as.numeric(E2_positive_Any))
  
  k_out_MN = kappa2(cbind(comp.1Any,comp.2Any), "unweighted")  # Cohen's kappa
  k_acc_MN = ( sum(comp.1Any==1 & comp.2Any==1) +  sum(comp.1Any==0 & comp.2Any==0) )/length(comp.1Any)  # Accuracy
  # k_outNL2 = kappa2(cbind(comp.1Any,comp.3Any), "unweighted")  # Cohen's kappa
  # k_accNL2 = ( sum(comp.1Any==1 & comp.3Any==1) +  sum(comp.1Any==0 & comp.3Any==0) )/length(comp.1Any)  # Accuracy
  k_out_ME = kappa2(cbind(comp.2Any,comp.3Any), "unweighted")  # Cohen's kappa
  k_acc_ME = ( sum(comp.2Any==1 & comp.3Any==1) +  sum(comp.2Any==0 & comp.3Any==0) )/length(comp.1Any)  # Accuracy
  
  # Create table of comparisons
  tab_kappa = rbind(
    c(paste(round(k_accEL,2)," (",round(k_outEL$value,2),")",sep=""),paste(round(k_acc_ME,2)," (",round(k_out_ME$value,2),")",sep=""),"-"), # ELISA
    c("-",paste(round(k_acc_MN,2)," (",round(k_out_MN$value,2),")",sep=""),"-"), # MIA any
    #c("-","-","-",paste(round(k_acc_MN,2)," (",round(k_out_MN$value,2),")",sep=""),"-"), # Neut Any
    c("-","-",paste(round(k_accNL3,2)," (",round(k_outNL3$value,2),")",sep="")) # MIA D3
    #c("-","-","-")  # Neut D3
  )
  
  colnames(tab_kappa) = c("MIA any DENV","Neut any DENV","Neut DENV-3")
  rownames(tab_kappa) = c("ELISA","MIA any DENV","MIA DENV-3")
  
  write.csv( tab_kappa,"plots/Kappa_compare.csv",row.names = T)
  
  
  
  # - - -
  # Plot combination of exposures
  if(plot_venn==T){
    d1 <- draw.quad.venn(sum(countEL$DENV1P==1), 
                         sum(countEL$DENV2P==1), 
                         sum(countEL$DENV3P==1), 
                         sum(countEL$DENV4P==1), 
                         sum(countEL$DENV1P==1 & countEL$DENV2P==1), #n12 
                         sum(countEL$DENV1P==1 & countEL$DENV3P==1), #n13 
                         sum(countEL$DENV1P==1 & countEL$DENV4P==1), #n14 
                         sum(countEL$DENV2P==1 & countEL$DENV3P==1), #n23 
                         sum(countEL$DENV2P==1 & countEL$DENV4P==1), #n24 
                         sum(countEL$DENV3P==1 & countEL$DENV4P==1), #n34 
                         sum(countEL$DENV1P==1 & countEL$DENV2P==1 & countEL$DENV3P==1), #n123
                         sum(countEL$DENV1P==1 & countEL$DENV2P==1 & countEL$DENV4P==1), #n124 
                         sum(countEL$DENV1P==1 & countEL$DENV3P==1 & countEL$DENV4P==1), #n134
                         sum(countEL$DENV2P==1 & countEL$DENV3P==1 & countEL$DENV4P==1), #n234
                         sum(countEL$DENV1P==1 & countEL$DENV2P==1 & countEL$DENV3P==1 & countEL$DENV4P==1),  #n1234
                         ind=F, category=c("DENV-1","DENV-2","DENV-3","DENV-4"),fontfamily="Arial",cat.fontfamily="Arial")
    
    tiff(filename = "plots/Figure_combination2013.tiff", compression = "lzw",width = 600, height = 700)
    print(grid.draw(d1))
    dev.off()
    
    d1 <- draw.quad.venn(sum(countEL$DENV1BP==1), 
                         sum(countEL$DENV2BP==1), 
                         sum(countEL$DENV3BP==1), 
                         sum(countEL$DENV4BP==1), 
                         sum(countEL$DENV1BP==1 & countEL$DENV2BP==1), #n12 
                         sum(countEL$DENV1BP==1 & countEL$DENV3BP==1), #n13 
                         sum(countEL$DENV1BP==1 & countEL$DENV4BP==1), #n14 
                         sum(countEL$DENV2BP==1 & countEL$DENV3BP==1), #n23 
                         sum(countEL$DENV2BP==1 & countEL$DENV4BP==1), #n24 
                         sum(countEL$DENV3BP==1 & countEL$DENV4BP==1), #n34 
                         sum(countEL$DENV1BP==1 & countEL$DENV2BP==1 & countEL$DENV3BP==1), #n123
                         sum(countEL$DENV1BP==1 & countEL$DENV2BP==1 & countEL$DENV4BP==1), #n124 
                         sum(countEL$DENV1BP==1 & countEL$DENV3BP==1 & countEL$DENV4BP==1), #n134
                         sum(countEL$DENV2BP==1 & countEL$DENV3BP==1 & countEL$DENV4BP==1), #n234
                         sum(countEL$DENV1BP==1 & countEL$DENV2BP==1 & countEL$DENV3BP==1 & countEL$DENV4BP==1),  #n1234
                         ind=F, category=c("DENV-1","DENV-2","DENV-3","DENV-4"),fontfamily="Arial",cat.fontfamily="Arial")
    
    tiff(filename = "plots/Figure_combination2015.tiff", compression = "lzw",width = 600, height = 700)
    print(grid.draw(d1))
    dev.off()
  }
  
  #- - 
  # Calculate proportion of each Luminex type
  
  outputLum = rbind( c(inputs$sero1 %>% sum() ,
                       inputs$sero2 %>% sum() ,length(inputs$sero2)),
                     c(L1_positive %>% sum() ,
                       L2_positive %>% sum() ,length(inputs$sero2)),
                        c(countEL[,"DENV1P"] %>% as.numeric() %>% sum() ,
                       countEL[,"DENV1BP"] %>% as.numeric() %>% sum() ,length(countEL[,"DENV1BP"])),
                     c(countEL[,"DENV2P"] %>% as.numeric() %>% sum() ,
                       countEL[,"DENV2BP"] %>% as.numeric() %>% sum() ,length(countEL[,"DENV1BP"])),
                     c(countEL[,"DENV3P"] %>% as.numeric() %>% sum() ,
                       countEL[,"DENV3BP"] %>% as.numeric() %>% sum() ,length(countEL[,"DENV1BP"])),
                     c(countEL[,"DENV4P"] %>% as.numeric() %>% sum() ,
                       countEL[,"DENV4BP"] %>% as.numeric() %>% sum() ,length(countEL[,"DENV1BP"]))
  )
  
  EL.2013 = c(binom.calc(outputLum[1,1],outputLum[1,3]),binom.calc(outputLum[2,1],outputLum[1,3]),binom.calc(outputLum[3,1],outputLum[1,3]),binom.calc(outputLum[4,1],outputLum[1,3]),binom.calc(outputLum[5,1],outputLum[5,3]),binom.calc(outputLum[6,1],outputLum[6,3]))
  EL.2015 = c(binom.calc(outputLum[1,2],outputLum[1,3]),binom.calc(outputLum[2,2],outputLum[1,3]),binom.calc(outputLum[3,2],outputLum[1,3]),binom.calc(outputLum[4,2],outputLum[1,3]),binom.calc(outputLum[5,2],outputLum[5,3]),binom.calc(outputLum[6,2],outputLum[6,3]))
  
  EL.difference = c(binom.calc(outputLum[1,2]-outputLum[1,1],outputLum[1,3]),binom.calc(outputLum[2,2]-outputLum[2,1],outputLum[1,3]),binom.calc(outputLum[3,2]-outputLum[3,1],outputLum[1,3]),binom.calc(outputLum[4,2]-outputLum[4,1],outputLum[1,3]),binom.calc(outputLum[5,2]-outputLum[5,1],outputLum[5,3]),binom.calc(outputLum[6,2]-outputLum[6,1],outputLum[6,3]))
  
  outputLum=cbind(outputLum,EL.2013,EL.2015,EL.difference
                  )
  
  colnames(outputLum) = c("2013","2015","N","2013pc","2015pc","increase"); rownames(outputLum) = c("ELISA","MIA any DENV","MIA DENV-1","MIA DENV-2","MIA DENV-3","MIA DENV-4")
  write.csv(outputLum,"plots/Table_2_seroprevalence.csv",row.names = T)
  
  list(tabinput = inputs,titrechange = titre.weights)
  
}

# - - - - - - - - - - - - - - - - - - - - 
# Functions for mixture model fitting

transform.lambda <- function(x){0.5/(1+exp(-x))}

fitmixture2<-function(param, val) {
  lambda1 = transform.lambda(param[["lambda1"]]) # constrain to be 0<.<1
  return(
    # likelihood function
    - sum(log(
      lambda1 * dgamma(val, shape = abs(param[["a1"]]), scale = abs(param[["b1"]])) +
        (1-lambda1) * dnorm(val, mean = 0, sd=abs(param[["b0"]]))
    ))
  )
}


# - - - - - - - - - - - - - - - - - - - - 
# Distributional bootstrap extraction - remove ELISA noise
# - - - - - - - - - - - - - - - - - - - - 

fit_ELISA_univariable <- function(inputs){ # use inputs from above fit_ELISA function
  
  # Load titre data
  loadinputs = fit_ELISA(ELISA.fit = T)
  inputsFit0 = loadinputs$tabinput
  
  titre.weights = loadinputs$titrechange
  
  inputsFit0$DconvertM = titre.weights %>% round()
  
  # merge urban and peri-urban (1= urban, 0=rural )
  inputsFit0[inputsFit0$GEOG==1 | inputsFit0$GEOG==2,"GEOG"] = 2
  inputsFit0$GEOG = as.numeric(inputsFit0$GEOG)
  inputsFit0$GEOG = 3-inputsFit0$GEOG
  inputsFit01 = inputsFit0
  inputsFit01[inputsFit01$ETHNIC!=1,"ETHNIC"]=0 # Any non-iTaukei individuals = 0
  
  # Seroconversion only
  inputsFit01=inputsFit01[inputsFit01$ELISA1<=9,] # ONLY PICK FULLY SERONEGATIVE
  inputsFit01$AGE_U_20 = as.numeric(inputsFit01$AGE_2015<20) # Age under 18
  inputsFit01$SEX = 1 - as.numeric(inputsFit01$SEX ) # Sex -  1=Female
  
  modelB.1 <- glm(Dconvert ~ AGE_U_20 , data = inputsFit01,family = "binomial") # Age in 2015
  modelB.1a <- glm(Dconvert ~ SEX , data = inputsFit01,family = "binomial") # Sex (Female = 1)
  modelB.1b <- glm(Dconvert ~ ETHNIC , data = inputsFit01,family = "binomial") # Ethnicity (iTaukei = 1; Other = 0)
  modelB.2 <- glm(Dconvert ~ I_MOS   , data = inputsFit01,family = "binomial") # Presence of mosquitoes
  modelB.3 <- glm(Dconvert ~ I_TIR    , data = inputsFit01,family = "binomial") # Presence of used car tires
  modelB.4 <- glm(Dconvert ~ I_WAT    , data = inputsFit01,family = "binomial") # Presence of open water containers (e.g. vase, bucket, basin, oil drum)
  modelB.5 <- glm(Dconvert ~ I_AC    , data = inputsFit01,family = "binomial") # Presence of long lasting puddles of water (i.e. do not dry up within a day after rain stops)
  modelB.6 <- glm(Dconvert ~ I_BLK   , data = inputsFit01,family = "binomial") # Presence of air conditioining
  modelB.7 <- glm(Dconvert ~ GEOG   , data = inputsFit01,family = "binomial") # Setting (urban/peri-urban = 1; rural = 0)
  modelB.8 <- glm(Dconvert ~ FEVER_2YR    , data = inputsFit01,family = "binomial") # Have you had a fever in the last 2 years? 
  modelB.9 <- glm(Dconvert ~ DOC_2YR   , data = inputsFit01,family = "binomial") # Have you visited a doctor with a fever in the last 2 years?
  modelB.10 <- glm(Dconvert ~ HH_D   , data = inputsFit01,family = "binomial") # Have any household members visited a doctor with a fever in the last 2 years?
  
  model.list=list(modelB.1,modelB.1a,modelB.1b,modelB.2,modelB.3,modelB.4,modelB.5,modelB.6,modelB.7,modelB.8,modelB.9,modelB.10)
  
  names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")
  
  data.tally = inputsFit01[,names.model]
  
  # Calculate odds ratios - iterate over different models
  store.oddsTable = NULL
  for(jj in 1:length(model.list)){
    modelT = model.list[[jj]]
    estm.D = exp(coef(modelT))
    conf.D = exp(confint(modelT)  )
    
    store.odds=names.model[jj]
    
    for(ii in 2){
      
      store.odds = cbind(store.odds,
                         sum(data.tally[jj]==1),
                         paste(round(estm.D[ii],2)," (",paste(round(conf.D[ii,],2),collapse="-"),")",sep="") ,
                         round(coef(summary(modelT))[,4][2],2))
    }
    
    store.oddsTable = rbind(store.oddsTable,store.odds )
    
  }
  
  
  store.oddsTable
  
  write.csv(store.oddsTable,"plots/Table_3_logistic_outputs.csv")
  
  
  # - - - - - 
  # Calculate rise in titre -- adjusted for baseline

  modelB.1 <- glm(RiseT ~ FEVER_2YR + ELISA1  , data = inputsFit0,family = "gaussian", weights = titre.weights)
  modelB.2 <- glm(RiseT ~ DOC_2YR + ELISA1 , data = inputsFit0,family = "gaussian", weights = titre.weights)

  model.list=list(modelB.1,modelB.2)
  
  names.model=c("FEV","DOC")
  
  # Calculate parameters - iterate over different models
  store.paramTable = NULL
  for(jj in 1:length(model.list)){
    modelT = model.list[[jj]]
    estm.D = coef(modelT)
    conf.D = confint(modelT)  
    
    store.param=names.model[jj]
    
    for(ii in 2){
      
      store.param = cbind(store.param,paste(signif(estm.D[ii],3)," (",paste(signif(conf.D[ii,],3),collapse="-"),")",sep=""),signif(coef(summary(modelT))[,4][2],3))
    }
    
    store.paramTable = rbind(store.paramTable,store.param )
    
  }

  store.paramTable
  colnames(store.paramTable) = c("variable","odds","p")
  write.csv(store.paramTable,"plots/Result_rise_vs_fever.csv")

  # Summary of seroconversion and fever visits
  store.fever = matrix(nrow=4,ncol=4) %>% data.frame()
  names(store.fever) = c("report","percent","report_inf","percent_inf")
                         
  total_part = length(inputsFit01$FEVER_2YR)
  store.fever$report = c(
    sum((inputsFit01$FEVER_2YR==0 )),
    sum((inputsFit01$FEVER_2YR==1 & inputsFit01$DOC_2YR==0)),
    sum((inputsFit01$DOC_2YR==1 )),
    total_part
  )
  
  store.fever$percent = sapply(store.fever$report, function(x){binom_output(x,total_part) } )
  
  # Adjusting for infection:
  
  total_inf = sum(inputsFit01$Dconvert) %>% round()
  store.fever$report_inf = c(
    sum((inputsFit01$FEVER_2YR==0 & inputsFit01$Dconvert)),
    sum((inputsFit01$FEVER_2YR==1 & inputsFit01$DOC_2YR==0 & inputsFit01$Dconvert)),
    sum((inputsFit01$DOC_2YR==1 & inputsFit01$Dconvert)),
    total_inf
  ) %>% round()
  
  store.fever$percent_inf = sapply(store.fever$report_inf, function(x){binom_output(x,total_inf) } )

  write.csv(store.fever,"plots/Table_S1_fever_breakdown.csv")
  
  
}


# - - - - - - - - - - - - - - - - - - - - 
# Distributional bootstrap extraction - remove ELISA noise
# - - - - - - - - - - - - - - - - - - - - 

remove_noise <- function(ELISA.fit=T){ # use inputs from above fit_ELISA function
  
  # Load serology data
  loadinputs = fit_ELISA(ELISA.fit)
  titre.weights = loadinputs$titrechange
  
  dcutoff = read.csv("data/serology_MIA_cutoff.csv") # Load MIA cutoffs
  
  inputs0 = loadinputs$tabinput
  
  if(ELISA.fit ==T){
    change.titre = inputs0$ELISA2-inputs0$ELISA1
    change.titre00 = change.titre[change.titre > -9] # Remove outliers
  }else{
    denv3.mia.1 = inputs0$DENV3/1000 #+dcutoff$DENV3
    denv3.mia.2 = inputs0$DENV3B/1000 #+dcutoff$DENV3
    change.titre = denv3.mia.2 - denv3.mia.1
    change.titre00 = change.titre[!is.na(change.titre)]
    change.titre00 = change.titre00[change.titre00>-5]
  }
  

  param = c(a0 = 0, b0=1, a1=2, b1=3, lambda1 = 1)
  result2 = optim(param, fitmixture2, method="L-BFGS-B", val=change.titre00,lower=c(rep(0,4),-10), hessian=FALSE) #, control=list(trace=1))
  parfit = result2$par #param #
  
  # - - - 
  # Bootstrap values according to the fitted distribution
  
  # Define probabilities of signal vs noise
  x.titre = seq(-10,30,1)
  y.error = (1-transform.lambda(parfit[["lambda1"]]))*dnorm(x.titre, mean = 0, sd = parfit[["b0"]])
  y.signl1 = transform.lambda(parfit[["lambda1"]])*dgamma(x.titre, shape = parfit[["a1"]], scale = parfit[["b1"]]) 
  y.signl2 = 0 
  p.signal = (y.signl1+y.signl2)/(y.error+y.signl1+y.signl2)
  
  # Plot titre dependent rise
  btstrap=1e2
  l.sample = length(inputs0$AGE_2015)
  x.rise.gam = seq(min(inputs0$ELISA1),max(inputs0$ELISA1),0.5)
  store.Y = NULL;   store.Y.lowess = NULL

  # - - - 
  # PLOT FIGURE WITH RISES
  
  par(mfrow=c(1,3),mar=c(4,4,2,2),mgp=c(1.7,0.5,0))
  
  #Plot Fig A
  hist(inputs0$ELISA1,freq=F,breaks=seq(-0.5,35.5,1),col=rgb(1,0.5,0,1),border="white",main="",xlab="DENV IgG ELISA value",ylab="proportion of samples",ylim=c(0,0.1),xlim=c(0,35),xaxs="i",yaxs="i")
  hist(inputs0$ELISA2,freq=F,breaks=seq(-0.5,35.5,1),add=T,border="white",col=rgb(0,0,1,0.5))
  txtSize=1
  text(x=25,y=0.09,adj=0,labels="2013",col=rgb(1,0.5,0),cex=txtSize)
  text(x=25,y=0.08,adj=0,labels="2015",col=rgb(0,0,1),cex=txtSize)
  
  lines(c(9,9),c(0,1),lty=2)
  lines(c(11,11),c(0,1),lty=2)

  title(adj=0,main=LETTERS[1])
  
  #Plot full distribution Fig B
  xx <- seq(-10,20,0.1)
  change.titreH = change.titre00
  hist(change.titreH,freq=F,breaks=seq(-6.5,20.5,1),col="light grey",main="",border="white",
       xlab="change in DENV IgG ELISA value",ylab="proportion of samples",ylim=c(0,0.25),xlim=c(-5,21),xaxs="i",yaxs="i")
  
  all.model = transform.lambda(parfit[["lambda1"]])*dgamma(xx, shape = parfit[["a1"]], scale = parfit[["b1"]]) + (1-transform.lambda(parfit[["lambda1"]]))*dnorm(xx, mean = 0, sd = parfit[["b0"]])
  lines(xx,all.model, col="black",lwd=1) # Simulation testing
  lines(xx,(1-transform.lambda(parfit[["lambda1"]]))*dnorm(xx, mean = 0, sd = parfit[["b0"]]), col=rgb(0.4,0.4,0.4),lwd=2) # Simulation testing
  lines(xx,transform.lambda(parfit[["lambda1"]])*dgamma(xx, shape = parfit[["a1"]], scale = parfit[["b1"]]), col=rgb(0.2,0.2,1),lwd=2) # Simulation testing
  title(adj=0,main=LETTERS[2])
  
  par(new=TRUE)
  plot(x.titre-0.5,p.signal,ylim=c(0,1),type="l",lty=2,xlim=c(-5.5,20.5),col="blue",xlab="",frame.plot=F,ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i") #,yaxt="n",ylab="",xlab="",pch=19)
  #lines(c(1:length(table.attack0$N))-shiftA,table.attack0$cases_P,col=rgb(1,0,0,0.2))
  
  axis(4,col="blue",col.axis="blue")
  mtext("probabilty infected", side=4, line=1.5,col="blue",cex=0.7) # Label for 2nd axis
  
  # Figure C
  plot(inputs0$ELISA1,inputs0$RiseT,pch=19,cex=0.7,ylim=c(0,21),xlim=c(0,30),
       col="white",xaxs="i",yaxs="i",xlab="DENV IgG ELISA value 2013",ylab="change in value")

  for(ii in 1:btstrap){
    # Bootstrap sample from data
    rand.pick = sample(x=c(1:l.sample),size=l.sample,replace = T)
    inputs01 = inputs0[rand.pick,]
    
    # Remove noise
    change.titre0 = inputs01$ELISA2-inputs01$ELISA1
    titre.weights0 = sapply(change.titre0,function(x){ p.signal[match(x,x.titre)]   })
    inputs1 = inputs01 #[randpick,]
    
    # Fit model and plot
    modelB.P <- gam(RiseT ~ s(ELISA1) , data = inputs1,family = "gaussian",weights = titre.weights0) 
    
    y.predict = predict(modelB.P, list(ELISA1=x.rise.gam), type = "link", se.fit = TRUE)
    
    store.Y = rbind(store.Y,as.numeric(y.predict$fit))
    
    randpick = runif(length(titre.weights0)) < titre.weights0 
    inputsRM = inputs01[randpick,]
    if(ii<100){
      points(inputsRM$ELISA1,inputsRM$RiseT,pch=19,cex=0.7,col=rgb(0,0,0,0.01))
    }
    
  }
  
  change.titre0 = inputs0$ELISA2-inputs0$ELISA1
  titre.weights0 = sapply(change.titre0,function(x){ p.signal[match(x,x.titre)]   })
  modelB.P <- gam(RiseT ~ s(ELISA1) , data = inputs0,family = "gaussian",weights = titre.weights0) 
  
  preds <- predict(modelB.P, newdata = list(ELISA1=x.rise.gam), type = "link", se.fit = TRUE)
  critval <- 1.96; upperCI <- preds$fit + (critval * preds$se.fit); lowerCI <- preds$fit - (critval * preds$se.fit)
  fit <- preds$fit
  fitPlotF <- modelB.P$family$linkinv(fit); CI1plotF <- modelB.P$family$linkinv(upperCI);  CI2plotF <- modelB.P$family$linkinv(lowerCI)
  
  polygon(c(x.rise.gam,rev(x.rise.gam)),c(CI1plotF,rev(CI2plotF)),col=rgb(0,0,1,0.2),lty=0)
  #polygon(c(x.rise.gam,rev(x.rise.gam)),c(apply(store.Y,2,c_50_1),rev(apply(store.Y,2,c_50_2))),col=rgb(0,0,1,0.2),lty=0)
  lines(x.rise.gam, fitPlotF ,col="blue",lwd=2)
  title(adj=0,main=LETTERS[3])

  dev.copy(pdf,paste("plots/Figure_2_ELISA_analysis.pdf",sep=""),width=7,height=2.5,useDingbats=FALSE)
  dev.off()
  
  
  # - - - - - - - - - - - - - - - - - - - - 
  # Calculate ELISA results by age, adjusting for primary and secondary

  change.titre = inputs0$ELISA2-inputs0$ELISA1
  titre.weightsAll = sapply(change.titre,function(x){ p.signal[match(x,x.titre)]   })

  age.LineListS = data.frame(read.csv("data/agedistribution_suspected.csv",stringsAsFactors = F))[,2]
  age.LineListP = data.frame(read.csv("data/agedistribution_positive.csv",stringsAsFactors = F))[,2]
  
  # Import population age distribution
  age.distn_central = read.csv("data/age_distn_central.csv",stringsAsFactors = F) %>% data.frame()
  
  under20report = (sum(age.LineListS<20)/length(age.LineListS))

  age.bins = c(0,10,20,30,50,100)
  age.bins = seq(0,80,10); age.bins[length(age.bins)]=90
  store.attack = NULL
  table.attack = NULL
  age.names = NULL
  inputs = inputs0
  
  for(ii in (1:(length(age.bins)-1))){
    
    age.ID = inputs0$AGE_2015>=age.bins[ii] & inputs0$AGE_2015<age.bins[ii+1]
    age.ID2013 = inputs0$AGE_2015>=(age.bins[ii]+2) & inputs0$AGE_2015<(age.bins[ii+1]+2)
    
    casesP.LL = sum(age.LineListP>=age.bins[ii] & age.LineListP<age.bins[ii+1])
    casesS.LL = sum(age.LineListS>=age.bins[ii] & age.LineListS<age.bins[ii+1])
    popn.tot = sum(age.distn_central[age.distn_central$age>=age.bins[ii] & age.distn_central$age<age.bins[ii+1],"n"] ) 
    
    # Baseline ELISA
    elisa1d = (inputs$sero1[age.ID]) %>% sum()
    elisa2d = (inputs$sero2[age.ID]) %>% sum()
    lumd3_1 = as.numeric(inputs$DENV3P[age.ID]) ; ltot = sum(!is.na(lumd3_1)); lumd3_1 = sum(lumd3_1[!is.na(lumd3_1)])
    lumd3_2 = as.numeric(inputs$DENV3BP[age.ID]) ; lumd3_2 = sum(lumd3_2[!is.na(lumd3_2)])
    
    # Use distributions to calculate
    pick.ages = inputs0[age.ID,"AGE_2015"]
    pick.signal.prob = titre.weightsAll[age.ID]
    pick.signal.Sndy = 0*titre.weightsAll[age.ID] #deprecated
    
    n.signl.All = sum(pick.signal.prob)
    n.signl.Sec = sum(pick.signal.Sndy)
    n.total = length(pick.signal.prob)
    
    # Naive attack rate using seroconversion
    sero.neg2013 = sum(inputs0[age.ID,"sero1"]==0)
    sero.convert = sum(inputs0[age.ID,"sero1"]==0 & inputs0[age.ID,"sero2"]==1)

    # Store data
    bc1=binom.calc(round((n.signl.All-n.signl.Sec)),n.total)
    bc2=binom.calc(round(n.signl.Sec),n.total)
    if(sero.neg2013>0){
      bc3=binom.calc(sero.convert,sero.neg2013)
    }else{
      bc3=0
    }

    table.attack = rbind(table.attack,c(age.bins[ii], age.bins[ii+1], n.total, sum(age.ID2013),popn.tot, 1000*casesP.LL/popn.tot, 1000*casesS.LL/popn.tot, n.signl.All/n.total, (n.signl.All-n.signl.Sec)/n.total, n.signl.Sec/n.total, sero.neg2013, sero.convert,elisa1d,elisa2d,lumd3_1,lumd3_2,ltot ))
    store.attack = rbind(store.attack,c(age.bins[ii], age.bins[ii+1], n.total, round(n.signl.All), bc1,sero.neg2013 ,sero.convert,bc3 ))

    age.names = c(age.names,paste(age.bins[ii], (age.bins[ii+1]-1),sep="-"))
  }
  
  table.attack = table.attack %>% data.frame(stringsAsFactors=F)
  names(table.attack) = c("age1","age2","N","N2013","pop","cases_P","cases_S","All_inf","Inf1","Inf2","seroneg","serocon","ELISA1","ELISA2","D3_1","D3_2","L_tot")
  table.attack0 = table.attack #[table.attack$N>1,]
  
  store.attack = rbind(store.attack,c("Total","Total", sum(table.attack$N),round(sum(table.attack$N*table.attack$Inf1)), binom.calc(round(sum(table.attack$N*table.attack$Inf1)),sum(table.attack$N)),sum(table.attack$seroneg),sum(table.attack$serocon),binom.calc(sum(table.attack$serocon),sum(table.attack$seroneg) )))
  
  store.attack = store.attack %>% data.frame(stringsAsFactors=F)
  names(store.attack) = c("age1","age2","N","Infect","Inf%","seroneg","serocon","sero%")
  
  write.csv(store.attack,"plots/Table_4_estimate_attack.csv")
  #write.csv(table.attack0,"plots/seroprevalance_data.csv")

  # - - - - - - - -
  # Plot ELISA and luminex distribution
  shiftA=0.07
  par(mfrow=c(2,1),mar=c(4,4,1,4),mgp=c(2,0.7,0))
  col1A =rgb(0.3,0.5,1);col1B=rgb(0,0.2,0.4); col2A= rgb(0.2,1,0.4); col2B= rgb(0,0.4,0.1)
  plot(table.attack0$ELISA1,col="white",ylim=c(0,1),ylab="proportion seropositive",xlab="age group",type="l",xaxt="n",xlim=c(0.8,length(table.attack0$N)+0.2))
  
  for(ii in 1:length(table.attack0$N)){
    htest1 = binom.test(round(table.attack0[ii,"ELISA1"] ),table.attack0[ii,"N"])$conf.int
    htest2 = binom.test(round(table.attack0[ii,"ELISA2"] ),table.attack0[ii,"N"])$conf.int
    mid1 = table.attack0[ii,"ELISA1"]/table.attack0[ii,"N"];    mid2 = table.attack0[ii,"ELISA2"]/table.attack0[ii,"N"]
    points(ii-shiftA,mid1,col=col1A,pch=15); lines(c(ii,ii)-shiftA,c(htest1[1],htest1[2]),col=col1A)
    points(ii+shiftA,mid2,col=col1B,pch=15); lines(c(ii,ii)+shiftA,c(htest2[1],htest2[2]),col=col1B )
    htest1 = binom.test(round(table.attack0[ii,"D3_1"] ),table.attack0[ii,"L_tot"])$conf.int
    htest2 = binom.test(round(table.attack0[ii,"D3_2"] ),table.attack0[ii,"L_tot"])$conf.int
    mid1 = table.attack0[ii,"D3_1"]/table.attack0[ii,"L_tot"];    mid2 = table.attack0[ii,"D3_2"]/table.attack0[ii,"N"]
    points(ii+2*shiftA,mid1,col=col2A,pch=19); lines(c(ii,ii)+2*shiftA,c(htest1[1],htest1[2]),col=col2A)
    points(ii+3*shiftA,mid2,col=col2B,pch=19); lines(c(ii,ii)+3*shiftA,c(htest2[1],htest2[2]),col=col2B )
  }
  axis(side = 1, at = c(1:length(table.attack0$N)),labels=age.names); title(main="A",adj=0)
  
  shiftY=0.07; maxY=1.05; txtSize=0.8
  points(0.8,maxY-shiftY,pch=15,col=col1A); text(x=0.9,y=maxY-shiftY,adj=0,labels="ELISA (2013)",col=col1A,cex=txtSize)
  points(0.8,maxY-2*shiftY,pch=15,col=col1B); text(x=0.9,y=maxY-2*shiftY,adj=0,labels="ELISA (2015)",col=col1B,cex=txtSize)
  points(0.8,maxY-3*shiftY,pch=19,col=col2A); text(x=0.9,y=maxY-3*shiftY,adj=0,labels="DENV-3 (2013)",col=col2A,cex=txtSize)
  points(0.8,maxY-4*shiftY,pch=19,col=col2B); text(x=0.9,y=maxY-4*shiftY,adj=0,labels="DENV-3 (2015)",col=col2B,cex=txtSize)

  # - - - - - - - -
  # Plot data 
  shiftA=0.05

  # Plot attack rates
  plot(table.attack0$Inf1,col="white",ylim=c(0,0.7),ylab="estimated proportion infected",xlab="age group",type="l",xaxt="n",xlim=c(0.8,length(table.attack0$N)+0.2))

  # Plot CIs
  for(ii in 1:length(table.attack0$N)){
    htest1 = binom.test(round(table.attack0[ii,"Inf1"] * table.attack0[ii,"N"]),table.attack0[ii,"N"])$conf.int
    htest2 = binom.test(round(table.attack0[ii,"Inf2"] * table.attack0[ii,"N"]),table.attack0[ii,"N"])$conf.int
    htestAll = binom.test(round(table.attack0[ii,"All_inf"] * table.attack0[ii,"N"]),table.attack0[ii,"N"])$conf.int
    mid1 = table.attack0[ii,"Inf1"];    mid2 = table.attack0[ii,"Inf2"];    midAll = table.attack0[ii,"All_inf"]
    points(ii+shiftA,mid1,col="black",pch=19); lines(c(ii,ii)+shiftA,c(htest1[1],htest1[2]),col="black")
    #points(ii+shiftA,mid2,col=rgb(0.4,0.7,1),pch=19); lines(c(ii,ii)+shiftA,c(htest2[1],htest2[2]),col=rgb(0.4,0.7,1) )
    #points(ii-shiftA,midAll,col="black",pch=19); lines(c(ii,ii)-shiftA,c(htestAll[1],htestAll[2]),col="black")
    
  }
  
  axis(side = 1, at = c(1:length(table.attack0$N)),labels=age.names); title(main="B",adj=0)

  par(new=TRUE)
  plot(c(1:length(table.attack0$N))-shiftA,table.attack0$cases_P,ylim=c(0,20),col="red",xlab="",ylab="",pch=19,xaxt="n",yaxt="n",xlim=c(0.8,length(table.attack0$N)+0.2)) #,yaxt="n",ylab="",xlab="",pch=19)
  #lines(c(1:length(table.attack0$N))-shiftA,table.attack0$cases_P,col=rgb(1,0,0,0.2))
  
  for(ii in 1:length(table.attack0$N)){
    htest1 = 1000*binom.test(round( table.attack0$pop[ii]*table.attack0$cases_P[ii]/1000) ,table.attack0$pop[ii] )$conf.int
    lines(c(ii,ii)-shiftA,c(htest1[1],htest1[2]),col="red")
  }
  
  axis(4,col="red",col.axis="red")
  mtext("cases per 1,000", side=4, line=2,col="red",cex=1) # Label for 2nd axis
  
  dev.copy(pdf,paste("plots/Figure_3_age_attack.pdf",sep=""),width=7,height=7,useDingbats=FALSE)
  dev.off()

}


# Function for binomial calculations

binom.calc <- function(x,n){

    htest <- binom.test(x,n, p = 1,conf.level=0.95)
    meanA=100*htest$estimate %>% as.numeric()  %>% signif(digits=3)
    conf1=100*htest$conf.int[1] %>% signif(digits=3)
    conf2=100*htest$conf.int[2] %>% signif(digits=3)
    paste(meanA,"% (",conf1,"-",conf2,"%)",sep="") 

}



# - - - - - - - - - - - - - - - - - - - - 
# Plot surveillance data
# - - - - - - - - - - - - - - - - - - - - 

plot_surveillance_data <- function(){


  time.seriesPOS = read.csv(paste("data/Central_positive.csv",sep=""), stringsAsFactors = F); time.seriesPOS$date = as.Date(time.seriesPOS$date) # Load confirmed data
  time.seriesD = read.csv(paste("data/Central_lab_tested.csv",sep=""), stringsAsFactors = F); time.seriesD$date = as.Date(time.seriesD$date) # Load lab tested dengue data
  time.series = read.csv(paste("data/Central_DLI.csv",sep=""), stringsAsFactors = F); time.series$date = as.Date(time.series$date) # Load DLI dengue data

  colourpick1 = c(rgb(1, 0.8, 0.1),rgb(0.3, 0.8, 1),rgb(0.1, 1, 0.3), # swap colours for West/North
                   rgb(1,0,0.2),rgb(0,0.2,1),rgb(0,0.6,0),rgb(1,0.5,0))
  
  div.names = c("Central","Western","North","East")
  minT=as.Date("2013-10-28"); maxT=as.Date("2014-08-31") 
  
  # PLOT 2013/14
  
  par(mfrow=c(2,2),mar=c(2,3,1,3),mgp=c(2,0.7,0))

  # PLOT lab confirmed data

  plot(time.series$date,rowSums(time.seriesPOS[,div.names]),col="white",type="l",lwd=1,xlim=c(minT,maxT),xaxs="i",yaxs="i",ylim=c(0,1100),xlab="",ylab="confirmed cases")
  
  for(ii in 1:4){
    lines(time.seriesD$date,time.seriesPOS[,div.names[ii]],col=colourpick1[ii],lwd=2)
    text(min(time.seriesD$date)+220,1100-100*ii,labels=div.names[ii],adj=0,col=colourpick1[ii])
  }
  text(min(time.seriesD$date)+220,1100-100*5,labels="All",adj=0)
  
  lines(time.series$date,rowSums(time.seriesPOS[,div.names]),col="black",lwd=1)
  title(main="A",adj=0)
  
  # PLOT Proportion positive
  
  plot(time.series$date,rowSums(time.seriesPOS[,div.names]),col="white",type="l",lwd=1,xlim=c(minT,maxT),xaxs="i",yaxs="i",ylim=c(0,1100),xlab="",ylab="lab cases in Central Division")
  
  lines(time.seriesD$date,time.seriesPOS[,div.names[1]],col=colourpick1[1],lwd=2)
  lines(time.seriesD$date,time.seriesD[,div.names[1]],col=colourpick1[1],lwd=2,lty=2)
  
  par(new=TRUE)
  plot(time.series$date, time.seriesPOS[,div.names[1]]/time.seriesD[,div.names[1]],col="grey",lwd=1,ylim=c(0,1.05),xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(minT,maxT),type="l",yaxs="i")
  axis(4,col="grey",col.axis="grey")
  mtext("proportion positive", side=4, line=1.5,col="grey",cex=0.8) # Label for 2nd axis
  
  title(main="B",adj=0)
  
  # PLOT DLI
  
  # Add non-confirmed
  #time.series[,div.names] =   time.series[,div.names] + time.seriesD[,div.names] -time.seriesPOS[,div.names]
  
  plot(time.series$date,rowSums(time.series[,div.names]),col="white",type="l",lwd=2,xlim=c(minT,maxT),xaxs="i",yaxs="i",ylim=c(0,1800),xlab="",ylab="DLI cases")
  
  for(ii in 1:4){
    lines(time.seriesD$date,time.series[,div.names[ii]],col=colourpick1[ii],lwd=2)
  }
  
  lines(time.series$date,rowSums(time.series[,div.names]),col="black",lwd=1)
  title(main="C",adj=0)
  
  # PLOT All suspected cases
  
  plot(time.series$date,rowSums(time.seriesD[,div.names])+rowSums(time.series[,div.names]),col="white",type="l",lwd=2,xlim=c(minT,maxT),xaxs="i",yaxs="i",ylim=c(0,2500),xlab="",ylab="all cases")

  # Tally cases
  print(c(colSums(time.seriesD[,div.names])+colSums(time.series[,div.names]),sum(time.seriesD[,div.names])+sum(time.series[,div.names])) )
  
  # Count up cases
  tallyC = NULL
  for(ii in 1:4){
    lines(time.seriesD$date,time.series[,div.names[ii]]+time.seriesD[,div.names[ii]],col=colourpick1[ii],lwd=2)

    time.seriesP = time.series[time.series$date<maxT,]
    tallyC = c(tallyC,sum(time.seriesP[,div.names[ii]]+time.seriesP[,div.names[ii]]))
  }
  tallyC
  
  lines(time.series$date,rowSums(time.series[,div.names])+rowSums(time.seriesD[,div.names]),col="black",lwd=1)
  title(main="D",adj=0)
  
  dev.copy(pdf,paste("plots/Figure_S1_timeseries_dengue_Geographic.pdf",sep=""),width=7,height=5)
  dev.off()
  
  # - - - - - - - - 
  # Plot ratio of reporting
  
  # par(mfrow=c(1,1),mar=c(2,3,1,3),mgp=c(2,0.7,0))
  # 
  # y.vals = time.seriesD[,div.names[1]] #time.seriesPOS[,div.names[1]]
  # y.vals2 = time.series[,div.names[1]] #+ time.series[,div.names[1]] - y.vals
  # 
  # plot(time.seriesD$date,y.vals/(y.vals+y.vals2),xlim=c(as.Date("2013-11-01"),as.Date("2014-07-01")),ylim=c(0,1),pch=19,ylab="confirmed cases/all suspected cases")
  # lines(time.seriesD$date,y.vals/(y.vals+y.vals2),xlim=c(as.Date("2013-11-01"),as.Date("2014-08-01")),ylim=c(0,1),pch=19)
  
  dev.copy(pdf,paste("plots/Figure_S1_timeseries_dengue_GeographicA.pdf",sep=""),width=7,height=5)
  dev.off()
  
  
}


