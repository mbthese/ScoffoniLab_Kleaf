require(likelihood)

#Define functions and parameter values
Exponential <-function (A,B,psi){A*exp(-B*psi)}
define_parsE <- function(input_df){
  parsE <- list(A=max(input_df$gs),B=1, sd=2)
  par_loE =list(A=0, B=0.1,sd=0.0005) 
  par_highE =list(A=max(input_df$gs)*2,B=10,sd=50)
  return(list(parsE=parsE, par_loE=par_loE,par_highE=par_highE))
} 
# The high values for all of these functions are sort of arbitrary. I set the original limits thinking that fitted gmax wouldn't be too much higher than measured gmax, and that P50 should be between 0 and -5 MPa, and then just ran the script a few times with different values to see if I could improve the number of curves it successfully fit. Using different par values for each model type helped a lot. 
##In the exponential function, the effect of changing Xo (the intercept) on gs is usually very small compared to that of any of the other parameters, so it is hard to fit with any kind of precision. This makes it hard to estimate error terms for Xo, and increases the error in the other parameters, since the model can't tell whether a given improvement in fit came from changing Xo or another parameter. I had very little success fitting Xo without NAs in the error terms, but the model fit really well when I removed Xo from the equation. I was concerned we wouldn't capture the tail of the function when gs is close to 0, but the graphs of the exponential fits look fine. I think this is a reasonable simplification and improves model fit, but if you think Xo is important, then you might have to fit it one curve at a time. 

Logistic<-function (A,B,Xo,psi){A/(1+((psi/Xo)^B))}
define_parsL <- function(input_df){
  parsL <- list(A=max(input_df$gs),B=5,Xo=2, sd=2)
  par_loL =list(A=0, B=0.1, Xo=-5, sd=0.0005) 
  par_highL=list(A=max(input_df$gs)*2, B=25, Xo=5, sd=50)
  return(list(parsL=parsL, par_loL=par_loL,par_highL=par_highL))
} 

Sigmoidal <-function (A,B,Xo,psi){A/(1+exp(-((psi-Xo)/B)))}
define_parsS <- function(input_df){
  parsS <- list(A=max(input_df$gs),B=-0.5,Xo=2, sd=2)
  par_loS =list(A=0, B=-1.25, Xo=-5, sd=0.0005) 
  par_highS=list(A=max(input_df$gs)*1.5, B=0, Xo=6, sd=50)
  return(list(parsS=parsS, par_loS=par_loS,par_highS=par_highS))
} 

Exponential2 <-function (A,B,C,psi){C+A*exp(-B*psi)}
define_parsE2 <- function(input_df){
  parsE2 <- list(A=max(input_df$gs),B=1,C=.1,sd=2)
  par_loE2 =list(A=0, B=0.1,C=.001, sd=0.0005) 
  par_highE2 =list(A=max(input_df$gs)*2,B=20,C=20,sd=50)
  return(list(parsE2=parsE2, par_loE2=par_loE2,par_highE2=par_highE2))
} 

Linear <-function (A,B,psi){B*(psi)+A}
#I changed the order of A and B so that A was still gmax

# Declare the data frame to store the results
modelfitting_results<-data.frame(Species= as.character(), data.type=as.character(),  A = numeric(), B = numeric(), C = numeric(), D = numeric(),loglikeli = numeric(), rsq = numeric(), slope = numeric(), AIC = numeric(), AICcorr = numeric(), sterror1 = numeric(), sterror2 = numeric(), sterror3 = numeric(), sterror4 = numeric(), N=numeric())


do_the_thing_nonlinear <- function(input_df, model_type, pars1, par_lo1, par_hi1) {
  model_type = model_type
  var <- list(psi="psi",
              x="gs",
              mean="predicted", 
              log=TRUE)
  #It looks weird that gs is the x variable here, but anneal calculates the slope and R2 of the fit using the predicted gs values as the y and the observed gs values as the x 			 
  pars = pars1
  par_lo = par_lo1
  par_hi= par_hi1
  
  res<-anneal(model = model_type, par= pars1, source_data = input_df, var = var, par_lo=par_lo1, par_hi= par_hi1, dep_var = "gs", pdf = dnorm, max_iter=50000, show_display=F, temp_red=0.05)
  #Setting the parameters to change slowly in the fitting procedure (the temp_red variable) helped a lot. You can watch the fitting proceed with show_display, but I've never found it very informative
  
  #AIC formula: -2LL + 2*parameters (incl nuisance, i.e.,sd)
  
  AIC<- res$aic
  
  #AICcorr formula: -2LL + (2*n*parameters (incl nuisance, i.e.,sd)/(n-parameters-1))
  
  AICcorr <- res$aic_corr
  
  slope <-sum(res$source_data$predicted*res$source_data$gs)/sum(res$source_data$predicted^2)
  
  rsq <- res$R2
  
  sterror<- res$std_errs
  
  N <- length(res$source_data$gs)
  
  
  parvecLog<-c(Species= paste(input_df[1,1]), data.type=paste(input_df[1,2]),A = res$best_pars[1], B = res$best_pars[2], C = res$best_pars[3], D = res$best_pars[4],loglikeli = res$max_likeli, rsq = rsq, slope = slope, AIC = AIC, AICcorr = AICcorr, sterror_1 = sterror[1], sterror = sterror[2], sterror = sterror[3], sterror = sterror[4], N =N)
  
  #Plot the fit
  plot(res$source_data$psi, res$source_data$gs, xlab = "Water Potential (-MPa)", ylab = "Stomatal Conductance (mmol m-2 s-1)")
  cbind(res$source_data$psi, res$source_data$predicted)-> for_plotting
  for_plotting[order(for_plotting[,1]),]-> for_plotting
  lines(for_plotting[,1], for_plotting[,2], col="blue")
  title(paste(input_df[1,1], input_df[1,2]))
  
  return(parvecLog)
  
}

do_the_thing_linear <- function(input_df, model_type) {
  model_type = model_type
  var <- list(psi="psi",
              x="gs",
              mean="predicted", 
              log=TRUE)
  
  lm(input_df$gs ~ input_df$psi)-> fita #a normal linear regression gave the best starting parameters
  pars = list(A=summary(fita)$coeff[1,1], B =summary(fita)$coeff[2,1], sd = 1)
  par_lo1= list(A=summary(fita)$coeff[1,1]*0.05, B =summary(fita)$coeff[2,1]*2, sd =0.005)
  par_hi1= list(A=summary(fita)$coeff[1,1]*4, B =summary(fita)$coeff[2,1]*0.1, sd =20)
  #########Update from Megan: I decreased the sd of the par_lo1, linear fits now match the lm results
  
  res<-anneal(model = model_type, par= pars, source_data = input_df, var = var, par_lo=par_lo1, par_hi= par_hi1, dep_var = "gs", pdf = dnorm, max_iter=50000, show_display=F, temp_red=0.05)
  
  
  #AIC formula: -2LL + 2*parameters (incl nuisance, i.e.,sd)
  
  AIC<- res$aic
  
  #AICcorr formula: -2LL + (2*n*parameters (incl nuisance, i.e.,sd)/(n-parameters-1))
  
  AICcorr <- res$aic_corr
  
  slope <-sum(res$source_data$predicted*res$source_data$gs)/sum(res$source_data$predicted^2)
  
  rsq <- res$R2
  
  sterror<- res$std_errs
  
  N <- length(res$source_data$gs)

  
  parvecLog<-c(Species= paste(input_df[1,1]), data.type=paste(input_df[1,2]),A = res$best_pars[1], B = res$best_pars[2], C = res$best_pars[3], D = res$best_pars[4],loglikeli = res$max_likeli, rsq = rsq, slope = slope, AIC = AIC, AICcorr = AICcorr, sterror_1 = sterror[1], sterror = sterror[2], sterror = sterror[3], sterror = sterror[4], N =N)
  
  plot(res$source_data$psi, res$source_data$gs, xlab = "Water Potential (-MPa)", ylab = "Stomatal Conductance (mmol m-2 s-1)")
  cbind(res$source_data$psi, res$source_data$predicted)-> for_plotting
  for_plotting[order(for_plotting[,1]),]-> for_plotting
  lines(for_plotting[,1], for_plotting[,2], col="blue")
  title(paste(input_df[1,1], input_df[1,2]))
  return(parvecLog)
  
}

