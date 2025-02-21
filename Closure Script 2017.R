source("R function MB 2017.R")
dd = read.csv("result/AVSA/AVSA_curve.csv",header=TRUE)

dd$to_split = paste(dd$species, dd$data.type)
aggregate(dd, by =list(dd$to_split), FUN=length) -> speciesN
speciesN[which(speciesN$to_split > 5),1] -> species_list
#Trying to fit curves with less than 6 points returns -inf for the AICcor. Are those the binned data? Would it be possible to make the bins smaller so there are more points?

# first, just to make sure that things work as expected, let's print out
# the kleaf column for each species
for (ii in 1:length(species_list)){
  subset(dd, dd$to_split == species_list[ii], select =c(1:5)) -> data_by_sp
  print(data_by_sp)
}
#I couldn't get the list format to save as a usable .csv file, so I changed everything to dataframes

# use these kleaf's as inputs to the likelihood function
#Linear
pdf(file = 'result/AVSA/Linear_fits.pdf', width =8.5, height=5, onefile=T)
#Open a pdf file to save fitted curve plots in
for (ii in 1:length(species_list)){
  subset(dd, dd$to_split == species_list[ii], select =c(1:5)) -> data_by_sp
  Linear_fits = do_the_thing_linear(data_by_sp, Linear)
  Linear_fits$D.NA<- NA # placeholder to keep all of the data frames the same size
  Linear_fits$sterror.NA <- NA
  rbind(modelfitting_results, as.data.frame(Linear_fits))-> modelfitting_results
  print(ii)
}
#For every species, fit the curve and record the results
dev.off()
modelfitting_results$Kmax<- modelfitting_results$A.A
modelfitting_results$psi_kleaf20<- ((0.8-1)*(-modelfitting_results$A.A)/(modelfitting_results$B.B))*-1
modelfitting_results$psi_kleaf50<- ((0.5-1)*(-modelfitting_results$A.A)/(modelfitting_results$B.B))*-1
modelfitting_results$psi_kleaf80<- ((0.2-1)*(-modelfitting_results$A.A)/(modelfitting_results$B.B))*-1
modelfitting_results$psi_kleaf95<- ((0.05-1)*(-modelfitting_results$A.A)/(modelfitting_results$B.B))*-1
#Every model has different formulas for kleaf50 and kleaf95
which(is.na(as.numeric(modelfitting_results[,12]))=='TRUE' | is.na(as.numeric(modelfitting_results[,13]))=='TRUE'| is.na(as.numeric(modelfitting_results[,14]))=='TRUE') -> problem_rows
#Look for curves that returns NAs for the parameter error values
write.csv(modelfitting_results, file = "result/AVSA/Linear_fits.csv")
cat("Linear has ", length(problem_rows), "curve(s) with NAs: ", problem_rows)


##The most common problem you will encounter is unsuitable limits on the parameter estimates. You can tell this is happening if the fitting function returns one of the limits as the best-fit parameter estimate. For example, for the logistic function in 'R function MB' (line 14 - 20), if the lower limit of A = 0 (par_loL(A = 0)), and modelfitting_results shows the fitted value of A = 0 (modelfitting_results$A.A = 0), the actual best-fit value of A is likely below this lower limit, and the fitting function is not returning the best model fit. Unfortunately, it is impossible to set apriori parameter limits that will fit every dataset. You will need to look through your results and make sure that the parameter estimates are within the limits, and that these limits produce reasonable results. If not, adjust the limits in 'R function MB' and call the source() command again.


#Logistic
#Declare the pdf file to hold the plots
pdf(file = 'result/AVSA/Logistic_fits.pdf', width =8.5, height=5, onefile=T)
#Empty the data frame of the previous results
modelfitting_results[FALSE,] -> modelfitting_results
for (ii in 1:length(species_list)){
  subset(dd, dd$to_split == species_list[ii], select =c(1:5)) -> data_by_sp
  define_parsL(data_by_sp) -> par_estimates
  parsL=par_estimates[[1]]
  par_loL=par_estimates[[2]]
  par_highL=par_estimates[[3]]
  Logistic_fits = do_the_thing_nonlinear(data_by_sp, Logistic, parsL, par_loL, par_highL)
  rbind(modelfitting_results, as.data.frame(Logistic_fits))-> modelfitting_results
  print(ii)
}
dev.off()
modelfitting_results$Kmax<- modelfitting_results$A.A
modelfitting_results$psi_kleaf20<- modelfitting_results$C.Xo*(1/0.8-1)^(1/modelfitting_results$B.B)
modelfitting_results$psi_kleaf50<- modelfitting_results$C.Xo*(1/0.5-1)^(1/modelfitting_results$B.B)
modelfitting_results$psi_kleaf80<- modelfitting_results$C.Xo*(1/0.2-1)^(1/modelfitting_results$B.B)
modelfitting_results$psi_kleaf95<- modelfitting_results$C.Xo*(1/0.05-1)^(1/modelfitting_results$B.B)
which(is.na(as.numeric(modelfitting_results[,12]))=='TRUE' | is.na(as.numeric(modelfitting_results[,13]))=='TRUE'| is.na(as.numeric(modelfitting_results[,14]))=='TRUE'| is.na(as.numeric(modelfitting_results[,15]))=='TRUE') -> problem_rows
write.csv(modelfitting_results, file = "result/AVSA/Logistic_fits.csv")
cat("Logistic has ", length(problem_rows), "curve(s) with NAs: ", problem_rows)


#Sigmoidal
#Declare the pdf file to hold the plots
pdf(file = 'result/AVSA/Sigmoidal_fits.pdf', width =8.5, height=5, onefile=T)
#Empty the data frame of results
modelfitting_results[FALSE,] -> modelfitting_results
for (ii in 1:length(species_list)){
  subset(dd, dd$to_split == species_list[ii], select =c(1:5)) -> data_by_sp
  define_parsS(data_by_sp) -> par_estimates
  parsS=par_estimates[[1]]
  par_loS=par_estimates[[2]]
  par_highS=par_estimates[[3]]
  Sigmoidal_fits = do_the_thing_nonlinear(data_by_sp, Sigmoidal, parsS, par_loS, par_highS)
  rbind(modelfitting_results, as.data.frame(Sigmoidal_fits))-> modelfitting_results
  print(ii)
}
dev.off()
modelfitting_results$Kmax<- modelfitting_results$A.A/(1+exp(modelfitting_results$C.Xo/modelfitting_results$B.B))
modelfitting_results$psi_kleaf20<- -modelfitting_results$B.B*(log(1+exp((modelfitting_results$C.Xo/modelfitting_results$B.B)-0.8)/0.8))+modelfitting_results$C.Xo
modelfitting_results$psi_kleaf50<- -modelfitting_results$B.B*(log(1+exp((modelfitting_results$C.Xo/modelfitting_results$B.B)-0.5)/0.5))+modelfitting_results$C.Xo
modelfitting_results$psi_kleaf80<- -modelfitting_results$B.B*(log(1+exp((modelfitting_results$C.Xo/modelfitting_results$B.B)-0.2)/0.2))+modelfitting_results$C.Xo
modelfitting_results$psi_kleaf95<- -modelfitting_results$B.B*(log(1+exp((modelfitting_results$C.Xo/modelfitting_results$B.B)-0.05)/0.05))+modelfitting_results$C.Xo
which(is.na(as.numeric(modelfitting_results[,12]))=='TRUE' | is.na(as.numeric(modelfitting_results[,13]))=='TRUE'| is.na(as.numeric(modelfitting_results[,14]))=='TRUE'| is.na(as.numeric(modelfitting_results[,15]))=='TRUE') -> problem_rows
write.csv(modelfitting_results, file = "result/AVSA/Sigmoidal_fits.csv")
cat("Sigmoidal has ", length(problem_rows), "curve(s) with NAs: ", problem_rows)

#Exponential
pdf(file = 'result/AVSA/Exponential_fits.pdf', width =8.5, height=5, onefile=T)
modelfitting_results[FALSE,] -> modelfitting_results
for (ii in 1:length(species_list)){
  subset(dd, dd$to_split == species_list[ii], select =c(1:5)) -> data_by_sp
  define_parsE(data_by_sp) -> par_estimates
  parsE=par_estimates[[1]]
  par_loE=par_estimates[[2]]
  par_highE=par_estimates[[3]]
  Exponential_fits = do_the_thing_nonlinear(data_by_sp, Exponential, parsE, par_loE, par_highE)
  Exponential_fits$D.NA<- NA
  Exponential_fits$sterror.NA <- NA
  rbind(modelfitting_results, as.data.frame(Exponential_fits))-> modelfitting_results
  print(ii)
}
dev.off()
modelfitting_results$Kmax<-modelfitting_results$A.A
modelfitting_results$psi_kleaf20<- log(0.8)/(-modelfitting_results$B.B)
modelfitting_results$psi_kleaf50<- log(0.5)/(-modelfitting_results$B.B)
modelfitting_results$psi_kleaf80<- log(0.2)/(-modelfitting_results$B.B)
modelfitting_results$psi_kleaf95<- log(0.05)/(-modelfitting_results$B.B)
which(is.na(as.numeric(modelfitting_results[,12]))=='TRUE' | is.na(as.numeric(modelfitting_results[,13]))=='TRUE'| is.na(as.numeric(modelfitting_results[,14]))=='TRUE') -> problem_rows
write.csv(modelfitting_results, file = "result/AVSA/Exponential_fits.csv")
cat("Exponential has ", length(problem_rows), "curve(s) with NAs: ", problem_rows)

#Exponential2
pdf(file = 'result/AVSA/Exponential2_fits.pdf', width =8.5, height=5, onefile=T)
modelfitting_results[FALSE,] -> modelfitting_results
for (ii in 1:length(species_list)){
  subset(dd, dd$to_split == species_list[ii], select =c(1:5)) -> data_by_sp
  define_parsE2(data_by_sp) -> par_estimates
  parsE2=par_estimates[[1]]
  par_loE2=par_estimates[[2]]
  par_highE2=par_estimates[[3]]
  Exponential2_fits = do_the_thing_nonlinear(data_by_sp, Exponential2, parsE2, par_loE2, par_highE2)
  rbind(modelfitting_results, as.data.frame(Exponential2_fits))-> modelfitting_results
  print(ii)
}
dev.off()
modelfitting_results$Kmax<- modelfitting_results$C.C+modelfitting_results$A.A
modelfitting_results$psi_kleaf20<- log(((modelfitting_results$C.C*0.8)+(0.8*modelfitting_results$A.A)-modelfitting_results$C.C)/(modelfitting_results$A.A))/(-modelfitting_results$B.B)
modelfitting_results$psi_kleaf50<- log(((modelfitting_results$C.C*0.5)+(0.5*modelfitting_results$A.A)-modelfitting_results$C.C)/(modelfitting_results$A.A))/(-modelfitting_results$B.B)
modelfitting_results$psi_kleaf80<- log(((modelfitting_results$C.C*0.2)+(0.2*modelfitting_results$A.A)-modelfitting_results$C.C)/(modelfitting_results$A.A))/(-modelfitting_results$B.B)
modelfitting_results$psi_kleaf95<- log(((modelfitting_results$C.C*0.05)+(0.05*modelfitting_results$A.A)-modelfitting_results$C.C)/(modelfitting_results$A.A))/(-modelfitting_results$B.B)
which(is.na(as.numeric(modelfitting_results[,12]))=='TRUE' | is.na(as.numeric(modelfitting_results[,13]))=='TRUE'| is.na(as.numeric(modelfitting_results[,14]))=='TRUE') -> problem_rows
write.csv(modelfitting_results, file = "result/AVSA/Exponential2_fits.csv")
cat("Exponential2 has ", length(problem_rows), "curve(s) with NAs: ", problem_rows)


