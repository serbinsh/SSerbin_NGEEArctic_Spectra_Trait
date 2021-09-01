####################################################################################################
#          
#   Invert leaf reflectance data to estimate transmittance and calculate absorption
#   
#   Approach: BayesianTools DEzs MCMC algorithm  (http://dream.r-forge.r-project.org/)
#   Using fewer-depends rrtm implementation of PROSPECT (https://github.com/ashiklom/rrtm)
#
#   Author: Shawn P. Serbin
#
#
#
#   Project: NGEE-Arctic
#
#    --- Last updated:  08.31.2021 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# load libraries
ok <- require(rrtm) ; if (! ok) {
  devtools::install_github("ashiklom/rrtm")
} else {
  print("*** Package found: rrtm ***")
}

# also requires distributions3 but loading creates namespace issues
# also uses coda
# also currently uses PEcAnRTM - but we dont load it until later.
# how to remove the PEcAn depends for autoburnin and RTM/spec response?  need
# to move spec response function to rrtm
list.of.packages <- c("here", "dplyr", "BayesianTools", "rrtm")
invisible(lapply(list.of.packages, library, character.only = TRUE))

# not in
`%notin%` <- Negate(`%in%`)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get data
inputdir <- file.path(here::here(),"data/compiled_data")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get spectra data: Fit by site and date
load(file.path(inputdir,"NGEEArctic_Leaf_and_Canopy_Reflectance.RData"))
head(NGEEArctic_Reflectance$Leaf_Reflectance)[,1:6]
unique(NGEEArctic_Reflectance$Leaf_Reflectance$Location)
unique(NGEEArctic_Reflectance$Leaf_Reflectance$Location)
unique(NGEEArctic_Reflectance$Leaf_Reflectance$Sample_Date)

# subset options
site_select <- "Utqiagvik" # Utqiagvik, Seward_Peninsula
samp_date <- "20140715" # 20130725 20150712 20140714 20160709 20140715 20160710 20140716 20150717 20160711
# 20150720 20150722 20150723 20150715 20160712 20160713
#

dataset <- NGEEArctic_Reflectance$Leaf_Reflectance %>%
  filter(Location == site_select) %>%
  filter(Sample_Date == samp_date)
head(dataset)[,1:15]

Start.wave <- 350
End.wave <- 2500
spec_waves <- names(dataset)[match(paste0("Wave_",seq(Start.wave,End.wave,1)),names(dataset))]
refl_samp_info <- dataset[,"Sample_ID"]
spectra <- droplevels(dataset[,spec_waves])
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Output directory
out.dir <- file.path(here::here(),"R_Output","PROSPECT5b",site_select,samp_date,'Range_400_700nm')
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create spec file + info to be inverted
output_sample_info_all <- data.frame(Spectra_Name=refl_samp_info,Location=dataset[,"Location"],
                                     Instrument=dataset[,"Instrument"],
                                     Species_Code=dataset[,"USDA_Species_Code"],
                                     Measurement_Date=dataset[,"Sample_Date"])
head(output_sample_info_all)
refl_spec_info2 <- output_sample_info_all
names_output_sample_info <- names(output_sample_info_all)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Subset spectra
Start.wave <- 415 # 470 480; to avoid the "ski jump" < 450 nm with some spectra
End.wave <- 2500
subset_waves <- seq(Start.wave,End.wave,1)

abs.Start.wave <- 400  # start abs calc wavelength
abs.End.wave <- 700    # end abs calc wavelength

sub_refl_data <- spectra[,which(names(spectra) %in% paste0("Wave_",seq(Start.wave,End.wave,1)))]
sub_refl_data <- sub_refl_data*0.01 # scale to 0-1 before proceeding

waves <- seq(Start.wave,End.wave,1)
prospect_waves <- seq(400,2500,1)
#plot(seq(400,500,1),sub_refl_data[124,1:101], type= "l", ylim=c(0,0.1))

## Quick diagnostics plots
plot(seq(Start.wave,End.wave,1),sub_refl_data[1,], type= "l", ylim=c(0,1))
lines(seq(Start.wave,End.wave,1),sub_refl_data[2,])
lines(seq(Start.wave,End.wave,1),sub_refl_data[4,])
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Run inversions - whole directory
output.LRT <- list(Spec.Info=array(NA,dim=c(dim(output_sample_info_all)[1],dim(output_sample_info_all)[2])),
                   obs.Reflectance=array(NA,dim=c(dim(sub_refl_data)[1],length(waves))),
                   mod.Reflectance=array(NA,dim=c(dim(sub_refl_data)[1],length(prospect_waves))),
                   mod.Transmittance=array(NA,dim=c(dim(sub_refl_data)[1],length(prospect_waves))))

mod.params <- array(NA,dim=c(dim(sub_refl_data)[1],18)) # PD param output
inv.samples <- NA
# names: N.mu, N.q25, N.q975, Cab.mu, Cab.q25, Cab.q75, Car.mu, Car.q25, Car.q75, Canth.mu,
# Canth.q25, Canth.q75, Cbrown.mu, Cbrown.q25, Cbrown.q75,
# Cw.mu, Cw.q25, Cw.q75, Cm.mu, Cm.q25, Cm.q75, gelman.diag
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Run inversion

# setup prospect error envelope list
p.refl.stats <- list(lower = array(data=NA,c(dim(sub_refl_data)[1],2101)),
                     upper = array(data=NA,c(dim(sub_refl_data)[1],2101)))

# setup likelihood function
likelihood <- function(params) {
  #N, Cab, Car, Cbrown, Cw, Cm
  mod <- prospect5(params[1], params[2], params[3], params[4], params[5], 
                   params[6])
  sum(dnorm(obs, mod[["reflectance"]], params[7], log = TRUE))
}
# !! in here the likelihood assumes the spectra to invert is labeled "obs"


## --------------------- Define priors --------------------- ##
# Define prior distribution objects here
# N density
Nd <- distributions3::Normal(1, 5)
curve(dnorm(x, 1, 5), 1, 5)

# Cab density
#Cabd <- distributions3::Normal(40, 15)
Cabd <- distributions3::LogNormal(4.1, 0.35)
#Cabd <- distributions3::Normal(65, 20)
#curve(dnorm(x, 40, 15), 0, 120)
curve(dlnorm(x, 4.1, 0.35), 0, 140)
#curve(dnorm(x, 65, 20), 0, 140)

# Car density
Card <- distributions3::LogNormal(2.1, 0.7)
curve(dlnorm(x, 2.1, 0.7), 0, 40)
#curve(dgamma(x, 2.1, 0.2), 0, 40)

# Cbrown density
Cbrownd <- distributions3::Normal(0.05, 0.03)
curve(dnorm(x, 0.03, 0.01), 0, 1)

# Cw density  
#Cwd <- distributions3::LogNormal(-4.456, 1.216)
#curve(dlnorm(x, -4.456, 1.216), 0, 0.2)
Cwd <- distributions3::LogNormal(-4.456, 1.3)
curve(dlnorm(x, -4.456, 1.3), 0, 0.2)

# Cm density
Cmd <- distributions3::LogNormal(-5.15, 1.328)
curve(dlnorm(x, -5.15, 1.328), 0, 0.15)
#curve(dlnorm(x, -3.15, 1.328), 0, 0.15)

# resid density
rsdd <- distributions3::Exponential(10)
curve(dexp(x, 10))
dev.off()
## --------------------- END Define priors --------------------- ##


## --------------------- Create priors --------------------- ##
prior <- createPrior(
  density = function(x) {
    # Note: Truncated at zero
    distributions3::log_pdf(Nd, x[1]) +
      distributions3::log_pdf(Cabd, x[2]) +
      distributions3::log_pdf(Card, x[3]) +
      distributions3::log_pdf(Cbrownd, x[4]) +
      distributions3::log_pdf(Cwd, x[5]) +
      distributions3::log_pdf(Cmd, x[6]) +
      distributions3::log_pdf(rsdd, x[7])
  },
  sampler = function(n = 1) {
    N <- distributions3::random(Nd, n)
    Cab <- distributions3::random(Cabd, n)
    Car <- distributions3::random(Card, n)
    Cbrown <- distributions3::random(Cbrownd, n)
    Cw <- distributions3::random(Cwd, n)
    Cm <- distributions3::random(Cmd, n)
    rsd <- distributions3::random(rsdd, n)
    cbind(N, Cab, Car, Cbrown, Cw, Cm, rsd)
  },
  lower <- c(1, 5, 0, 0, 0, 0, 0)
)
## --------------------- END Create priors --------------------- ##

# setup the inversion
setup <- createBayesianSetup(
  likelihood = likelihood,
  prior = prior
)

# Plot title variable
title_var <- "Spectra_Name"

# Set progress bar
# ?see txtProgressBar
print("Starting Inversion:")
print(paste0("Inverting: ", dim(sub_refl_data)[1]))
print(" ")
pb <- txtProgressBar(min = 0, max = dim(sub_refl_data)[1], width= 50,style = 3)
system.time(for (i in seq_along(1:dim(sub_refl_data)[1]) ) {
#system.time(for (i in seq_along(1:6) ) {
  print(" ")
  print(paste0("Inverting: ",unlist(refl_spec_info2[i,title_var])))
  obs <- unlist(sub_refl_data[i,])
  settings <- list(iterations = 45000)
  samples <- BayesianTools::runMCMC(setup, sampler = "DEzs", settings = settings)
  #samples <- BayesianTools::runMCMC(setup)

  ## !! HERE WE COULD ADD A NON PECAN AUTOBURNIN !!
  #temp <- coda::gelman.plot(BayesianTools::getSample(samples, coda = TRUE), bin.width = 10, max.bins = 50,
  #            confidence = 0.95, transform = FALSE, autoburnin=TRUE, auto.layout = TRUE)
  
  # !! replace this with a non PEcAn dependent approach !!
  samples_burned <- PEcAn.assim.batch::autoburnin(BayesianTools::getSample(samples, coda = TRUE), method = 'gelman.plot')
  coda::varnames(samples_burned) <- c("N", "Cab", "Car", "Cbrown", "Cw", "Cm", "rsd")
  mean_estimates <- do.call(cbind, summary(samples_burned)[c('statistics', 'quantiles')])
  row.names(mean_estimates) <- c("N", "Cab", "Car", "Cbrown", "Cw", "Cm", "rsd")
  
  grDevices::pdf(file = file.path(out.dir,paste0(unlist(refl_spec_info2[i,title_var]),'_MCMC_trace_diag.pdf')), 
      width = 8, height = 6, onefile=T)
  par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  plot(samples_burned)
  dev.off()
  
  # include process error in error envelop (i.e. generate prediction interval)
  #n_target <- 1000
  spec.length <- 2101
  param.samples <- do.call(rbind, samples_burned)
  if(nrow(param.samples)<1000) n_target <- nrow(param.samples) else n_target <- 1000 
  param.samples <- param.samples[sample(nrow(param.samples), n_target), ]
  RT_pred <- array(data=NA,c(n_target,2101))
  print("*** Calculating error stats ***")
  for (r in seq_len(n_target)) {
    perturbed.prospect.ref <- rnorm(spec.length, rrtm::prospect5(param.samples[r,1], param.samples[r,2],
                                                                param.samples[r,3], param.samples[r,4],
                                                                param.samples[r,5], param.samples[r,6])$reflectance,
                                    param.samples[r,7])
    RT_pred[r,] <- perturbed.prospect.ref
  }
  
  # stats
  p.refl.stats$lower[i,] <- apply(RT_pred,2,quantile,probs=c(0.05), na.rm=T)
  p.refl.stats$upper[i,] <- apply(RT_pred,2,quantile,probs=c(0.95), na.rm=T)

  # Generate modeled spectra
  num_params <- 6 # PROSPECT-5b
  input.params <- as.vector(unlist(mean_estimates[,1]))[1:num_params]
  LRT <- prospect5(input.params[1], input.params[2], input.params[3], input.params[4],
                   input.params[5], input.params[6])
  output_sample_info <- droplevels(output_sample_info_all[i,])
  
  output.LRT$Spec.Info[i,] <- unlist(lapply(output_sample_info, as.character))
  output.LRT$obs.Reflectance[i,] <- as.vector(unlist(sub_refl_data[i,]))
  output.LRT$mod.Reflectance[i,] <- LRT$reflectance
  output.LRT$mod.Transmittance[i,] <- LRT$transmittance
  
  # Extract results
  mod.params[i,] <- c(mean_estimates[row.names(mean_estimates)=="N",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="N",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="N",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Cab",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Cab",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Cab",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Car",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Car",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Car",colnames(mean_estimates)=="97.5%"],
                      #mean_estimates[row.names(mean_estimates)=="Canth",colnames(mean_estimates)=="Mean"],
                      #mean_estimates[row.names(mean_estimates)=="Canth",colnames(mean_estimates)=="25%"],
                      #mean_estimates[row.names(mean_estimates)=="Canth",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Cbrown",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Cbrown",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Cbrown",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Cw",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Cw",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Cw",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Cm",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Cm",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Cm",colnames(mean_estimates)=="97.5%"])
  
  setTxtProgressBar(pb, i)
  rm(samples,samples_burned,input.params,LRT,mean_estimates,
     param.samples,RT_pred,perturbed.prospect.ref, output_sample_info)
  print(" ")
  print(" ")
  print("*** Starting new inversion ***")
  print(" ")
  print(" ")
  
  flush.console()
  
}) ## End of inversion loop
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Clean up
# names: N.mu, N.q25, N.q975, Cab.mu, Cab.q25, Cab.q75, Car.mu, Car.q25, Car.q75,
# Cw.mu, Cw.q25, Cw.q75, Cm.mu, Cm.q25, Cm.q75
mod.params <- as.data.frame(mod.params)
#names(mod.params) <- c("N.mu", "N.q25", "N.q975", "Cab.mu", "Cab.q25", "Cab.q975", "Car.mu", "Car.q25",
#                       "Car.q975", "Canth.mu", "Canth.q25", "Canth.q975", "Cbrown.mu", "Cbrown.q25", 
#                       "Cbrown.q975","Cw.mu", "Cw.q25", "Cw.q975", "Cm.mu", "Cm.q25", "Cm.q975")
names(mod.params) <- c("N.mu", "N.q25", "N.q975", "Cab.mu", "Cab.q25", "Cab.q975", "Car.mu", "Car.q25",
                       "Car.q975", "Cbrown.mu", "Cbrown.q25", "Cbrown.q975","Cw.mu", "Cw.q25", "Cw.q975", 
                       "Cm.mu", "Cm.q25", "Cm.q975")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Plot comparison
grDevices::pdf(file=file.path(out.dir,'PROSPECTD_Inversion_Diagnostics.pdf'),height=8,width=9)
par(mfrow=c(1,1), mar=c(4.3,4.3,1.0,4.3), oma=c(0.1,0.1,0.1,0.1)) # B L T R
for (i in seq_along(1:dim(sub_refl_data)[1] )) {
#for (i in seq_along(1:6) ) {
  plot(waves,output.LRT$obs.Reflectance[i,], type="l", col="black",xlab="Wavelength (nm)",ylab="Reflectance (0-1)",
       lwd=3,main=paste0(output.LRT$Spec.Info[i,1]," ", output.LRT$Spec.Info[i,3]),
       ylim=c(min(p.refl.stats$lower[i,]), max(p.refl.stats$upper[i,])+0.1) )
  lines(prospect_waves, rrtm::prospect5(mod.params[i,"N.mu"],mod.params[i,"Cab.mu"],mod.params[i,"Car.mu"],
                                        mod.params[i,"Cbrown.mu"],mod.params[i,"Cw.mu"],
                                        mod.params[i,"Cm.mu"])$reflectance,col="red",lwd=2)
  polygon(c(prospect_waves ,rev(prospect_waves)),c(p.refl.stats$upper[i,], rev(p.refl.stats$lower[i,])),
          col="grey70",border=NA)
  lines(waves,output.LRT$obs.Reflectance[i,],col="black",lwd=2)
  lines(prospect_waves, rrtm::prospect5(mod.params[i,"N.mu"],mod.params[i,"Cab.mu"],mod.params[i,"Car.mu"],
                                        mod.params[i,"Cbrown.mu"],mod.params[i,"Cw.mu"],
                                        mod.params[i,"Cm.mu"])$reflectance,col="red",lwd=2)
  legend("topright",legend=c("Observed","Modeled","Modeled 95% PI"),lty=1,col=c("black","red","grey70"),
         lwd=c(3,3,8),bty = "n")
  box(lwd=2.2)
  plot(waves,output.LRT$obs.Reflectance[i,], type="l", col="black",xlab="Wavelength (nm)",ylab="Reflectance (0-1)",
       lwd=3,ylim=c(0,1),main=paste0(output.LRT$Spec.Info[i,1]," ", output.LRT$Spec.Info[i,3]),cex.lab=1.7)
  lines(prospect_waves, rrtm::prospect5(mod.params[i,"N.mu"],mod.params[i,"Cab.mu"],mod.params[i,"Car.mu"],
                                        mod.params[i,"Cbrown.mu"],mod.params[i,"Cw.mu"],
                                        mod.params[i,"Cm.mu"])$reflectance,col="red",lwd=2)
  lines(prospect_waves,1-rrtm::prospect5(mod.params[i,"N.mu"],mod.params[i,"Cab.mu"],mod.params[i,"Car.mu"],
                                         mod.params[i,"Cbrown.mu"],mod.params[i,"Cw.mu"],
                                         mod.params[i,"Cm.mu"])$transmittance,col="grey70",lwd=3) 
  axis(4,labels = rev(seq(0.0, 1.0, 0.2)), at = seq(0.0, 1.0, 0.2))
  mtext("Transmittance (0-1)",side=4,line=3,cex=1.7)
  legend("right",legend=c("Meas. Reflectance","Mod. Reflectance","Mod. Transmittance"),
         lty=1,col=c("black","red","grey70"),
         lwd=c(3,3,8),bty = "n")
  box(lwd=2.2)
}
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Calculate absorption - use measured refl, modeled trans
Start.wave <- 400; End.wave <- 2500
waves <- seq(Start.wave,End.wave,1)
#full_spectrum_absorption <- data.frame(1 - output.LRT$mod.Transmittance - output.LRT$obs.Reflectance)
full_spectrum_absorption <- data.frame(1 - output.LRT$mod.Transmittance - 
                                         output.LRT$mod.Reflectance) # based on modeled results only
names(full_spectrum_absorption) <- paste0("Wave_",seq(Start.wave,End.wave,1))
plot(waves, full_spectrum_absorption[1,],type="l")
lines(waves, full_spectrum_absorption[2,])
lines(waves, full_spectrum_absorption[3,])

sub_spec_abs <- full_spectrum_absorption[,which(names(full_spectrum_absorption) %in% 
                                                  paste0("Wave_",seq(abs.Start.wave,abs.End.wave,1)))]
leaf_vis_absorption <- rowMeans(sub_spec_abs,na.rm = T)
leaf_vis_absorption_final <- data.frame(output.LRT$Spec.Info, leaf_vis_absorption)
names(leaf_vis_absorption_final) <- c(names_output_sample_info,"Leaf_VIS_Spectral_Absorption")

## Licor Fs absorption
library(PEcAnRTM)
data(sensor.rsr)
licor.abs <- array(NA,dim=c(dim(mod.params)[1],3))
for (i in 1:dim(mod.params)[1] ) {
  #for (i in 1:6 ) {
  # Using observed refl data instead of modeled
  #refl <- spectral.response(unlist(sub_refl_data[i,]), 'licor')
  # Using modeled reflectance spectra
  
  # !!!! HERE ALSO NEED TO REMOVE DEPENDS ON PECAN !!!!
  refl <- PEcAnRTM::spectral.response(as.vector(rrtm::prospect5(mod.params[i,"N.mu"],mod.params[i,"Cab.mu"],
                                                                mod.params[i,"Car.mu"],
                                                                mod.params[i,"Cbrown.mu"],mod.params[i,"Cw.mu"],
                                                                mod.params[i,"Cm.mu"])$reflectance), 'licor')
  trans <- PEcAnRTM::spectral.response(as.vector(rrtm::prospect5(mod.params[i,"N.mu"],mod.params[i,"Cab.mu"],
                                                                mod.params[i,"Car.mu"],
                                                                mod.params[i,"Cbrown.mu"],mod.params[i,"Cw.mu"],
                                                                mod.params[i,"Cm.mu"])$transmittance), 'licor')
  licor.abs[i,] <- 1-trans-refl
  rm(refl,trans)
}
licor.abs <- data.frame(licor.abs)
names(licor.abs) <- c("Blue","Red","Far_Red")
licor.abs <- data.frame(licor.abs$Red,licor.abs$Blue,licor.abs$Far_Red)
names(licor.abs) <- c("Red","Blue","Far_Red")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Output data

# Modeled Reflectance
output.refl <- data.frame(output.LRT$Spec.Info, output.LRT$mod.Reflectance)
names(output.refl) <- c(names_output_sample_info, paste0("Wave_",seq(Start.wave,End.wave,1)))
write.csv(output.refl,file=file.path(out.dir,'PROSPECTD_Modeled_Reflectance.csv'),
          row.names=FALSE)

prospect_refl <- output.LRT$mod.Reflectance
names(prospect_refl) <- paste0("Wave_",prospect_waves)
keep <- which(names(prospect_refl) %in% paste0("Wave_",subset_waves))
prospect_refl <- prospect_refl[,keep]
#plot(subset_waves,prospect_refl[1,], type="l")
resids <- as.matrix(prospect_refl)-as.matrix(output.LRT$obs.Reflectance)
output.resids <- data.frame(output.LRT$Spec.Info,resids)
#plot(subset_waves,resids[1,], type="l")
names(output.resids) <- c(names_output_sample_info, paste0("Wave_",subset_waves))
write.csv(output.resids,file=file.path(out.dir,'PROSPECTD_Reflectance_Residuals.csv'),
          row.names=FALSE)

# Modeled Transmittance
output.trans <- data.frame(output.LRT$Spec.Info, output.LRT$mod.Transmittance)
names(output.trans) <- c(names_output_sample_info, paste0("Wave_",seq(Start.wave,End.wave,1)))
write.csv(output.trans,file=file.path(out.dir,'PROSPECTD_Modeled_Transmittance.csv'),
          row.names=FALSE)

## Modeled Absorption
# full spectra
output.abs <- data.frame(output.LRT$Spec.Info, full_spectrum_absorption)
names(output.abs) <- c(names_output_sample_info, paste0("Wave_",seq(Start.wave,End.wave,1)))
write.csv(output.abs,file=file.path(out.dir,'PROSPECTD_Modeled_Full_Spec_Absorption.csv'),
          row.names=FALSE)

# vis
write.csv(leaf_vis_absorption_final,file=file.path(out.dir,'PROSPECTD_Estimated_Leaf_VIS_Absorption.csv'),
          row.names=FALSE)
# Licor Fs
output.licor.abs <- data.frame(output.LRT$Spec.Info,licor.abs)
names(output.licor.abs) <- c(names_output_sample_info,"Red","Blue","Far_Red")
write.csv(output.licor.abs,file=file.path(out.dir,'PROSPECTD_Modeled_LiCor_Absorption.csv'),
          row.names=FALSE)

# Modeled params
output.params <- data.frame(output.LRT$Spec.Info,mod.params)
names(output.params) <- c(names_output_sample_info, names(mod.params))
write.csv(output.params,file=file.path(out.dir,'PROSPECTD_Modeled_Parameters.csv'),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#

### EOF