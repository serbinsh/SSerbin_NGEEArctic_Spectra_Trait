####################################################################################################
#          
#   Invert leaf reflectance data to estimate transmittance and calculate absorption
#   
#   Approach: BayesianTools DEzs MCMC algorithm  (http://dream.r-forge.r-project.org/)
#
#   Author: Shawn P. Serbin
#
#
#
#   Project: NGEE-Arctic
#
#    --- Last updated:  06.21.2021 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# load libraries
ok <- require(PEcAn.assim.batch) ; if (! ok) {
  devtools::install_github("PEcAnProject/pecan", subdir="modules/assim.batch", 
                           ref = "develop") # use development version of PEcAn
} else {
  print("*** Package found: PEcAn.assim.batch ***")
}

ok <- require(PEcAnRTM) ; if (! ok) {
  devtools::install_github("PEcAnProject/pecan", subdir="modules/rtm",
                           ref = "develop") # use development version of PEcAn
} else {
  print("*** Package found: PEcAnRTM ***")
}

ok <- require(rrtm) ; if (! ok) {
  devtools::install_github("ashiklom/rrtm")
} else {
  print("*** Package found: rrtm ***")
}

#library(PEcAnRTM)
#library(rrtm)

list.of.packages <- c("here","PEcAnRTM","dplyr")
invisible(lapply(list.of.packages, library, character.only = TRUE))

# not in
`%notin%` <- Negate(`%in%`)

custom_prior <- TRUE # TRUE/FALSE
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

# subset options
site_select <- "Utqiagvik" # Utqiagvik, Seward_Peninsula
samp_date <- "20130725"

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
out.dir <- file.path(here::here(),"R_Output","PROSPECTD",site_select,samp_date,'Range_400_700nm')
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

mod.params <- array(NA,dim=c(dim(sub_refl_data)[1],21)) # PD param output
inv.samples <- NA
# names: N.mu, N.q25, N.q975, Cab.mu, Cab.q25, Cab.q75, Car.mu, Car.q25, Car.q75, Canth.mu,
# Canth.q25, Canth.q75, Cbrown.mu, Cbrown.q25, Cbrown.q75,
# Cw.mu, Cw.q25, Cw.q75, Cm.mu, Cm.q25, Cm.q75, gelman.diag
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Run inversion
prospect_ver <- "D"

# setup prospect error envelope list
p.refl.stats <- list(lower = array(data=NA,c(dim(sub_refl_data)[1],2101)),
                     upper = array(data=NA,c(dim(sub_refl_data)[1],2101)))

# setup inversion model
model <- function(x) prospect(x, prospect_ver)[min(which(prospect_waves %in% waves, 
                                                         arr.ind = TRUE)):2101,1] # subset

if (custom_prior) {
  prior <- prospect_bt_prior("D")
  prior$lower <- c(1,5,1,1,0,0,0,0) # "N"      "Cab"    "Car"    "Canth"  "Cbrown" "Cw"     "Cm" , resid
  #prior$upper <- c(1,Inf,80,80,Inf,Inf,Inf,Inf) # "N"      "Cab"    "Car"    "Canth"  "Cbrown" "Cw"     "Cm" , resid
  prior$upper <- c(1,Inf,25,5,1,Inf,Inf,Inf) # "N"      "Cab"    "Car"    "Canth"  "Cbrown" "Cw"     "Cm" , resid
  prior$best <- c(2,55,8,2.5,0.003,0.014,0.001,0.001)
  
  #cust_prior <- list(Cab = list('Cab', 'lnorm', log(3.6683201), 0.1952258, 1))
  #prior <- prospect_bt_prior("D", custom_prior=cust_prior)
} else {
  prior <- prospect_bt_prior("D")
}

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
  #https://cran.r-project.org/web/packages/BayesianTools/vignettes/BayesianTools.html#the-different-mcmc-samplers
  #Metropolis, DE, DEzs, DREAM, DREAMzs
  samples <- invert_bt(observed = t(sub_refl_data[i,]), model = model, prior = prior,
                       custom_settings = list(init = list(iterations = 2000),
                                              loop = list(iterations = 1000),
                                              other = list(sampler = "DEzs",
                                                           max_iter = 1e6, min_samp = 5000, threshold = 1.1)))
  samples_burned <- PEcAn.assim.batch::autoburnin(BayesianTools::getSample(samples, coda = TRUE), method = 'gelman.plot')
  mean_estimates <- do.call(cbind, summary(samples_burned)[c('statistics', 'quantiles')])
  
  # include process error in error envelop (i.e. generate prediction interval)
  n_target <- 1000
  spec.length <- 2101
  param.samples <- do.call(rbind, samples_burned)
  param.samples <- param.samples[sample(nrow(param.samples), n_target), ]
  RT_pred <- array(data=NA,c(n_target,2101))
  print("*** Calculating error stats ***")
  for (r in seq_len(n_target)) {
    perturbed.prospect.ref <- rnorm(spec.length,prospect(param.samples[r,1:7], prospect_ver)[,1], param.samples[r,6])
    RT_pred[r,] <- perturbed.prospect.ref
  }
  
  # stats
  p.refl.stats$lower[i,] <- apply(RT_pred,2,quantile,probs=c(0.05), na.rm=T)
  p.refl.stats$upper[i,] <- apply(RT_pred,2,quantile,probs=c(0.95), na.rm=T)
  
  pdf(file = file.path(out.dir,paste0(unlist(refl_spec_info2[i,title_var]),'_MCMC_trace_diag.pdf')), 
      width = 8, height = 6, onefile=T)
  par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  plot(samples_burned)
  dev.off()
  
  # Generate modeled spectra
  if (prospect_ver=="D") {
    num_params <- 7
  }
  input.params <- as.vector(unlist(mean_estimates[,1]))[1:num_params]
  LRT <- prospect(param = input.params, version=prospect_ver)
  output_sample_info <- droplevels(output_sample_info_all[i,])
  
  # for testing
  #output_sample_info
  #str(output_sample_info)
  #unlist(lapply(output_sample_info, as.character))
  #output.LRT2 <- output.LRT
  #output.LRT2$Spec.Info[i,] <- unlist(lapply(output_sample_info, as.character))
  
  output.LRT$Spec.Info[i,] <- unlist(lapply(output_sample_info, as.character))
  output.LRT$obs.Reflectance[i,] <- as.vector(unlist(sub_refl_data[i,]))
  output.LRT$mod.Reflectance[i,] <- LRT[,1]
  output.LRT$mod.Transmittance[i,] <- LRT[,2]
  
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
                      mean_estimates[row.names(mean_estimates)=="Canth",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Canth",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Canth",colnames(mean_estimates)=="97.5%"],
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
names(mod.params) <- c("N.mu", "N.q25", "N.q975", "Cab.mu", "Cab.q25", "Cab.q975", "Car.mu", "Car.q25",
                       "Car.q975", "Canth.mu", "Canth.q25", "Canth.q975", "Cbrown.mu", "Cbrown.q25", 
                       "Cbrown.q975","Cw.mu", "Cw.q25", "Cw.q975", "Cm.mu", "Cm.q25", "Cm.q975")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Plot comparison
pdf(file=file.path(out.dir,'PROSPECTD_Inversion_Diagnostics.pdf'),height=8,width=9)
par(mfrow=c(1,1), mar=c(4.3,4.3,1.0,4.3), oma=c(0.1,0.1,0.1,0.1)) # B L T R
for (i in seq_along(1:dim(sub_refl_data)[1] )) {
#for (i in seq_along(1:6) ) {
  plot(waves,output.LRT$obs.Reflectance[i,], type="l", col="black",xlab="Wavelength (nm)",ylab="Reflectance (0-1)",
       lwd=3,main=paste0(output.LRT$Spec.Info[i,1]," ", output.LRT$Spec.Info[i,3]) )
  lines(prospect_waves,prospect(param = mod.params[i,c(1,4,7,10,13,16,19)], version=prospect_ver)[,1],col="red",lwd=2)
  polygon(c(prospect_waves ,rev(prospect_waves)),c(p.refl.stats$upper[i,], rev(p.refl.stats$lower[i,])),
          col="grey70",border=NA)
  lines(waves,output.LRT$obs.Reflectance[i,],col="black",lwd=2)
  lines(prospect_waves,prospect(param = mod.params[i,c(1,4,7,10,13,16,19)], version=prospect_ver)[,1],col="red",lwd=2)
  legend("topright",legend=c("Observed","Modeled","Modeled 95% PI"),lty=1,col=c("black","red","grey70"),
         lwd=c(3,3,8),bty = "n")
  box(lwd=2.2)
  plot(waves,output.LRT$obs.Reflectance[i,], type="l", col="black",xlab="Wavelength (nm)",ylab="Reflectance (0-1)",
       lwd=3,ylim=c(0,1),main=paste0(output.LRT$Spec.Info[i,1]," ", output.LRT$Spec.Info[i,3]),cex.lab=1.7)
  lines(prospect_waves,prospect(param = mod.params[i,c(1,4,7,10,13,16,19)], version=prospect_ver)[,1],col="red",lwd=3)
  lines(prospect_waves,1-prospect(param = mod.params[i,c(1,4,7,10,13,16,19)], version=prospect_ver)[,2],col="grey70",lwd=3) 
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
data(sensor.rsr)
licor.abs <- array(NA,dim=c(dim(mod.params)[1],3))
for (i in 1:dim(mod.params)[1] ) {
#for (i in 1:6 ) {
  # Using observed refl data instead of modeled
  #refl <- spectral.response(unlist(sub_refl_data[i,]), 'licor')
  # Using modeled reflectance spectra
  refl <- spectral.response(as.vector(prospect(param = mod.params[i,c(1,4,7,10,13,16,19)], 
                                               version=prospect_ver)[,1]), 'licor')
  trans <- spectral.response(as.vector(prospect(param = mod.params[i,c(1,4,7,10,13,16,19)], 
                                                version=prospect_ver)[,2]), 'licor')
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