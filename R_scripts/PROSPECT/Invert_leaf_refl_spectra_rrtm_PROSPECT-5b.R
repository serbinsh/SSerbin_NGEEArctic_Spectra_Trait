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

# setup likelihood function
likelihood <- function(params) {
  #N, Cab, Car, Cbrown, Cw, Cm
  mod <- prospect5(params[1], params[2], params[3], params[4], params[5], 
                   params[6])
  sum(dnorm(obs, mod[["reflectance"]], log = TRUE))
  #mod <- rrtm::prospect5(params[1], params[2], params[3], params[4], params[5],
  #                       params[6])$reflectance[min(which(prospect_waves %in% waves, 
  #                                                 arr.ind = TRUE)):2101]
  #sum(dnorm(obs, mod, log = TRUE))
}
# !! in here the likelihood assumes the spectra to invert is labeled "obs"

# setup prospect error envelope list
p.refl.stats <- list(lower = array(data=NA,c(dim(sub_refl_data)[1],2101)),
                     upper = array(data=NA,c(dim(sub_refl_data)[1],2101)))

# create prior
# N 
curve(dnorm(x, 1, 5), 1, 5)
# Cab density
#curve(dnorm(x, 40, 15), 0, 120)
curve(dlnorm(x, 4.0, 0.28), 0, 120)
# Car density
curve(dlnorm(x, 2.1, 0.7), 0, 40)
#curve(dgamma(x, 2.1, 0.2), 0, 40)
# Cbrown density
curve(dnorm(x, 0.07, 0.04), 0, 1)
#curve(dlnorm(x, 0.01, 3), 0, 1)
#curve(dgamma(x, 2, 12), 0, 1)
# Cw density
curve(dlnorm(x, -4.456, 1.216), 0, 0.2)
# Cm density
#curve(dlnorm(x, -5.15, 1.328), 0, 0.15)
curve(dlnorm(x, -4.15, 1.328), 0, 0.15)
# resid
curve(dexp(x, 10))

prior <- createPrior(
  density = function(x) {
    # Note: Truncated at zero
    N <- dnorm(x[1], 1, 5, log = TRUE)
    Cab <- dlnorm(x[2], 4.0, 0.28, log = TRUE)
    Car <- dlnorm(x[3], 2.1, 0.7, log = TRUE)
    Cbrown <- dnorm(x[4], 0.07, 0.04, log = TRUE)
    Cw <- dlnorm(x[5], -4.456, 1.216, log = TRUE)
    Cm <- dlnorm(x[6], -4.15, 1.328, log = TRUE)
    rsd <- dexp(x[7], 10, log = TRUE)
    N + Cab + Car + Cbrown + Cw + Cm + rsd
  },
  sampler = function(n = 1) {
    N <- rnorm(n, 1, 5)
    Cab <- rlnorm(n, 4.0, 0.28)
    Car <- rlnorm(n, 2.1, 0.7)
    Cbrown <- rnorm(n, 0.07, 0.04)
    Cw <- rlnorm(n, -4.456, 1.216)
    Cm <- rlnorm(n, -4.15, 1.328)
    rsd <- rexp(n, 10)
    cbind(N, Cab, Car, Cbrown, Cw, Cm, rsd)
  },
  lower = c(1, 0, 0, 0, 0, 0, 0)
)

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
#system.time(for (i in seq_along(1:dim(sub_refl_data)[1]) ) {
system.time(for (i in seq_along(1:3) ) {
  print(" ")
  print(paste0("Inverting: ",unlist(refl_spec_info2[i,title_var])))
  obs <- unlist(sub_refl_data[i,])
  settings <- list(iterations = 45000)
  samples <- BayesianTools::runMCMC(setup, settings = settings)
  #samples <- BayesianTools::runMCMC(setup)
  
  samples_burned <- PEcAn.assim.batch::autoburnin(BayesianTools::getSample(samples, coda = TRUE), method = 'gelman.plot')
  #samples_burned <- BayesianTools::getSample(samples, coda = TRUE)
  mean_estimates <- do.call(cbind, summary(samples_burned)[c('statistics', 'quantiles')])
  row.names(mean_estimates) <- c("N", "Cab", "Car", "Cbrown", "Cw", "Cm", "resid")
  
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
  
  pdf(file = file.path(out.dir,paste0(unlist(refl_spec_info2[i,title_var]),'_MCMC_trace_diag.pdf')), 
      width = 8, height = 6, onefile=T)
  par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  plot(samples_burned)
  dev.off()
  
  
  # Generate modeled spectra
  num_params <- 6 # PROSPECT-5b
  input.params <- as.vector(unlist(mean_estimates[,1]))[1:num_params]
  LRT <- prospect5(input.params[1], input.params[2], input.params[3], input.params[4],
                   input.params[5], input.params[6])
  output_sample_info <- droplevels(output_sample_info_all[i,])
  
  # for testing
  # output_sample_info
  # str(output_sample_info)
  # unlist(lapply(output_sample_info, as.character))
  # output.LRT2 <- output.LRT
  # output.LRT2$Spec.Info[i,] <- unlist(lapply(output_sample_info, as.character))
  
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
pdf(file=file.path(out.dir,'PROSPECTD_Inversion_Diagnostics.pdf'),height=8,width=9)
par(mfrow=c(1,1), mar=c(4.3,4.3,1.0,4.3), oma=c(0.1,0.1,0.1,0.1)) # B L T R
#for (i in seq_along(1:dim(sub_refl_data)[1] )) {
for (i in seq_along(1:3) ) {
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




#plot(samps)
