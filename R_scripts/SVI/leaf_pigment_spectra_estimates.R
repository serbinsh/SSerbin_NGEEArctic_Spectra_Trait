####################################################################################################
#
#
#
#    --- Last updated: 05.07.2021 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Load libraries
list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra",
                      "spectratrait")
invisible(lapply(list.of.packages, library, character.only = TRUE))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup other functions and options
`%notin%` <- Negate(`%in%`)
output_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/R_Output/SVI/leaf/Chl_ab")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Load datasets
data_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
load(file.path(data_dir,"NGEEArctic_Leaf_Pigments.RData"))
head(Utqiagvik_2013_2015_leaf_pigment_data)

load(file.path(data_dir,"NGEEArctic_Leaf_and_Canopy_Reflectance.RData"))
head(NGEEArctic_Reflectance$Leaf_Reflectance)[,1:6]

# What is the target variable?
#inVar <- "Narea_g_m2"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory
if (output_dir=="tempdir") {
  outdir <- tempdir()
} else {
  if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
  outdir <- file.path(path.expand(output_dir))
}
setwd(outdir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create dataset
Utqiagvik_2013_2015_leaf_pigment_data <- Utqiagvik_2013_2015_leaf_pigment_data %>%
  mutate(Chl_ab_area_L_mgcm2=Chl_ab_area_L*0.0001) %>%
  select(Location,Sample_ID,Sample_Date,USDA_Species_Code,Chl_ab_area_L,
         Chl_ab_area_L_mgcm2)
head(Utqiagvik_2013_2015_leaf_pigment_data)

# merge data
lr <- NGEEArctic_Reflectance$Leaf_Reflectance %>%
  select(Sample_ID,Instrument,starts_with("Wave_"))
head(lr)[,1:6]

# remove spec outliers
lr <- lr %>%
  filter(Wave_830>35)

merged_data <- merge(x = Utqiagvik_2013_2015_leaf_pigment_data, y = lr, by = "Sample_ID")
head(merged_data)[,1:15]
unique(merged_data$Sample_Date)

Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
Spectra <- as.matrix(merged_data[,names(merged_data) %in% paste0("Wave_",wv)])
head(Spectra)[1:6,1:10]

sample_info <- merged_data[,names(merged_data) %notin% paste0("Wave_",seq(350,2500,1))]
head(sample_info)

sample_info2 <- sample_info %>%
  select(Location, Sample_ID, Sample_Date, USDA_Species_Code, Spec_Instrument=Instrument, Chl_ab_area_L, 
         Chl_ab_area_L_mgcm2)

analysis_data <- data.frame(sample_info2,Spectra)
rm(sample_info,sample_info2,Spectra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
analysis_data <- analysis_data %>%
  mutate(ChlNDI=(Wave_750-Wave_705)/(Wave_750+Wave_705)) %>%
  mutate(Ratio=Wave_736/Wave_751) %>%
  mutate(pred_tot_chlab=1.81*10^-4+(4.6*10^-2*ChlNDI)+(5.12*10^-2*ChlNDI^2)) %>%
  select(!starts_with("Wave_"))
head(analysis_data)

plot(analysis_data$Chl_ab_area_L_mgcm2*1000, analysis_data$pred_tot_chlab*1000,
     xlim=c(0,75),ylim = c(0,75))
abline(0,1,lty=2)


plot(analysis_data$Ratio,analysis_data$Chl_ab_area_L_mgcm2*1000)
#--------------------------------------------------------------------------------------------------#


