####################################################################################################
#
#    --- Last updated: 12.04.2020 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
list.of.packages <- c("devtools","readr","RCurl","httr","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))

# Load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
devtools::source_url("https://raw.githubusercontent.com/TESTgroup-BNL/PLSR_for_plant_trait_prediction/master/R_Scripts/functions.R")

# not in
`%notin%` <- Negate(`%in%`)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### output location
output_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Barrow datasets

#NGEE Arctic Leaf Spectral Reflectance Utqiagvik (Barrow) Alaska 2013
ecosis_id <- "224c23b4-9e85-4275-b5f4-308db02547b3"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:45]

# clean up
dat_raw[dat_raw==-9999]=NA
head(dat_raw)

dat_raw$Sample_ID <- lapply(dat_raw$Sample_ID, function(x) paste0('BNL', x))
dat_raw$Sample_ID <- as.character(dat_raw$Sample_ID)
head(dat_raw)
names(dat_raw)[1:45]

nga_2013_barrow_leaf_spec <- dat_raw %>%
  select(Sample_ID,Sample_Date=`Measurement Date`,USDA_Species_Code=`USDA Symbol`,Instrument=`Instrument Model`,
         `350`:`2500`)
head(nga_2013_barrow_leaf_spec)
spectra <- nga_2013_barrow_leaf_spec %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))
spec_info <- nga_2013_barrow_leaf_spec[,names(nga_2013_barrow_leaf_spec) %notin% seq(350,2500,1)]

if (max(spectra)<1) {
  spectra <- spectra*100
}
nga_2013_barrow_leaf_spec_out <- data.frame(spec_info,spectra)
head(nga_2013_barrow_leaf_spec_out)[,1:5]
#plot(seq(350,2500,1), nga_2013_barrow_leaf_spec_out[1,5:dim(nga_2013_barrow_leaf_spec_out)[2]])

write.csv(nga_2013_barrow_leaf_spec_out, 
          file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_2013_Leaf_Spectral_Reflectance.csv"), 
          row.names = F)


# NGEE Arctic Leaf Spectral Reflectance and Transmittance Data 2014 to 2016 Utqiagvik (Barrow) Alaska
ecosis_id <- "bf41fff2-8571-4f34-bd7d-a3240a8f7dc8"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:40]

# clean up
dat_raw[dat_raw==-9999]=NA
names(dat_raw)

nga_20142016_barrow_leaf_spec <- dat_raw %>%
  filter(Measurement=="Reflectance") %>%
  select(Sample_ID=BNL_Barcode,Sample_Date=`Measurement Date`,USDA_Species_Code=`USDA Symbol`,
         Instrument=`Instrument Model`,`350`:`2500`)
head(nga_20142016_barrow_leaf_spec)
spectra <- nga_20142016_barrow_leaf_spec %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))
spec_info <- nga_20142016_barrow_leaf_spec[,names(nga_20142016_barrow_leaf_spec) %notin% seq(350,2500,1)]
if (mean(spectra$Wave_850, na.rm=T)<1) {
  spectra <- spectra*100
}
nga_20142016_barrow_leaf_spec_out <- data.frame(spec_info,spectra)
head(nga_20142016_barrow_leaf_spec_out)[,1:5]

write.csv(nga_20142016_barrow_leaf_spec_out, 
          file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_2014_2016_Leaf_Spectral_Reflectance.csv"), 
          row.names = F)

# NGEE Arctic BNL Canopy Spectral Reflectance Utqiagvik (Barrow) Alaska 2014 to 2016
ecosis_id <- "ae1463cc-9984-4a0e-8277-0ba784eda5fd"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:40]

# clean up
dat_raw[dat_raw==-9999]=NA
names(dat_raw)

nga_20142016_barrow_canopy_spec <- dat_raw %>%
  select(Sample_ID,Sample_Date=`Measurement Date`,USDA_Species_Code=`USDA Symbol`,
         Instrument=Instrument,`350`:`2500`)
head(nga_20142016_barrow_canopy_spec)
spectra <- nga_20142016_barrow_canopy_spec %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))
spec_info <- nga_20142016_barrow_canopy_spec[,names(nga_20142016_barrow_canopy_spec) %notin% seq(350,2500,1)]
if (mean(spectra$Wave_850, na.rm=T)<1) {
  spectra <- spectra*100
}
nga_20142016_barrow_canopy_spec_out <- data.frame(spec_info,spectra)
head(nga_20142016_barrow_canopy_spec_out)[,1:5]

write.csv(nga_20142016_barrow_canopy_spec_out, 
          file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_2014_2016_Canopy_Spectral_Reflectance.csv"), 
          row.names = F)

### SewPen datasets
#NGEE Arctic 2016 Leaf Spectral Reflectance Kougarok Road Seward Peninsula Alaska 2016
ecosis_id <- "960dbb0c-144e-4563-8117-9e23d14f4aa9"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)

# clean up
dat_raw[dat_raw==-9999]=NA
names(dat_raw)[1:40]

nga_2016_sewpen_leaf_spec <- dat_raw %>%
  select(Sample_ID=BNL_Barcode,Sample_Date=`Sample Collection Date`,USDA_Species_Code=`USDA Symbol`,
         Instrument=`Instrument Model`,`350`:`2500`)
head(nga_2016_sewpen_leaf_spec)[,1:6]
spectra <- nga_2016_sewpen_leaf_spec %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))
spec_info <- nga_2016_sewpen_leaf_spec[,names(nga_2016_sewpen_leaf_spec) %notin% seq(350,2500,1)]
if (mean(spectra$Wave_850, na.rm=T)<1) {
  spectra <- spectra*100
}
nga_2016_sewpen_leaf_spec_out <- data.frame(spec_info,spectra)
head(nga_2016_sewpen_leaf_spec_out)[,1:5]

write.csv(nga_2016_sewpen_leaf_spec_out, 
          file = file.path(output_dir,"NGEE-Arctic_SewardPeninsula_2016_Leaf_Spectral_Reflectance.csv"), 
          row.names = F)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Combine datasets
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
Utqiagvik_2013_Leaf_Reflectance <- read.csv(file = file.path(input_dir,"NGEE-Arctic_Utqiagvik_2013_Leaf_Spectral_Reflectance.csv"))
orig_names <- names(Utqiagvik_2013_Leaf_Reflectance)
Utqiagvik_2013_Leaf_Reflectance$Location <- rep("Utqiagvik",times=dim(Utqiagvik_2013_Leaf_Reflectance)[1])
Utqiagvik_2013_Leaf_Reflectance <- Utqiagvik_2013_Leaf_Reflectance %>%
  select(Location,paste(orig_names))
head(Utqiagvik_2013_Leaf_Reflectance)[,1:10]

Utqiagvik_2014_2016_Leaf_Reflectance <- read.csv(file = file.path(input_dir,"NGEE-Arctic_Utqiagvik_2014_2016_Leaf_Spectral_Reflectance.csv"))
orig_names <- names(Utqiagvik_2014_2016_Leaf_Reflectance)
Utqiagvik_2014_2016_Leaf_Reflectance$Location <- rep("Utqiagvik",times=dim(Utqiagvik_2014_2016_Leaf_Reflectance)[1])
Utqiagvik_2014_2016_Leaf_Reflectance <- Utqiagvik_2014_2016_Leaf_Reflectance %>%
  select(Location,paste(orig_names))
head(Utqiagvik_2014_2016_Leaf_Reflectance)[,1:10]

SewardPeninsual_2016_Leaf_Reflectance <- read.csv(file = file.path(input_dir,"NGEE-Arctic_SewardPeninsula_2016_Leaf_Spectral_Reflectance.csv"))
orig_names <- names(SewardPeninsual_2016_Leaf_Reflectance)
SewardPeninsual_2016_Leaf_Reflectance$Location <- rep("Seward_Peninsula",times=dim(SewardPeninsual_2016_Leaf_Reflectance)[1])
SewardPeninsual_2016_Leaf_Reflectance <- SewardPeninsual_2016_Leaf_Reflectance %>%
  select(Location,paste(orig_names))
head(SewardPeninsual_2016_Leaf_Reflectance)[,1:10]

Utqiagvik_2014_2016_Canopy_Reflectance <- read.csv(file = file.path(input_dir,"NGEE-Arctic_Utqiagvik_2014_2016_Canopy_Spectral_Reflectance.csv"))
orig_names <- names(Utqiagvik_2014_2016_Canopy_Reflectance)
Utqiagvik_2014_2016_Canopy_Reflectance$Location <- rep("Utqiagvik",times=dim(Utqiagvik_2014_2016_Canopy_Reflectance)[1])
Utqiagvik_2014_2016_Canopy_Reflectance <- Utqiagvik_2014_2016_Canopy_Reflectance %>%
  select(Location,paste(orig_names))
head(Utqiagvik_2014_2016_Canopy_Reflectance)[,1:10]

NGEEArctic_leaf_reflectance <- rbind(Utqiagvik_2013_Leaf_Reflectance,Utqiagvik_2014_2016_Leaf_Reflectance,
                                     SewardPeninsual_2016_Leaf_Reflectance)
head(NGEEArctic_leaf_reflectance)[,1:10]
unique(NGEEArctic_leaf_reflectance$Location)
unique(NGEEArctic_leaf_reflectance$USDA_Species_Code)
unique(NGEEArctic_leaf_reflectance$Instrument)

NGEEArctic_canopy_reflectance <- rbind(Utqiagvik_2014_2016_Canopy_Reflectance)

NGEEArctic_Reflectance <- list(Leaf_Reflectance=NGEEArctic_leaf_reflectance,
                               Canopy_Reflectance=NGEEArctic_canopy_reflectance)
head(NGEEArctic_Reflectance$Leaf_Reflectance)[,1:5]
head(NGEEArctic_Reflectance$Canopy_Reflectance)[,1:5]

save(NGEEArctic_Reflectance, file = file.path(input_dir,"NGEEArctic_Leaf_and_Canopy_Reflectance.RData"))
#--------------------------------------------------------------------------------------------------#




