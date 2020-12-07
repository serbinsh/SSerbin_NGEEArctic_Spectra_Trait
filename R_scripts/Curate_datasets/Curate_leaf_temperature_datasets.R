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
library(dplyr)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### output location
output_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
temp_data_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/raw_data/")
temp_data <- read.csv(file = file.path(temp_data_dir,"NGEE-Arctic_Barrow_2015_2016_spectra_leaf_temperatures.csv"), header=T)
head(temp_data)

temp_data[temp_data==-9999]=NA
temp_data <- temp_data %>%
  select(Sample_ID,Sample_Date=Measurement_Date,USDA_Species_Code,Air_Temperature_degC,
         Leaf_Temperature_PreSpec_degC,Leaf_Temperature_PostSpec_degC)
head(temp_data)  

write.csv(temp_data, 
          file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_spectra_leaf_temperatures_2015_2016.csv"), 
          row.names = F)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Combine into a single dataset
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
barrow_2015_2016_leaf_temp_data <- read.csv(file = file.path(input_dir,
                                                             "NGEE-Arctic_Utqiagvik_spectra_leaf_temperatures_2015_2016.csv"))
orig_names <- names(barrow_2015_2016_leaf_temp_data)
barrow_2015_2016_leaf_temp_data$Location <- rep("Utqiagvik",times=dim(barrow_2015_2016_leaf_temp_data)[1])
barrow_2015_2016_leaf_temp_data <- barrow_2015_2016_leaf_temp_data %>%
  select(Location,paste(orig_names))
head(barrow_2015_2016_leaf_temp_data)

# currently missing SewPen leaf temp data
leaf_temperature_data <- data.frame(barrow_2015_2016_leaf_temp_data)
save(leaf_temperature_data, file = file.path(input_dir,"NGEEArctic_Leaf_Temperature.RData"))
#--------------------------------------------------------------------------------------------------#

