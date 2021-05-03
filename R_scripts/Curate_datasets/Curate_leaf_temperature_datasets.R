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
# Barrow data
temp_data_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/raw_data/")
temp_data <- read.csv(file = file.path(temp_data_dir,"NGEE-Arctic_Barrow_2015_2016_spectra_leaf_temperatures.csv"), header=T)
head(temp_data)

# cleanup
temp_data[temp_data==-9999]=NA
#orig_names <- names(barrow_2015_2016_leaf_temp_data)
#temp_data$Location <- rep("Utqiagvik",times=dim(temp_data)[1])
#temp_data <- temp_data %>%
#  select(Location,paste(orig_names))
#head(temp_data)

temp_data <- temp_data %>%
  mutate(Location=rep("Utqiagvik",times=dim(temp_data)[1])) %>%
  select(Location,Sample_ID,Sample_Date=Measurement_Date,USDA_Species_Code,Air_Temperature_degC,
         Leaf_Temperature_PreSpec_degC,Leaf_Temperature_PostSpec_degC)
head(temp_data)  

write.csv(temp_data, 
          file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_spectra_leaf_temperatures_2015_2016.csv"), 
          row.names = F)
rm(temp_data,temp_data_dir)

# SewPen 2019 data - to be used with shrub Vcmax data
temp_data_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/raw_data/")
temp_data <- read.csv(file = file.path(temp_data_dir,"Seward_2019_VegTemperature.csv"), header=T)
head(temp_data)

temp_data[temp_data==-9999]=NA
temp_data <- temp_data %>%
  mutate(Location = rep("Seward_Peninsula")) %>%
  mutate(Tcanopy_minus_Tair = Tcanopy-TairCanopy) %>%
  mutate(Tleaf_minus_Tair = Tleaf-TairLeaf) %>%
  select(Location,Sample_ID=SampleID,Sample_Date=Date,USDA_Species_Code=Species,Tcanopy_degC=Tcanopy,
         TairCanopy_degC=TairCanopy, Tleaf_degC=Tleaf,TairLeaf_degC=TairLeaf,
         Tleaf_minus_Tair,Tcanopy_minus_Tair)
head(temp_data)  

write.csv(temp_data, 
          file = file.path(output_dir,"NGEE-Arctic_sewpen_spectra_leaf_temperatures_2019.csv"), 
          row.names = F)
rm(temp_data,temp_data_dir)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Combine into a single dataset
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
barrow_2015_2016_leaf_temp_data <- read.csv(file = file.path(input_dir,
                                                             "NGEE-Arctic_Utqiagvik_spectra_leaf_temperatures_2015_2016.csv"))
head(barrow_2015_2016_leaf_temp_data)

input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
sewpen_2019_temp_data <- read.csv(file = file.path(input_dir,
                                                             "NGEE-Arctic_sewpen_spectra_leaf_temperatures_2019.csv"))
head(sewpen_2019_temp_data)

# combine into a single Rdata file
temperature_data <- list(Utqiagvik_data=barrow_2015_2016_leaf_temp_data,
                         SewPen_data = sewpen_2019_temp_data)
save(temperature_data, file = file.path(input_dir,"NGEEArctic_Leaf_and_Canopy_Temperature.RData"))
#--------------------------------------------------------------------------------------------------#

