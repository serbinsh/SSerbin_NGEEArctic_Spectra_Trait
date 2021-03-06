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
library(RCurl)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### output location
output_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Pull down NGEE-Arctic leaf pigment data

# Barrow
url <- "ftp://ngee.ornl.gov/data/outgoing/lma_leaf_carbon_nitrogen/data/NGEEArctic_BNL_leaf_C_N_LMA_2012-2016.csv"
dataset <- getURL(url, ssl.verifypeer = FALSE)
data_header <- read.csv(textConnection(dataset), skip = 8, nrows = 1, header=F)
head(data_header)
chn_data <- read.csv(textConnection(dataset), skip = 10, header=F)
head(chn_data)
names(chn_data) <- data_header
head(chn_data)

# clean up
chn_data[chn_data==-9999]=NA
head(chn_data)

chn_data <- chn_data %>%
  mutate(lwc = rep(NA,dim(chn_data)[1])) %>%
  select(Sample_ID,Sample_Date=Date,USDA_Species_Code,Cmass_mg_g=C_mass,Nmass_mg_g=N_mass,LMA_gDW_m2=LMA,LWC=lwc,
         Carea_g_m2=C_area,Narea_g_m2=N_area,CN_ratio)
head(chn_data)

write.csv(chn_data, 
          file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_leaf_C_N_LMA_2012-2016.csv"), 
          row.names = F)
rm(chn_data)

# SewPen
url <- "ftp://ngee.ornl.gov/data/outgoing/NGA107_BNL_lma_leaf_carbon_nitrogen_kougarok_teller/data/NGEEArctic_BNL_Seward_LMA_CN_2016.csv"
dataset <- getURL(url, ssl.verifypeer = FALSE)
data_header <- read.csv(textConnection(dataset), skip = 6, nrows = 1, header=F)
head(data_header)
chn_data <- read.csv(textConnection(dataset), skip = 8, header=F)
head(chn_data)
names(chn_data) <- data_header
head(chn_data)

# clean up
chn_data[chn_data==-9999]=NA
head(chn_data)

chn_data$Sample_ID <- lapply(chn_data$Sample_ID, function(x) paste0('BNL', x))
chn_data$Sample_ID <- as.character(chn_data$Sample_ID)
head(chn_data)

chn_data <- chn_data %>%
  mutate(lwc = rep(NA,dim(chn_data)[1])) %>%
  select(Sample_ID,Sample_Date=Date,USDA_Species_Code=Species_USDA,Cmass_mg_g=C_mass,Nmass_mg_g=N_mass,LMA_gDW_m2=LMA,
         LWC=lwc,Carea_g_m2=C_area,Narea_g_m2=N_area,CN_ratio=`C:N`)
head(chn_data)

write.csv(chn_data, 
          file = file.path(output_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2016.csv"), 
          row.names = F)
rm(chn_data)

# 2017 SewPen CHN/LMA - temporary approach
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/raw_data/")
chn_data <- read.csv(file = file.path(input_dir,"Seward_2017_CN_LMA.csv"),header = T, sep = ",")
head(chn_data)

# clean up
chn_data[chn_data==-9999]=NA
head(chn_data)

names(chn_data)[3] <- "Sample_ID"

#chn_data$Sample_ID <- lapply(chn_data$Sample_ID, function(x) paste0('BNL', x))
#chn_data$Sample_ID <- as.character(chn_data$Sample_ID)
#head(chn_data)


chn_data <- chn_data %>%
  mutate(lwc = rep(NA,dim(chn_data)[1])) %>%
  select(Sample_ID,Sample_Date=Date,USDA_Species_Code=Species,Cmass_mg_g=Cmass,Nmass_mg_g=Nmass,LMA_gDW_m2=LMA,
         LWC=lwc,Carea_g_m2=Carea,Narea_g_m2=Narea,CN_ratio=CNratio)
head(chn_data)

write.csv(chn_data, 
          file = file.path(output_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2017.csv"), 
          row.names = F)
rm(chn_data)

# 2018 SewPen CHN/LMA - temporary approach
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/raw_data/")
chn_data <- read.csv(file = file.path(input_dir,"Seward_2018_LMA_LWC_CN.csv"),header = T, sep = ",")
head(chn_data)

# clean up
chn_data[chn_data==-9999]=NA
head(chn_data)

chn_data <- chn_data %>%
  select(Sample_ID=SampleID,Sample_Date=Date,USDA_Species_Code=Species,Cmass_mg_g=Cmass,
         Nmass_mg_g=Nmass,LMA_gDW_m2=LMA,LWC,
         Carea_g_m2=Carea,Narea_g_m2=Narea,CN_ratio=CNratio)
head(chn_data)

write.csv(chn_data, 
          file = file.path(output_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2018.csv"), 
          row.names = F)
rm(chn_data)


# 2019 SewPen CHN/LMA - temporary approach
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/raw_data/")
chn_data <- read.csv(file = file.path(input_dir,"Seward_2019_LMA_LWC_CN.csv"),header = T, sep = ",")
head(chn_data)

# clean up
chn_data[chn_data==-9999]=NA
head(chn_data)

chn_data <- chn_data %>%
  mutate(CNratio=Cmass/Nmass) %>%
  select(Sample_ID=SampleID,Sample_Date=Date,USDA_Species_Code=Species,Cmass_mg_g=Cmass,
         Nmass_mg_g=Nmass,LMA_gDW_m2=LMA,LWC,
         Carea_g_m2=Carea,Narea_g_m2=Narea,CN_ratio=CNratio)
head(chn_data)

write.csv(chn_data, 
          file = file.path(output_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2019.csv"), 
          row.names = F)
rm(chn_data)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Combine into a single dataset
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
barrow_2012_2016_chn_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_Utqiagvik_leaf_C_N_LMA_2012-2016.csv"))
orig_names <- names(barrow_2012_2016_chn_data)
barrow_2012_2016_chn_data$Location <- rep("Utqiagvik",times=dim(barrow_2012_2016_chn_data)[1])
barrow_2012_2016_chn_data <- barrow_2012_2016_chn_data %>%
  select(Location,paste(orig_names))
head(barrow_2012_2016_chn_data)

sewpen_2016_chn_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2016.csv"))
orig_names <- names(sewpen_2016_chn_data)
sewpen_2016_chn_data$Location <- rep("Seward_Peninsula",times=dim(sewpen_2016_chn_data)[1])
sewpen_2016_chn_data <- sewpen_2016_chn_data %>%
  select(Location,paste(orig_names))
head(sewpen_2016_chn_data)

sewpen_2017_chn_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2017.csv"))
orig_names <- names(sewpen_2017_chn_data)
sewpen_2017_chn_data$Location <- rep("Seward_Peninsula",times=dim(sewpen_2017_chn_data)[1])
sewpen_2017_chn_data <- sewpen_2017_chn_data %>%
  select(Location,paste(orig_names))
head(sewpen_2017_chn_data)

sewpen_2018_chn_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2018.csv"))
orig_names <- names(sewpen_2018_chn_data)
sewpen_2018_chn_data$Location <- rep("Seward_Peninsula",times=dim(sewpen_2018_chn_data)[1])
sewpen_2018_chn_data <- sewpen_2018_chn_data %>%
  select(Location,paste(orig_names))
head(sewpen_2018_chn_data)

sewpen_2019_chn_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_SewardPeninsula_C_N_LMA_2019.csv"))
orig_names <- names(sewpen_2018_chn_data)
sewpen_2019_chn_data$Location <- rep("Seward_Peninsula",times=dim(sewpen_2019_chn_data)[1])
sewpen_2019_chn_data <- sewpen_2019_chn_data %>%
  select(Location,paste(orig_names))
head(sewpen_2019_chn_data)


chn_lma_data <- rbind(barrow_2012_2016_chn_data,sewpen_2016_chn_data,sewpen_2017_chn_data,
                      sewpen_2018_chn_data,sewpen_2019_chn_data)
head(chn_lma_data)

save(chn_lma_data, file = file.path(input_dir,"NGEEArctic_CHN_LMA.RData"))
#--------------------------------------------------------------------------------------------------#


