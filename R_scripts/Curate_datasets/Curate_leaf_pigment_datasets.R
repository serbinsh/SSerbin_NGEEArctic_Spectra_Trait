####################################################################################################
#
#    --- Last updated: 12.03.2020 By Shawn P. Serbin <sserbin@bnl.gov>
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
url <- "ftp://ngee.ornl.gov/data/outgoing/leaf_chlorophyll_carotenoids_barrow/data/NGEE_Arctic_chlorophyll_and_carotenoids_Barrow_2013_2015.csv"
dataset <- getURL(url, ssl.verifypeer = FALSE)
data_header <- read.csv(textConnection(dataset), skip = 7, nrows = 1, header=F)
head(data_header)
pigment_data <- read.csv(textConnection(dataset), skip = 9, header=F)
head(pigment_data)
names(pigment_data) <- data_header
head(pigment_data)

# clean up
pigment_data[pigment_data==-9999]=NA
head(pigment_data)

pigment_data <- pigment_data %>%
  select(Sample_ID,Sample_Date=Date,USDA_Species_Code,Chl_a_L,Chl_b_L,Carot_tot_L,Chl_a_area_L,
         Chl_b_area_L,Chl_ab_area_L,Chl_a_b_ratio_L=`Chl a_b_ratio_L`,Carot_tot_area_L,Chl_a_P,Chl_b_P,Chl_a_area_P,
         Chl_b_area_P,Chl_ab_area_P,Chl_a_b_ratio_P,Chl_a_mol_P)
head(pigment_data)


write.csv(pigment_data, 
          file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_leaf_pigments_2013_and_2015.csv"), 
          row.names = F)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Combine datasets
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
barrow_2013_2015_pigment_data <- read.csv(file = file.path(input_dir,
                                                           "NGEE-Arctic_Utqiagvik_leaf_pigments_2013_and_2015.csv"))
orig_names <- names(barrow_2013_2015_pigment_data)
head(barrow_2013_2015_pigment_data)
barrow_2013_2015_pigment_data$Location <- rep("Utqiagvik",times=dim(barrow_2013_2015_pigment_data)[1])
barrow_2013_2015_pigment_data <- barrow_2013_2015_pigment_data %>%
  select(Location,paste(orig_names))
head(barrow_2013_2015_pigment_data)

Utqiagvik_2013_2015_leaf_pigment_data=barrow_2013_2015_pigment_data
save(Utqiagvik_2013_2015_leaf_pigment_data, file = file.path(input_dir,"NGEEArctic_Leaf_Pigments.RData"))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#