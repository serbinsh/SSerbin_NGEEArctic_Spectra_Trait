####################################################################################################
#
#    --- Last updated: 12.02.2020 By Shawn P. Serbin <sserbin@bnl.gov>
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
## for scaling to spec temp
arrhenius <- function(Tleaf.1,Tleaf.2,VJ,Ev) {
  R <- 0.008314472
  VJ.T2 <- VJ*exp((Ev*((Tleaf.2+273.15)-(Tleaf.1+273.15)))/((Tleaf.1+273.15)*R*(Tleaf.2+273.15)))
  return(VJ.T2)
} 
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### output location
output_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Pull down NGEE-Arctic gas exchange datasets. Fitted and 1pt data

# Utqiagvik 2012-2015 data
url <- "ftp://ngee.ornl.gov/data/outgoing/leaf_photosynthetic_parameters_Vcmax_Jmax/data/NGEE-Arctic_Vcmax_Jmax_2012-2015.csv"
dataset <- getURL(url, ssl.verifypeer = FALSE)
data_header <- read.csv(textConnection(dataset), skip = 8, nrows = 1, header=F)
head(data_header)
gasex_data <- read.csv(textConnection(dataset), skip = 10, header=F)
head(gasex_data)
names(gasex_data) <- data_header
head(gasex_data)

# clean up
gasex_data[gasex_data==-9999]=NA
head(gasex_data)
gasex_data$Sample_ID <- lapply(gasex_data$Sample_Barcode, function(x) paste0('BNL', x))
gasex_data$Sample_ID <- as.character(gasex_data$Sample_ID)
head(gasex_data)

gasex_data <- gasex_data %>%
  mutate(Vcmax25_Rog=arrhenius(Mean_Tleaf,25,Vcmax_Tleaf,54.26)) %>%
  mutate(Jmax25_Rog=arrhenius(Mean_Tleaf,25,Jmax_Tleaf,36.21)) %>%
  mutate(Vcmax25_Bern=arrhenius(Mean_Tleaf,25,Vcmax_Tleaf,65.33)) %>%
  mutate(Jmax25_Bern=arrhenius(Mean_Tleaf,25,Jmax_Tleaf,43.9)) %>%
  mutate(Location=rep("Utqiagvik",times=dim(gasex_data)[1])) %>%
  select(Location,Sample_ID,Sample_Date,USDA_Species_Code,Temp_Response_flagged_with_1,Replicate,Tair=Mean_Tair,
         Tleaf=Mean_Tleaf,VpdL=Mean_VpdL,Vcmax_Tleaf,Jmax_Tleaf,Rd_Tleaf,Vcmax25_Rog,Jmax25_Rog,
         Vcmax25_Bern,Jmax25_Bern)
head(gasex_data)
hist(gasex_data$Vcmax25_Rog)
plot(gasex_data$Tleaf,gasex_data$Vcmax_Tleaf)

write.csv(gasex_data, file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_Vcmax_Jmax_2012-2015.csv"), row.names = F)

rm(gasex_data)

# 2016 data
url <- "ftp://ngee.ornl.gov/data/outgoing/leaf_photosynthetic_parameters_Vcmax_Jmax/data/NGEE_Arctic_Vcmax_1point_2016.csv"
dataset <- getURL(url, ssl.verifypeer = FALSE)
data_header <- read.csv(textConnection(dataset), skip = 8, nrows = 1, header=F)
head(data_header)
gasex_data <- read.csv(textConnection(dataset), skip = 10, header=F)
head(gasex_data)
names(gasex_data) <- data_header
head(gasex_data)

# cleanup all of the empty rows and columns
gasex_data <- gasex_data[, colMeans(is.na(gasex_data)) != 1]
gasex_data <- gasex_data[rowSums(is.na(gasex_data)) == 0, ]

# other clean up
gasex_data[gasex_data==-9999]=NA
gasex_data$Sample_ID <- lapply(gasex_data$Sample_ID, function(x) paste0('BNL', x))
gasex_data$Sample_ID <- as.character(gasex_data$Sample_ID)
head(gasex_data)

gasex_data <- gasex_data %>%
  mutate(Vcmax25_Rog=arrhenius(Tleaf,25,VcmaxT_one_point,54.26)) %>%
  mutate(Vcmax25_Bern=arrhenius(Tleaf,25,VcmaxT_one_point,65.33)) %>%
  mutate(Location=rep("Utqiagvik",times=dim(gasex_data)[1])) %>%
  select(Location, Sample_ID,Sample_Date=Date,USDA_Species_Code,Asat,Tleaf,Ci,VcmaxT_one_point,Vcmax25_Rog,
         Vcmax25_Bern)
head(gasex_data)
hist(gasex_data$Vcmax25_Rog)
hist(gasex_data$Vcmax25_Bern)
plot(gasex_data$Tleaf,gasex_data$VcmaxT_one_point)

write.csv(gasex_data, file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_Vcmax_1point_2016.csv"), row.names = F)

rm(data_header,gasex_data)

#url <- "ftp://ngee.ornl.gov/data/outgoing/NGA175_leaf_photosynthetic_param_barrow_2016/data/NGA175_AQfittedParams_Barrow_2016.csv"
#dataset <- getURL(url, ssl.verifypeer = FALSE)
#data_header <- read.csv(textConnection(dataset), skip = 6, nrows = 1, header=F)
#head(data_header)
#gasex_data <- read.csv(textConnection(dataset), skip = 8, header=F)
#head(gasex_data)
#names(gasex_data) <- data_header
#head(gasex_data)

url <- "https://ngee.ornl.gov/ngeedata/NGA175/data/NGA175_AQfittedParams_Barrow_2016_v2.csv"
dataset <- getURL(url, ssl.verifypeer = FALSE)
gasex_data <- read.csv(textConnection(dataset), header=T)
head(gasex_data)
rm(dataset,url)

# clean up
gasex_data[gasex_data==-9999]=NA
gasex_data$Sample_ID <- lapply(gasex_data$Sample_ID, function(x) paste0('BNL', x))
gasex_data$Sample_ID <- as.character(gasex_data$Sample_ID)
head(gasex_data)

gasex_data <- gasex_data %>%
#  select(Sample_ID,Sample_Date,USDA_Species_Code,Tleaf=Mean_Tleaf,VpdL=Mean_VpdL,CO2S=Mean_CO2S,
#         PhiCO2=PhiCO2.i,Convex,Asat=Asat.g,Rd)
  mutate(Location=rep("Utqiagvik",times=dim(gasex_data)[1])) %>%
  select(Location,Sample_ID,Sample_Date,USDA_Species_Code,Tleaf=Mean_Tleaf,RHs=Mean_RHs,CO2s=Mean_CO2s,Ci=Mean_Ci,
         PhiCO2i,AsatG,Theta,Rday)
head(gasex_data)

plot(gasex_data$Tleaf,gasex_data$Asat)
plot(gasex_data$PhiCO2,gasex_data$Asat)

write.csv(gasex_data, file = file.path(output_dir,"NGEE-Arctic_Utqiagvik_AQfittedParams_2016.csv"), row.names = F)

rm('gasex_data')

# SewPen 2019 data
url <- "https://ngee.ornl.gov/ngeedata/NGA215_plant_phys_SP_2019/data/Seward_2019_Vcmax_Jmax.xlsx"
dataset <- tempfile()
getBinaryURL(url, ftp.use.epsv = FALSE,crlf =TRUE)%>%writeBin(con=dataset)
gasex_data <- readxl::read_xlsx(dataset, sheet = 2)
head(gasex_data)

# clean up
gasex_data[gasex_data==-9999]=NA
gasex_data <- gasex_data %>%
  mutate(Vcmax25_Rog=arrhenius(Mean_Tleaf,25,Vcmax_Tleaf,54.26)) %>%
  mutate(Jmax25_Rog=arrhenius(Mean_Tleaf,25,Jmax_Tleaf,36.21)) %>%
  mutate(Vcmax25_Bern=arrhenius(Mean_Tleaf,25,Vcmax_Tleaf,65.33)) %>%
  mutate(Jmax25_Bern=arrhenius(Mean_Tleaf,25,Jmax_Tleaf,43.9)) %>%
  mutate(Location=rep("Seward_Peninsula",times=dim(gasex_data)[1])) %>%
  select(Location,Sample_ID,Sample_Date,USDA_Species_Code=Species,Tair=Mean_Tair,Tleaf=Mean_Tleaf,
         VpdL=Mean_VpdL,Vcmax_Tleaf,Jmax_Tleaf,Rd_Tleaf,Vcmax25_Rog,Jmax25_Rog,
         Vcmax25_Bern,Jmax25_Bern)
head(gasex_data)
hist(gasex_data$Vcmax_Tleaf)
hist(gasex_data$Vcmax25_Rog)
plot(gasex_data$Tleaf,gasex_data$Vcmax_Tleaf)

write.csv(gasex_data, file = file.path(output_dir,"NGEE-Arctic_SewardPeninsula_Vcmax_Jmax_2019.csv"), 
          row.names = F)

rm(gasex_data)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Combine datasets
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
barrow_2016_1pt_Vcmax_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_Utqiagvik_Vcmax_1point_2016.csv"))
head(barrow_2016_1pt_Vcmax_data)

barrow_2016_Aq_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_Utqiagvik_AQfittedParams_2016.csv"))
head(barrow_2016_Aq_data)

barrow_2012_2015_VJ_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_Utqiagvik_Vcmax_Jmax_2012-2015.csv"))
head(barrow_2012_2015_VJ_data)

sewpen_2019_VJ_data <- read.csv(file = file.path(input_dir,"NGEE-Arctic_SewardPeninsula_Vcmax_Jmax_2019.csv"))
head(sewpen_2019_VJ_data)

nga_gasex_data <- list(Utqiagvik_2016_1pt_Vcmax=barrow_2016_1pt_Vcmax_data,
                       Utqiagvik_2016_Aq_data=barrow_2016_Aq_data,
                       Utqiagvik_2012_2015_VcmaxJmax_data=barrow_2012_2015_VJ_data,
                       SewPen_2019_VcmaxJmax_data=sewpen_2019_VJ_data)
names(nga_gasex_data)

save(nga_gasex_data, file = file.path(input_dir,"NGEEArctic_GasEx.RData"))
#--------------------------------------------------------------------------------------------------#


### EOF