####################################################################################################
#
#    --- Last updated: 12.29.2020 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
list.of.packages <- c("devtools","dplyr","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))

# Load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))

# not in
`%notin%` <- Negate(`%in%`)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set ouput directory
out.dir <- file.path('~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/R_Output/Data_exploration/')
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Trait data
load(file.path('~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/NGEEArctic_GasEx.RData'))
names(nga_gasex_data)
#Utqiagvik_2016_Aq_data
aq_data <- nga_gasex_data$Utqiagvik_2016_Aq_data
head(aq_data)
unique(aq_data$Sample_Date)

### spectra
load(file.path('~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/NGEEArctic_Leaf_and_Canopy_Reflectance.RData'))
spectra <- rbind(NGEEArctic_Reflectance$Leaf_Reflectance)
head(spectra)[,1:10]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Merge datasets
#spectra_trait_data <- merge(aq_data, spectra,by="Sample_ID")
spectra_trait_data <- merge(aq_data, spectra)
head(spectra_trait_data)[1:15]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
pri <- (spectra_trait_data$Wave_571-spectra_trait_data$Wave_530) / (spectra_trait_data$Wave_571+spectra_trait_data$Wave_530)
hist(pri)
pri <- (spectra_trait_data$Wave_531-spectra_trait_data$Wave_570) / (spectra_trait_data$Wave_531+spectra_trait_data$Wave_570)
hist(pri)


RGR <- (spectra_trait_data$Wave_600-spectra_trait_data$Wave_699) / 
  (spectra_trait_data$Wave_500+spectra_trait_data$Wave_599)
hist(RGR)


plot(pri+1/2,spectra_trait_data$PhiCO2i)
plot(pri+1/2,spectra_trait_data$AsatG)
plot(pri+1/2,spectra_trait_data$Theta)
plot(pri+1/2,spectra_trait_data$Rday)


plot(spectra_trait_data$Wave_760,spectra_trait_data$PhiCO2i)
plot(spectra_trait_data$Wave_720,spectra_trait_data$aQY)


spec <- spectra_trait_data[,which(names(spectra_trait_data) %in% paste0("Wave_", seq(350,2500,1)))]
corr <- cor(spec,spectra_trait_data$PhiCO2i, use = "pairwise.complete.obs")
plot(seq(350,2500,1),corr)
row.names(corr)[which(corr==min(corr))]

plot(spectra_trait_data$Wave_1387,spectra_trait_data$PhiCO2i)

corr <- cor(spec,spectra_trait_data$AsatG, use = "pairwise.complete.obs")
plot(seq(350,2500,1),corr)
row.names(corr)[which(corr==max(corr))]

plot(spectra_trait_data$Wave_1956,spectra_trait_data$AsatG)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
temp <- spectra_trait_data %>%
  group_by(USDA_Species_Code) %>%
  mutate(pri= (Wave_571-Wave_530) / (Wave_571+Wave_530)) %>%
  summarise_at(c("PhiCO2i","AsatG","Theta","Rday","pri"), c(mean=mean,sd=sd), na.rm = TRUE)

temp

plot(temp$pri_mean,temp$AsatG_mean)
#--------------------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------------------#
### EOF