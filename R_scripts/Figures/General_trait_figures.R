####################################################################################################
#
#    --- Last updated: 12.07.2020 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
list.of.packages <- c("devtools","remotes","readr","RCurl","httr","dplyr","reshape2","here",
                      "ggplot2","gridExtra", "lubridate")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
                                                                       
# Load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))

# not in
`%notin%` <- Negate(`%in%`)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### output location
output_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/figures_and_tables/leaf_traits")
if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
data_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/")
load(file.path(data_dir,'NGEEArctic_CHN_LMA.RData'))
head(chn_lma_data)

load(file.path(data_dir,'NGEEArctic_GasEx.RData'))
head(nga_gasex_data$Utqiagvik_2012_2015_VcmaxJmax_data)

load(file.path(data_dir,'NGEEArctic_Leaf_Pigments.RData'))
head(Utqiagvik_2013_2015_leaf_pigment_data)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Figures
# years in data file
unique(lubridate::year(as.Date(as.character(chn_lma_data$Sample_Date), format="%Y%m%d")))
# 2012 2013 2014 2015 2016

chn_lma_data <- chn_lma_data %>%
  filter(lubridate::year(as.Date(as.character(chn_lma_data$Sample_Date), format="%Y%m%d"))>2012) %>%
  filter(USDA_Species_Code %notin% c("ARRU","CHLA13","SAALA"))
  
LMA.bp <-  ggplot(data=chn_lma_data, aes(x = USDA_Species_Code, y = LMA_gDW_m2, fill="white")) + theme_bw()
LMA.bp <- LMA.bp + geom_boxplot(alpha=I(.7)) + labs(x="", y=expression(paste("LMA"," (",g~m^{-2},")")))
LMA.bp <- LMA.bp + theme(axis.text=element_text(size=18),legend.text=element_text(size=18),
                         legend.title=element_text(size=0),
                         axis.title.y=element_text(size=18,face="bold"), 
                         axis.title.x=element_blank(),
                         axis.text.x=element_text(size=18,face="bold", angle = 90),
                         axis.ticks.x = element_blank(),legend.position="none")
LMA.bp

Nmass.bp <-  ggplot(data=chn_lma_data, aes(x = USDA_Species_Code, y = Nmass_mg_g*0.1, fill="white")) + 
  theme_bw() 
Nmass.bp <- Nmass.bp + geom_boxplot(alpha=I(.7)) + labs(x="", y=expression(paste("N[mass]")))
Nmass.bp <- Nmass.bp + theme(axis.text=element_text(size=18),legend.text=element_text(size=18),
                         legend.title=element_text(size=0),
                         axis.title.y=element_text(size=18,face="bold"), 
                         axis.title.x=element_blank(),
                         axis.text.x=element_text(size=18,face="bold", angle = 90),
                         axis.ticks.x = element_blank(),legend.position="none")
Nmass.bp


#--------------------------------------------------------------------------------------------------#


