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
load(file.path('~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/NGEEArctic_Leaf_Pigments.RData'))
names(Utqiagvik_2013_2015_leaf_pigment_data)
#Utqiagvik_2016_Aq_data
pigment_data <- Utqiagvik_2013_2015_leaf_pigment_data
head(pigment_data)
unique(pigment_data$Sample_Date)

### spectra
load(file.path('~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/data/compiled_data/NGEEArctic_Leaf_and_Canopy_Reflectance.RData'))
spectra <- rbind(NGEEArctic_Reflectance$Leaf_Reflectance)
head(spectra)[,1:10]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### pigment plots

Chl_ab_area_L_bp <- ggplot(data=pigment_data, aes(x = USDA_Species_Code, y = Chl_ab_area_L)) + geom_boxplot(alpha=I(.7)) + theme_bw()
Chl_ab_area_L_bp <- Chl_ab_area_L_bp + theme(axis.text=element_text(size=18),legend.text=element_text(size=18),legend.title=element_text(size=0),
                         axis.title.y=element_text(size=18,face="bold"), axis.title.x=element_blank(),
                         axis.text.x=element_text(size=18,face="bold"), axis.ticks.x = element_blank(),
                         legend.position="none")
Chl_ab_area_L_bp

Chl_ab_area_P_bp <- ggplot(data=pigment_data, aes(x = USDA_Species_Code, y = Chl_ab_area_P)) + geom_boxplot(alpha=I(.7)) + theme_bw()
Chl_ab_area_P_bp <- Chl_ab_area_P_bp + theme(axis.text=element_text(size=18),legend.text=element_text(size=18),legend.title=element_text(size=0),
                                             axis.title.y=element_text(size=18,face="bold"), axis.title.x=element_blank(),
                                             axis.text.x=element_text(size=18,face="bold"), axis.ticks.x = element_blank(),
                                             legend.position="none")
Chl_ab_area_P_bp



#--------------------------------------------------------------------------------------------------#






