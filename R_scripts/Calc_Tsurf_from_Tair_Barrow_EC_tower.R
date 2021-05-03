####################################################################################################
#
#    --- Last updated: 05.03.2021 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

library(lubridate)
library(readxl)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set ouput directory
out.dir <- '~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/R_Output/Tsurf_from_Tair/'
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Read in EC dataset
data_source <- '~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/SSerbin_NGEEArctic_Spectra_Trait/Data/EC_data/Barrow-2015/'
ec_data_all <- read_xlsx(path = file.path(data_source,'NGEE-Barrow-2015-b1.xlsm'),sheet = "Data",
                         skip = 4)
names(ec_data_all)
names(ec_data_all)[91]
#ec_temp_data <- ec_data_all[,c("X__2","X__4","(umol m-2 s-1)__4","(W m-2)__8","( C )__1","( C )__2")]
#ec_temp_data <- ec_data_all[,c("...2","...4","(umol m-2 s-1)..4","(W m-2)..8","( C )..1","( C )..2")]
ec_temp_data <- ec_data_all[,c(2,4,85,73,115,124)]
names(ec_temp_data) <- c("Timestamp","Frac_DOY","PAR_umol_m2_s1","Solar_W_m2","Tair_degC","Tsurf_degC")

ec_temp_data[ec_temp_data==-9999] <- NA
ec_temp_data$Tair_degC[ec_temp_data$Tair_degC>40] <- NA
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Calculate relationship
reg <- lm(formula = ec_temp_data$Tsurf_degC~ec_temp_data$Tair_degC)
summary(reg)

png(file=file.path(out.dir,'2015_Barrow_EC_Tsurf_vs_Tair_comparison.png'),height=2500,
    width=3500, res=340)
par(mfrow=c(1,1), mar=c(4.5,5.7,0.3,0.4), oma=c(0.3,0.9,0.3,0.1)) # B, L, T, R
plot(ec_temp_data$Tair_degC,ec_temp_data$Tsurf_degC, pch=21,bg="grey80",
     xlab="Tair (degC)", ylab="Tsurf (degC)")
abline(0,1,lty=2, col="grey40", lwd=2)
#mtext(text = "Tsurf = 1.06*Tair + 0.99")
legend("topleft", legend = c("Tsurf = 1.06*Tair + 0.99"), bty="n",
       cex=1.5)
box(lwd=2.2)
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Update time, write out data for later use
ec_timestamp <- ec_temp_data$Timestamp
ec_date <- as.Date(as.numeric(substr(ec_timestamp,5,7))-1, origin="2015-01-01")
ec_time <- as.numeric(substr(ec_timestamp,8,11))
ec_dateime <- force_tz(seq.POSIXt(ISOdatetime(2015,1,1,0,0,0), 
                                  ISOdatetime(2015,12,31,23,59,59), 
                                  by=(60*30), tz="UTC"), tzone="UTC")
output_data <- data.frame(Date=as_date(ec_dateime), Date_Time = ec_dateime,
                          Frac_DOY=ec_temp_data$Frac_DOY-1, 
                          PAR_umol_m2_s1=ec_temp_data$PAR_umol_m2_s1,
                          Solar_W_m2=ec_temp_data$Solar_W_m2,
                          Tair_degC=ec_temp_data$Tair_degC, 
                          Tsurf_degC=ec_temp_data$Tsurf_degC)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Subset to peak season
subset_data <- output_data[(as_date(output_data$Date_Time) >= "2015-06-01" & 
                              (as_date(output_data$Date_Time)) <= "2015-09-05"),]

plot(subset_data$Date_Time[subset_data$Date>="2015-07-20" & subset_data$Date<= "2015-07-22"],
     subset_data$PAR_umol_m2_s1[subset_data$Date>="2015-07-20" & subset_data$Date<= "2015-07-22"])
# def looks like UTC 

reg <- lm(formula = subset_data$Tsurf_degC~subset_data$Tair_degC)
summary(reg)

png(file=file.path(out.dir,'2015_Barrow_EC_Tsurf_vs_Tair_comparison_June-August.png'),height=2500,
    width=3500, res=340)
par(mfrow=c(1,1), mar=c(4.5,5.7,0.3,0.4), oma=c(0.3,0.9,0.3,0.1)) # B, L, T, R
plot(subset_data$Tair_degC,subset_data$Tsurf_degC, pch=21,bg="grey80",xlim=c(-5,25),ylim=c(-5,25),
     xlab="Tair (degC)", ylab="Tsurf (degC)")
abline(0,1,lty=2, col="grey40", lwd=2)
#mtext(text = "Tsurf = 1.06*Tair + 0.99")
legend("topleft", legend = c("Tsurf = 1.07*Tair + 1.26"), bty="n",
       cex=1.5)
box(lwd=2.2)
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#output data
write.csv(output_data,file = file.path(out.dir,"2015_Barrow_EC_tower_30min_Tair_Tsurf_data.csv"), 
          row.names = F)
write.csv(subset_data,file = file.path(out.dir,"2015_June-August_Barrow_EC_tower_30min_Tair_Tsurf_data.csv"), 
          row.names = F)
#--------------------------------------------------------------------------------------------------#
### EOF