# calc the CH4 RF -- relative to the PI
RF_calc <- function(ch4) {
  # use the IPCC tables to calculate the RF of methane since the PI
  # Data from Table 8.SM.1
  M_conv = 2.5e19 
  alpha = 0.036 
  N0 = 332 # from https://www.esrl.noaa.gov/gmd/webdata/ccgg/trends/n2o_trend_all_gl.png
  M0 = 760
  M = (ch4/M_conv)*1E9
  fMN = function(M, N) {
    0.47*log(1+ (2.01E-5*M*N*0.75) + (5.31E-15*M*(M*N)^1.52) )
  }
  RF = alpha*(sqrt(M) - sqrt(M0) ) #- fMN(M, N0) - fMN(M0, N0)
  return(RF)
}

par(mfrow=c(1,3) )
# plot the methane as a mixing ratio
plot(out_data$time/t_step/12, 
     (out_data$CH4/M)*1E6, 
     ylab = "CH4 /ppm", 
     xlab = "Time since 2020 /years",
     type="l" )

#plot the H2 as a mixing ratio
plot(out_data$time/t_step/12, 
     (out_data$H2/M)*1E6, 
     ylab = "H2 /ppm", 
     xlab = "Time since 2020 /years",
     type="l" )


# plot the methane RF
plot(out_data$time/t_step/12, 
     RF_calc(out_data$CH4), 
     ylim=c(0, 10), 
     ylab = "RF from CH4 /Wm-2", 
     xlab = "Time since 2020 /years",
     type="l" )
abline(h=0.48, col="red", lty=2)

