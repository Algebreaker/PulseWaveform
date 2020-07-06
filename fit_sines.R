
find_sd_sine <- function(p, wuvn, osndx, osndy, pw, plot=FALSE){
  x. = 1:(length(p$x))
  s_sine <- list()
  d_sine <- list()
  
  for(i in 1:length(osndx)){
    #Define the period of the sine curve as 3*V-U 
    period <- 3*(wuvn$v_x[i]-wuvn$u_x[i])
    #period = 25
    
    #Define phi
    phi <- ((2*pi)/period)  * (x.-(wuvn$w[i]+(period/4)))
    #The phi values that are above pi and below -pi are set to those values - this means we specify only the sine curve we want 
    phi[phi>(pi)] <- pi
    phi[phi<(-pi)] <- -pi
    
    # Finding first Sine
    y <-((osnd_y[[c(i, 2)]] -  osnd_y[[c(i, 1)]])/2) * cos(phi) + ((osnd_y[[c(i,1)]] + osnd_y[[c(i, 2)]])/2) 
    #The last term ((S+O)/2) can be simplified to S/2 if O is always set to 0
    
    #Defining phi for the D peak
    d_phi <- ((2*pi)/period)  * (x.-(osnd_x[[c(i,4)]]))
    d_phi[d_phi>(pi)] <- pi
    d_phi[d_phi<(-pi)] <- -pi
    
    #Finding the sine for the D peak
    d_y <- ((osnd_y[[c(i,4)]]-osnd_y[[c(i, 1)]])/2) * cos(d_phi) + ((osnd_y[[c(i,1)]] + osnd_y[[c(i,4)]])/2)
    
    
    
    # Plot back on poly_wave[[i]]:
    if(plot){
      plot(pw[[i]])
      lines(x., y, col = 'red')
      lines(x., d_y, col='purple')
    }
    
    
    
    s_sine[[i]] <- y 
    d_sine[[i]] <- d_y
  }
  return(c(s_sine, d_sine))
}



fit_n_sine <- function(p, ss, ds, osndx, osndy, wuvn, plot=FALSE){
  x. = 1:(length(p$x))
  n_sine <- list()
  resid_test <- list()
  resid_peaks_x <- c()
  resid_peaks_y <- c()
  const <- c()
  
  for(i in 2:ncol(p)){
    
    #The constant needs to be added back in every time we subtract a sine from the PPG trace
    #Create a residual of the PPG trace after subtracting the S and D peaks 
    const <- rep(ss[[i-1]][[length(ss[[i-1]])]],length(ss[[i-1]]))
    d_y_resid <- p[,i] - ss[[i-1]]
    d_y_resid <- d_y_resid + const
    d_y_resid <- d_y_resid - ds[[i-1]]
    d_y_resid <- d_y_resid + const
    
    if(plot){
      plot(1:(nrow(p)), p[,i], type = 'l', ylab = 'Plotting residual of S and D sines')
      lines(1:(nrow(p)),ss[[i-1]], col='red')
      lines(1:(nrow(p)),ds[[i-1]], col='purple')
      lines(1:(nrow(p)),d_y_resid, col='green')
    }
    
    #Define the period as 3*(x-V - x-U)
    period <- 3*(wuvn$v_x[i-1]-wuvn$u_x[i-1])
    #period <- 25
    
    #Find the peaks on the residual between the x values of the S and D peaks
    pks <- findpeaks(d_y_resid[osndx[[c(i-1,2)]]:osndx[[c(i-1,4)]]])
    #If no peaks are found, stick with the original N-values
    if(is.null(pks)){
      resid_peaks_x[i-1] <- osndx[[c(i-1,3)]]
      resid_peaks_y[i-1] <- osndy[[c(i-1,3)]]
    #Otherwise, find the peak that is closest to our current estimate of the x-N value, and choose that one as the peak of the N-sine
    }else{
    pk_loc <- which(abs(osndx[[c(i-1,3)]] - (osndx[[c(i-1,2)]] + pks[,2]))==min(abs(osndx[[c(i-1,3)]] - (osndx[[c(i-1,2)]] + pks[,2]))))
    resid_peaks_x[i-1] <- osndx[[c(i-1,2)]] + pks[pk_loc,2]
    resid_peaks_y[i-1] <- pks[pk_loc,1]
    }
    
    #Define phi for the n-sine
    n_phi <- ((2*pi)/period)  * (x.-(resid_peaks_x[i-1]))
    n_phi[n_phi>(pi)] <- pi
    n_phi[n_phi<(-pi)] <- -pi
    
    #create notch sine curve
    n_y <- ((resid_peaks_y[i-1]-osndy[[c(i-1, 1)]])/2) * cos(n_phi) + ((osndy[[c(i-1,1)]]+resid_peaks_y[i-1])/2)
    
    if(plot){
      plot(1:(nrow(p)), p[,i], type = 'l', ylab = 'S/N/D sines')
      lines(1:(nrow(p)),ss[[i-1]], col='red')
      lines(1:(nrow(p)),ds[[i-1]], col='purple')
      lines(1:(nrow(p)), n_y, col='green')
    }
    
    n_sine[[i-1]] <- n_y
    
    #n_resid[[i-1]] <- resid
    
    
  }
  return(n_sine)
}
