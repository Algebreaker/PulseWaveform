#p=pulse2
#ss = s_sines
#ds = d_sines
#wuvn = wuv
#osndx = osnd_x
#osndy = osnd_y

fit_n_sine <- function(p, ss, ds, osndx, osndy, wuvn, plot=FALSE){
  x. = 1:(length(p$x))
  n_sine <- list()
  resid_test <- list()
  resid_peaks_x <- c()
  resid_peaks_y <- c()
  const <- c()
  
  for(i in 2:length(p)){
    
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
    
    resid_test[[i-1]] <- d_y_resid
    

    period <- 3*(wuvn$v_x[i-1]-wuvn$u_x[i-1])
  
    pks <- findpeaks(d_y_resid[osndx[[c(i-1,2)]]:osndx[[c(i-1,4)]]])
    pk_loc <- which(abs(osnd_x[[c(i-1,3)]] - (osndx[[c(i-1,2)]] + pks[,2]))==min(abs(osnd_x[[c(i-1,3)]] - (osndx[[c(i-1,2)]] + pks[,2]))))
    resid_peaks_x[i-1] <- osndx[[c(i-1,2)]] + pks[pk_loc,2]
    resid_peaks_y[i-1] <- pks[pk_loc,1]
       
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
    
    resid <- resid_test[[i-1]] - n_y
    
    #n_resid[[i-1]] <- resid
    
    
  }
  return(n_sine)
}