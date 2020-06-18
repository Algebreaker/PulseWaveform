#p=pulse2?
#wuvn = wuv
#osnd x /y
#poly_wave

find_sd_sine <- function(p, wuvn, osndx, osndy, pw, plot=FALSE){
  x. = 1:(length(p$x))
  s_sine <- list()
  d_sine <- list()
  
  for(i in 1:length(osndx)){
    period <- 3*(wuvn$v_x[i]-wuvn$u_x[i])
    
    phi <- ((2*pi)/period)  * (x.-(wuvn$w[i]+(period/4)))
    phi[phi>(pi)] <- pi
    phi[phi<(-pi)] <- -pi
    
    # Finding first Sine
    y <-((osnd_y[[c(i, 2)]] -  osnd_y[[c(i, 1)]])/2) * cos(phi) + ((osnd_y[[c(i,1)]] + osnd_y[[c(i, 2)]])/2) 
    #the last bit (s/2) only works if O is always 0 - otherwise should be (O+S)/2 
    
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