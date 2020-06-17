p_w = poly_wave
p = pulse2
col_len = source_data_column_length
wuvn = wuv
ow = o_w_difference


find_osnd <- function(p, p_w, col_len, wuvn, ow, plot=FALSE){
  x_osnd <- list()
  osnd <- list()
  s <- c()
  s_yval <- c()
  
  for(i in 2:(ncol(p))){        
  
  sfunction <- splinefun(1:(col_len+4), p[, i], method = "natural")
  deriv1_wave <- sfunction(seq(1, (col_len+4)), deriv = 1)
  deriv1_wave_poly <- CubicInterpSplineAsPiecePoly(1:(col_len+4), deriv1_wave, "natural") 
  
  # Find inflexion points on deriv1_wave_poly
  inflexion_points_deriv1_wave_poly <- solve(deriv1_wave_poly, b = 0, deriv = 1)
  inflexion_points_deriv1_wave_poly_yval <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)
  #if(plot){
   # plot(deriv1_wave_poly)
  #  points(inflexion_points_deriv1_wave_poly, inflexion_points_deriv1_wave_poly_yval, pch = 19)
  #}

  
  # Find OSND
  inflexion_points_new <- solve(p_w[[i-1]], b = 0, deriv = 1)
  inflexion_points_new_yval <- predict(p_w[[i-1]], inflexion_points_new)
  o <- wuvn$w[i-1] - o_w_difference[i-1]
  o_yval <- predict(p_w[[i-1]], o)
  
  if(length(inflexion_points_new) >= 4){
    
    #Find the four inflexion points that are 1 to the left and 3 to the right of W
    old_o <- max(which(inflexion_points_new < wuvn$w[i-1]))
    # Find x coords of OSND
    x_osnd[[i-1]] <- inflexion_points_new[old_o:(old_o+3)]
    # Replace old o with new o                 # the new o is already confirmed as the correct one to use from a previous for loop          
    x_osnd[[c(i-1, 1)]] <- o
    
    # Find new S
    s[i-1] <- wuvn$w[i-1] + (2*(wuvn$v_x[i-1] - wuvn$w[i-1]))
    s_yval[i-1] <- predict(p_w[[i-1]], s[i-1])
    
    # Continue to define OSND (divergence here between waveforms based on if Paul type or not)
    # if s - w is greater than new s - w, aka if its a Paul type:
    if((inflexion_points_new[old_o+1] - wuvn$w[i-1]) > (s[i-1] - wuvn$w[i-1])){
      
      # Remove false N and D values 
      x_osnd_precursor <- x_osnd[[i-1]]
      if(length(which(complete.cases(x_osnd_precursor) ==1)) != 4){
        x_osnd_precursor <- x_osnd_precursor[-(which(is.na(x_osnd_precursor)))]
      }
      if(length(x_osnd_precursor) == 2){
        x_osnd[[i-1]] <- x_osnd_precursor
      }else{
        false_points <- which(x_osnd_precursor > wuvn$notch_x[i-1]) 
        x_osnd_precursor <- x_osnd_precursor[-false_points]
        x_osnd[[i-1]] <- x_osnd_precursor
      }
      # Find y coords of OSND
      osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
      osnd[[c(i-1, 1)]] <- o_yval
      
      # Remove false values again
      osnd_precursor <- osnd[[i-1]]
      if(length(which(complete.cases(osnd_precursor) ==1)) != 4){
        osnd_precursor <- osnd_precursor[-(which(is.na(osnd_precursor)))]
      }
      if(length(osnd_precursor) == 2){
        osnd[[i-1]] <- osnd_precursor
      }else{
        osnd_precursor <- osnd_precursor[-false_points]
        osnd[[i-1]] <- osnd_precursor
      }
      
      osnd[[c(i-1, 3)]] <- osnd[[c(i-1, 2)]]
      osnd[[c(i-1, 2)]] <- s_yval[i-1]
      x_osnd[[c(i-1, 3)]] <- x_osnd[[c(i-1, 2)]]
      x_osnd[[c(i-1, 2)]] <- s[i-1]
      osnd[[c(i-1, 4)]] <- wuvn$notch_y[i-1]
      x_osnd[[c(i-1, 4)]] <- wuvn$notch_x[i-1]
      
      
    }else{
      osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
      osnd[[c(i-1, 1)]] <- o_yval
      #osnd[[i-1]] <- osnd[[i-1]] - inflexion_points_new_yval[o]
      
      if((x_osnd[[c(i-1, 4)]] - x_osnd[[c(i-1, 3)]]) < 3 & (osnd[[c(i-1, 3)]] / (osnd[[c(i-1, 2)]] - osnd[[c(i-1, 1)]])) > 0.5){         # these lines are for cases when the canonical waveform has a prominent
        osnd[[c(i-1, 3)]] <- inflexion_points_new_yval[old_o + 4]                                                                         # renal wave such that o+3 and o+4 are inflexion points between the 
        x_osnd[[c(i-1, 3)]] <- inflexion_points_new[old_o + 4]                                                                            # systolic peak and the notch. This shifts them along by two. 
        osnd[[c(i-1, 4)]] <- inflexion_points_new_yval[old_o + 5]                                                                         # A fixed threshold isn't ideal, but a threshold relative to O-S could work. 
        x_osnd[[c(i-1, 4)]] <- inflexion_points_new[old_o + 5]
      }
      
      if((x_osnd[[c(i-1, 4)]] - x_osnd[[c(i-1, 3)]]) < 3 &  x_osnd[[c(i-1, 4)]] < 30){
        osnd[[c(i-1, 3)]] <- inflexion_points_new_yval[old_o + 4]                             # This is a special case for when there are two notches apparent
        x_osnd[[c(i-1, 3)]] <- inflexion_points_new[old_o + 4]                                # between systolic and diastolic waves
        osnd[[c(i-1, 4)]] <- inflexion_points_new_yval[old_o + 5]                             
        x_osnd[[c(i-1, 4)]] <- inflexion_points_new[old_o + 5]
      }
      
      
      
      if(osnd[[c(i-1, 3)]] < 0 & (osnd[[c(i-1, 4)]] - osnd[[c(i-1, 3)]]) > 1){        # If N is less than 0 and D-N >1, this implies that N has risen above D such 
        # find points either side of the notch                                        # that neither are inflection points and hence the inflection points
        d <- wuvn$notch_x[i-1] + 2                                                                # of next O and S are incorrectly assigned to N and D. 
        n <- wuvn$notch_y - 2                                                                # These lines correct for that by finding the inflection point instead, 
        d_yval <- predict(p_w[[i-1]], d)                                        # and taking values either side of it that approximate N and D. 
        n_yval <- predict(p_w[[i-1]], n) 
        # assign those points
        osnd[[c(i-1, 3)]] <- n_yval
        x_osnd[[c(i-1, 3)]] <- n
        osnd[[c(i-1, 4)]] <- d_yval
        x_osnd[[c(i-1, 4)]] <- d
      }
      
      if((x_osnd[[c(i-1, 3)]] - x_osnd[[c(i-1, 2)]]) < 5){     # despite this being on the canonical waveform side, these lines
        osnd[[c(i-1, 3)]] <- osnd[[c(i-1, 4)]]                 # are here because when an inflection point exists between systolic
        x_osnd[[c(i-1, 3)]] <- x_osnd[[c(i-1, 4)]]             # peaks on non-canonical waveforms, they fulfill the condition
        osnd[[c(i-1, 4)]] <- wuvn$notch_y[i-1]                  # to be processed on the canonical side. 
        x_osnd[[c(i-1, 4)]] <- wuvn$notch_x[i-1]
      }
    }    
  }else{
    
    #define old_o
    old_o <- max(which(inflexion_points_new < wuvn$w[i-1]))
    #define osnd
    x_osnd[[i-1]] <- inflexion_points_new[old_o:(old_o+3)]
    osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
    # Find new S
    s[i-1] <- wuvn$w[i-1] + (2*(wuvn$v_x[i-1] - wuvn$w[i-1]))
    s_yval[i-1] <- predict(p_w[[i-1]], s[i-1])
    
    # Distinguish between canonical vs non-canonical
    if((inflexion_points_new[old_o+1] - wuvn$w[i-1]) > (s[i-1] - wuvn$w[i-1])){
      
      # Remove false N and D values 
      x_osnd_precursor <- x_osnd[[i-1]]
      if(length(which(complete.cases(x_osnd_precursor) ==1)) != 4){
        x_osnd_precursor <- x_osnd_precursor[-(which(is.na(x_osnd_precursor)))]
      }
      if(length(x_osnd_precursor) == 2){
        x_osnd[[i-1]] <- x_osnd_precursor
      }else{
        false_points <- which(x_osnd_precursor > wuvn$notch_x[i-1]) 
        x_osnd_precursor <- x_osnd_precursor[-false_points]
        x_osnd[[i-1]] <- x_osnd_precursor
      }
      # Find y coords of OSND
      osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
      osnd[[c(i-1, 1)]] <- o_yval
      
      # Remove false values again
      osnd_precursor <- osnd[[i-1]]
      if(length(which(complete.cases(osnd_precursor) ==1)) != 4){
        osnd_precursor <- osnd_precursor[-(which(is.na(osnd_precursor)))]
      }
      if(length(osnd_precursor) == 2){
        osnd[[i-1]] <- osnd_precursor
      }else{
        osnd_precursor <- osnd_precursor[-false_points]
        osnd[[i-1]] <- osnd_precursor
      }
      
      osnd[[c(i-1, 3)]] <- osnd[[c(i-1, 2)]]
      osnd[[c(i-1, 2)]] <- s_yval[i-1]
      x_osnd[[c(i-1, 3)]] <- x_osnd[[c(i-1, 2)]]
      x_osnd[[c(i-1, 2)]] <- s[i-1]
      osnd[[c(i-1, 4)]] <- wuvn$notch_y[i-1]
      x_osnd[[c(i-1, 4)]] <- wuvn$notch_x[i-1]
      
    }else{
      x_osnd[[c(i-1, 1)]] <- o                                 # Cases where there are fewer than 4 inflexion points identified 
      x_osnd[[c(i-1, 2)]] <- inflexion_points_new[2]           # can occur if N and D converge and the spline isn't long enough for
      x_osnd[[c(i-1, 3)]] <- wuvn$notch_x[i-1]                             # the next S to register as an inflexion point.
      x_osnd[[c(i-1, 4)]] <- wuvn$notch_x[i-1]   
      osnd[[c(i-1, 1)]] <- o_yval
      osnd[[c(i-1, 2)]] <- inflexion_points_new_yval[2]
      osnd[[c(i-1, 3)]] <- wuvn$notch_y[i-1]   
      osnd[[c(i-1, 4)]] <- wuvn$notch_y[i-1]   
      
    }   
  }


  
  
  }
  
  return(c(osnd, x_osnd))
}