#p = pulse2
#col_len = source_data_column_length
#p_w = poly_wave

find_wuv <- function(p, col_len, p_w){ 
  u_x <- c()
  v_x <- c()
  u_y <- c()
  v_y <- c()
  w <- c()
  w_y <- c()
  notch_x <- c()
  notch_y <- c()
  
  for(i in 2:(ncol(p))){        
    
    sfunction <- splinefun(1:(col_len+4), p[, i], method = "natural")
    deriv1_wave <- sfunction(seq(1, (col_len+4)), deriv = 1)
    deriv1_wave_poly <- CubicInterpSplineAsPiecePoly(1:(col_len+4), deriv1_wave, "natural") 
    
    # Find inflexion points on deriv1_wave_poly
    inflexion_points_deriv1_wave_poly <- solve(deriv1_wave_poly, b = 0, deriv = 1)
    inflexion_points_deriv1_wave_poly_yval <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)

    # Find correct threshold using histogram
    # hdat_2<-hist(deriv1_wave)
    quantiles_2 <- quantile(deriv1_wave, probs=c(.025,.95))
    threshold_2 <- quantiles_2[2]
    
    # Identifying peaks of deriv1_poly:
    w_poly_peaks_wave <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)
    w_poly_peaks_wave <- which(w_poly_peaks_wave > threshold_2)     
    w_poly_peaks_wave <- inflexion_points_deriv1_wave_poly[w_poly_peaks_wave]
    
    
    # Finding the 'notch' (renal wave)    
    notch <- inflexion_points_deriv1_wave_poly[(which(inflexion_points_deriv1_wave_poly_yval == min(inflexion_points_deriv1_wave_poly_yval)))+1]
    notch_poly_yval <- predict(p_w[[i-1]], notch)
    
    
    #Find U and V
    
    # Find half the height of w (on derivative y-axis)
    w_half_height_wave <- predict(deriv1_wave_poly, w_poly_peaks_wave[1])/2
    # Find u and v for derivative:
    half_heights_wave_new <- solve(deriv1_wave_poly, b = w_half_height_wave[1])
    half_heights_wave_new_yval <- predict(deriv1_wave_poly, half_heights_wave_new)
    u <- half_heights_wave_new[1]
    v <- half_heights_wave_new[2]
    # Find u and v y-values for original wave:
    u_v_yval_wave <- predict(p_w[[i-1]], half_heights_wave_new)
    u_yval <- u_v_yval_wave[1]
    v_yval <- u_v_yval_wave[2]
    
    u_y[i-1] <- u_yval
    v_y[i-1] <- v_yval
    u_x[i-1] <- u
    v_x[i-1] <- v
    w[i-1] <- w_poly_peaks_wave[1]
    notch_x[i-1] <- notch
    notch_y[i-1] <- notch_poly_yval
    

  }
  wuv <- data.frame(w, u_x, u_y, v_x, v_y, notch_x, notch_y)
  return(wuv)
}
    