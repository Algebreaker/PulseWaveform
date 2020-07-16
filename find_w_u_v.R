find_w <- function(dat, d1, d1p, sp, plot = FALSE){
  
  #Plot the y values:
  if(plot){
    plot(spline_poly)
    points(inflexion_points, inflexion_points_yval, pch = 19)
  }
  ## Finding W  (note that the quantile threshold 0.95 needs adjusting for some datasets)
  # Find inflexion points on deriv1_poly
  inflexion_points_deriv1 <- solve(d1p, b = 0, deriv = 1)
  inflexion_points_deriv1_yval <- predict(d1p, inflexion_points_deriv1)
  if(plot){
    plot(d1p)
    points(inflexion_points_deriv1, inflexion_points_deriv1_yval, pch = 19)
  }
  # Find correct threshold using histogram
  # hdat<-hist(deriv1)
  quantiles<-quantile(d1,probs=c(.025,.95))      ## this still needs adjusting based on source data              
  threshold<-quantiles[2]
  # Identifying peaks of deriv1_poly:
  w_poly_peaks <- predict(d1p, inflexion_points_deriv1)
  w_poly_peaks <- which(w_poly_peaks > threshold)     
  w_poly_peaks <- inflexion_points_deriv1[w_poly_peaks]
  w_poly_deriv_yval <- predict(d1p, w_poly_peaks)
  w_poly_peaks_yval <- predict(sp, w_poly_peaks)
  
  
  # Plot on deriv1_poly
  if(plot){
    plot(d1p)
    points(w_poly_peaks, w_poly_deriv_yval, pch = 19)
    # Plot back on spline_poly:
    plot(sp)
    points(w_poly_peaks, w_poly_peaks_yval, pch = 19)
  }
  
  df <- data.frame(w_poly_peaks, w_poly_peaks_yval, w_poly_deriv_yval)
  return(df)
}


find_u_v <- function(dat, wx, wy, d1, d1p, plot = FALSE){
  # Find half the height of w (on derivative y-axis)
  w_half_height <- predict(d1p, wx)/2
  # Find u and v:
  half_heights <- c()
  half_heights_yval <- c()
  for(i in 1:length(w_half_height)){ 
    deriv1_poly_peak_subset <- CubicInterpSplineAsPiecePoly((wx[i]-9):(wx[i]+9), d1[(wx[i]-9):(wx[i]+9)], "natural") 
    half_heights_precursor <- solve(deriv1_poly_peak_subset, b = w_half_height[i])
    half_heights[c((2*(i)-1), (2*(i)))] <- half_heights_precursor
    half_heights_yval[c((2*(i)-1), (2*(i)))] <- predict(deriv1_poly_peak_subset, half_heights[c((2*(i)-1), (2*(i)))])
  }
  
  # Plot u's and v's on deriv1_poly
  if(plot){
    plot(d1p)
    points(half_heights, half_heights_yval, pch = 19)
  }
  # Find u and v 
  u <- half_heights[seq_along(half_heights) %%2 != 0] 
  v <- half_heights[seq_along(half_heights) %%2 == 0]   
  # Find u and v y-values for spline_poly
  u_v_yval <- c()
  for(i in 1:length(wx)){
    spline_poly_peak_subset <- CubicInterpSplineAsPiecePoly((wx[i]-10):(wx[i]+10), dat[(wx[i]-10):(wx[i]+10)], "natural") 
    u_v_yval[c((2*(i)-1), (2*(i)))] <- predict(spline_poly_peak_subset, half_heights[c((2*(i)-1), (2*(i)))])
  }
  u_yval <- u_v_yval[seq_along(u_v_yval) %%2 != 0] 
  
  v_yval <- u_v_yval[seq_along(u_v_yval) %%2 == 0]  
  # Plot u's and v's on spline_poly
  if(plot){
    plot(spline_poly)
    points(half_heights, u_v_yval, pch = 19)
  }  
  df <- data.frame(u, u_yval, v, v_yval)
  return(df)
}

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
    
    sfunction <- splinefun(p$x, p[, i], method = "natural")
    deriv1_wave <- sfunction(p$x, deriv = 1)
    deriv1_wave_poly <- CubicInterpSplineAsPiecePoly(p$x, deriv1_wave, "natural") 
    
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
    notch_range <- which(inflexion_points_deriv1_wave_poly < (0.25*max(inflexion_points_deriv1_wave_poly)))  ## Look for the notch in the first quarter of the range from W to the max inflexion point
    notch <- inflexion_points_deriv1_wave_poly[(which(inflexion_points_deriv1_wave_poly_yval[notch_range] == min(inflexion_points_deriv1_wave_poly_yval[notch_range])))+1]
    #Find corresponding y value on poly_wave
    notch_poly_yval <- predict(poly_wave[[i-1]], notch)
    
    
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

