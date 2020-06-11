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
  
