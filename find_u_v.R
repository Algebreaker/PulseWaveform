find_u_v <- function(dat, wx, wy, d1, d1p, plot = FALSE){
# Find half the height of w (on derivative y-axis)
  w_half_height <- predict(d1p, wx)/2
  # Find u and v:
  half_heights <- c()
  half_heights_yval <- c()
  for(i in 1:length(w_half_height)){ 
    deriv1_poly_peak_subset <- CubicInterpSplineAsPiecePoly((wx[i]-10):(wx[i]+10), d1[(wx[i]-10):(wx[i]+10)], "natural") 
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
