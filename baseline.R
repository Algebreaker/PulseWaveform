baseline <- function(plot = FALSE){
  # Making a (non-polynomial) spline to fit the baseline
  sfunction2 <- splinefun(inflexion_points[o], inflexion_points_yval[o], method = "natural")
  spline_base <- sfunction2(seq(1, length(undetrended_data$undetrended)), deriv = 0)

  # Plotting spline_base on spline_poly
  if(plot){
    plot(spline_poly)
    points(inflexion_points[o], inflexion_points_yval[o], pch = 19)
    lines(spline_base)
  }

  # Correcting for baseline:
  baseline_corrected <- undetrended_data$undetrended - spline_base
  if(plot){
    plot(baseline_corrected, type = "l")
  # Plot new baseline (y = 0)
    lines(1:length(undetrended_data$undetrended), seq(from = 0, to = 0, length.out = length(undetrended_data$undetrended)))
  }
  
  return(baseline_corrected)
}
