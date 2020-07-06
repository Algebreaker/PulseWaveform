##Takes in the current estimated values for the S, N and D curves and refits them using the trace derivatives 
#p=pulse2
#ss = s_sines
#ds = d_sines
#ns = n_sines
#period = average period

refit_peaks <- function(p, ss, ds, ns, period){
  
  max_shift=4
  sample_rate = 75
  s_x <- c()
  s_y <- c()
  d_x <- c()
  d_y <- c()
  n_x <- c()
  n_y  <- c()
  
  for(i in 2:(ncol(p))){
    
    const = c()
 
    fitted_s <- FALSE
    fitted_n <- FALSE
    fitted_d <- FALSE
    
    while(!fitted_s | !fitted_n | !fitted_d){
      
      ##define the peaks of the current waveform
      s_peak <- max(ss[[i-1]])
      n_peak <- max(ns[[i-1]])
      d_peak <- max(ds[[i-1]])
      d_peak_x <- match(d_peak, ds[[i-1]])
      n_peak_x <- match(n_peak, ns[[i-1]])
      s_peak_x <- match(s_peak, ss[[i-1]])
      
      #if we've refitted any of the peaks already, define them as -1
      if(fitted_s){
        s_peak = -1
      }
      if(fitted_n){
        n_peak = -1
      }
      if(fitted_d){
        d_peak = -1
      }
      
      #Create a constant with the final value of the S-sine (should be the same as the D-sine) - this needs to be added in every time
      #we subtract a sine from the PPG trace
      const <- rep(ss[[i-1]][[length(ss[[i-1]])]],length(ss[[i-1]]))
      
      if(s_peak > n_peak){
        #Subtract the D and N peaks from the original trace, adding back in the constant each time to create the residual
        d_resid <- p[,i] - ds[[i-1]]
        d_resid <- d_resid + const
        d_n_resid <- d_resid - ns[[i-1]]
        d_n_resid <- d_n_resid + const
        
        #Create a spline of the residual, and solve for the inflexion points
        d_n_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(p[, i]), d_n_resid, "natural")
        inflex_test_x <- solve(d_n_resid_spline, b=0, deriv=1)
        inflex_test_y <- predict(d_n_resid_spline, inflex_test_x)
        
        #find which peak is closest to the current x-value of the S peak
        pk_loc <- which(abs(inflex_test_x-s_peak_x)==min(abs(inflex_test_x - s_peak_x)))
        
        s_x[i-1] <- inflex_test_x[pk_loc]
        s_y[i-1] <- inflex_test_y[pk_loc]
        
        
        fitted_s <- TRUE
        
        
      } else if(n_peak > d_peak){
        
        #Subtract the S and D peaks from the original trace, adding back in the constant each time to create the residual 
        s_resid <- p[,i] - ss[[i-1]]
        s_resid <- s_resid + const
        s_d_resid <- s_resid - ds[[i-1]]
        s_d_resid <- s_d_resid + const
        
        #Create a spline, and estimate the y-value on the first derivative at the current x coordinate of the N-peak
        s_d_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(p[, i]), s_d_resid, "natural")
        current_deriv_n <- predict(s_d_resid_spline, n_peak_x, deriv=1)
        
        #Get the inflexion points on the first derivative
        inflex_test_x <- solve(s_d_resid_spline, b=0, deriv=1)
        inflex_test_y <- predict(s_d_resid_spline, inflex_test_x)
        
        #If the current_deriv_n is < 0, we look between the S and N peaks for the refitted N value
        if(current_deriv_n < 0){
          subset_x <- inflex_test_x[inflex_test_x > s_peak_x & inflex_test_x < n_peak_x]
          x_loc <- match(subset_x, inflex_test_x)
          subset_y <- inflex_test_y[x_loc]
          pk_loc <- which(abs(subset_x - n_peak_x)==min(abs(subset_x - n_peak_x)))
         #Otherwise we look between the N and D peaks for the refitted N value
        } else if(current_deriv_n > 0){
          subset_x <- inflex_test_x[inflex_test_x > n_peak_x & inflex_test_x < d_peak_x]
          x_loc <- match(subset_x, inflex_test_x)
          subset_y <- inflex_test_y[x_loc]
          pk_loc <- which(abs(subset_x - n_peak_x)==min(abs(subset_x - n_peak_x)))
          
        } 
        
        n_y[i-1] <- subset_y[pk_loc]
        n_x[i-1] <- subset_x[pk_loc]
        
        if(fitted_s & (!fitted_d)){
          fitted_s = FALSE;
        }
        
        fitted_n = TRUE;
        
        
      } else{
        #Subtract the S and N peaks from the original trace, adding back in the constant each time to create the residual
        s_resid <- p[,i] - ss[[i-1]]
        s_resid <- s_resid + const
        s_n_resid <- s_resid - ns[[i-1]]
        s_n_resid <- s_n_resid + const
        
        #Create a spline, and estimate the y-value on the first derivative at the current x coordinate of the D-peak
        s_n_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(p[, i]), s_n_resid, "natural")
        current_deriv_d <- predict(s_n_resid_spline, d_peak_x, deriv=1)
        
        #Get the inflexion points on the first derivative
        inflex_test_x <- solve(s_n_resid_spline, b=0, deriv=1)
        inflex_test_y <- predict(s_n_resid_spline, inflex_test_x)
       
        #If the current_deriv_d is < 0, we look between the N and D peaks for the refitted N value
        if(current_deriv_d < 0){
          subset_x <- inflex_test_x[inflex_test_x > n_peak_x & inflex_test_x < d_peak_x]
          x_loc <- match(subset_x, inflex_test_x)
          subset_y <- inflex_test_y[x_loc]
          pk_loc <- which(abs(subset_x - d_peak_x)==min(abs(subset_x - d_peak_x)))
        #Otherwise, look between the D peak and the end of the trace
        } else if(current_deriv_d > 0){
          subset_x <- inflex_test_x[inflex_test_x > d_peak]
          x_loc <- match(subset_x, inflex_test_x)
          subset_y <- inflex_test_y[x_loc]
          pk_loc <- which(abs(subset_x - d_peak_x)==min(abs(subset_x - d_peak_x)))
          
        }
        
        d_x_est <- subset_x[pk_loc]
        d_y_est <- subset_y[pk_loc]
        
        # if(d_x_est > (d_peak_x + max_shift)){
        #d_y_est <- inflex_test_y[which.max(inflex_test_y)]
        #d_x_est <- inflex_test_x[which.max(inflex_test_y)]
        #d_y_est = 2 * d_y_est * sample_rate * avg_period
        #d_x_est = d_x_est + 0.5 * pi * avg_period
        #  }
        
        d_y[i-1] <- d_y_est
        d_x[i-1] <- d_x_est
        
        
        fitted_d = TRUE
        
      }
    }
  }
  tot <- data.frame(s_x, s_y, n_x, n_y, d_x, d_y)
  return(tot)

}
