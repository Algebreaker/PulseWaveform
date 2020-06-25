
preproc <- function(dat){
  dat<-dat[!(dat$PPG.PulseOx1=='NaN'),]
  
  #Downsample
  # The BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second, 
  # so we typically have repeated values in a pattern of 3-3-4 repeating. DownSample tries to retrieve the 
  # unique values, with an effort to be robust against variation in the repeat pattern and also against 
  # genuine repeated values.
  
  list<-rle(dat$PPG.PulseOx1)
  ID <- rep(1:length(list$values), times = list$lengths)
  data2 <- cbind(data, ID)
  data_downsampled <-c()
  
  for (i in 1:max(ID)){
    sub.data <- filter(data2, ID == i)
    if(nrow(sub.data) <= 4){
      data_downsampled <- rbind(data_downsampled, sub.data[1,])
    }else if(nrow(sub.data) > 4 ){data_downsampled <- rbind(data_downsampled, sub.data[1,], sub.data[5,])}
  }
  
  
  #Undetrend
  # Analysis of device output indicates that the PPG signal is detrended by application of the following
  # formula: OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1]), where the constance 0.96875 is
  # an approximation fitted to the data.
  # Individual pulse events are more comprehensible if the detrending is not used, so this function 
  #removes it by inverting the above function. 
  undetrended <-replicate(length(data_downsampled$PPG.PulseOx1)-1, 0) 
  undetrended<-c(data_downsampled$PPG.PulseOx1[1],undetrended) #add first detrended value to vector
  for (i in 2:length(data_downsampled$PPG.PulseOx1)){
    undetrended[i]<-((data_downsampled$PPG.PulseOx1[i]-80) - ((data_downsampled$PPG.PulseOx1[i-1]-80) * 0.96875) + (undetrended[i-1]))
  }
  undetrended_data<-cbind(data_downsampled,undetrended)
  
  return(undetrended_data)
}


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


