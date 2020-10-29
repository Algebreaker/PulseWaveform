# Find average:

# p = pulse
# ao = after_o

find_average <- function(p, ao){
  
  # Find the last 10 values of each wave:
  last10 <- list()
  for(i in 2:(ncol(p))){
    last10[[i-1]] <- p[, i][!is.na(p[, i])][(length(p[, i][!is.na(p[, i])])-10):(length(p[, i][!is.na(p[, i])]))]
  }
  
  # Find the average last 10 values
  mean_last10 <- c()
  for(i in 1:10){   
    row_vector <- c()
    for(j in 1:length(last10)){      
      row_vector[j] <- last10[[c(j, i)]] 
    }
    mean_last10[i] <- mean(row_vector[!is.na(row_vector)])   
  }
  
  # Find the consecutive y-axis differences (gradient essentially) of last 10 values:
  mean_diff_last10 <- c()
  for(i in 1:9){                    # 9 here since 10 values will give 9 values for differences between them
    mean_diff_last10[i] <- mean_last10[i+1] - mean_last10[i]
  }
  
  # Simply calculating the average wave by averaging each row of the dataframe doesn't work at the end of the wave, 
  # since as waves end, they no longer factor in to the calculation of the average. The average would be erroneously 
  # drawn out until the end of the last wave. Therefore, a clone dataframe of p is used to continue the trajectories
  # of waves after they have ended (by using the average gradient above) and thus make the end of the average wave
  # a more accurate approximation of the true average. 
  
  # Replace NA values with continuing downward gradient in a clone dataframe:
  p_for_finding_average <- p
  for(i in 2:(ncol(p_for_finding_average))){
    for(j in 1:length(p_for_finding_average[, i][ao[[(i-1)]][-1]])){
      p_for_finding_average[, i][ao[[(i-1)]][-1]][j] <- p_for_finding_average[, i][ao[[(i-1)]][1]] + j*mean(mean_diff_last10[1:5])
    }
  }
  
  # Calculate the average wave by row:
  average_wave <- c()
  sd_wave <- c()
  median_wave <- c()
  for(i in 1:nrow(p_for_finding_average)){
    row_vector <- c()
    for(j in 2:(ncol(p_for_finding_average))){
      row_vector[j-1] <- p_for_finding_average[i, j] 
    }
    average_wave[i] <- mean(row_vector[!is.na(row_vector)])   
    sd_wave[i] <- sd(row_vector[!is.na(row_vector)]) 
    median_wave[i] <- median(row_vector[!is.na(row_vector)])
  }
  
  # Find where the end of the average wave should end: 
  # This is done by finding the average (mode) y-value of the final value of each wave (before trajectory continuation)
  end <- c()
  for(i in 2:ncol(p)){
    end[i-1] <- p_for_finding_average[, i][ao[[(i-1)]][1]]
  }
  # define mode function
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  # Made this into a for loop because there can be waves that have no values after o, and if all waves are like that, end will be null
  if(length(end) > 1){     
    average_end <- mode(round(end[!is.na(end)], digits = 2))
    # Remove elements of average wave after where its end should be:
    below_o <- which(average_wave < average_end)
    # Find the first point in the last quarter of the wave that goes below baseline
    b. <- below_o[min(which(below_o > (length(average_wave)/(4/3))))]
    # If the wave doesn't go below baseline, don't remove any 
    if(is.na(b.) == FALSE){
      # otherwise, remove those elements from the average:
      average_wave[b.:length(average_wave)] <- NA
    }
  }
  
  # Find the median of the first x-values - make that the start of average:
  first_element <- c()
  for(i in 2:ncol(p)){
    stpqw <- p[, i][1:200]
    stq <- stpqw[is.na(stpqw)]
    first_element[i-1] <- length(stq) + 1
  }
  start <- median(first_element)-2
  average_wave[1:start] <- NA
  
  
  return(as.vector(average_wave))
}
