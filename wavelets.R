# Wavelet Transformation of the time series for pre-processing / visualiation / feature extraction:

library(WaveletComp)

# Complete step 1 of main script before running this:

# Create the wavelet object (currently arguments are tuned for bioradio sampling rate):
my.w <- analyze.wavelet(undetrended_data, "undetrended",
                        loess.span = 0,
                        dt = 1/75, dj = 1/1000,
                        lowerPeriod = 0.03,
                        upperPeriod = round(length(undetrended_data$undetrended)/75),
                        make.pval = FALSE, n.sim = 10)

# Gerenerate time/frequency spectrogram (also tuned for bioradio with regard to axis labeling):
wt.image(my.w, periodlab = "period length (seconds)", timelab = "time (seconds)", main = "test",         
         legend.params = list(lab = "wavelet power levels", mar = 4.7, label.digits = 2),  
         label.time.axis = TRUE, spec.time.axis = list(at = seq(0, length(undetrended_data$undetrended), by = 75), 
              labels = seq(0, length(seq(0, length(undetrended_data$undetrended), by = 75))-1, by = 1)))  

# Reconstruct the time series from all constituent wavelets:
reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2),
            legend.coords = "bottomleft") 

# Reconstruct time series based on just one frequency (/period) (or a range of frequencies):
# In this case have simply set lower limit (similar to low-pass filter)
x. <- reconstruct(my.w, sel.lower = 0.16, plot.waves = FALSE, lwd = c(1,2))

#Alternative arguments for reconstruct():
# only.ridge = TRUE   - only reconstruct from ridge frequencies
# sel.period =        - specify as a single period value, or a concatenation of values 

# Save as a new variable
post_wavelet_time_series <- x.$series$undetrended.r
plot(post_wavelet_time_series, type = "l")

## Make the reconstructed wave into a form that can be fed back into step 2 of the main script (undetrended_data$undetrended)
undetrended_data <- data_frame(1:length(post_wavelet_time_series))   
nrow(undetrended_data)
undetrended_data <- cbind(undetrended_data, post_wavelet_time_series)
colnames(undetrended_data)[2] <- "undetrended"

