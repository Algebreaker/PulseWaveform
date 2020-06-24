spectrum <- function(baseline_corrected){
  #filtered <- filter.fft(spline_poly,fc=0.0925,BW=0.0525,n=50)
  #plot.fft(filtered)
  #spectral<-spec.fft(spline_poly)
  #plot(data_downsampled$PPG.PulseOx1)
  #plot(baseline_corrected)
  
  powerspectrum<-spectrum(data_downsampled$PPG.PulseOx1)
  powerspectrum<-spectrum(baseline_corrected)
  
  LF<-ffilter(data_downsampled$PPG.PulseOx1, f=75, from = 0.04, to = 0.145, bandpass = TRUE) #low bandpass, sampling frequency of 75 Hz (75 times per second)
  powerspectrum<-spectrum(LF)
  HF<-ffilter(data_downsampled$PPG.PulseOx1, f=75, from = 0.145, to = 0.45, bandpass = TRUE) #high bandpass
  powerspectrum<-spectrum(HF)
  
  spectralratio<-sum(LF)/sum(HF) #ratio of LF/HF
  
  Y1 <- fft(baseline_corrected[1:1000])
  
  plot(abs(Y1), type="h")
  
  library(DescTools)
  bands<-c(0.04, 0.145, 0.45)
  frebands<-Freq(baseline_corrected, breaks = bands)
  Y1 <- fft(frebands[1])
  
  library('signal')
  low<-0.04
  high<-0.145
  bf <- butter(2, c(low, high), type = "pass")
  signal.filtered <- filtfilt(bf, baseline_corrected)
  fourier <- fft(signal.filtered)
  plot(abs(fourier), type="h")
  
  low<-0.145
  high<-0.45
  bf <- butter(2, c(low, high), type = "pass")
  signal.filtered <- filtfilt(bf, baseline_corrected)
  fourier2 <- fft(signal.filtered)
  plot(abs(fourier2), type="h")
  
  spectralratio<-sum(fourier)/sum(fourier2) #ratio of LF/HF
  spectralratio
  return(spectralratio)
  
}