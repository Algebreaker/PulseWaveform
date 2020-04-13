setwd("/home/johanna/Documents/Ninja theory/PulseAnalysis/Data/Craig")
data <- read.table("Source.csv", header=T, sep=",") #first line of the csv file needs to be deleted
library()

min_beats_per_minute <- 30
max_beats_per_minute <- 240
min_interval <- 0.25
max_interval <- 4
min_amplitude <- 0.5
min_frac_amplitude <- 0.2
sample_rate <- 75
peak_interval <- 9
peak_width <- 2.7   # Should peak width be an array of differing widths per peak?
skipped_beat_fraction <- 1.75

#main function
FNTBioRadioPulse<-function(input)
{
}

#attaching log events
AddLogEvent<-function(constFString, LogMessage)
{
}

#Debug draw 
Draw<-function()
{

}

#Downsample
# The BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second, 
# so we typically have repeated values in a pattern of 3-3-4 repeating. DownSample tries to retrieve the 
# unique values, with an effort to be robust against variation in the repeat pattern and also against 
# genuine repeated values.

DownSample <- function(input)   #Clarify TraceCount and SampleModulo
{
  
}

#Undetrend
# Analysis of device output indicates that the PPG signal is detrended by application of the following
# formula: OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1]), where the constance 0.96875 is
# an approximation fitted to the data.
# Individual pulse events are more comprehensible if the detrending is not used, so this function 
# removes it by inverting the above function. 

Undetrend <- function(Input, PrevInput, PrevOutput)
{
  
}

#VerifySignal
# Detect when the input signal indicates device removal.  If the un-detrended data is a constant,
# the PPG device has probably been removed.

VerifySignal <- function()
{
  
}


###more functions...



list<-rle(data$PPG.PulseOx1)
realrepeats<-rep(list$lengths > 4,times = list$lengths)

cbind(data,realrepeats)
