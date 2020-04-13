setwd("/home/johanna/Documents/Ninja theory/PulseAnalysis/Data/Craig")
data <- read.table("Source.csv", header=T, sep=",") #first line of the csv file needs to be deleted
library()

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

###more functions...


#The BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second, so we typically have repeated values in a pattern of 3-3-4 repeating.  DownSample tries to retrieve the unique values, with an effort to be robust against variation in the repeat pattern and also against genuine repeated values.
list<-rle(data$PPG.PulseOx1)
realrepeats<-rep(list$lengths > 4,times = list$lengths)

cbind(data,realrepeats)