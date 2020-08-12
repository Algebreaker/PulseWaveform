setwd("/home/johanna/Documents/Ninja theory/Khalsa/ISO_3.0/ISO_3.0/Iso_3_PhysioData/AB765_ISOonly/scan_20161111/physiological_files")
setwd("/home/johanna/Documents/Ninja theory/Other PPG/jarchin_casson2017_wrist-ppg-during-exercise-1.0.0/wrist-ppg-during-exercise-1.0.0")
library('e1071')
library('swdft')
library('signal')

samplfreq=40

#useful functions:
plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

# Plot the i-th harmonic
# Xk: the frequencies computed by the FFt
#  i: which harmonic
# ts: the sampling time points
# acq.freq: the acquisition rate
plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
  Xk.h <- rep(0,length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
  points(ts, harmonic.trajectory, type="l", col=color)
}


#########
data <- read.table("ECG_10_ISO_R1.1D", header=T)
library("R.matlab")
data <- read.table("s1_high_resistance_bike.atr", sep=",")
data<-data$V1
data<-data[1,]
plot(data[1:90], type="h")

plot(data[,1:90])
plot(data[,1], type="h")

#pick PPG column
ppgdata<-data$X.256
plot(ppgdata, type='l')
#just some plotting to see how long a beat is
shortdata<-data$X.256[1:200] #50 data points is roughly one beat
plot(shortdata, type='l')
#######Splitting Data into Windows#######
# Size of sliding window
win<- 200 #40 to make each window 1sec long
# Size of overlap
o <- 3/4*win

# Define start and end point (s and e)
s <- 1
e <- win

# Loop to create fragments
for(i in 1:(length(ppgdata)/o)){
  
  assign(paste0("x", i), ppgdata[s:e])
  s <- s + o
  e <- (s + win) - 1
  
}

# Call fragments  
x1
x2
x3
#etc.


#####


#out<-stft(ppgdata, win=50, inc=1/4*win, coef=40, wtype="hamming.window") #coefficient numbers need to be the same as window size, for 1s it's 40
#vals<-out$values
#multvec<-c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#filtered<-sweep(vals, MARGIN=2, multvec, `*`)
#inverse<-  fft(vals, inverse=TRUE) / length(vals)
#library(seewave)
#filtered<-as.complex(filtered)
#inv<-ifft(filtered)
#i<-fft(filtered, inverse=TRUE) / length(filtered)
#plot(abs(i[1:100]), type='l')
#plot(Re(i[1:40]), type='l', col="#ff0000")
#par(new=TRUE)
#plot(shortdata, type="l")



#vals[vals > 3000] <- 0 #set vaues to 0

#plot(Mod(vals), log='y', main='FFT of ecg vs index')
#inverse<-  fft(vals, inverse=TRUE) / length(vals)
#plot(inverse, type='l')

###
#filtlength<-win-1
#filtvector<-c(0.02, 0.2)
#filtoutput<-fir1(n=filtlength, w=filtvector, type = c("pass"),window = hamming(filtlength + 1), scale = TRUE) #  Nyquist frequency must be 20 Hz since sampling frequency is 40Hz, so 0.4 and 4 scaled to a range fo 0 to 1 must be 0.02 and 0.2 Hz respectively
#plot(filtoutput)

########Apply Hamming Window#######
hamlength<-win
plot(hamming.window(hamlength))
ham<-hamming.window(hamlength)
w<-ham*x1
plot(w, type='l', col='green')
#y<-stft(x1, wtype="hamming.window")
X.k <- fft(w)
#coefficients are 1/T where T=win/40Hz
T=win/samplfreq
#calculate Hz for each coefficients to design filter between 0.4 and 4Hz
coef<-vector("numeric", length(X.k)/2)
for (j in 1:(length(X.k)/2)){
  coef<-j/T
}
  

X.k[6:36]<-0 #set all coefficients apart from first 4 to zero
X.k[6:36]<-0 #set all coefficients apart from first 4 to zero
inverse<-  fft(X.k, inverse=TRUE) / length(X.k)

#c <- 2*Re(inverse)/win
#c[1] <- Re(inverse[1])/win
#s <- -2*Im(inverse)/win

plot(Re(inverse), type='l', col="red")

#c <- 2*Re(X.k)/win
#c[1] <- Re(X.k[1])/win
#s <- -2*Im(X.k)/win


lines(1:win,c[1]+c[2]*cos(2*pi*0:39/win)+s[2]*sin(2*pi*0:39/win),col="#ff0000")
lines(1:win,c[1]+c[2]*cos(2*pi*0:39/win)+s[2]*sin(2*pi*0:39/win)+c[3]*cos(4*pi*0:39/win)+s[3]*sin(4*pi*0:39/win)+c[4]*cos(6*pi*0:39/win)+s[4]*sin(6*pi*0:39/win)+c[5]*cos(8*pi*0:39/win)+s[5]*sin(8*pi*0:39/win),col="#ff8000",lty=2)

y <- 1.5+sin(2*pi*(1:N/N))+cos(4*pi*(1:N/N))

# Generate basis functions
dt=0:1/60:3;
df=[3:3:12];
basis1=exp(1j*2*pi*df(1)*dt);
basis2=exp(1j*2*pi*df(2)*dt);
basis3=exp(1j*2*pi*df(3)*dt);
basis4=exp(1j*2*pi*df(4)*dt);

% Reconstruct var
var_recon=basis1*f_useful(1)+...
basis2*f_useful(2)+...
basis3*f_useful(3)+...
basis4*f_useful(4);
var_recon=real(var_recon);


plot(Re(inverse), type='l', col="red")
par(new=TRUE)
plot(x1, type='l')

plot(y$values[,2], type='l', col="red")

#bartlett window
plot(bartlett(hamlength))
bart<-bartlett(win)
t<-bart*x1
plot(t, type='l', col='green')
y<-stft(x1, wtype="bartlett.window")

plot(x1, type='l')
par(new=TRUE)
plot(t, type='l', col="red")
######FFT each window#######

library(stats)
X.k <- fft(x1)
plot(abs(X.k)[1:16], type="h")
frequencies<-abs(X.k)
frequencies2<-Re(X.k)
plot(frequencies2)
plot(Arg(X.k), type="h") #Plots the phase of each value of X.k

acq.freq <- 40 #sampling rate was 40Hz
plot.frequency.spectrum(X.k,xlimits=c(0,acq.freq/2))
plot(X.k, type='l')

plot(Mod(X.k), type='l', log='y', main='FFT of ecg vs index')

plot(Re(fft(x1))^2)

#####bandpass filter#########




#########IFFT#########
inverse<-  fft(X.k, inverse=TRUE) / length(X.k)
plot(inverse, type='l')




#trying stuff
src<-x1
fft<-fft(x1)
fft[2:16]<-0
inv<-fft((fft), inverse=TRUE)
plot(Re(inv), type='l')
par(new=TRUE)
plot(src, type='l')


ex<-
  fft(fft(src), inverse=TRUE)
output<-Re(fft(fft(src), inverse=TRUE))
plot(output, type='l')

hamming.window(X.k)

y<-stft(X.k, wtype="hamming.window")
plot(y$values)
