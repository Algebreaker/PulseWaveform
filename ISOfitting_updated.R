# Model2 (Updated for across time series fitting) for ISO data:

# Cleaned version:

setwd("Desktop/Scripts")

library(tidyverse)                                 
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils) 
library(spectral)
library(DescTools)

source("~/Desktop/WorkshopCS/Workshop2/model2.R")
source("~/Desktop/WorkshopCS/Workshop2/simplex.R")
source("~/Desktop/Scripts/PPG_funcs.R")
source("~/Desktop/Scripts/iso_funcs.R")

samplingRate <- 40
beats_in <- 10  
batch_number <- 10

# Load time series:
ppg <- read.csv("~/Desktop/Khalsa_Fletcher_collaboration copy/Physio Dial/ISO_3.0/Iso_3_PhysioData/AS495/scan_20170130/physiological_files/ECG_7_ISO_R1.1D", sep = "")   
ppg <- data.frame(
  time = (0:(nrow(ppg)-1)) / 40,
  ppg = ppg[,1]
)
names(ppg)[1] <- "time (s)"
names(ppg)[2] <- "Detrended"

# Find beats:
n <- dim(ppg)[1]
vpg <- ppg[2:n,2] - ppg[1:(n-1),2]
beat <- ppg[which(vpg[1:(n-1)] < 300 & vpg[2:n] >= 300),1]
nBeat <- length(beat)
beat <- data.frame(
  beat = beat,
  dt = (1:nBeat)*0.0
)
rm(vpg)

# Check time series:
plot(ppg$`time (s)`[1:100], ppg$Detrended[1:100], type = "l")

# Adjust factor:
factor_value <- FactorAdjust(ppg, beat, fs = model2.FindSegment, gs = model2.GetSegment, u = UnDetrend, factorCutoff = -40, plot = T)
ppg3 <- data.frame(ppg[,1],UnDetrend(ppg,factor=factor_value,offset=1))

# Adjust offset:
offset_value <- OffsetAdjust(ppg3, ppg, u = UnDetrend, factor_value, plot = F)
ppg[,2] = UnDetrend(ppg,factor=factor_value,offset=offset_value)  

# Baseline correct:
undetrended <- ppg$Detrended
sfunction <- splinefun(1:length(undetrended), undetrended[1:length(undetrended)], method = "natural")
deriv1 <- sfunction(seq(1, length(undetrended)), deriv = 1)
spline1 <-  sfunction(seq(1, length(undetrended)), deriv = 0)
splinePoly <- CubicInterpSplineAsPiecePoly(1:length(undetrended), undetrended[1:length(undetrended)], "natural")
deriv1Poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended), deriv1, "natural") 
inflexX <- solve(splinePoly, b = 0, deriv = 1)
inflexY <- predict(splinePoly, inflexX)
w <- find_w(d1p = deriv1Poly, deriv1 = deriv1, sp = splinePoly, sr = samplingRate)
uv <- find_u_v(dat = undetrended, wx = w$wX, wy = w$wY, d1 = deriv1, d1p = deriv1Poly, spline = splinePoly, spline_o = spline1, plot=FALSE)
o <- find_o(wx = w$wX, inx = inflexX, iny = inflexY, d1p = deriv1Poly, sp = splinePoly)
tmp <- preclean_wuv(w=w, uv=uv, o=o, samp = samplingRate, sp = spline1, q = F)
w <- tmp[[1]]
uv <- tmp[[2]]
rm(tmp)
baseCor <- baseline(inx = inflexX, iny = inflexY, o = o, dat = undetrended, sp = splinePoly, plot=F)
ppg[, 2] <- baseCor


# Check time series:
plot(ppg[1:1000,1],ppg[1:1000,2],t='l')

# Add output columns to ppg and beat:
ppg$Baseline = 1:nrow(ppg) * 0
ppg$Excess   = 1:nrow(ppg) * 0
ppg$Residue  = 1:nrow(ppg) * 0
beat <- AddOutput(beat)

# Fill beat and ppg with parameters derived from the excess:
temp <- FindStartParams(batch_number, beats_in, beat, ppg, fs = model2.FindSegment, gs = model2.GetSegment, e = model2.Excess, sep = model2.SubtractExcessPeak)
beat <- temp[[1]]
ppg <- temp[[2]]

# Make and run simplex for each batch:
beat_orig <- beat
fit_check <- list()
for(k in 1:batch_number){         
  
  beat <- beat_orig[((k*beats_in)-(beats_in-1)):(k*beats_in), ]
  renal_param <- median(beat$NTime)
  dias_param <- median(beat$DTime)
  par <- as.numeric(beat[1,5:16])   
  beat_start <- beat[, 3]
  beat_end <- beat[, 4]
  beat_vector <- list(beats_in, beat_start, beat_end)
  
  # Refine parameters:
  for(i in 1:4){
    if(i == 1){new_beat <- beat}
    within_params <- FindWithinParams(beats_in, ppg, beat = new_beat, gs = model2.GetSegment, fp = model2.FixParams3, ms = simplex.MakeSimplex3, m2 = model2.ChiSq)
    across_params <- simplex.MakeSimplex2(data=ppg, param = par, f = model2.ChiSq3, inScale = 0.1, beat_vector = beat_vector, beat = new_beat, renal_param = renal_param, dias_param = dias_param)
    mat <- make_matrix(across_params, within_params)
    sim <- simplex.Run2(data = ppg, simplexParam = mat, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, renal_param = renal_param, dias_param = dias_param, run = c("run", i))
    output <- extractOutput(beats_in, sim)
    fixed <- FixOutput(beats_in, beat = new_beat, ppg, gs = model2.GetSegment, fp = model2.FixParams3, across = output[[1]], within = output[[2]])
    new_beat <- UpdateBeat(beats_in, beat, fixed)
  }
  
  # Assess fit:
  fit_check[[k]] <- model2.ChiSq3(data = ppg, params = NULL, beats = beat_vector, beat = new_beat, a = sim[1, ], plot = FALSE, renal_param = renal_param, dias_param = dias_param)
  
  # Finalise:
  beat2 <- new_beat       
  colnames(beat2) <- colnames(beat)
  beat2 <- beat2[, -c(1:4)]
  
  # Plot
  PlotFits(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2)
}

