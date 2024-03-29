---
title: "ReadMe"
output: md_document
---

```{r setup, include=FALSE}
```

Please see https://github.com/stw32/PulseWaveform for an up to date version of this repository. 

## PPG Model 2

This is a ReadMe for the photopletysmography (PPG) analysis published in []. It is supposed to both provide an overview of the workings of the code, any issues we encountered and some further explanations that did not make it into the paper. The basic assumption is that the signal is decaying towards a baseline and that the three relevant peaks can be modelled more consistently by discounting the decay element first. General_Purpose_script can be used an adapted for independent datasets. The ISO_study scripts were used for our analyses. 


# Prerequisites
The script makes use of the following packages:  
- library(splines2)  
- library(pracma)  
- library(SplinesUtils)    
- library(spectral)
- library(zoo)

# ISOFitting
This is the main script from which the other functions are called. First the working directory needs to be changed so the other functions can be called and the sampling rate need to be adjusted, i.e. `samplingRate <- 40`. The `beats_in` and `batch_number` parameters determine how many beats the parameters are estimated over and how many batches are generated for analysis and optimisation. The code uses a chi-square goodness of fit test to improve model fit and will do so by minimising chi-square over all the beats within a batch. The default values are 10 for both.
If you would like to select a specific subsection of your dataset rather than analyse all beats change the parameter `all_beats <- TRUE` to `FALSE`.
Then the path of the to-be-read-in files needs to specified. 


In the first instance we are trying to estimate the various detrending algorithms that are often applied to the raw data by the hardware: As a first step in the preprocessing pipeline the factor value is adjusted to reverse-engineer an assumed positive gradient at the tail-end of each beat. The values are changed until they reach a plausible-looking threshold for each individual beat: 
![](factorvalue.png)


The next step involves calculating and applying an offset such that the main trendline of the data is approximately horizontal:
![](offset1.png)
Then the main fun starts!

We interpolate a cubic spline of the provided data points 
`sfunction <- splinefun(1:length(undetrended), undetrended, method = "natural")`
and take its first derivative
`deriv1 <- sfunction(seq(1, length(undetrended)), deriv = 1)`. 
Then the maximum of the first derivative, called the W point(s) are found:  
`w <- find_w(d1p = deriv1Poly, deriv1 = deriv1, sp = splinePoly, sr = samplingRate)`
![](Wpoint.png)  
Beat segmentation is based on the correct identification of the W points.  
Next U, V and O are detected.
Based on the O points we remove a baseline in order to filter out low frequency modulation of the signal without altering the morphology of individual waveforms. 
![](baseline.png)  
 
The `sep_beats` function then takes care of some rudimentary cleaning procedures - plotting options of the outliers can be added by setting `q = TRUE`:   
- long waves (where a distance of more than two waves (Os) was counted as one) are removed  provided the waves to either side of them are not similarly long (if so there is a plotting option to make a manual decision)
- abnormally high amplitude
- waves that have a second systolic upstroke are removed  
- waves that fall below the O-O threshold (currently hardcoded to 4?) are excluded  
- waves that deviate more than 2 SDs from the average in their residuals
- waves that deviate more than 5 standard deviatiosn from the average wave (calculated after all of the above are already excluded) are removed  
- when `subset == T` is only necessary for the identification of a subset of data with an autonomic arousal manipulation (based on inter.beat intervals)  

In the next step, in the `FindStartParams` function each beat is modelled as a composition of three peaks (S, D and N) with parameters determining the amplitude, timing and width of each within the dataframe `beat`: `STime`, `SAmplitude` `SWidth` and so forth. Initial parameters are estimated using the excess of each beat which is what is left once the assumed decay towards a variable baseline is factored out:     
`excess[1] = data[1,2] - (baseline + config.rate*(yPrev-baseline))`

The width of each peak is initialised at 0.25 i.e. `par[4] <- 0.25`. 


## Improving fit using downhill simplex

To estimate the timing of the first and second reflectance peaks, we take the median of the time values generated by the initial parameter estimation.

`  renal_param <- median(beat$NTime)`

`  dias_param <- median(beat$DTime)`

Then we begin to refine the parameters originally generated by `FindStartParams`. We extract both within and across beat parameters. 
The across beat parameters extracted are those that are expected to be held constant across beats. These are:

- The width of S, R1 and R2
- The timing for R1 and R2
- The config rate

These parameters are fed into the `simplex.MakeSimplex2` function. This function slightly improves each of the across beat parameters, and outputs a suitable matrix to feed into the downhill simplex algorithm. For further explanation about the simplex creation and algorithm, please refer to Numerical Recipes (Press, Teukolsky, Vetterling & Flannery, 2007).

The within beat parameters produced are those which are able to vary from beat to beat. These are:

- The timing of S
- The amplitude of S, R1 and R2
- The two baseline values 

These parameters are fed into the `FindWithinParams` function, which subsequently feeds into a simplex. This improves each of the within beat parameters and outputs a suitable matrix for the downhill simplex algorithm. The distinction between this and the across beat parameters is that the within beat parameter values are calculated for each individual beat, whereas the across beat parameters are calculated for each batch.

The within and across beat parameters are combined into a matrix, and then fed into the downhill simplex algorithm, which iteratively improves the parameter estimation by minimizing the chi-sq value. The improved parameter values are extracted using `extractOutput`. Finally, `fixOutput` is utilised to ensure that there are no excessive deviations of the new parameters from our original estimations.

This process is then iterated upon four times, to ensure a the simplex does not get stuck in local minima. 


## Penalties applied to the chi-sq value
Penalties are applied to the chi-sq value if any of the below conditions are met: 

- If any peak amplitudes are estimated as negative 
- If the peak width is estimated as <0.05s or >0.5s (units?) 
- If the width of the R1 peak is <0.1 or >0.25
- If the diastolic width is >0.45
- If the renal peak timing strays too far from initial parameter estimation
- If the S-peak deviates too far from the maxima of the beat
- If the timing of the diastolic is < 0.2
- If there is less than 0,1s between the peaks
- If there is a large shift between baselines (Currently not applied!)
- If the configuration rate is above 0.95
- The renal peak is gradually penalised as the amplitude increase

# Contributing 
We would love to hear from people who would like to contribute or have ideas for developing out model further. 
