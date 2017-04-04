

JTPcalibrateToGlucose <- function(realSpectra, ppm){
JTP=list()

# Locate the target region (n ppm either side of 5.233)
regionMask = (ppm < 5.733) & (ppm > 4.98)
# nosaltres tenim la part positiva de l'espectre a l'esquerra
maskOffset = which(ppm > 5.733)
maskOffset = maskOffset[length(maskOffset)]

# Take the approximate second derivative
realSpectra2 = diff(realSpectra[regionMask], differences = 2)

# Find the two lowest points, corresponding to the top of the two sharpest
# peaks.
min1Index = which.min(realSpectra2) ##ok
min1Index = maskOffset + min1Index
while (realSpectra[min1Index+1]>realSpectra[min1Index]) {
  min1Index = min1Index + 1
}
while (realSpectra[min1Index-1]>realSpectra[min1Index]) {
  min1Index = min1Index - 1
}

peak1PPM = ppm[min1Index];

# Having found the first peak, flatten it so we can locate the second.
peak1Mask = (ppm[regionMask] > (peak1PPM - 0.004)) & (ppm[regionMask] < (peak1PPM + 0.004));
realSpectra2[peak1Mask] = 0

min2Index = which.min(realSpectra2); ##ok
min2Index = maskOffset + min2Index
while (realSpectra[min2Index+1]>realSpectra[min2Index]) {
  min2Index = min2Index + 1
}
peak2PPM = ppm[min2Index]

# Reference to the midpoint of peak 1 and 2.

JTP$deltaPPM = mean(c(peak1PPM, peak2PPM)) - 5.233
JTP$ppm = ppm - JTP$deltaPPM
return(JTP)
}




JTPcalibrateToTSP <- function(realSpectra, ppm){
JTP=list()
# Locate the target region (n ppm either side of 5.233)
regionMask = (ppm < 0.1) & (ppm > -0.1)
# nosaltres tenim la part positiva de l'espectre a l'esquerra
maskOffset = which(ppm > 0.1)
maskOffset = maskOffset[length(maskOffset)]

# Take the approximate second derivative
TSP_max_position = which.max(realSpectra[regionMask])
deltaPPM = round(TSP_max_position + maskOffset)
JTP$deltaPPM = ppm[deltaPPM]

JTP$ppm = ppm - JTP$deltaPPM;

return(JTP)
}



JTPcalibrateToPeak2 <- function(realSpectra, ppm, PeakPos, winsize){
JTP=list()
# Locate the target region (n ppm either side of PeakPos)
regionMask = (ppm < PeakPos+winsize) & (ppm > PeakPos-winsize)
# nosaltres tenim la part positiva de l'espectre a l'esquerra
maskOffset = which(ppm > PeakPos+winsize)
maskOffset = maskOffset[length(maskOffset)]

Peak_max_position = which.max(realSpectra[regionMask])
deltaPPM = round(Peak_max_position + maskOffset)
JTP$deltaPPM = ppm[deltaPPM] - PeakPos

JTP$ppm = ppm - JTP$deltaPPM

return(JTP)
}




JTPcalibrateToPeak <- function(realSpectra, ppm, PeakPos){
JTP=list()
# Locate the target region (n ppm either side of PeakPos)
regionMask = (ppm < PeakPos+0.1) & (ppm > PeakPos-0.1)
#regionMask = (ppm < PeakPos+0.05) & (ppm > PeakPos-0.05);
# nosaltres tenim la part positiva de l'espectre a l'esquerra
maskOffset = which(ppm > PeakPos+0.1)
maskOffset = maskOffset[length(maskOffset)]

# Take the approximate second derivative
realSpectra2 = diff(realSpectra[regionMask], differences = 2)

# Find the lowest point, corresponding to the top of the peak.
min1Index = which.min(realSpectra2) ##ok
min1Index = maskOffset + min1Index
while (realSpectra[min1Index+1]>realSpectra[min1Index]) {
  min1Index = min1Index + 1
}
while (realSpectra[min1Index+1]>realSpectra[min1Index]) {
  min1Index = min1Index - 1
}

JTP$deltaPPM = ppm[min1Index] - PeakPos

JTP$ppm = ppm - JTP$deltaPPM

return(JTP)
}