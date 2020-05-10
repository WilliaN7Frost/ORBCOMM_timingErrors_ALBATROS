# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 02:48:03 2020

@author: wilia


    Once cross-correlations and their respective noises are obtained for a given range of observation times, a Monte Carlo approach can be used to 
estimate the timing error for each observation time of a certain satellite signal. 
    For simplicity, this will be defined as the standard deviation of many trials where the noise function is randomly shifted and added in to its 
respective cross-correlation, and then a least-squares quadratic fit to the cross-correlation peak is performed. Again, for simplicity, this peak 
is defined to be the point with maximal distance to the time axis (where cross-cor = 0).
   
In addition, data regarding the location of the inferred xcorr peak as well as the SNR is recorded as a function of observation time.

"""

import numpy as np
import matplotlib.pyplot as plt
from os import path
import monteCarloFuncs as mc



"""
A short setup phase before performing some Monte Carlo's on the data
"""



    # sat_A is was the satellite we were interested in cross-correlating
            # The numbers here are based on the prior labelling of satellite signals. For example, to cross-correlate the 2nd channel interval
            # found in 'freq_ranges_good', we would designate sat_A = 2 and satTypeA='Good'
    # sat_B was the satellite we used to determine the (over-estimated) noise by aligning the different polarizations of the A and B signals 
    # in Fourier space and cross-correlating them
sat_A = 2;    satTypeA='Good'
sat_B = 3;    satTypeB='Good'

    # What were the min and max observation times?
minObsTime = 5.;    maxObsTime = 125.



    # ZERO-PADDING
    # This controls how many additional zeros we had added to the end of a Fourier Transform to increase interpolation 
    # between points when inverting back to time-stream space. It is based of a multiple of the original FT length. 
    # Therefore, a 'zpadCoeff' of 1 adds ~100% more points (as zeros) to the end of the FT
zpadCoeff = 1
if zpadCoeff==0: zpad_str = '_noZPad'
else:            zpad_str = '_zPadCoeff'+str(int(zpadCoeff))+'p'+str(int((zpadCoeff % int(zpadCoeff))*10))



    # The length in the x-axis of our cross-cor. It increases as we zero-pad in the cross-correlation process, but the
    # actual time intervals represented remain with the same bounds
chunk_len = 1000
cc_len = int( chunk_len*(zpadCoeff+1) ) + int(2*zpadCoeff)



    # The relative path from the current directory to where the raw data had been stored
relPath = "dataFromEarlyFeb"
    # Name of the directory where we want to funnel the cross-cor and noise data generated
    # This directory will be created in the same directory the raw data is in
#output_directory = "xcorr_"+satTypeA+"Sat"+str(sat_A)+"_wNoiseUsing_"+satTypeB+"Sat"+str(sat_B)
output_directory = "xcorr_"+satTypeA+"Sat"+str(sat_A)+"_withNoiseUsing_"+satTypeB+"Sat"+str(sat_B)+"_obsTime"+str(int(minObsTime))+"to"+str(int(maxObsTime))+"secs"
#output_directory = ""
output_directory = path.join(relPath, output_directory)

    

    # Boolean value to check if we want to save the timing error data generated for future use
save=False
    # Boolean value to check if  we want to plot various things as we do the MC
doPlot=False
doPlotInMC=False
    # If we only want to plot the cross-cors and noise obtained, set doMC = False and doPlot=True
doMC = True


    # These parameters control plot label sizes
plt.rc('font', size=14)
plt.rc('axes', titlesize=30);     plt.rc('axes', labelsize=25)    
plt.rc('xtick', labelsize=18);    plt.rc('ytick', labelsize=18);    plt.rc('legend', fontsize=18)   






"""
The setup phase is now complete. Let's get crackin'
"""



obsTimes = np.load( output_directory+r"\obsTimes"+".npy" )
numberOfObs = len(obsTimes)

SNR = np.empty(numberOfObs)
mean_Xpeak = np.empty(numberOfObs)
timing_errors = np.empty(numberOfObs)

for i in range(numberOfObs):
    
    
    obsTime_str = str(int(obsTimes[i]))+'p'+str(int((obsTimes[i] % int(obsTimes[i]))*1000))

        # Getting our time-lag, cross-cor and noise values back
    time_lag = np.load( output_directory+r"\ccLen"+str(cc_len)+zpad_str+".npy" )
    crossCor = np.load( output_directory+r"\ccLen"+str(cc_len)+"_crossCor_"+satTypeA+"Sat"+str(sat_A)+"_obsTime"+obsTime_str+"secs"+zpad_str+".npy" )
    if satTypeA==satTypeB:
        noise = np.load(output_directory+r"\ccLen"+str(cc_len)+"_noiseCC_"+satTypeA+"Sats"+str(sat_A)+"and"+str(sat_B)+"_obsTime"+obsTime_str+"secs"+zpad_str+".npy" )
    else:
        noise = np.load(output_directory+r"\ccLen"+str(cc_len)+"_noiseCC_"+satTypeA+'Sat'+str(sat_A)+satTypeB+'Sat'+str(sat_B)+"_obsTime"+obsTime_str+"secs"+zpad_str+".npy" )
   
     
    SNR[i] = round(  np.max(np.abs(crossCor))/np.std(noise) , 2  )
    print("SNR at "+str(obsTimes[i])+" secs of obsTime = "+str(SNR[i]))
    
    
    if doPlot:
        
        plt.figure()
        plt.plot(time_lag , crossCor)
        plt.grid()
        plt.title("Cross-Correlation Of An ORBCOMM signal")
        plt.xlabel("Time Lag (ns)")
        plt.ylabel("Cross-Correlation")
        plt.xlim(-20000,20000)
        
        plt.figure()
        plt.plot(time_lag , noise)
        plt.grid()
        plt.title("Cross-Correlation Of 2 ORBCOMM signals (noise)")
        plt.xlabel("Time Lag (ns)")
        plt.ylabel("Cross-Correlation")
        plt.xlim(-20000,20000)
        
    
        # NOTE: The following commented code was originally used for another timing error estimation technique, where instead of fitting a quadratic to a peak,
        #       a linear fit was performed on the 2 points to the left and right of a zero-crossing close to a max peak. This technique is not used anymore,
        #       but its corpse remains for education purposes
    if doMC:
        
            # Finding the location of the peak in the cross-correlation
        peakIndex = mc.findPeak( time_lag , crossCor )
        #near0x , near0y = mc.getZeroCrossingCoords( peakIndex , time_lag , crossCor , returnIndex=False)
        
            # returns the indices to the left and right of our peak
        near0leftOfPeak , near0rightOfPeak = mc.getZeroCrossingCoords( peakIndex , time_lag , crossCor , returnIndex=True , doPlot=doPlotInMC)
        
        #zeroScat , median , quantiles = f.runZeroCrossingMC( near0leftOfPeak , time_lag , crossCor , noise , numTrials=10000)
        #print("Difference in quantiles = "+str((quantiles[1]-quantiles[0])))
        #print("Mean = "+str(np.mean(zeroScat)))
        #print("STD = "+str(np.std(zeroScat)))
        
        
            # Runs the quadratic least-squares Monte Carlo about the peak of our cross-correlation
        peakScatter , meds , quantiles = mc.runParabolaFitMC( near0leftOfPeak[0] , near0rightOfPeak[1] , time_lag , crossCor , noise ,
                                                             numTrials=15000 , doPlot=doPlotInMC , doPrint=False)
        peakScatter = np.transpose(peakScatter)
        mean = np.mean(peakScatter[1])
        std = np.std(peakScatter[1])
        print("Mean xPeak = "+ str(mean) +".   STD xPeak = "+ str(std))
        mean_Xpeak[i] = mean
        timing_errors[i] = std
    


if save and doMC:
    np.save(output_directory + r"\obsTimes_timeErr_xPeak_SNR" + "_crossCor_"+satTypeA+"Sat"+str(sat_A)+zpad_str+".npy" , [obsTimes , timing_errors , mean_Xpeak , SNR] )







