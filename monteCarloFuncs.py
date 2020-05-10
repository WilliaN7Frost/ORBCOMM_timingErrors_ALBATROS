# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 22:49:31 2020

@author: wilia
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import corner




def findPeak(timelag , cc , returnIndex=True):
    
    maxPeak = 0.
    maxPeak_distToZeroTimeLag = timelag[-1]
    index=0
    for i in range(len(timelag)):
        if (np.abs(cc[i]) > np.abs(maxPeak)): maxPeak = cc[i];  maxPeak_distToZeroTimeLag = timelag[i];  index=i
    if returnIndex: return index
    else: return maxPeak_distToZeroTimeLag , maxPeak
    
    
   
    
    
def getZeroCrossingCoords( peakIndex , timelag , cc , returnIndex=True , doPlot=False):
    
    notPastZero_left=True;    notPastZero_right=True
    i_left = peakIndex-1;    i_right = peakIndex+1
    if cc[peakIndex] > 0:
        peakSign = 1
    else:
        peakSign = -1
    
    while (  notPastZero_left and notPastZero_right ):
        
        try:
            if (peakSign*cc[i_left] > 0 and notPastZero_left):
                i_left -= 1
            else:
                notPastZero_left = False
            if (peakSign*cc[i_right] > 0 and notPastZero_right):
                i_right += 1
            else:
                notPastZero_right = False
        
        except:
            print("getZeroCrossingsAround() went out of bounds of your input data array")
            return [i_left , i_right]
        
    near0_x = [ timelag[i_left],timelag[i_left+1],timelag[i_right-1],timelag[i_right]  ]
    near0_y = [ cc[i_left],cc[i_left+1],cc[i_right-1],cc[i_right]  ]
    
    if doPlot:
        plt.figure()
        plt.plot(timelag,cc , c='b')
        plt.scatter( timelag[peakIndex] , cc[peakIndex] , c='r')
        plt.scatter(  near0_x , near0_y , c='g')
        
    if returnIndex: return [i_left , i_left+1] , [i_right-1 , i_right]
    else: return near0_x , near0_y






def runZeroCrossingMC( near0index , timelag , cc , noise , numTrials=1000 , quantiles=[0.16,0.84] , doPrint=True , doPlot=False):
    
    zeroPointScatter = np.empty(numTrials)
    xVals = [ timelag[near0index[0]] , timelag[near0index[1]] ]
    cc_len = len(cc)
    
    for i in range(numTrials):
        
        noise = np.roll( noise , np.random.randint(0,cc_len))
        yVals = [ cc[near0index[0]]+noise[near0index[0]] , cc[near0index[1]]+noise[near0index[1]] ]
        
        slope = (yVals[1]-yVals[0])/(xVals[1]-xVals[0])
        intercept = yVals[0] - xVals[0]*slope
        
        zeroPointScatter[i] = -intercept/slope
    
    median = np.median(zeroPointScatter)
    quantiles = [ np.quantile(zeroPointScatter,quantiles[0]) , np.quantile(zeroPointScatter,quantiles[1]) ]
    
    if doPlot:
        plt.figure()
        plt.plot(timelag,cc , c='b')
        plt.scatter(  median , 0 , c='g' , label='median zero-crossing value')
        plt.scatter( xVals , [ cc[near0index[0]] , cc[near0index[1]] ] , c='r', label='Zero-crossing border points')
        plt.legend()
        
        plt.figure()
        plt.hist(zeroPointScatter , bins=25)
        plt.axvline(median , c='r')
        plt.axvline(quantiles[0] , c='black')
        plt.axvline(quantiles[1] , c='black')
    if doPrint:
        print("\nMedian For Linear Fit of 0-crossing = "+str(median))
        print("1-sigma Quandtiles = "+str( (quantiles[0]-median,quantiles[1]-median) ))
        print("Mean of 0-crossings = " + str(np.mean(zeroPointScatter)))
        print("STD of 0-crossings = " + str(np.std(zeroPointScatter)))
    
    return zeroPointScatter , median , quantiles







def quadratic(x,a,h,k):
    return a*(x-h)**2 + k

def runParabolaFitMC( startIndex , stopIndex , timelag , cc , noise , numTrials=1000 , 
                                                   quantiles=[0.16,0.84] , doPrint=True , doPlot=False ):
    
    peakScatter = np.empty((numTrials,3))
    cc_len = len(cc)
    
    for i in range(numTrials):
        
        noise = np.roll( noise , np.random.randint(0,cc_len))
        xVals = timelag[ startIndex : stopIndex+1 ]
        yVals = cc[ startIndex : stopIndex+1 ] + noise[ startIndex : stopIndex+1 ]
        guessPeakY = np.max(yVals)
        indexGuess = np.where(yVals==guessPeakY)[0][0]
        guessPeakX = xVals[indexGuess]
        
        params , pcov = curve_fit( quadratic , xVals , yVals , p0=[1. , guessPeakX , guessPeakY] , maxfev=5000)
        peakScatter[i] = params
    
        
    if doPlot:
        corner.corner(peakScatter , labels=["a","peak_X","peak_Y"] , quantiles=[0.16,0.5,0.84] , show_titles=True)
    
    meds = np.median(peakScatter , axis=0)
    quantiles = [ np.quantile(peakScatter , quantiles[0] , axis=0) , np.quantile(peakScatter , quantiles[1] , axis=0) ]
    
    if doPlot:
        plt.figure()
        plt.plot(timelag,cc , c='b')
        plt.scatter(  meds[1] , meds[2] , c='g')
        
        inds = np.random.randint(numTrials, size=100)
        x = np.linspace(xVals[0] , xVals[-1] , 100 , endpoint=True)
        for ind in inds:
            sample = peakScatter[ind]
            plt.plot(x, sample[0] * (x-sample[1])**2 + sample[2], alpha=0.05, color='red')
            
        plt.title("Monte Carlo Fit of Peak Using Quadratic Least Squares")
        plt.xlabel("Time Lag (ns)")
        plt.ylabel("Cross-Correlation")
    
    if doPrint:
        peakScatter = np.transpose(peakScatter)
        print("\nMedians for peak XY coord = " + str( (meds[1],meds[2]) ))
        print("1-sigma Quantiles" + str( (quantiles[0],quantiles[1]) ))
        print("Mean xPeak = "+ str(np.mean(peakScatter[1])) + ". Mean yPeak = " + str(np.mean(peakScatter[2])))
        print("STD xPeak = "+ str(np.std(peakScatter[1])) + ". STD yPeak = " + str(np.std(peakScatter[2])))
        peakScatter = np.transpose(peakScatter)
    
    return peakScatter , meds , quantiles