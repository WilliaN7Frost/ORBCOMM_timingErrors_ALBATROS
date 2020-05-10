# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 22:28:34 2020

@author: wilia

A collection of functions used to:
    
    - Invert PFB data
    - Auto-correlate rechannelized PFB data
    - Perform cross-correlation on time streams

"""

import numpy as np
import matplotlib.pyplot as plt
import read_4bit
import pfb
import time as t
import gc


    # not used but could be useful if necessary
def db(x):
    """ Convert linear value to dB value """
    return 10*np.log10(x)


    # Used to remove sharp risig edges due to unwanted code-generated spikes
    # Can also act as a general smoothing function for any 2D plots
def getRidOfSpikes(porcupine , fracOfMax_dy , runs , doPlot=False , doPrint=False):
    """
    porcupine : 1D array
        Fourier Space data to be un-spiked, hence the name 'porcupine'
    fracOfMax_dy : float
        Fraction of the max change in y of any consecutive points. Used to determine if a spike should be removed or not
    runs : int
        Total number of getRidOfSpikes trials
    """
    if doPlot:
        plt.figure()
        plt.plot(porcupine , label='Initial FT')
        plt.xlabel("Channels")
        plt.ylabel("PSD")
        plt.title("View the evolution of getRidOfSpikes")
        plt.legend()
    
    for r in range(runs):
              
        max_dy = 0
        for i in range(1,len(porcupine)):
            dy = porcupine[i]-porcupine[i-1]
            if(dy > max_dy):
               max_dy = dy
        
        count = 0
        for i in range(2,len(porcupine)-2):
            prev2 = porcupine[i-2]
            prev1 = porcupine[i-1]
            next1 = porcupine[i+1]
            next2 = porcupine[i+2]
            if((next1-porcupine[i]) < -1.0*fracOfMax_dy*max_dy or (porcupine[i]-prev1) > fracOfMax_dy*max_dy):
                count+=1
                porcupine[i] = (prev2+next2)/2
                
        if doPrint: print("getRidOfSpikes Run #"+str(r+1)+". Max dy was " + str(max_dy) + " and # of spikes removed = " + str(count))
        if doPlot:
            if(r==runs-1):
                plt.plot(porcupine , label='Final FT')
            else:
                plt.plot(porcupine)
        
    return porcupine




    # Used to get an array of regions of interest corresponding to known satellite channel intervals
def getSatsFromFT(ft , satRegions , totRefChans):
    """
    ft : 1D or 2D array
        Fourier Space data for which we want to isolate certain regions (zero-pad the other channels)
    satRegions : 2D array
        Each element of this array contains lower and upper bounds that constrain channel regions of interest
    totRefChans : int
        Total number of channels in the fourier space used by the user to define the subset channel regions of interest 
    """
    sats = []
    ft_dims= np.ndim(ft)
    if    ft_dims==1 : ft_len = len(ft);  rangeAdjust = (ft_len)/totRefChans
    elif  ft_dims==2 : ft_len , inside_ft_len = ft.shape;  rangeAdjust = (inside_ft_len)/totRefChans
    else: print("Fourier Transform entered in 'getSatsFromFT' was not 1D or 2D"); return [ft]
    
    for i in range(len(satRegions)):
        ft_copy = ft.copy()
        l_bound = int(satRegions[i][0]*rangeAdjust)
        u_bound = int(satRegions[i][1]*rangeAdjust)
        for q in range(ft_len):
            if ft_dims==1:
                if(q > u_bound or q < l_bound):
                    ft_copy[q] = 0
            else:
                for u in range(inside_ft_len):
                    if(u > u_bound or u < l_bound):
                        ft_copy[q][u] = 0
        sats.append(ft_copy)
    return np.array(sats)








    # Used to align 2 signals that are not found in the same channels
    # Used in noise (over)estimation
def alignSats( ft1,ft2 , sat1 , sat2 , totRefChans , doPlot=False):
    """
    ft1, ft2 : 1D or 2D array
        Fourier Space data for which we want to shift one such as to align it with a channel region in the other
    sat1, sat2 : 1D array
        Each element of these arrays contains lower and upper bounds that constrain channel regions of interest to be aligned with each other
    totRefChans : int
        Total number of channels in the fourier space used by the user to define the subset channel regions of interest 
    """
    
    #print("Initially, sat regions are " + str( (sat1[0],sat1[1]) ) +" "+ str( (sat2[0],sat2[1]) ) )
    
    ft_dims= np.ndim(ft1)
    if    ft_dims==1 : ft_len = len(ft1);  rangeAdjust = (ft_len)/totRefChans
    elif  ft_dims==2 : ft_len , inside_ft_len = ft1.shape;  rangeAdjust = (inside_ft_len)/totRefChans
    else: print("Fourier Transform entered in 'alignSats' was not 1D or 2D"); return 
    
    if doPlot:
        plt.figure()
        plt.title('Before sat shift')
        if ft_dims==1:
            plt.plot( np.abs(ft1[::500])**2 , c='blue' )
            plt.plot( np.abs(ft2[::500])**2 , c='red'  )
        else:
            plt.plot( np.mean( np.abs(ft1)**2 , axis=0), c='blue' )
            plt.plot( np.mean( np.abs(ft2)**2 , axis=0), c='red' )
    
    
    sat1_width = sat1[1]-sat1[0]
    sat2_width = sat2[1]-sat2[0]
    
        # This conditional statement ensures that the center of the bigger channel interval is shifted towards the center of the smaller one.
        # It also sets the bigger channel interval equal to the smaller one, such that once the shift is performed, the isolated channel
        # intervals contain at most one complete satellite signal, and do not include parts of another.
    if ( sat1_width > sat2_width ):
        shiftBy = int( sat2[0]-sat1[0] + (sat2_width-sat1_width)/2 )
        sat1 = sat2.copy()
        if ft_dims==1:
            ft1 = np.roll( ft1 , int(shiftBy*rangeAdjust) )
        else:
            ft1 = np.roll( ft1 , int(shiftBy*rangeAdjust) , axis=1 )
    
    else:
        shiftBy = int( sat1[0]-sat2[0] + (sat1_width-sat2_width)/2 )
        sat2 = sat1.copy()
        if ft_dims==1:
            ft2 = np.roll( ft2 , int(shiftBy*rangeAdjust) )
        else:
            ft2 = np.roll( ft2 , int(shiftBy*rangeAdjust) , axis=1 )
          
    #print("After shifting, sat regions are " + str( (sat1[0],sat1[1]) ) +" "+ str( (sat2[0],sat2[1]) ) )   
    if doPlot:
        plt.figure()
        plt.title('After sat shift')
        if ft_dims==1:
            plt.plot( np.abs(ft1[::500])**2 , c='blue' )
            plt.plot( np.abs(ft2[::500])**2 , c='red'  )
        else:
            plt.plot( np.mean( np.abs(ft1)**2 , axis=0), c='blue' )
            plt.plot( np.mean( np.abs(ft2)**2 , axis=0), c='red' )
            
    return ft1 , ft2 , sat1 , sat2     









    # Used to produce a time-lag axis centered at x=0
def x_fftshift(x):
    x_shift = np.fft.fftshift(x)
    i = 0
    xLen = len(x)
    while(i < xLen/2):
        x_shift[i] -= xLen
        i+=1
    return x_shift

    
def doCrossCor(x,y):
        return np.fft.irfft( np.fft.rfft(x) * np.conj(np.fft.rfft(y)) )

    # Used to perform cross-cor with zero-padding in the FTs
def doCrossCor_zpad(x,y,zpadCoeff):
    xft = np.fft.rfft(x)
    yft = np.fft.rfft(y)
    #print("xft,yft have shape " + str( (xft.shape,yft.shape) ))
    if np.ndim(x)==1:
        xft = np.array( np.append( xft , np.zeros( int(len(xft)*zpadCoeff)) ) )
        yft = np.array( np.append( yft , np.zeros( int(len(yft)*zpadCoeff)) ) )
    elif np.ndim(x)==2:
        xft = np.array( np.append( xft , np.zeros( (len(xft) , int(len(xft[0])*zpadCoeff)) ) , axis=1) )
        yft = np.array( np.append( yft , np.zeros( (len(yft) , int(len(yft[0])*zpadCoeff)) ) , axis=1) )
    #print("xft,yft (zeropad) have new shape " + str( (xft.shape,yft.shape) ))
    return np.fft.irfft( xft * np.conj(yft) )















    # Used to perform cross correlation on satellite signals found in Fourier Space
    # This function is mostly used as a helper in the 'performInv_Rechan_AC_CC' function to perform cross-cor
    # on many satellites. For a specific satellite cross-correlation only,  the 'performCC_on2Sats' function is more appropriate
def performCC_analysisOnSats( ft1,ft2 , satRegions , totRefChans , dt , newSubLen=500 ,
                              dt_units='unspecified' , satDesignation='' , zpadCoeff=0 , useChunksOfRTS=True):
    """
    ft1 , ft2 : 1D or 2D array
            the fourier transforms containing intervals of interest (in this case satellite signals)
    satRegions : 2D array 
            Each element of this array contains lower and upper bounds that constrain channel regions of interest
    totRefChans , dt , newSubLen , zpadCoeff : int
            - 'totRefChans' is the total number of channels in the fourier space used to determine the 'satRegions'
            - 'dt' is the time increment in the recovered timestreams
            - 'newSubLen' is the desired length of our cross-cor data. Shouldn't be too big since we expect satellite cross-cors to be about 0
            - 'zpadCoeff' is a coefficient that, when multiplied with a FT length, corresponds to the amount of zeros appended at the end of that FT.
    """
    if np.ndim(satRegions)==1: satRegions = [satRegions];
    
    sats_ft1 = getSatsFromFT(ft1 , satRegions , totRefChans)
    sats_ft2 = getSatsFromFT(ft2 , satRegions , totRefChans)
    
    sats_rts1 = np.fft.irfft(sats_ft1[0]).ravel()
    lenCC = len(sats_rts1)
    
    sats_ft1_0pad = []
    sats_ft2_0pad = []
    for i in range(len(satRegions)):
        sats_ft1_0pad.append( np.append( sats_ft1[i], np.zeros(zpadCoeff*len(sats_ft1[i])) ) )
        sats_ft2_0pad.append( np.append( sats_ft2[i], np.zeros(zpadCoeff*len(sats_ft1[i])) ) )
    sats_ft1 = np.array(sats_ft1_0pad)
    sats_ft2 = np.array(sats_ft2_0pad)
    gc.collect()


    sats_CC = []
    dummyBool=True
    
    for i in range(len(satRegions)):
        
        sats_rts1 = np.fft.irfft(sats_ft1[i]).ravel()
        sats_rts2 = np.fft.irfft(sats_ft2[i]).ravel()
        lenCC_0pad = len(sats_rts1)
        if dummyBool:
            dt = dt * lenCC/lenCC_0pad
            print("new dt is = " + str(int(dt)))
            dummyBool=False
            
        # Reshaping the signal before passing it through the cross-cor. Allows to discard regions far from the peaks we want to look at
        newLength = int(lenCC_0pad/newSubLen)*newSubLen
        sats_rts1 = sats_rts1[:newLength].reshape(-1,newSubLen)
        sats_rts2 = sats_rts2[:newLength].reshape(-1,newSubLen)
        
        sats_CC.append( doCrossCor( sats_rts1 , sats_rts2 ) )
    
    sats_CC = np.sum( sats_CC , axis=1 )
    
    plt.figure()
    plt.title("Cross-Correlation of " +satDesignation+ " Satellite Sources")
    plt.xlabel("Time Lag ("+dt_units+")")
    plt.ylabel("Cross-Correlation")
    for i in range(len(satRegions)):
        plt.plot( x_fftshift( np.arange(newSubLen) )*dt , np.fft.fftshift( sats_CC[i] ) , label=("Sat %s" % str(i+1)) )
    plt.legend()
    plt.grid()
    






def performCC_on2Sats( ft1,ft2 , satRegion1,satRegion2 , totRefChans , dt , newSubLen=500 ,
                        dt_units='' , zpadCoeff=0 , plotTitle='' , doPlot=True , doPrint=False , useChunksOfRTS=True):
    """
    ft1 , ft2 : 1D or 2D array
            the fourier transforms containing intervals of interest (in this case satellite signals)
    satRegion1,2 : 2D arrays
            Different fourier space regions containing lower and upper bounds that constrain channel regions of interest in ft1, ft2 respectively
    totRefChans , dt , newSubLen , zpadCoeff : int
            - 'totRefChans' is the total number of channels in the fourier space used to determine the 'satRegions'
            - 'dt' is the time increment in the recovered timestreams
            - 'newSubLen' is the desired length of our cross-cor data. Shouldn't be too big since we expect satellite cross-cors to be about 0
            - 'zpadCoeff' is a coefficient that, when multiplied with a FT length, corresponds to the amount of zeros appended at the end of that FT.
    """
   
        # if both satellite regions are not the same, align them by shifting one towards the other in Fourier space
    if np.ndim(satRegion1)==1 and np.ndim(satRegion2)==1 and satRegion1[0]!=satRegion2[0] and satRegion1[1]!=satRegion2[1]:
        ft1 , ft2 , satRegion1 , satRegion2 = alignSats( ft1 , ft2 , satRegion1 , satRegion2 , totRefChans)
        # Safety check to make sure the correct data type is input to 'getSatsFromFT'
    if np.ndim(satRegion1)==1: satRegion1 = [satRegion1];
    if np.ndim(satRegion2)==1: satRegion2 = [satRegion2];
        # Getting sat regions in the FTs
    sat_ft1 = getSatsFromFT(ft1 , satRegion1 , totRefChans)
    sat_ft1 = sat_ft1[0]
    sat_ft2 = getSatsFromFT(ft2 , satRegion2 , totRefChans)
    sat_ft2 = sat_ft2[0]
        # Adjusting 'dt' based on the zero-padding done to the FTs
    dt = dt / (1+zpadCoeff)

        
    if useChunksOfRTS:
        sat_rts1 = np.fft.irfft(sat_ft1).ravel()
        sat_rts2 = np.fft.irfft(sat_ft2).ravel()
        lenCC = len(sat_rts1)
            # Reshaping the signal before passing it through the cross-cor. Allows to limit the amount of points plotted far from x=0
        newLength = int(lenCC/newSubLen)*newSubLen
        sat_rts1 = sat_rts1[:newLength].reshape(-1,newSubLen)
        sat_rts2 = sat_rts2[:newLength].reshape(-1,newSubLen)
        if doPrint: print("New rts shape is " + str(sat_rts1.shape))
        crossCor = np.fft.fftshift(  np.sum( doCrossCor_zpad(sat_rts1,sat_rts2 , zpadCoeff) , axis=0)  )
        time_lag = x_fftshift( np.arange(len(crossCor)) )*dt
    else:
        sat_ft1 = np.array( np.append( sat_ft1, np.zeros(int(len(sat_ft1)*zpadCoeff))) )
        sat_ft2 = np.array( np.append( sat_ft2, np.zeros(int(len(sat_ft2)*zpadCoeff))) )
        crossCor = np.fft.fftshift( np.fft.irfft( sat_ft1 * np.conj(sat_ft2) ) ) 
        time_lag = x_fftshift( np.arange(len(crossCor)) )*dt 
        lowerPlotBound = int(len(crossCor)/2)-int(newSubLen/2)
        upperPlotBound = int(len(crossCor)/2)+int(newSubLen/2)
    
    gc.collect()
    
    
    if doPlot:
        plt.figure()
        plt.title("Cross-Cor of " + plotTitle)
        plt.xlabel("Time Lag ("+dt_units+")")
        plt.ylabel("Cross-Correlation")
        plt.grid()
        if useChunksOfRTS: plt.plot( time_lag , crossCor , marker='.')
        else: plt.plot( time_lag[lowerPlotBound:upperPlotBound] , crossCor[lowerPlotBound:upperPlotBound] )
        #plt.xlim(-20000,20000)
    
    if useChunksOfRTS: return time_lag , crossCor
    else: return time_lag[lowerPlotBound:upperPlotBound] , crossCor[lowerPlotBound:upperPlotBound]









    # Function which allows to invert PFB, rechannelize PFB-inverse data, with options of plotting auto-cor and cross-cor
    # Its base functionality is returning the inverse of PFB data called the recovered time-stream (RTS), although function arguments 
    # can be changed to return or plot data in different forms
def performInv_Rechan_AC_CC(fileID , fileNum , reChan=1025 , originalNumChans=2048 , originalSampRate = 2*(125e6) , saveFigs = False ,
                    trim=True , showRTS=False , performAutoCor=False , performCrossCor=False , doAllSats=False , zpadCoeff=0 ,
                    newSubLen=500 , useChunksOfRTS=True , returnRTS=True , returnFFT=False , freq_ranges=[] , freq_range_types=[]):
    
    start = t.time()
    
        # Extracting PFB data from the files (both polarizations and the channels used)
    pfb_raw1 , pfb_raw2 = read_4bit.read_4bit_new(fileID)
    rows , nchan = pfb_raw1.shape
    print('\nRaw data dims in file #'+str(fileNum)+' = (' + str(rows) + ',' + str(nchan) + ')')
    
        # Calculate the effective time per sample in the recovered timestream
    dt = originalNumChans / originalSampRate / nchan
    
        # Performing the inverse PFB on both polarizations. Assume nTaps=4
    rts_pol1 = pfb.inverse_pfb(pfb_raw1 , 4)
    rts_pol2 = pfb.inverse_pfb(pfb_raw2 , 4)
    print("Inverted PFB shape in file #"+str(fileNum)+" = "+str(rts_pol1.shape))
    
        # Ravel the output of the inverse PFB (make it one big timestream) and remove edges due to high residuals
    rts_pol1 = rts_pol1.ravel()[1000:-1000]
    rts_pol2 = rts_pol2.ravel()[1000:-1000]
    lenRTS = len(rts_pol1)
    print("Recovered timestream length (removed edges) in file #"+str(fileNum)+" = " + str(lenRTS))
    #print("Effective 'dt' in recovered timestream is " + str(int(dt*10**9)) + " nanoseconds")
    print("Total time of recovered timestream (removed edges) in file #"+str(fileNum)+" = "+ str(round(lenRTS*dt, 3)) + "s")
    nanosecsPerSamp = dt * 10**9
    
    
    
    if (returnFFT or performCrossCor):
        # Do rechannelization using FFT
        fft1 = np.fft.rfft(rts_pol1[1000:-1000])
        fft2 = np.fft.rfft(rts_pol2[1000:-1000])
        print("Rechannelized FFT (max channels) shape = " + str(fft1.shape))
    if (performAutoCor):
        # Do rechannelization using PFB
        pfb_rechan1 = pfb.pfb(rts_pol1[1000:-1000] , reChan)
        pfb_rechan2 = pfb.pfb(rts_pol2[1000:-1000] , reChan)
        print("Rechannelized PFB ("+str(reChan)+" channels) shape = " + str(pfb_rechan1.shape))
    
        # These parameters control plot label size
    plt.rc('font', size=14)          # controls default text sizes
    plt.rc('axes', titlesize=30);    plt.rc('axes', labelsize=25)    
    plt.rc('xtick', labelsize=18);   plt.rc('ytick', labelsize=18);   plt.rc('legend', fontsize=18)    
    
    
    
    if performAutoCor:
        
            # Performing auto-correlation on the original raw PFB data to see what it looks like
        pfb_raw1_AC = np.mean(np.abs(pfb_raw1)**2 , axis=0)
        pfb_raw2_AC = np.mean(np.abs(pfb_raw2)**2 , axis=0)
            # Performing auto-correlation Power on the recovered and rechannelized signals. Option to trim unwanted spikes or not
        pfb_rechan1_AC_PSD = np.mean( np.abs(pfb_rechan1)**2 , axis=0)
        pfb_rechan2_AC_PSD = np.mean( np.abs(pfb_rechan2)**2 , axis=0)
        if trim:    # Get rid of artificial spikes
            pfb_rechan1_AC_PSD = getRidOfSpikes(pfb_rechan1_AC_PSD , 0.2 , 6 , doPlot=False)
            pfb_rechan2_AC_PSD = getRidOfSpikes(pfb_rechan2_AC_PSD , 0.2 , 30 , doPlot=False)
        
        
        plt.figure()
        plt.plot(pfb_raw1_AC , c='b' , label='polar_1')
        plt.plot(pfb_raw2_AC , c='r' , label='polar_2')
        plt.title("Auto-Correlation of Raw PFB Data for Both Polarizations")
        plt.xlabel("Frequency Channels"); plt.ylabel("Power Spectrum Density")
        plt.legend()
        if saveFigs:
            plt.savefig(r"dataFromEarlyFeb\Figs\f"+str(fileNum)+"_rawPFB_26chan_AC.png")
        
        if showRTS:
            plt.figure()
            plt.subplot(1,3,1);  plt.plot(np.arange(1000) , rts_pol1[:1000])
            plt.subplot(1,3,2)
            plt.plot(np.arange( int(lenRTS/2) -500 , int(lenRTS/2) +500 ) , rts_pol1[int(lenRTS/2) -500 : int(lenRTS/2) +500])
            plt.title("Recovered Timestream");  plt.xlabel("Timestep = " + str(int(nanosecsPerSamp)) + " ns")
            plt.subplot(1,3,3);  plt.plot(np.arange( lenRTS-1000,lenRTS ) , rts_pol1[-1000:])
            plt.ticklabel_format(useOffset=False)
            if saveFigs:
                plt.savefig(r"dataFromEarlyFeb\Figs\f"+str(fileNum)+"_recoveredTimestream.png")
        
        plt.figure()
        plt.plot(pfb_rechan1_AC_PSD, c='b', label='polar_1')
        plt.plot(pfb_rechan2_AC_PSD, c='r', label='polar_2')
        plt.xlabel("Frequency Channels"); plt.ylabel("Power Spectrum Density")
        plt.title("Rechannelized ("+str(reChan)+" channels) Auto-Cor For Both Polarizations")
        plt.legend()
        if saveFigs:
            plt.savefig(r"dataFromEarlyFeb\Figs\f"+str(fileNum)+"_"+str(reChan)+"chan_AutoCor_TimeSqrAv.png")
    
    
    
    
    if performCrossCor:
        
        if len(freq_ranges) == 0:
            """ Freq ranges for satellites seen by eye using rechannelized auto-cor graph
                Good Signals:
                        - [75,214]  out of 1025
                        - [350,490] out of 1025
                        - [738,900] out of 1025
                Ok Signals:
                        - [12,75]   out of 1025
                        - [215,340] out of 1025
                        - [490,615] out of 1025
                What are those:
                        - [620,740] out of 1025
            `"""
            
            freq_ranges_goodSats = np.array([[75,214],[350,490],[738,900]])
            if doAllSats:
                freq_ranges_okSats = np.array([[12,75],[215,340],[490,615]])
                freq_ranges_poorSats = np.array([[620,740]])
                freq_ranges = [ freq_ranges_goodSats , freq_ranges_okSats , freq_ranges_poorSats  ]
            else:
                freq_ranges = [freq_ranges_goodSats]
            
            freq_range_types = ['Good', 'Ok', 'Poor']
            
            
        
        performCC_analysisOnSats( fft1,fft2 , freq_ranges[0] , reChan , newSubLen=newSubLen , useChunksOfRTS=useChunksOfRTS ,
                                  dt=nanosecsPerSamp , dt_units='ns' , satDesignation=freq_range_types[0] , zpadCoeff=zpadCoeff)
        if doAllSats:
            for i in range(1,len(freq_ranges)):
                performCC_analysisOnSats( fft1,fft2 , freq_ranges[i] , reChan , newSubLen=newSubLen , useChunksOfRTS=useChunksOfRTS ,
                                      dt=nanosecsPerSamp , dt_units='ns' , satDesignation=freq_range_types[i] , zpadCoeff=zpadCoeff)
        
    
    
    end = t.time()
    print("Total execution time in 'performInv_Rechan_AC_CC' was "+str(round(end-start, 3))+" seconds")
    if returnRTS: return rts_pol1 , rts_pol2
    elif returnFFT: return fft1 , fft2
    else: return pfb_rechan1 , pfb_rechan2














def extractSatSignalsFromTimeStream( rts1 , rts2 , satRegion1 , satRegion2 , totRefChans ):
    ft1 = np.fft.rfft( rts1 )
    ft2 = np.fft.rfft( rts2 )
    
    if np.ndim(satRegion1)==1 and np.ndim(satRegion2)==1 and satRegion1[0]!=satRegion2[0] and satRegion1[1]!=satRegion2[1]:
        ft1 , ft2 , satRegion1 , satRegion2 = alignSats( ft1 , ft2 , satRegion1 , satRegion2 , totRefChans)
    
    if np.ndim(satRegion1)==1: satRegion1 = [satRegion1];
    if np.ndim(satRegion2)==1: satRegion2 = [satRegion2];
    
    rts1 = np.fft.irfft(  getSatsFromFT(ft1, satRegion1, totRefChans)[0]  ).ravel()
    rts2 = np.fft.irfft(  getSatsFromFT(ft2, satRegion1, totRefChans)[0]  ).ravel()

    return rts1 , rts2