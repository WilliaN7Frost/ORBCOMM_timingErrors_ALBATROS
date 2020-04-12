# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 22:28:34 2020

@author: wilia
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


    # Useful to integrate only portions of the signal to test out different integration time
def integrate(signal , intTime):
    signal = signal[:np.round(signal.shape[0]//intTime)*intTime]
    signal = signal.reshape(signal.shape[0]//intTime, intTime, signal.shape[1])
    return np.mean(signal , axis=1)


    # Used to remove sharp risig edges due to unwanted code-generated spikes
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



    # Used to perform cross correlation on satellite signals found in Fourier Space
def performCC_analysisOnSats( ft1,ft2 , satRegions , totRefChans , dt , ft_isPFB=False , newSubLen=250 ,
                              dt_units='unspecified' , satDesignation='' , zeropad=0):
    """
    ft1 , ft2 : 1D or 2D array
            the fourier transforms containing intervals of interest (in this case satellite signals)
    satRegions : 2D array 
            Each element of this array contains lower and upper bounds that constrain channel regions of interest
    totRefChans , dt , newSubLen , zeropad : int
            - 'totRefChans' is the total number of channels in the fourier space used to determine the 'satRegions'
            - 'dt' is the time increment in the recovered timestreams
            - 'newSubLen' is the desired length of our cross-cor data. Shouldn't be too big since we expect satellite cross-cors to be about 0
            - 'zeropad' is the amout of zeros we want to add to our fourier space
    ft_isPFB : boolean
            True if ft1, ft2 have passed thru a PFB frontend before being FT'd . False if they come from straight FFT
    """
    if np.ndim(satRegions)==1: satRegions = [satRegions];
    if zeropad != 0: newSubLen = 500
    
    sats_ft1 = getSatsFromFT(ft1 , satRegions , totRefChans)
    sats_ft2 = getSatsFromFT(ft2 , satRegions , totRefChans)
    #print("\n Preparing to zeropad")
    #print("Shape of initial ft = " + str(sats_ft1.shape))
    if ft_isPFB:
        sats_rts1 = pfb.inverse_pfb(sats_ft1[0] , 4).ravel()
    else:
        sats_rts1 = np.fft.irfft(sats_ft1[0]).ravel()
    lenCC = len(sats_rts1)
    
    sats_ft1_0pad = []
    sats_ft2_0pad = []
    for i in range(len(satRegions)):
        sats_ft1_0pad.append( np.append( sats_ft1[i], np.zeros(zeropad)) )
        sats_ft2_0pad.append( np.append( sats_ft2[i], np.zeros(zeropad)) )
    sats_ft1 = np.array(sats_ft1_0pad)
    sats_ft2 = np.array(sats_ft2_0pad)
    #print("Length of zeropadded ft = " + str(sats_ft1.shape))
    gc.collect()

    sats_CC = []
    dummyBool=True
    for i in range(len(satRegions)):
        if ft_isPFB:
            sats_rts1 = pfb.inverse_pfb(sats_ft1[i] , 4).ravel()
            sats_rts2 = pfb.inverse_pfb(sats_ft2[i] , 4).ravel()
            lenCC_0pad = len(sats_rts1)
            if dummyBool:
                dt = dt * lenCC/lenCC_0pad
                #print("new dt is = " + str(int(dt)))
                dummyBool=False
            #print("In PFB rts length is " + str(lenCC))
        else:
            sats_rts1 = np.fft.irfft(sats_ft1[i]).ravel()
            sats_rts2 = np.fft.irfft(sats_ft2[i]).ravel()
            lenCC_0pad = len(sats_rts1)
            #print("lenCC = "+str(lenCC)+" and lenCC_0pad = "+str(lenCC_0pad))
            if dummyBool:
                dt = dt * lenCC/lenCC_0pad
                print("new dt is = " + str(int(dt)))
                dummyBool=False
            #print("In FFT rts length is " + str(lenCC))
        
            # Reshaping the signal before passing it through the cross-cor. Allows to discard regions far from the peaks we want to look at
        newLength = int(lenCC_0pad/newSubLen)*newSubLen
        sats_rts1 = sats_rts1[:newLength].reshape(-1,newSubLen)
        sats_rts2 = sats_rts2[:newLength].reshape(-1,newSubLen)
        #print("New rts shape is " + str(sats_rts1.shape))
        
        sats_CC.append( doCrossCor( sats_rts1 , sats_rts2 ) )
    
    sats_CC = np.mean( sats_CC , axis=1 )
    #print("CrossCor of sats after mean is " + str(sats_CC.shape))
    
    plt.figure()
    plt.title("Cross-Correlation of " +satDesignation+ " Satellite Sources")
    plt.xlabel("Time Lag ("+dt_units+")")
    plt.ylabel("Cross-Correlation")
    for i in range(len(satRegions)):
        plt.plot( x_fftshift( np.arange(newSubLen) )*dt , np.fft.fftshift( sats_CC[i] ) , label=("Sat %s" % str(i+1)) )
    plt.legend()
    plt.grid()





def performCC_diffSats( ft1,ft2 , satRegion1,satRegion2 , totRefChans , dt , ft_isPFB=False , newSubLen=250 ,
                        dt_units='' , samePol=True , zeropad=0):
    """
    ft1 , ft2 : 1D or 2D array
            the fourier transforms containing intervals of interest (in this case satellite signals)
    satRegion1,2 : 2D arrays
            Different fourier space regions containing lower and upper bounds that constrain channel regions of interest in ft1, ft2 respectively
    totRefChans , dt , newSubLen , zeropad : int
            - 'totRefChans' is the total number of channels in the fourier space used to determine the 'satRegions'
            - 'dt' is the time increment in the recovered timestreams
            - 'newSubLen' is the desired length of our cross-cor data. Shouldn't be too big since we expect satellite cross-cors to be about 0
            - 'zeropad' is the amout of zeros we want to add to our fourier space
    ft_isPFB , samePol : boolean
            - 'ft_isPFB' is true if ft1, ft2 have passed thru a PFB frontend before being FT'd . False if they come from straight FFT
            - 'samePol' is simply used to mention in the graph title if the signals we cross-cor are from  the same polarisation or not
    """
    if np.ndim(satRegion1)==1: satRegion1 = [satRegion1];
    if np.ndim(satRegion2)==1: satRegion2 = [satRegion2];
    if zeropad != 0: newSubLen=500
    
    sat_ft1 = getSatsFromFT(ft1 , satRegion1 , totRefChans)
    sat_ft2 = getSatsFromFT(ft2 , satRegion2 , totRefChans)
    print("\n Preparing to zeropad")
    print("Shape of initial ft = " + str(sat_ft1.shape))
    if ft_isPFB:
        sat_rts1 = pfb.inverse_pfb( sat_ft1[0] , 4 ).ravel()
    else:
        sat_rts1 = np.fft.irfft(sat_ft1[0]).ravel()
    lenCC = len(sat_rts1)
    
    sat_ft1 = np.array( np.append( sat_ft1, np.zeros(zeropad)) )
    sat_ft2 = np.array( np.append( sat_ft2, np.zeros(zeropad)) )
    print("Length of zeropadded ft = " + str(sat_ft1.shape))
    gc.collect()
    
    if ft_isPFB:
        sat_rts1 = pfb.inverse_pfb(sat_ft1 , 4).ravel()
        sat_rts2 = pfb.inverse_pfb(sat_ft2 , 4).ravel()
        lenCC_0pad = len(sat_rts1)
        dt = dt *(lenCC/lenCC_0pad)
        #print("In PFB rts length is " + str(lenCC))
    else:
        sat_rts1 = np.fft.irfft(sat_ft1).ravel()
        sat_rts2 = np.fft.irfft(sat_ft2).ravel()
        lenCC_0pad = len(sat_rts1)
        dt = dt *(lenCC/lenCC_0pad)
        #print("In FFT rts length is " + str(lenCC))
    
        # Reshaping the signal before passing it through the cross-cor. Allows to discard regions far from the peaks we want to look at
    newLength = int(lenCC/newSubLen)*newSubLen
    sat_rts1 = sat_rts1[:newLength].reshape(-1,newSubLen)
    sat_rts2 = sat_rts2[:newLength].reshape(-1,newSubLen)
    #print("New rts shape is " + str(sats_rts1.shape))
    
    time_lag = x_fftshift( np.arange(newSubLen) )*dt 
    crossCor = np.fft.fftshift(  np.mean( doCrossCor(sat_rts1,sat_rts2) , axis=0)  )
    
    if samePol==True: s = '(Same Polars)'
    else: s = '(Diff Polars)'
    plt.figure()
    plt.title("Cross-Correlation of 2 Different Satellite Sources "+s)
    plt.xlabel("Time Lag ("+dt_units+")")
    plt.ylabel("Cross-Correlation")
    plt.plot( time_lag , crossCor )
    plt.grid()
    
    return time_lag , crossCor
    
    


    # The master function which allows to rechannelize PFB data, with options of plotting auto-cor and cross-cor plots
def performRechan_AC_CC(fileID , fileNum , reChan=1025 , originalNumChans=2048 , originalSampRate = (125*10**6) , saveFigs = False , returnFFT=True ,
                    trim=True , showRTS=False , performAutoCor=True , performCrossCor=True , satRegions=[] , doAllSats=False , ccUsingPFB=False ,
                    addZerosFT=0 , newSubLen=300):
    
    start = t.time()
        # Extracting PFB data from the files (both polarizations and the channels used)
    pfb_raw1 , pfb_raw2 = read_4bit.read_4bit_new(fileID)
    rows , nchan = pfb_raw1.shape
    print('\nRaw data dims (rows,nchan) = (' + str(rows) + ',' + str(nchan) + ')')
    
    # calculate the effective time per sample in the recovered timestream
    # divide by 2 because we assume sampled at Nyquist rate
    dt = originalNumChans / originalSampRate / nchan / 2
    
        # Performing the inverse PFB on both polarizations. Assume nTaps=4
    rts_pol1 = pfb.inverse_pfb(pfb_raw1 , 4)
    rts_pol2 = pfb.inverse_pfb(pfb_raw2 , 4)
    print("Inverted PFB shape = "+str(rts_pol1.shape))
    
        # Ravel the output of the inverse PFB (make it one big timestream)
    ravelRTS_pol1 = rts_pol1.ravel()
    ravelRTS_pol2 = rts_pol2.ravel()
    lenRTS = len(ravelRTS_pol1)
    print("Recovered timestream length = " + str(lenRTS))
    print("Effective 'dt' in recovered timestream is " + str(int(dt*10**9)) + " nanoseconds")
    print("Total time length of recovered timestream = "+ str(lenRTS*dt) + " seconds")
    nanosecsPerSamp = dt * 10**9
    
        # Reperform the PFB (both polarizations) with larger number of channels, while excluding the
        # edges of the recovered timestream due to probable high residual difference with original timestream
    pfb_rechan1 = pfb.pfb(ravelRTS_pol1[1000:-1000] , reChan)
    pfb_rechan2 = pfb.pfb(ravelRTS_pol2[1000:-1000] , reChan)
    print("Rechannelized PFB ("+str(reChan)+" channels) shape = " + str(pfb_rechan1.shape))
        # Do rechannelization for FFT as well since useful down the line
    fft1 = np.fft.rfft(ravelRTS_pol1[1000:-1000])
    fft2 = np.fft.rfft(ravelRTS_pol2[1000:-1000])
    print("Rechannelized FFT (max channels) shape = " + str(fft1.shape))
    
        # These parameters control plot label size
    plt.rc('font', size=14)          # controls default text sizes
    plt.rc('axes', titlesize=30)     # fontsize of the axes title
    plt.rc('axes', labelsize=25)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    plt.rc('legend', fontsize=18)    # legend fontsize
    
    
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
            plt.subplot(1,3,1);  plt.plot(np.arange(1000) , ravelRTS_pol1[:1000])
            plt.subplot(1,3,2)
            plt.plot(np.arange( int(lenRTS/2) -500 , int(lenRTS/2) +500 ) , ravelRTS_pol1[int(lenRTS/2) -500 : int(lenRTS/2) +500])
            plt.title("Recovered Timestream");  plt.xlabel("Timestep = " + str(int(nanosecsPerSamp)) + " ns")
            plt.subplot(1,3,3);  plt.plot(np.arange( lenRTS-1000,lenRTS ) , ravelRTS_pol1[-1000:])
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
        
            # The number of zeros we want to add to the end of our fourier space data. Used to improve cross-cor resolution
        zeropad = addZerosFT
        
            # The chosen FT to perform cross-cor on is PFB data
        if ccUsingPFB:
            performCC_analysisOnSats( pfb_rechan1,pfb_rechan2 , freq_ranges_goodSats , reChan , newSubLen=newSubLen ,
                                      dt=nanosecsPerSamp , ft_isPFB=True , dt_units='ns' , satDesignation='Good' , zeropad=zeropad)
            if doAllSats:
                performCC_analysisOnSats( pfb_rechan1,pfb_rechan2 , freq_ranges_okSats , reChan , newSubLen=newSubLen ,
                                          dt=nanosecsPerSamp , ft_isPFB=True , dt_units='ns' , satDesignation='Mediocre' , zeropad=zeropad)
                performCC_analysisOnSats( pfb_rechan1,pfb_rechan2 , freq_ranges_poorSats , reChan , newSubLen=newSubLen ,
                                          dt=nanosecsPerSamp , ft_isPFB=True , dt_units='ns' , satDesignation='Poor' , zeropad=zeropad)
        
        # The chosen FT to perform cross-cor on is FFT data
        else:
            performCC_analysisOnSats( fft1,fft2 , freq_ranges_goodSats , reChan , newSubLen=newSubLen ,
                                      dt=nanosecsPerSamp , ft_isPFB=False , dt_units='ns' , satDesignation='Good' , zeropad=zeropad)
            if doAllSats:
                performCC_analysisOnSats( fft1,fft2 , freq_ranges_okSats , reChan , newSubLen=newSubLen ,
                                          dt=nanosecsPerSamp , ft_isPFB=False , dt_units='ns' , satDesignation='Mediocre' , zeropad=zeropad)
                performCC_analysisOnSats( fft1,fft2 , freq_ranges_poorSats , reChan , newSubLen=newSubLen ,
                                          dt=nanosecsPerSamp , ft_isPFB=False , dt_units='ns' , satDesignation='Poor' , zeropad=zeropad)
        
    
    end = t.time()
    print("Total execution time was "+str(end-start)+" seconds")
    if returnFFT: return fft1 , fft2
    else: return pfb_rechan1 , pfb_rechan2







