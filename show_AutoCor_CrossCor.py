# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 04:21:13 2020

@author: wilia
"""

import numpy as np
from glob import glob
from os import path
import timingErrorFuncs as f





"""============================================================================"""

    # Location of the files you want to unpack
absPath = "/Users/wilia/OneDrive/Documents/McGill/Academics/Endgame/PHYS489_AntennaTimeAccuracy/pfb_final/dataFromEarlyFeb/*.raw"
    # The relative path from where this code is located to where the data is
relPath = "dataFromEarlyFeb"

"""^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"""
    # The names of the files to read. The amount actually used will depend on observation time required
files = [path.join(relPath, path.basename(x)) for x in glob(absPath)]
    
    # The file which you want to use in the cross-cor and auto cor. I reccomend a small file size since we only want a good
    # idea of what things look like, and this does not require a large observation time
fileNum = 5
fileChosen = files[fileNum-1]

    # The number of channels we want to use to use during rechannelization
rechan = 1025
    # Boolean value to determine if we perform cross-cor on the deduced satellite signal regions
doCC = True


ADC_chans = 2048;   sampRate = 250e6;   numChansTaken = 26
dtSecs = (ADC_chans/sampRate/numChansTaken)
dtNanoSecs = dtSecs*10**9




    # Perform rechannelized auto-cor to view and select channel regions of interest. Return FFTs of the recovered time-streams
ft1 , ft2 = f.performInv_Rechan_AC_CC( fileChosen , fileNum , reChan=rechan , performAutoCor=True , returnFFT = True , returnRTS=False , 
                                       originalNumChans=ADC_chans , originalSampRate=sampRate )


    # An initial cross-cor of the satellite signals found in the auto-cor graph can be peformed here
if doCC:
        # From the (rechannelized) auto-cor data of any file, these are the regions in fourier space showing satellite signals
        # ThESE MUST BE SET MANUALLY ONCE YOU HAVE SEEN THE AUTO-COR GRAPH
    freq_ranges_good = np.array([[75,214],[350,490],[738,900]])
    freq_ranges_ok = np.array([[12,75],[215,340],[490,615]])
    freq_ranges_poor = np.array([[620,740]])
        # Number of channels in the Fourier space plot used to determine satellite frequency ranges/channels
    totRefChans = 1025
    
    
    freq_ranges = [freq_ranges_good , freq_ranges_ok , freq_ranges_poor]
    freq_range_types = ['Good' , 'Ok' , 'Poor']
    doAllSats = True
    
    zpadCoeff = 1

    """
        # Cross-cor of Good Sat #2
    f.performCC_on2Sats( ft1,ft2 , freq_ranges_good[1],freq_ranges_good[1] , totRefChans , dtNanoSecs , newSubLen=500 ,
                        dt_units='' , zpadCoeff=zpadCoeff , plotTitle='' , doPlot=True , doPrint=False , useChunksOfRTS=True)
        # Noise Cross-cor of Good Sats #2 and #3
    f.performCC_on2Sats( ft1,ft2 , freq_ranges_good[1],freq_ranges_good[2] , totRefChans , dtNanoSecs , newSubLen=500 ,
                        dt_units='' , zpadCoeff=zpadCoeff , plotTitle='' , doPlot=True , doPrint=False , useChunksOfRTS=True)
    """
    
    f.performCC_analysisOnSats( ft1,ft2 , freq_ranges_good , totRefChans , dtNanoSecs , newSubLen=500 ,
                              dt_units='ns' , satDesignation=freq_range_types[0] , zpadCoeff=zpadCoeff , useChunksOfRTS=True)
    if doAllSats:
        f.performCC_analysisOnSats( ft1,ft2 , freq_ranges_ok , totRefChans , dtNanoSecs , newSubLen=500 ,
                              dt_units='ns' , satDesignation=freq_range_types[1] , zpadCoeff=zpadCoeff , useChunksOfRTS=True)
        
        f.performCC_analysisOnSats( ft1,ft2 , freq_ranges_poor , totRefChans , dtNanoSecs , newSubLen=500 ,
                              dt_units='ns' , satDesignation=freq_range_types[2] , zpadCoeff=zpadCoeff , useChunksOfRTS=True)
    
    



