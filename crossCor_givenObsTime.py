# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 00:21:10 2020

@author: William Frost


In the beginning of this code, the important variables that govern how the cross-cor are run are highlighted as such:
    
    =================================================
    
    importantVariable = this
   
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are the ones you should be modifying most of the time
"""

import numpy as np
import matplotlib.pyplot as plt
import time as t
from glob import glob
import os
import gc   # Used to try and free up space as much as possible if we have to deal with big arrays. I might not quite understand
            # the 'gc' behaviour, but I figured it would not hurt to have it there after certain functions are done running
import timingErrorFuncs as f






"""======================================================================================================"""

    # sat_A is the satellite we are interested in cross-correlating (using its orthogonal polarizations)
    # sat_B is the satellite we use to determine the (over-estimated) noise by aligning the different polarizations of the A and B signals 
    # in Fourier space and cross-correlating them
    # The numbers and type are based on the prior labelling of satellite signals. They are classified first by their type (Good, Ok and Poor) based on
    # the power difference their peaks have with the ambient noise. Then the numbers indicate their order from left to right in an auto-cor plot.
    # For example, to cross-correlate the 2nd signal found in the family of signals with 'good' quality, we would designate sat_A = 2 and satTypeA='Good'
sat_A = 2;    satTypeA='Good'
sat_B = 3;    satTypeB='Good'

    # amount of different observation times we want to perform
numberOfObs = 25
    # min and max observation times in seconds
minObsTime = 5.;    maxObsTime = 125.
obsTimes = np.linspace( minObsTime , maxObsTime , numberOfObs , endpoint=True )

    # ZERO-PADDING
    # This controls how many additional zeros we add to the end of a Fourier Transform to increase interpolation between points when 
    # inverting back to time-stream space. It is based of a multiple of the original FT length. 
    # Therefore, a 'zpadCoeff' of 1 adds ~100% more points (as zeros) to the end of the FT
zpadCoeff = 1
if zpadCoeff==0: zpad_str = '_noZPad'
else:            zpad_str = '_zPadCoeff'+str(int(zpadCoeff))+'p'+str(int((zpadCoeff % int(zpadCoeff))*10))

    # The length in the x-axis of our cross-cor. It increases as we zero-pad in the cross-correlation process, but the
    # actual time intervals represented remain with the same bounds
chunk_len = 1000

    # Boolean value to check if we want to save the cross-cor data generated for future use
save=True
    # Boolean value to check if  we want to plot the cross-cor and noise outputs as we create them
doPlot=False
    # Boolean value to decide if we observe the raw files in order or in reverse
readFilesInOrder = True
    

"""^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"""





"""
This next section of code is just a setup phase before we get to observing what's inside our files
"""


    # Location of the files you want to unpack, as well as the file type we wish to analyse
absPathOfData = "/Users/wilia/OneDrive/Documents/McGill/Academics/Endgame/PHYS489_AntennaTimeAccuracy/pfb_final/dataFromEarlyFeb/*.raw"
    # The relative path from where this code is located to where the data is
relPathOfData = "dataFromEarlyFeb"
    # Name of the directory where we want to funnel the cross-cor and noise data generated. It can be named "" if the output is to be recorded in 
    # the same directory as the raw data. This output directory, if needed, will be created in the same directory the raw data is in.
output_directory = "xcorr_"+satTypeA+"Sat"+str(sat_A)+"_withNoiseUsing_"+satTypeB+"Sat"+str(sat_B)+"_obsTime"+str(int(minObsTime))+"to"+str(int(maxObsTime))+"secs"
#output_directory = ""

    # Creating the output_directory if it does not exist
cwd = os.getcwd()
output_directory = os.path.join(relPathOfData, output_directory)
dirr = os.path.join(cwd, output_directory)
if not os.path.exists(dirr):
    print("Creating a new output directory")
    os.mkdir(dirr)

    # The names of the files to read. The amount actually used will depend on observation time required
files = [os.path.join(relPathOfData, os.path.basename(x)) for x in glob(absPathOfData)]
amountOfFiles = len(files);     filesUsed = 0
    # 'fileNumOrdering' can be customized to open files in forward or backwards order
            # Ex: Have a total of 5 files. You only want 5 seconds of observations, but file #1 has 30 seconds and file #5 has 5 seconds.
            # Simply open the files backwards and waste less computation time opening a file too large for your needs.
if readFilesInOrder:
    fileNumOrdering = np.linspace(1, len(files), len(files), endpoint=True, dtype=np.int8 )
else:
    fileNumOrdering = np.linspace(len(files), 1, len(files), endpoint=True, dtype=np.int8 )







    # From the (rechannelized) auto-cor data of any file, these are the regions in fourier space showing satellite signals
    # THESE HAVE TO BE SET MANUALLY. UNLESS SOMEONE SMART DEVELOPS A GOOD ALGO TO DO IT.
    # For now, one can use the 'performInv_Rechan_AC_CC' found in the 'timingAccuracyFuncs.py' to perform this auto-correlation
freq_ranges_good = np.array([[75,214],[350,490],[738,900]])
freq_ranges_ok = np.array([[12,75],[215,340],[490,615]])
freq_ranges_poor = np.array([[620,740]])
    # Number of channels in the Fourier space plot used to determine satellite frequency ranges/channels
totRefChans = 1025


    # Used to refer to which channel region in our Fourier space we are refering to. The labeling is meant
    # to indicate which signals are the strongest and thus will have better Signal to Noise Ratio (SNR)
if satTypeA=='Good' : freqsA = freq_ranges_good
elif satTypeA=='Ok' : freqsA = freq_ranges_ok
elif satTypeA=='Poor' : freqsA = freq_ranges_poor
if satTypeB=='Good' : freqsB = freq_ranges_good
elif satTypeB=='Ok' : freqsB = freq_ranges_ok
elif satTypeB=='Poor' : freqsB = freq_ranges_poor




    # units of time in our recovered signals. 
    # It is a function of the ADC channels, sampling rate (based on Nyquist) and the number of channels taken from the total available
    # The 'ADC_chans' and 'sampRate' values here are based of the Casper SNAP board specs used for ALBATROS
ADC_chans = 2048;   sampRate = 250e6;   numChansTaken = 26
dtSecs = (ADC_chans/sampRate/numChansTaken)
dtNanoSecs = dtSecs*10**9





cc_len = int( chunk_len*(zpadCoeff+1) ) + int(2*zpadCoeff)
    # Initializing the crossCor and noise arrays
crossCor = np.zeros( cc_len )
noise = np.zeros( cc_len )
    # This variable is used to adjust the 'dt' value for accurate x-axis values regardless of the value of 'zpadCoeff'
zpad_adjust = dtNanoSecs / ( 1+zpadCoeff + ( 2*zpadCoeff/(2*zpadCoeff + chunk_len*(zpadCoeff+1)) ) )
    # Creating the x-axis values for our cross-cor and noise graohs
time_lag = f.x_fftshift( np.arange( cc_len ) ) * zpad_adjust




if save:
        # Saving the x-axis of the cross-correlation and the observation times for this run
    np.save( str(output_directory)+r"\ccLen"+str(cc_len)+zpad_str+".npy" , time_lag )
    np.save( str(output_directory)+r"\obsTimes"+".npy" , obsTimes )
    

    # These parameters control plot label sizes
plt.rc('font', size=14)
plt.rc('axes', titlesize=30);     plt.rc('axes', labelsize=25)    
plt.rc('xtick', labelsize=18);    plt.rc('ytick', labelsize=18);    plt.rc('legend', fontsize=18)    

 
     

"""
Now that the setup process is done, the calculations may begin
"""




  
t1 = t.time()
    # 'current_rta_len' will keep track of the total length (unitless) in time-stream space that has been observed up until now
current_rts_len = 0
    # Boolean value used to check if we need to change to the next file because the current one has been fully explored
doFileChange = True


    # Go through all the observation times you desire
for i in range(numberOfObs):
    
    
        # The desired unitless length in time-stream space for a given observation time
    wanted_rts_len = int(obsTimes[i]/dtSecs)
    
        # so long as we haven't reached our deisred observation time or that we haven't run out of files to look at, we continue 
    while ( current_rts_len < wanted_rts_len   and   filesUsed < amountOfFiles ):
        
        
        
            # If we have looked at everything in our current file, or haven't started observing yet, 
            # we switch to the next one and unpack its contents
        if doFileChange:
            
            doFileChange = False
            currentFile = fileNumOrdering[filesUsed]
            print("\n---------------\nABOUT TO LOAD IN DATA FROM FILE #" + str(currentFile) + "\n---------------")
                  # Getting the timestreams for both polarizations in a given file
            rts1 , rts2 = f.performInv_Rechan_AC_CC( files[currentFile-1] , currentFile , returnFFT=False , returnRTS=True , 
                                                                      originalNumChans=ADC_chans , originalSampRate = sampRate )
            gc.collect()
            rts_len = len(rts1)
                # Getting the time-streams in preparation for cross-correlation
            rts1_cc , rts2_cc = f.extractSatSignalsFromTimeStream( rts1 , rts2 , freqsA[sat_A-1] , freqsA[sat_A-1] , totRefChans )
            rts1_noise , rts2_noise = f.extractSatSignalsFromTimeStream( rts1 , rts2 , freqsA[sat_A-1] , freqsB[sat_B-1] , totRefChans )
            gc.collect()
                # Setting the index pointer of a file to 0, as we are exploring a new file
            indexInRTS = 0
        
        
        
            # As long as we have not reached our desired integration time and the file we are in is not fully explored, continue to
            # cross-correlate more chunks of the satellite signals
        while ( current_rts_len < wanted_rts_len   and   indexInRTS+chunk_len < rts_len ):
            
            crossCor = crossCor + np.fft.fftshift(   f.doCrossCor_zpad( rts1_cc[ indexInRTS:indexInRTS+chunk_len ] , rts2_cc[ indexInRTS:indexInRTS+chunk_len ] , zpadCoeff ) )
            noise = noise + np.fft.fftshift(   f.doCrossCor_zpad( rts1_noise[ indexInRTS:indexInRTS+chunk_len ] , rts2_noise[ indexInRTS:indexInRTS+chunk_len ] , zpadCoeff ) )
            indexInRTS += chunk_len
            current_rts_len += chunk_len
        
        
        
        # If we have completely finished looking at a file, switch to the next one
        if (indexInRTS+chunk_len >= rts_len   and   current_rts_len < wanted_rts_len):
            print("\nFile #" +str(fileNumOrdering[filesUsed])+ " did not contain enough obsTime (had "+str(round(rts_len*dtSecs, 5))
                   +"s) to reach a total of "+str(obsTimes[i])+"s. Switching to the next file")
            filesUsed += 1
            doFileChange = True
            continue
        # If the file we were looking at had enough observation time, we save the cross-cor and noise data we have for now
        # and continue observing based on the next desired observation time.
        elif (current_rts_len >= wanted_rts_len):
            print("\nFile #" +str(fileNumOrdering[filesUsed])+ " contained enough obsTime."+
                  "\nObsTime needed was "+str(round(wanted_rts_len*dtSecs,5))+ "s and what was used was "+str(round((current_rts_len-chunk_len)*dtSecs,5)) + "s" )
            break
        
    
    
    
    t2 = t.time()
    
    if filesUsed >= amountOfFiles:
        print("\nYou do not have enough files to observe " +str(obsTimes[i])+ "s.  You only got to see " + str(round(current_rts_len*dtSecs, 5))+"s")
        break
    else:    
        print("\nObserving "+str(obsTimes[i])+" seconds to get cross-cor AND noise of satellite with "
               +str(zpadCoeff)+"x zero-padding took " +str(round(t2-t1,5))+ " secs to complete")
    
    
    
        # Saving cross-cor and noise data for a given observation time
    if save:
            # Some technicalities regarding the naming of files
        if filesUsed >= amountOfFiles:
            obsTimes[i] = round(current_rts_len*dtSecs, 3)
        obsTime_str = str(int(obsTimes[i]))+'p'+str(int((obsTimes[i] % int(obsTimes[i]))*1000))
            # Saving Cross-cor data
        np.save(output_directory+r"\ccLen"+str(cc_len)+"_crossCor_"+satTypeA+"Sat"+str(sat_A)+"_obsTime"+obsTime_str+"secs"+zpad_str+".npy" , crossCor )
            # Saving Noise data
        if satTypeA==satTypeB:
            np.save(output_directory+r"\ccLen"+str(cc_len)+"_noiseCC_"+satTypeA+"Sats"+str(sat_A)+"and"+str(sat_B)+"_obsTime"+obsTime_str+"secs"+zpad_str+".npy" , noise )
        else:
            np.save(output_directory+r"\ccLen"+str(cc_len)+"_noiseCC_"+satTypeA+'Sat'+str(sat_A)+satTypeB+'Sat'+str(sat_B)+"_obsTime"+obsTime_str+"secs"+zpad_str+".npy" , noise )
    
     
        
        
    if doPlot:
        plt.figure()
        plt.title("Cross-Cor of " + satTypeA +'Sat'+ str(sat_A) + ' ('+zpad_str+') (obsTime = '+str(obsTimes[i])+'secs)')
        plt.xlabel("Time Lag (ns)")
        plt.ylabel("Cross-Correlation")
        plt.grid()
        plt.plot( time_lag , crossCor )
        plt.xlim(-20000,20000)
        
        plt.figure()
        plt.title('Noise Using '+satTypeA+'Sat'+str(sat_A)+' and '+satTypeB+'Sat'+str(sat_B)+ ' ('+zpad_str+') (obsTime = '+str(obsTimes[i])+'secs)' )
        plt.xlabel("Time Lag (ns)")
        plt.ylabel("Cross-Correlation")
        plt.grid()
        plt.plot( time_lag , noise )
        plt.xlim(-20000,20000)
        
    
    


    
        
