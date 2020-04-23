# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 03:39:16 2020

@author: wilia

"""

import matplotlib.pyplot as plt
import numpy as np
from os import path



sat_A = 2;   satTypeA='Good'
sat_B = 3;   satTypeB='Good'



save=True
    # The relative path from the current directory to where the raw data had been stored
relPath = "dataFromEarlyFeb"
    # Name of the directory where we want to funnel the cross-cor and noise data generated
    # This directory will be created in the same directory the raw data is in
output_directory = "xcorr_"+satTypeA+"Sat"+str(sat_A)+"_wNoiseUsing_"+satTypeB+"Sat"+str(sat_B)
output_directory = path.join(relPath, output_directory)


    # The different amount of curves you want to show for each graph, based on the zero-padding behind each timing error result
zpadCoeffs = [0]

plotCrossCorPeaks = True
plotSNR = True




    # These parameters control plot label sizes
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=30)     # fontsize of the axes title
plt.rc('legend', fontsize=18)    # legend fontsize


fig1 , ax1 = plt.subplots()
ax1.set_title('Timing Error vs. Observation Time')
ax1.set_xlabel('Observation time (s)')
ax1.set_ylabel('Timing Error (ns)')
ax1.xaxis.get_label().set_fontsize(25)    
ax1.yaxis.get_label().set_fontsize(25)    
ax1.tick_params(labelsize=18)
ax1.grid()
#ax1.set_ylim(0,1.1)


if plotCrossCorPeaks:
    fig2 , ax2 = plt.subplots()
    ax2.set_title('Time-Lag (@ maximum peak) Found vs. Observation Time')
    ax2.set_xlabel('Observation time (s)')
    ax2.set_ylabel('Time-Lag Between Polarizations of Same Signal (ns)')
    ax2.xaxis.get_label().set_fontsize(25)    
    ax2.yaxis.get_label().set_fontsize(25)    
    ax2.tick_params(labelsize=18)
    ax2.grid()


if plotSNR:
    fig3 , ax3 = plt.subplots()
    ax3.set_title('Signal-to-Noise vs. Observation Time')
    ax3.set_xlabel('Observation time (s)')
    ax3.set_ylabel('SNR')
    ax3.xaxis.get_label().set_fontsize(25)    
    ax3.yaxis.get_label().set_fontsize(25)    
    ax3.tick_params(labelsize=18)
    ax3.grid()
    


colors = ['black','blue','red']

for zpadCoeff in zpadCoeffs:
    
    if zpadCoeff==0:
        zpad_str = '_noZPad'
        zpad_str_forGraph = "No Zero-Pad"
    else:
        zpad_str = '_zPadWith'+str(int(zpadCoeff))+'p'+str(int((zpadCoeff % int(zpadCoeff))*10))+'xLenOfFT'
        zpad_str_forGraph = "Zero-pad x"+str(zpadCoeff)
                    
    
    data = np.load(output_directory + r"\obsTimes_timeErr_xPeak_SNR" + "_crossCor_"+satTypeA+"Sat"+str(sat_A)+zpad_str+".npy") 
    
    ax1.plot( data[0] , data[1] , label=zpad_str_forGraph , color=colors[zpadCoeff] , marker='.')
    if plotCrossCorPeaks: ax2.plot( data[0] , data[2] , label=zpad_str_forGraph , color=colors[zpadCoeff] , marker='.')
    if plotSNR:           ax3.plot( data[0] , data[3] , label=zpad_str_forGraph , color=colors[zpadCoeff])

ax1.legend()
ax2.legend()
ax3.legend()