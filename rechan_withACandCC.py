# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:52:43 2020

@author: wilia
"""

import numpy as np
import matplotlib.pyplot as plt
import read_4bit
import pfb
import projectFunctions as f
import time as t
import gc


    # The names of the files to read
file1 = r'dataFromEarlyFeb\935278822.raw'
file2 = r'dataFromEarlyFeb\935278854.raw'
file3 = r'dataFromEarlyFeb\935278887.raw'
file4 = r'dataFromEarlyFeb\935278920.raw'
file5 = r'dataFromEarlyFeb\935278953.raw'

    # From the data, these are the regions in fourier space showing satellite signals
freq_ranges_good = np.array([[75,214],[350,490],[738,900]])
freq_ranges_ok = np.array([[12,75],[215,340],[490,615]])
freq_ranges_poor = np.array([[620,740]])

dt = (2048/125e6/26/2)*10**9


#f.performRechan_AC_CC( file4 , 4 , ccUsingPFB=False , showRTS=True , doAllSats=True , returnFFT=True)

ft1 , ft2 = f.performRechan_AC_CC( file5 , 5 , performAutoCor=False , performCrossCor=True , returnFFT=True ,
                                    doAllSats=True , addZerosFT=3*10**7 )

f.performCC_diffSats( ft1,ft2 , [freq_ranges_good[1]],[freq_ranges_good[2]] , 1025 , dt , dt_units='ns' , samePol=False , zeropad=3*10**7 )
f.performCC_diffSats( ft1,ft1 , [freq_ranges_good[1]],[freq_ranges_good[2]] , 1025 , dt , dt_units='ns' , samePol=True , zeropad=3*10**7)

gc.collect()

