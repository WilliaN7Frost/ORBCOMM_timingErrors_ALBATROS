# PHYS489_ProjectCode

- 'show_AutoCor_CrossCor.py' is used to determine what channels in our FT correspond to satellite signals. It can also be used to plot preliminary cross-correlations of the deduced satellite signals.

- 'crossCor_givenObsTime.py' outputs cross-correlation and noise data for a set amount of different observation times. Used to save the cross-correlation and noise data to disk.

- 'findTimingError_givenObsTime.py' outputs timing errors in the cross-correlation of a given observation time. Usesthe output of 'crossCor_givenObsTime.py' to perform it's calculations.

- 'plotTimingErrors_vs_ObsTime.py' uses the timing errors output by 'findTimingError_givenObsTime.py' and plots them as a function of observation time.

- 'timingErrorFuncs.py' contains the functions used in 'show_AutoCor_CrossCor.py' and 'crossCor_givenObsTime.py' .

- 'monteCarloFuncs.py' contains the functions used in 'findTimingError_givenObsTime.py' .
