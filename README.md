# PHYS489_ProjectCode

Short Description:
  - 'show_AutoCor_CrossCor.py' is used to determine what channels in our FT correspond to satellite signals. It can also be used to plot        preliminary cross-correlations of the deduced satellite signals.

  - 'crossCor_givenObsTime.py' outputs cross-correlation and noise data for a set amount of different observation times. Used to save the     cross-correlation and noise data to disk.

  - 'findTimingError_givenObsTime.py' outputs timing errors in the cross-correlation of a given observation time. Usesthe output of           'crossCor_givenObsTime.py' to perform it's calculations.

  - 'plotTimingErrors_vs_ObsTime.py' uses the timing errors output by 'findTimingError_givenObsTime.py' and plots them as a function of        observation time.

  - 'timingErrorFuncs.py' contains the functions used in 'show_AutoCor_CrossCor.py' and 'crossCor_givenObsTime.py' .

  - 'monteCarloFuncs.py' contains the functions used in 'findTimingError_givenObsTime.py' .




Detailed Description:

What this code repository accomplishes is taking in .raw files which contains Fourier space data of 2 orthogonal polarizations of a signal and performs cross-correlation between certain channel regions in that Fourier space. These can then be used to quantify the timing error in a satellite signal with respect to integration/observation time. This can be used to estimate the synchronization capabilites of that satellite signal when used as a clock mechanism in radio interferometry. To obtain this timing error, the following protocol is followed using the pieces of code in this repo:
    
    
    1) Auto-correlate the Fourier Space data unpacked from the .RAW files to visualize which channel regions correspond to satellite            signals. You might need to rechannelize this unpacked data to better resolve signals of interest. 
       NOTE: This can be accomplished using the   'show_AutoCor_CrossCor.py'   code.


  NOTE: Steps 2-3-4 are performed by  'crossCor_givenObsTime.py'
    
    2) Once satellite channel regions have been identified (manually I'm afraid), isolate those channels (zero-out the others) and              invert back to a time-stream to obtain what should be a 'pure' satellite signal.
    
    3) With the relevant signals now obtained, cross-correlation can be performed on both polarizations of a satellite signal of                interest. To quantify the noise in this cross-correlation, an over-estimate of it is produced by cross-correlating different            polarizitions of different satellite signals.
    
    4) For each observation time, different supplemental techniques can be applied before cross-correlation. One of them is zero-padding        the fourier transform, where the end result is a better interpolation between the points of the cross-correlations. To note that        this simply smooths out the data representing our cross-cor, and does really improve resolution. The option to zeropad or not is        given as an option in this code.

      
    5) Once cross-correlations and their respective noises are obtained/saved for each observation time, a Monte Carlo approach can be          used to estimate the timing error. For simplicity, the timing error is defined as the standard deviation of many least-squares          quadratic fits performed on the cross-correlation peak. Again, for simplicity, this peak is defined to be the point with maximal        distance to the time axis (where cross-cor = 0). The cross-correlation data used for each least-squares fit is changed by                randomly shifting the noise function and adding it in to the cross-cor.
       NOTE: This MC approach is accomplished using the   'findTimingErrors_givenObsTime.py'   code to estimate the phase error in a            satellite signal.
       
    
    6) To plot timing error as a function of observation time, simply use   'plotTimeErrors_vs_ObsTime.py'   , importing the relevant          saved outputs generated from   'findTimingErrors_givenObsTime.py' .
