# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:23:23 2020

@author: wilia
"""
import numpy as np
np.fft.restore_all()
import scipy
from scipy.signal import firwin, freqz, lfilter
import matplotlib.pyplot as plt
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, r"\Users\wilia\OneDrive\Documents\McGill\Academics\Endgame\pfb_tryout")
import pfb



def naiveFFT(data,sampRate):
    """
    data: 2D array containing time values in data[0], y-vals in data[1] 
    """
    N = len(data[0])
    yFTraw = np.fft.fft(data[1])
    yFT = 2.0*np.abs(yFTraw)/N
    #yInvFT = np.fft.ifft(yFTraw)
    freqs = np.fft.fftfreq(len(yFTraw))*sampRate
    mask = freqs > 0
    
    plt.figure()
    plt.suptitle("Naive FFT")
    plt.subplot(3,1,1); plt.plot( data[0],data[1] ); plt.xlabel('Time'); plt.ylabel('Signal')
    plt.subplot(3,1,2); plt.plot(freqs[mask], yFT[mask] ); plt.xlabel('Frequency Bins'); plt.ylabel('Magnitude')
    plt.subplot(3,1,3); plt.plot( freqs[mask],np.angle(yFTraw)[mask]); plt.xlabel('Frequency Bins'); plt.ylabel('Phase?')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    return yFTraw


def singleWindowPFB(data , window , mTaps , sampRate):
    """
    data:   2D array containing time values in data[0], y-vals in data[1] 
    window: 1D array containing the weight coefficients to apply to the data
    mTaps:  Number of taps
    """
    N = len(data[0])
    pBranches = int(N/mTaps)
    windedData = data[1]*window
    """
    Taking windowed data and making it a 2D array with rows corresponding to
    the data elements to be added together
    """
    windedDataT = windedData.reshape(mTaps,pBranches).T
    
    """
    Performing the polyphase filtering
    """
    sumArr = np.zeros(pBranches)
    xx = range(0,pBranches)
    for i in range(0,pBranches):
        sumArr[i] = windedDataT[i].sum()
    
    yFTraw = np.fft.fft(sumArr)
    yFT = 2.0*np.abs(yFTraw/pBranches)
    freqs = np.fft.fftfreq(len(yFTraw))*sampRate
    mask = freqs > 0
    
    plt.figure()
    plt.suptitle("PFB with only 1 window")
    plt.subplot(3,2,1); plt.plot( data[0],data[1] ); plt.ylabel('Signal')
    plt.subplot(3,2,3); plt.plot( data[0],window ); plt.xlabel('Window Func')
    plt.subplot(3,2,5); plt.plot( data[0],windedData ); plt.xlabel('Time'); plt.ylabel('Windowed Signal')
    plt.subplot(3,2,2); plt.plot( xx,sumArr ); plt.ylabel('Summed Parsed Data')
    plt.subplot(3,2,4); plt.plot( freqs[mask],yFT[mask] ); plt.xlabel('Freq Bins'); plt.ylabel('Magnitude')
    plt.subplot(3,2,6); plt.plot( freqs[mask],np.angle(yFTraw[mask]) ); plt.xlabel('Freq Bins'); plt.ylabel('Phase amplitude?')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
    return yFTraw




def multiWindowPFB(data , win_coeffs , mTaps , P  ,sampRate):
    """ 
    From the CASPER jupyter notebook 'pfb_introduction' @ https://github.com/telegraphic/pfb_introduction
    
    data:     2D array containing time values in data[0], y-vals in data[1] 
    window:   1D array containing the weight coefficients to apply to the data
    mTaps:    Number of taps
    P:        Number of Branches
    SampRate: Trivial meaning
    """
    
    p_FitsHowManyTimes = int(len(data[1])/(mTaps*P))
    sumArr = pfb_fir_frontend(data[1] , win_coeffs , mTaps , P)
    """
    this boolean value is True if mTaps*P == len(data[1]), False otherwise
    """
    bool = False
    if(p_FitsHowManyTimes == 1):
            bool = True
    
    if(bool):
        yFTraw = np.fft.fft(sumArr , P)
        yFT = 2.0*np.abs(yFTraw)/pBranches
    else:
        yFTraw = np.fft.fft(sumArr , P, axis=1)
        yFT = (2.0*np.abs(yFTraw)/pBranches).mean(axis=0)

    freqs = np.fft.fftfreq(len(yFT))*sampRate
    windedData = data[1][:mTaps*P]*win_coeffs
    mask = freqs > 0
    xx = range(0,P)
    
    print("Window fits "+str(p_FitsHowManyTimes)+" times in data")
    print("sumArr shape = " + str(sumArr.shape))
    print("yFTnotrue shape = " + str(yFTraw.shape))
    print("freqs shape = " + str(freqs.shape))
    
    plt.figure()
    plt.suptitle("PFB with "+str(int(len(data[1])/(mTaps*P)))+" windows")
    plt.subplot(3,2,1); plt.plot( data[0],data[1] ); plt.ylabel('Signal')
    plt.subplot(3,2,3); plt.plot( data[0][:mTaps*P],win_coeffs ); plt.xlabel('Window Func')
    plt.subplot(3,2,5); plt.plot( data[0][:mTaps*P],windedData ); plt.xlabel('Time'); plt.ylabel('Windowed Signal')
    if(bool):
        plt.subplot(3,2,2); plt.plot( xx,sumArr ); plt.ylabel('Summed Parsed Data')
        plt.subplot(3,2,6); plt.plot( freqs[mask],np.angle(yFTraw[mask]) ); plt.xlabel('Freq Bins'); plt.ylabel('Phase amplitude?')
    else:
        plt.subplot(3,2,2); plt.plot( xx,sumArr[0] ); plt.ylabel('Summed Parsed Data')
        plt.subplot(3,2,6); plt.plot( freqs[mask],np.angle(yFTraw.mean(axis=0)[mask]) ); plt.xlabel('Freq Bins'); plt.ylabel('Phase amplitude?')
    plt.subplot(3,2,4); plt.plot( freqs[mask],yFT[mask] ); plt.xlabel('Freq Bins'); plt.ylabel('Power?')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    return yFTraw
      
# The following 2 functions are taken from the CASPER Intro to PFBs @ https://github.com/telegraphic/pfb_introduction
def pfb_fir_frontend(x, win_coeffs, M, P):
    W = int(x.shape[0] / M / P)
    x_p = x.reshape((W*M, P)).T
    h_p = win_coeffs.reshape((M, P)).T
    if(W==1):
        return (x_p*h_p).sum(axis=1)
    else:
        x_summed = np.zeros((P, M * W - M))
        for t in range(0, M*W-M):
            x_weighted = x_p[:, t:t+M] * h_p
            x_summed[:, t] = x_weighted.sum(axis=1)
        return x_summed.T

def generate_win_coeffs(M, P, window_fn="hamming"):
    win_coeffs = scipy.signal.get_window(window_fn, M*P)
    sinc       = scipy.signal.firwin(M * P, cutoff=1.0/P, window="rectangular")
    win_coeffs *= sinc
    return win_coeffs




def pfb_Inv(rawPFB , mTaps , pBranches , sampRate):
    
    print("I tried but its now 4am and I need to sleep")


N = 4*1024                    # Number of samples
sampRate = 20.5               # Sampling Rate >= 2*freq_max
dt = 1/(sampRate)             # Time increments
T = N*dt                      # Total time of observation
mTaps = 4                     # Number of taps in PFB
pBranches = int(N/mTaps)      # Number of branches in each tap
pBsubdiv = 4                  # Used to dvide pBranches for use in multiWindowPFB implementation

w = 2*np.pi
x = np.arange(0,T,dt)
y1 = 0.2*np.sin(10.1*w*x)
y2 = 30.0*np.sin(10.2*w*x + np.pi/3)
y3 = 40.0*np.sin(10.0*w*x - np.pi/4)
noise = np.random.randint(-2,2,N)
ytot = y1+y2+y3+noise
data = np.array([x,ytot])
window1 = generate_win_coeffs( mTaps , pBranches , window_fn="hamming")
window2 = generate_win_coeffs( mTaps , int(pBranches/pBsubdiv) , window_fn="hamming")


naiveFFT(data, sampRate)
singleWindowPFB(data , window1 , mTaps , sampRate)
multiWindowPFB(data , window2 , mTaps , int(pBranches/pBsubdiv) , sampRate)

# This next portion of code compares the CHIME PFB implementation to a naive FFT
frec = np.fft.rfftfreq(pBranches)*sampRate
Nfreqs = len(frec[frec>0])
chimePFB = pfb.pfb(ytot , Nfreqs , mTaps)

plt.figure()
plt.title("Chime PFB output")
plt.plot(frec[frec>0] , (2.0*np.abs(chimePFB[0]/pBranches)))
plt.figure()
N = len(data[0])
yFTraw = np.fft.fft(data[1])
yFT = 2.0*np.abs(yFTraw)/N
freqs = np.fft.fftfreq(len(yFT))*sampRate
mask = freqs > 0
plt.title("Naive FFT output")
plt.plot( freqs[mask],yFT[mask] ); plt.xlabel('Frequency Bins'); plt.ylabel('Magnitude')





















    