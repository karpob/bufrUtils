#!/usr/bin/env python3
import h5py
import numpy as np
from matplotlib import pyplot as plt
import argparse, glob, os

def PlanckFreqInv(wn,Radiance):
    """
    Convert to Brightness temperature. For a given radiance.
    Radiance -> Radiance in mW/(m2-sr-cm-1)
    wn->wavenumbers in cm^-1
    <-T in Kelvins
    """
    C2_inv_cm = 1.4387752
    C1_inv_cm = 1.191042e-5
    B=np.abs(Radiance)
    print(B)
    T=(C2_inv_cm*wn)/(np.log(1.+(np.power(wn,3)*C1_inv_cm)/B))
    return T

def go(cntlPath, expPath):
    cntList = glob.glob(os.path.join(cntlPath,'SCRIF*.h5'))
    expList = glob.glob(os.path.join(expPath, 'SCRIF*.h5'))
    lwAve = []
    mwAve = []
    swAve = []
    cntList.sort()
    expList.sort() 
    cntTotal = 0
    for i,f in enumerate(cntList):
        print(i, len(cntList), f)
        if(os.path.basename(f)[0:46] != os.path.basename(expList[i])[0:46]): continue 
        cLw,cMw,cSw = readSdr(f)
        eLw,eMw,eSw = readSdr( expList[i] )
        lwAveTmp = eLw - cLw
        mwAveTmp = eMw - cMw
        swAveTmp = eSw - cSw
        cnt = 4*30*9 
        print (lwAveTmp.shape, mwAveTmp.shape, swAveTmp.shape)
        if (len(lwAve) == 0) and (len(mwAve) == 0.0) and (len(swAve) == 0 ):
            lwAve = lwAveTmp.mean(axis=0).mean(axis=0).mean(axis=0)
            mwAve = mwAveTmp.mean(axis=0).mean(axis=0).mean(axis=0)
            swAve = swAveTmp.mean(axis=0).mean(axis=0).mean(axis=0)
        else:
            lwAve = (cntTotal*lwAve+cnt*lwAveTmp.mean(axis=0).mean(axis=0).mean(axis=0))/(cnt+cntTotal)
            mwAve = (cntTotal*mwAve+cnt*mwAveTmp.mean(axis=0).mean(axis=0).mean(axis=0))/(cnt+cntTotal)
            swAve = (cntTotal*swAve+cnt*swAveTmp.mean(axis=0).mean(axis=0).mean(axis=0))/(cnt+cntTotal)

        cntTotal+= cnt
        if(i==10): break
    wnLw = np.arange(648.75,1096.75,0.625)
    wnMw = np.arange(1208.75,1751.25+0.625,0.625)
    wnSw = np.arange(2153.75,2551.25+0.625,0.625)
    plt.figure()        
    plt.plot(wnLw[1:-1], lwAve)
    plt.savefig('lw.png')
    plt.figure()
    plt.plot(wnSw[1:-1], swAve)
    plt.savefig('sw.png')
    plt.figure()
    plt.plot(wnMw[1:-1], mwAve)
    plt.savefig('mw.png')
def smooth(x,window_len=3,window='hamming'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),x,mode='valid')
    return y
    
def readSdr( fname ):
    h5 = h5py.File(fname,'r')
    wnLw = np.arange(648.75,1096.75,0.625)
    wnMw = np.arange(1208.75,1751.25+0.625,0.625)
    wnSw = np.arange(2153.75,2551.25+0.625,0.625)
    lwData = np.asarray( h5['All_Data/CrIS-FS-SDR_All/ES_RealLW'])
    mwData = np.asarray( h5['All_Data/CrIS-FS-SDR_All/ES_RealMW'])
    swData = np.asarray( h5['All_Data/CrIS-FS-SDR_All/ES_RealSW'])
    lwShapeOut = lwData.shape
    mwShapeOut = mwData.shape
    swShapeOut = swData.shape
    lwShapeOut = lwShapeOut[0],lwShapeOut[1],lwShapeOut[2], lwShapeOut[3]-2
    mwShapeOut = mwShapeOut[0],mwShapeOut[1],mwShapeOut[2], mwShapeOut[3]-2
    swShapeOut = swShapeOut[0],swShapeOut[1],swShapeOut[2], swShapeOut[3]-2
    lwOut = np.zeros(lwShapeOut)
    mwOut = np.zeros(mwShapeOut)
    swOut = np.zeros(swShapeOut)

    for i in range(lwShapeOut[0]):
        for ii in range(lwShapeOut[1]):
            for iii in range(lwShapeOut[2]):
                lwOut[i,ii,iii,:] = smooth(lwData[i,ii,iii,:])
     

    for i in range(mwShapeOut[0]):
        for ii in range(mwShapeOut[1]):
            for iii in range(mwShapeOut[2]):
                mwOut[i,ii,iii,:] = smooth(mwData[i,ii,iii,:])

    for i in range(swShapeOut[0]):
        for ii in range(swShapeOut[1]):
            for iii in range(swShapeOut[2]):
                swOut[i,ii,iii,:] = smooth(swData[i,ii,iii,:])

    LW = PlanckFreqInv(wnLw[1:-1], lwOut )
    MW = PlanckFreqInv(wnMw[1:-1], mwOut )
    SW = PlanckFreqInv(wnSw[1:-1], swOut )
  
    return LW, MW, SW

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = ".")
    parser.add_argument('--cntl', help = "input path where h5 are stored.", required = True, dest = 'cntl' )
    parser.add_argument('--exp',help = " where you want to put the output files.", required = True, dest= 'exp')
    a = parser.parse_args()
    go (a.cntl, a.exp)  
