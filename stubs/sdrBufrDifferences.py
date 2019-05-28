#!/usr/bin/env python3
import h5py
import numpy as np
from matplotlib import pyplot as plt
import argparse, glob, os
import traceback
import sys
from eccodes import *
import numpy as np

def bufr_decode(input_file):
    f = open(input_file, 'rb')
    radAll = []
    chanAll= []
    m = 0
    while True:
        chan = []
        rad = []
        ibufr = codes_bufr_new_from_file(f)

        if(ibufr is None): break

        codes_set(ibufr, 'unpack', 1)
        iVals = codes_get_array(ibufr, 'delayedDescriptorReplicationFactor')
        iVals = codes_get_array(ibufr, 'extendedDelayedDescriptorReplicationFactor')
        iVal = codes_get(ibufr, 'edition')
        iVal = codes_get(ibufr, 'masterTableNumber')
        iVal = codes_get(ibufr, 'bufrHeaderCentre')
        iVal = codes_get(ibufr, 'bufrHeaderSubCentre')
        iVal = codes_get(ibufr, 'updateSequenceNumber')
        iVal = codes_get(ibufr, 'dataCategory')
        iVal = codes_get(ibufr, 'internationalDataSubCategory')
        iVal = codes_get(ibufr, 'dataSubCategory')
        iVal = codes_get(ibufr, 'masterTablesVersionNumber')
        iVal = codes_get(ibufr, 'localTablesVersionNumber')
        iVal = codes_get(ibufr, 'typicalYear')
        iVal = codes_get(ibufr, 'typicalMonth')
        iVal = codes_get(ibufr, 'typicalDay')
        iVal = codes_get(ibufr, 'typicalHour')
        iVal = codes_get(ibufr, 'typicalMinute')
        iVal = codes_get(ibufr, 'typicalSecond')
        iVal = codes_get(ibufr, 'numberOfSubsets')
        iVal = codes_get(ibufr, 'observedData')
        iVal = codes_get(ibufr, 'compressedData')
        iValues = codes_get_array(ibufr, 'unexpandedDescriptors')
        iVal = codes_get(ibufr, 'satelliteIdentifier')
        iVal = codes_get(ibufr, 'centre')
        iVal = codes_get(ibufr, '#1#satelliteInstruments')
        iVal = codes_get(ibufr, 'satelliteClassification')
        iVal = codes_get(ibufr, 'year')
        iVal = codes_get(ibufr, 'month')
        iVal = codes_get(ibufr, 'day')
        iVal = codes_get(ibufr, 'hour')
#        iVal = codes_get(ibufr, 'minute')
        dVals = codes_get_array(ibufr, 'second')
        dVals = codes_get_array(ibufr, 'latitude')
        dVals = codes_get_array(ibufr, 'longitude')
        nscans = np.asarray(dVals).shape[0]
        dVals = codes_get_array(ibufr, 'satelliteZenithAngle')
        dVals = codes_get_array(ibufr, 'bearingOrAzimuth')
        dVals = codes_get_array(ibufr, 'solarZenithAngle')
        dVals = codes_get_array(ibufr, 'solarAzimuth')
        iVal = codes_get(ibufr, 'orbitQualifier')
        iValues = codes_get_array(ibufr, 'fieldOfRegardNumber')
        iValues = codes_get_array(ibufr, 'fieldOfViewNumber')
        iVal = codes_get(ibufr, 'orbitNumber')
        iValues = codes_get_array(ibufr, 'heightOfLandSurface')
        dVals = codes_get_array(ibufr, 'height')
        dVals = codes_get_array(ibufr, 'landFraction')
        iValues = codes_get_array(ibufr, 'landOrSeaQualifier')
        iVal = codes_get(ibufr, 'radianceTypeFlags')
        iVal = codes_get(ibufr, 'scanLevelQualityFlags')
        iVal = codes_get(ibufr, '#1#band')
        dVal = codes_get(ibufr, '#1#waveNumber')
        dVal = codes_get(ibufr, '#2#waveNumber')
        iVal = codes_get(ibufr, '#1#startChannel')
        iVal = codes_get(ibufr, '#1#endChannel')
        iVal = codes_get(ibufr, '#1#calibrationQualityFlags')
        iVal = codes_get(ibufr, '#1#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, '#2#band')
        dVal = codes_get(ibufr, '#3#waveNumber')
        dVal = codes_get(ibufr, '#4#waveNumber')
        iVal = codes_get(ibufr, '#2#startChannel')
        iVal = codes_get(ibufr, '#2#endChannel')
        iVal = codes_get(ibufr, '#2#calibrationQualityFlags')
        iVal = codes_get(ibufr, '#2#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, '#3#band')
        dVal = codes_get(ibufr, '#5#waveNumber')
        dVal = codes_get(ibufr, '#6#waveNumber')
        iVal = codes_get(ibufr, '#3#startChannel')
        iVal = codes_get(ibufr, '#3#endChannel')
        iVal = codes_get(ibufr, '#3#calibrationQualityFlags')
        iVal = codes_get(ibufr, '#3#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, 'geolocationQuality')
        iVal = codes_get(ibufr, 'qualityInformation')
        rad = np.zeros([2223,nscans])
        for ii,ch in enumerate(range(1,2224)):
            rad[ii,:] = np.asarray(codes_get_array(ibufr, '#{}#channelRadiance'.format(ch)))  
        codes_release(ibufr)
        radAll.append(rad)
        m+=1

    f.close()
    msgCnt = len(radAll)
    obCountList = []
    obCounts = 0
    for i in range(msgCnt):
        obCountList.append(radAll[i].shape[1])
        obCounts += radAll[i].shape[1]
    outRad = np.zeros([2223,obCounts])
    pos = 0 
    for i in range(msgCnt):
        outRad[:,pos:pos+obCountList[i]] = radAll[i][:]
        pos += obCountList[i]
    return outRad

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
    T=(C2_inv_cm*wn)/(np.log(1.+(np.power(wn,3)*C1_inv_cm)/B))
    return T

def go(cntlPath, expPath):
    cntList = glob.glob(os.path.join(cntlPath,'*.bufr'))
    expList = glob.glob(os.path.join(expPath, '*.bufr'))
    lwAve = []
    mwAve = []
    swAve = []
    cntList.sort()
    expList.sort() 
    cntTotal = 0
    for i,f in enumerate(cntList):
        print(i, len(cntList), f)
        #if(os.path.basename(f)[0:46] != os.path.basename(expList[i])[0:46]): continue 
        cLw,cMw,cSw = readSdrBufr(f)
        eLw,eMw,eSw = readSdrBufr( expList[i] )
        lwAveTmp = eLw - cLw
        mwAveTmp = eMw - cMw
        swAveTmp = eSw - cSw
        cnt = swAveTmp.shape[1] 
        print (lwAveTmp.shape, mwAveTmp.shape, swAveTmp.shape)
        if (len(lwAve) == 0) and (len(mwAve) == 0.0) and (len(swAve) == 0 ):
            lwAve = lwAveTmp.mean(axis=1)
            mwAve = mwAveTmp.mean(axis=1)
            swAve = swAveTmp.mean(axis=1)
        else:
            lwAve = (cntTotal*lwAve+cnt*lwAveTmp.mean(axis=1))/(cnt+cntTotal)
            mwAve = (cntTotal*mwAve+cnt*mwAveTmp.mean(axis=1))/(cnt+cntTotal)
            swAve = (cntTotal*swAve+cnt*swAveTmp.mean(axis=1))/(cnt+cntTotal)

        cntTotal+= cnt
        if(i == 5): print(lwAveTmp,mwAveTmp,swAveTmp)

    wnLw = np.arange(648.75,1096.75,0.625)
    wnMw = np.arange(1208.75,1751.25+0.625,0.625)
    wnSw = np.arange(2153.75,2551.25+0.625,0.625)
    plt.figure()        
    plt.plot(wnLw[1:-1], lwAve)
    plt.savefig('lw_bufr.png')
    plt.figure()
    plt.plot(wnSw[1:-1], swAve)
    plt.savefig('sw_bufr.png')
    plt.figure()
    plt.plot(wnMw[1:-1], mwAve)
    plt.savefig('mw_bufr.png')
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
def readSdrBufr ( fname ):
    rad = bufr_decode( fname )
    lw = rad[0:717,:]
    mw = rad[717:717+869,:]
    sw = rad[717+869:717+869+637,:]
    lwOut = np.zeros([lw.shape[0]-2,lw.shape[1]])
    mwOut = np.zeros([mw.shape[0]-2,lw.shape[1]])
    swOut = np.zeros([sw.shape[0]-2,lw.shape[1]])

    for i in range(0,lw.shape[1]):
        lwOut[:,i] = 1000.0*smooth(lw[:,i])
    for i in range(0,mw.shape[1]):
        mwOut[:,i] = 1000.0*smooth(mw[:,i])
    for i in range(0,sw.shape[1]):
        swOut[:,i] = 1000.0*smooth(sw[:,i])

    wnLw = np.arange(648.75,1096.75,0.625)
    wnMw = np.arange(1208.75,1751.25+0.625,0.625)
    wnSw = np.arange(2153.75,2551.25+0.625,0.625)
    LW = np.zeros(lwOut.shape)
    MW = np.zeros(mwOut.shape)
    SW = np.zeros(swOut.shape)
    for i in range(lwOut.shape[1]):
        LW[:,i] = PlanckFreqInv(wnLw[1:-1], lwOut[:,i] )
        MW[:,i] = PlanckFreqInv(wnMw[1:-1], mwOut[:,i] )
        SW[:,i] = PlanckFreqInv(wnSw[1:-1], swOut[:,i] )
    return LW, MW, SW 
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = ".")
    parser.add_argument('--cntl', help = "input path where h5 are stored.", required = True, dest = 'cntl' )
    parser.add_argument('--exp',help = " where you want to put the output files.", required = True, dest= 'exp')
    a = parser.parse_args()
    go (a.cntl, a.exp)  
