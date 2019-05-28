#!/usr/bin/env python3
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse, glob, os
import traceback
import sys
from eccodes import *
import numpy as np
from datetime import datetime
from multiprocessing import Pool
from functools import partial
from maps import plotMapHist
#from Ipython import embed
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def isinFun(eTimes, cTimes):
    idxOut, = np.where( np.isin(np.asarray(eTimes), cTimes, invert=True) )
    return idxOut
 
def bufr_decode(input_file):
    f = open(input_file, 'rb')
    radAll = []
    timesAll = []
    latsAll = []
    lonsAll = []
    chanAll= []
    m = 0
    while True:
        chan = []
        rad = []
        try:ibufr = codes_bufr_new_from_file(f)
        except: continue
        if(ibufr is None): break

        codes_set(ibufr, 'unpack', 1)
        #iVals = codes_get_array(ibufr, 'delayedDescriptorReplicationFactor')
        #iVals = codes_get_array(ibufr, 'extendedDelayedDescriptorReplicationFactor')
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
        year = codes_get(ibufr, 'year')
        month = codes_get(ibufr, 'month')
        day = codes_get(ibufr, 'day')
        hour = codes_get(ibufr, 'hour')
        minutes = codes_get_array(ibufr, 'minute')
        seconds = codes_get_array(ibufr, 'second')
        lats = codes_get_array(ibufr, 'latitude')
        lons = codes_get_array(ibufr, 'longitude')
        nscans = np.asarray(lons).shape[0]
        dVals = codes_get_array(ibufr, 'satelliteZenithAngle')
        dVals = codes_get_array(ibufr, 'bearingOrAzimuth')
        dVals = codes_get_array(ibufr, 'solarZenithAngle')
        dVals = codes_get_array(ibufr, 'solarAzimuth')
        #iVal = codes_get(ibufr, 'orbitQualifier')
        forNumber = codes_get_array(ibufr, 'fieldOfRegardNumber')
        
        fovNumber = codes_get_array(ibufr, 'fieldOfViewNumber')
        orbit = codes_get(ibufr, 'orbitNumber')
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
        #iVal = codes_get(ibufr, '#1#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, '#2#band')
        dVal = codes_get(ibufr, '#3#waveNumber')
        dVal = codes_get(ibufr, '#4#waveNumber')
        iVal = codes_get(ibufr, '#2#startChannel')
        iVal = codes_get(ibufr, '#2#endChannel')
        iVal = codes_get(ibufr, '#2#calibrationQualityFlags')
        #iVal = codes_get(ibufr, '#2#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, '#3#band')
        dVal = codes_get(ibufr, '#5#waveNumber')
        dVal = codes_get(ibufr, '#6#waveNumber')
        iVal = codes_get(ibufr, '#3#startChannel')
        iVal = codes_get(ibufr, '#3#endChannel')
        iVal = codes_get(ibufr, '#3#calibrationQualityFlags')
        #iVal = codes_get(ibufr, '#3#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, 'geolocationQuality')
        iVal = codes_get(ibufr, 'qualityInformation')
        rad = np.zeros([399,nscans])
        for ii,ch in enumerate(range(1,400)):
            #print(ch)
            rad[ii,:] = np.asarray(codes_get_array(ibufr, '#{}#channelRadiance'.format(ch)))  
        codes_release(ibufr)
        dates = []
        for i,lat in enumerate(lats):
            if(len(minutes)>1): m = minutes[i]
            else: m = minutes[0]
            dates.append( int( "{:03d}".format(orbit)+ str(int( datetime(year ,month, day, hour, m, int(seconds[0]),0 ).timestamp() ))+"{:04d}".format ( int((lons[i]%360)*10000)) +"{:04d}".format (int((lat%180)*10000)) ) )   
        timesAll.extend(dates)
        latsAll.extend(lats)
        lonsAll.extend(lons)  
        radAll.append(rad)
        m+=1

    f.close()
    msgCnt = len(radAll)
    obCountList = []
    obCounts = 0
    for i in range(msgCnt):
        obCountList.append(radAll[i].shape[1])
        obCounts += radAll[i].shape[1]
    outRad = np.zeros([399,obCounts])
    pos = 0 
    for i in range(msgCnt):
        outRad[:,pos:pos+obCountList[i]] = radAll[i][:]
        pos += obCountList[i]
    return outRad, timesAll,latsAll, lonsAll

def PlanckFreqInv(wn,Radiance):
    """
    Convert to Brightness temperature. For a given radiance.
    Radiance -> Radiance in mW/(m2-sr-cm-1)
    wn->wavenumbers in cm^-1
    <-T in Kelvins
    """
    C2_inv_cm = 1.4387752
    C1_inv_cm = 1.191042e-5
    B=Radiance
    T=(C2_inv_cm*wn)/(np.log(1.+(np.power(wn,3)*C1_inv_cm)/B))
    return T

def go(cntlPath, expPath):
    p = Pool(28)
    cntList = glob.glob(os.path.join(cntlPath,'*.bufr_d'))
    expList = glob.glob(os.path.join(expPath, '*.bufr_d'))
    lwAve = []
    mwAve = []
    swAve = []
    cntList.sort()
    expList.sort() 
    cntTotalLw = 0
    cntTotalMw = 0
    cntTotalSw = 0
    h5 = h5py.File('cris_wavenumbers.h5','r')
    idx = np.asarray(h5['idxBufrSubset']) - 1
    wn = np.asarray(h5['wavenumbers'])[idx]
    idxAssimilated = np.asarray(h5['geosAssimilated'])-1
    wnAssim = np.asarray(h5['wavenumbers'])[idxAssimilated]
    wnLw = wn[np.where(wn<=1095)]
    wnMw = wn[np.where( (wn>1200) & (wn<=1750) )] 
    wnSw = wn[np.where( (wn>2100)  )] 
 
    for i,f in enumerate(cntList):
        print(i, len(cntList), f)
        #if(os.path.basename(f)[0:46] != os.path.basename(expList[i])[0:46]): continue 
        cLw,cMw,cSw,cTimes,cLat,cLon = readSdrBufr(f)
        eLw,eMw,eSw,eTimes,eLat,eLon = readSdrBufr( expList[i] )
        print('done reading {}'.format(f))
        print('etimes,ctimes',len(eTimes),len(cTimes))
        print('eLw shapes,cLw.shape', eLw.shape, cLw.shape)
        print('cLat, cLon shapes', cLat.shape, cLon.shape)
        print('eLat, eLon shapes', eLat.shape, eLon.shape)
        cntSw = eSw.shape[1]
        cntMw = eMw.shape[1]
        cntLw = eLw.shape[1]
        print('etimes,ctimes',len(eTimes),len(cTimes))
        """
        eTimesChunks = list( chunks(eTimes.tolist(),28) )
        chunkSizes = []
        for chunk in eTimesChunks:
            chunkSizes.append(len(chunk))
        idxOutOfControlList = p.map( partial(isinFun, cTimes=cTimes), eTimesChunks )
        idxOffset = 0
        idxOutOfControl = []
        print('idxOutOfControlList', idxOutOfControlList) 
        for i,idxChunk in enumerate(idxOutOfControlList):
            ii = np.asarray(idxChunk)
            if(len(idxChunk)>0):
                offsetChunk = ii + idxOffset
                idxOutOfControl.extend(offsetChunk.tolist())
            idxOffset += chunkSizes[i]
        print('idxOutOfControl:',idxOutOfControl)
        with h5py.File('idxOutOfControl.h5','w') as f:
            dset = f.create_dataset('idxOutOfControl',data = np.asarray(idxOutOfControl) )
        

        cTimesChunks = list( chunks(cTimes.tolist(),28) )
        chunkSizes = []
        for chunk in cTimesChunks:
            chunkSizes.append(len(chunk))

        idxOutOfExperimentList = p.map( partial(isinFun, cTimes=eTimes), cTimesChunks )
        idxOffset = 0
        idxOutOfExperiment = []
        for i,idxChunk in enumerate(idxOutOfExperimentList):
            ii = np.asarray(idxChunk)
            if(len(idxChunk)>0):
                offsetChunk = ii + idxOffset
                idxOutOfExperiment.extend(offsetChunk.tolist())
            idxOffset += chunkSizes[i]
        with h5py.File('idxOutOfExperiment.h5','w') as f:
            dset = f.create_dataset('idxOutOfExperiment',data = np.asarray(idxOutOfExperiment) )
        """
        h5 = h5py.File('idxOutOfExperiment.h5','r')
        idxOutOfExperiment = np.asarray(h5['idxOutOfExperiment'])
       
        h5 = h5py.File('idxOutOfControl.h5','r')
        idxOutOfControl = np.asarray(h5['idxOutOfControl'])
 
        eTimes = np.delete(eTimes, idxOutOfControl)
        cTimes = np.delete(cTimes, idxOutOfExperiment)
        eTimes = np.asarray(eTimes)
        cTimes = np.asarray(cTimes)
        idxExpSort = np.argsort(eTimes)
        idxCtlSort = np.argsort(cTimes)

        eLw = np.delete(eLw, idxOutOfControl, axis=1)
        cLw = np.delete(cLw, idxOutOfExperiment, axis=1)
        print('lw shapes',eLw.shape,cLw.shape)
        print('lon shapes',eLon.shape,cLon.shape)
        eLat = np.delete(np.asarray(eLat), idxOutOfControl)
        eLon = np.delete(np.asarray(eLon), idxOutOfControl)
        cLat = np.delete(np.asarray(cLat), idxOutOfExperiment)
        cLon = np.delete(np.asarray(cLon), idxOutOfExperiment)

        cLonO = np.asarray(cLon[idxCtlSort])
        eLonO = np.asarray(eLon[idxExpSort])

        cLatO = np.asarray(cLat[idxCtlSort])
        eLatO = np.asarray(eLat[idxExpSort])
        eLwO = np.asarray(eLw[:,idxExpSort])
        cLwO = np.asarray(cLw[:,idxCtlSort])
 
        print('lon shapes',eLon.shape,cLon.shape)
        for i in range(len(wnLw)):
            if wnLw[i] in wnAssim:
                plt.figure()
                #plt.title('Control Channel Wavenumber {}'.format(wnLw[i]))
                #plt.hist(cLw[i,:],bins='auto')
                #plt.ylabel('Number of Observations [count]')
                #plt.xlabel('Brightness Temperature [K]')
                #plt.savefig( 'histControl_{}.png'.format(wnLw[i]) )
                wtf = []
                wtf2 = []
                #plotMapHist(cLat, cLon,\
                #    cLw[i,:],'Brightness Temperature Wavenumber {}'.format(wnLw[i]), 'map_control{}'.format(wnLw[i]) )
                #plotMapHist(eLat, eLon,\
                #    eLw[i,:],'Brightness Temperature Wavenumber {}'.format(wnLw[i]), 'map_experiment{}'.format(wnLw[i]) )
                #plotMapHist(eLat[idxOutOfControl], eLon[idxOutOfControl],\
                #    eLw[i,idxOutOfControl],'Brightness Temperature Wavenumber {}'.format(wnLw[i]), 'map_out_of_control_{}.png'.format(wnLw[i]) )
                
                plotMapHist(cLatO, cLonO,\
                    eLonO-cLonO,'Difference Lon {}'.format(wnLw[i]), 'lon_difference{}'.format(wnLw[i]) )

                plotMapHist(cLatO, cLonO,\
                    eLatO-cLatO,'Difference Lat {}'.format(wnLw[i]), 'lat_difference{}'.format(wnLw[i]) )

                plotMapHist(cLatO, cLonO,\
                    eLwO[i,:]-cLwO[i,:],'Brightness Temperature Difference Wavenumber {}'.format(wnLw[i]), 'map_difference{}'.format(wnLw[i]) )
                #plt.title('Experiment Channel Wavenumber {}'.format(wnLw[i])+'cm${-1}$')
                #plt.hist(eLw[i,:],bins='auto')
                #plt.ylabel('Number of Observations [count]')
                #plt.xlabel('Brightness Temperature [K]')
                #plt.savefig( 'histExperiment_{}.png'.format(wnLw[i]) )
"""
        for i in range(len(wnMw)):
            if wnMw[i] in wnAssim:
                #plt.figure()
                #plt.title('Control Channel Wavenumber {}'.format(wnMw[i])+'cm${-1}$')
                #plt.hist(cMw[i,:],bins='auto')
                #plt.ylabel('Number of Observations [count]')
                #plt.xlabel('Brightness Temperature [K]')
                #plt.savefig( 'histControl_{}.png'.format(wnMw[i]) )

                plt.figure()
                plt.title('Experiment Channel Wavenumber {}'.format(wnMw[i])+'cm${-1}$')
                plt.hist(eMw[i,idxOutOfControl],bins='auto')
                plt.ylabel('Number of Observations [count]')
                plt.xlabel('Brightness Temperature [K]')
                plt.savefig( 'histExperiment_{}.png'.format(wnMw[i]) )

        for i in range(len(wnSw)):
            if wnSw[i] in wnAssim:
                #plt.figure()
                #plt.title('Control Channel Wavenumber {}'.format(wnSw[i])+'cm${-1}$')
                #plt.hist(cSw[i,:],bins='auto')
                #plt.xlabel('Number of Observations [count]')
                #plt.ylabel('Brightness Temperature [K]')
                #plt.savefig( 'histControl_{}.png'.format(wnSw[i]) )

                plt.figure()
                plt.title('Experiment Channel Wavenumber {}'.format(wnSw[i])+'cm${-1}$')
                plt.hist(eSw[i,idxOutOfControl],bins='auto')
                plt.xlabel('Number of Observations [count]')
                plt.ylabel('Brightness Temperature [K]')
                plt.savefig( 'histExperiment_{}.png'.format(wnSw[i]) )
 
 
        lwAveTmp = np.delete(eLw, idxOutOfControl, axis=1).mean(axis=1) - cLw.mean(axis=1)
        mwAveTmp = np.delete(eMw, idxOutOfControl, axis=1).mean(axis=1) - cMw.mean(axis=1)
        swAveTmp = np.delete(eSw, idxOutOfControl, axis=1).mean(axis=1) - cSw.mean(axis=1)
        print ('experiment,control',eLw.shape,cLw.shape)
        if (len(lwAve) == 0) and (len(mwAve) == 0.0) and (len(swAve) == 0 ):
            lwAve = lwAveTmp
            mwAve = mwAveTmp
            swAve = swAveTmp
        else:
            lwAve = (cntTotalLw*lwAve+cntLw*lwAveTmp)/(cntLw+cntTotalLw)
            mwAve = (cntTotalMw*mwAve+cntMw*mwAveTmp)/(cntMw+cntTotalMw)
            swAve = (cntTotalSw*swAve+cntMw*swAveTmp)/(cntSw+cntTotalMw)

        cntTotalLw += cntLw
        cntTotalSw += cntSw
        cntTotalMw += cntMw
        if(i == 10): break

 
    plt.figure()
    plt.title('Difference between Average Temperature Experiment - Control')        
    plt.plot(wnLw, lwAve)
    plt.ylabel('Brightness Temperature Difference [K]')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    plt.savefig('lw_bufr.png')

    plt.figure()
    plt.title('Difference between Average Temperature Experiment - Control')        
    plt.plot(wnSw, swAve)
    plt.ylabel('Brightness Temperature Difference [K]')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    plt.savefig('sw_bufr.png')

    plt.figure()
    plt.title('Difference between Average Temperature Experiment - Control')        
    plt.plot(wnMw, mwAve)
    plt.ylabel('Brightness Temperature Difference [K]')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    plt.savefig('mw_bufr.png')

    assimilatedWnLw = []
    assimilatedDiffLw = []
    for i,wn in enumerate(wnLw):
        if wn in wnAssim:
            assimilatedWnLw.append(wn)
            assimilatedDiffLw.append(lwAve[i])

    assimilatedWnMw = []
    assimilatedDiffMw = []
    for i,wn in enumerate(wnMw):
        if wn in wnAssim:
            assimilatedWnMw.append(wn)
            assimilatedDiffMw.append(mwAve[i])


    assimilatedWnSw = []
    assimilatedDiffSw = []
    for i,wn in enumerate(wnSw):
        if wn in wnAssim:
            assimilatedWnSw.append(wn)
            assimilatedDiffSw.append(mwAve[i])


    plt.figure()
    plt.title('Difference between Average Temperature Experiment - Control')        
    plt.plot(assimilatedWnLw, assimilatedDiffLw)
    plt.ylabel('Brightness Temperature Difference [K]')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    plt.savefig('lw_bufr_assim.png')

    plt.figure()
    plt.title('Difference between Average Temperature Experiment - Control')        
    plt.plot(assimilatedWnMw, assimilatedDiffMw)
    plt.ylabel('Brightness Temperature Difference [K]')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    plt.savefig('mw_bufr_assim.png')


    plt.figure()
    plt.title('Difference between Average Temperature Experiment - Control')        
    plt.plot(assimilatedWnSw, assimilatedDiffSw)
    plt.ylabel('Brightness Temperature Difference [K]')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    plt.savefig('sw_bufr_assim.png')

"""
        


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
    rad, times, lat, lon = bufr_decode( fname )
    
    h5 = h5py.File('cris_wavenumbers.h5','r')
    idx = np.asarray(h5['idxBufrSubset']) - 1
    wn = np.asarray(h5['wavenumbers'])[idx] 
    wnLw = wn[np.where(wn<=1095)]
    wnMw = wn[np.where( (wn>1200) & (wn<=1750) )] 
    wnSw = wn[np.where( wn>2100 )]
    
    lw = rad[np.where(wn<=1095)]
    mw = rad[np.where( (wn>1200) & (wn<1800))] 
    sw = rad[np.where( (wn>2100))]

    idxLw = np.where 
    LW = np.zeros(lw.shape)
    MW = np.zeros(mw.shape)
    SW = np.zeros(sw.shape)
    for i in range(lw.shape[1]):
        LW[:,i] = PlanckFreqInv(wnLw, lw[:,i]*1000.0 )
        MW[:,i] = PlanckFreqInv(wnMw, mw[:,i]*1000.0 )
        SW[:,i] = PlanckFreqInv(wnSw, sw[:,i]*1000.0 )
    return LW, MW, SW, np.asarray(times), np.asarray(lat), np.asarray(lon) 
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = ".")
    parser.add_argument('--cntl', help = "input path where h5 are stored.", required = True, dest = 'cntl' )
    parser.add_argument('--exp',help = " where you want to put the output files.", required = True, dest= 'exp')
    a = parser.parse_args()
    go (a.cntl, a.exp)  
