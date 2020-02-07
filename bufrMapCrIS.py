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
from maps import plotMapHist,plotMap
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
    satIdents = []
    m = 0
    msgCnt = 0 
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
        #iVal = codes_get(ibufr, 'internationalDataSubCategory')
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
        try: satIdent = codes_get(ibufr, 'satelliteIdentifier')
        except: continue
        if(satIdent != 225 and satIdent != 224): continue
        if (satIdent not in satIdents): satIdents.append(satIdent)

        iVal = codes_get(ibufr, 'centre')
        #iVal = codes_get(ibufr, '#1#satelliteInstruments')
        #iVal = codes_get(ibufr, 'satelliteClassification')
        year = codes_get(ibufr, 'year')
#        try: year = codes_get(ibufr, 'year')
#        except: break
        month = codes_get(ibufr, 'month')
        day = codes_get(ibufr, 'day')
        hour = codes_get_array(ibufr, 'hour')
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
        #iVal = codes_get(ibufr, '#1#calibrationQualityFlags')
        #iVal = codes_get(ibufr, '#1#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, '#2#band')
        dVal = codes_get(ibufr, '#3#waveNumber')
        dVal = codes_get(ibufr, '#4#waveNumber')
        iVal = codes_get(ibufr, '#2#startChannel')
        iVal = codes_get(ibufr, '#2#endChannel')
        #iVal = codes_get(ibufr, '#2#calibrationQualityFlags')
        #iVal = codes_get(ibufr, '#2#fieldOfViewQualityFlags')
        iVal = codes_get(ibufr, '#3#band')
        dVal = codes_get(ibufr, '#5#waveNumber')
        dVal = codes_get(ibufr, '#6#waveNumber')
        iVal = codes_get(ibufr, '#3#startChannel')
        iVal = codes_get(ibufr, '#3#endChannel')
        #iVal = codes_get(ibufr, '#3#calibrationQualityFlags')
        #iVal = codes_get(ibufr, '#3#fieldOfViewQualityFlags')
        #iVal = codes_get(ibufr, 'geolocationQuality')
        #iVal = codes_get(ibufr, 'qualityInformation')
        rad = []
        ch = 1
        while True:
            try: rad.append(np.asarray(codes_get_array(ibufr, '#{}#channelRadiance'.format(ch))))
            except: break
            ch += 1
        rad = np.asarray(rad)
        nchan = rad.shape[0]  
        codes_release(ibufr)
        dates = []
        for i,lat in enumerate(lats):
            if(len(minutes)>1): m = minutes[i]
            else: m = minutes[0]
            if(len(hour)>1): h = hour[i]
            else: h =hour[0]
            if(len(seconds)>1): s = int(seconds[i])
            else: s = int(seconds[0])
            # this line is admittedly bizarre, but without thinking too much, these seems easiest...
            #dateStr = "{}/{}/{} {}:{}:{:.4f}".format(month, day, year, h, m, s)
            #dates.append( int( datetime.strptime(dateStr,'%m/%d/%Y %H:%M:%S.%f').timestamp() ) )
   
            #dates.append( int( "{:03d}".format(orbit)+ str(int( datetime(year ,month, day, hour, m, int(seconds[0]),0 ).timestamp() ))+"{:04d}".format ( int((lons[i]%360)*10000)) +"{:04d}".format (int((lat%180)*10000)) ) )   
            dates.append( int( "{:03d}".format(orbit)+ str(int( datetime(year ,month, day, h, m, s,0 ).timestamp() ))+"{:04d}".format ( int((lons[i]%360)*10000)) +"{:04d}".format (int((lat%180)*10000)) ) )   
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
    outRad = np.zeros([nchan,obCounts])
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

def go(cntlPath):
    cntList = glob.glob(os.path.join(cntlPath,'*.bufr_d'))
    h5 = h5py.File('cris_wavenumbers.h5','r')
    idx = np.asarray(h5['idxBufrSubset']) - 1
    wn = np.asarray(h5['wavenumbers'])[idx]
    idxAssimilated = np.asarray(h5['geosAssimilated'])-1
    wnAssim = np.asarray(h5['wavenumbers'])[idxAssimilated]


    for i,f in enumerate(cntList):
        #if(os.path.basename(f)[0:46] != os.path.basename(expList[i])[0:46]): continue 
        cLw, cMw, cSw, cTimes, cLat, cLon, wnLw, wnMw, wnSw = readSdrBufr(f)
        for i in range(len(wnLw)):
            if wnLw[i] in wnAssim:
                plotMapHist(cLat, cLon,\
                    cLw[i,:],'Brightness Temperature Wavenumber {} '.format(wnLw[i])+'cm$^{-1}$', 'map_control_{0:08.3f}'.format(wnLw[i]) )
        for i in range(len(wnMw)):
            if wnMw[i] in wnAssim:
                plotMapHist(cLat, cLon,\
                    cMw[i,:],'Brightness Temperature Wavenumber {} '.format(wnMw[i])+'cm$^{-1}$', 'map_control_{0:08.3f}'.format(wnMw[i]) )
        for i in range(len(wnSw)):
            if wnSw[i] in wnAssim:
                plotMapHist(cLat, cLon,\
                    cSw[i,:],'Brightness Temperature Wavenumber {} '.format(wnSw[i])+'cm$^{-1}$', 'map_control_{0:08.3f}'.format(wnSw[i]) )
        


def readSdrBufr ( fname ):
    rad, times, lat, lon = bufr_decode( fname )
    if (rad.shape[0] == 1305):
        lw = rad[0:713,:]
        mw = rad[713:713+433,:]
        sw = rad[713+433:713+433+159,:]

        wnLw = np.arange(650,1095+0.625,0.625)
        wnMw = np.arange(1210,1750+1.25, 1.25)
        wnSw = np.arange(2155,2550+2.5, 2.5)
        print('wavenumbers',len(wnLw)+len(wnMw)+len(wnSw))
    else:
        print('warning, this is not NSR 1305. This is probably not going to work.') 

    LW = np.zeros(lw.shape)
    MW = np.zeros(mw.shape)
    SW = np.zeros(sw.shape)
    for i in range(lw.shape[1]):
        LW[:,i] = PlanckFreqInv(wnLw, lw[:,i]*1000.0 )
        MW[:,i] = PlanckFreqInv(wnMw, mw[:,i]*1000.0 )
        SW[:,i] = PlanckFreqInv(wnSw, sw[:,i]*1000.0 )
    return LW, MW, SW, np.asarray(times), np.asarray(lat), np.asarray(lon), wnLw, wnMw, wnSw 
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = ".")
    parser.add_argument('--cntl', help = "input path where h5 are stored.", required = True, dest = 'cntl' )
    a = parser.parse_args()
    go (a.cntl)  
