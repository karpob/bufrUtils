#!/usr/bin/env python3
from dateutil import rrule
from datetime import datetime, timedelta
import glob, argparse, os


def go ( a ):
    """
    Main to this script: where most of the happens.
    Input --> a arguments for argparse
    Output --> potentially lots of files, but fewer than came in.
    """

    #start period of our data
    now = datetime.strptime(a.start,'%Y%m%d%H')

    #end day of our data
    theEnd = datetime.strptime(a.end,'%Y%m%d%H')

    # calculate lengths. used later to point strptime to right point in filename
    pLen = len(a.prefix)
    bLen = len(a.between)

    # all source bufr files in the directory 
    allBufrFiles = glob.glob( os.path.join(a.inpath, a.prefix+"*.bufr") )
    allBufrFiles.sort()
    print(a.inpath) 

    # Make list of files to generate and calculate time windows.
    filesToMake = []
    startTimes = []
    endTimes = []
 
    for dt in rrule.rrule(rrule.HOURLY, interval = 6, dtstart = now, until = theEnd):
        #file to make (destination of cat for 6 hour window) what the GSI currently requires.
        filesToMake.append( dt.strftime('gdas1.%Y%m%d.t%Hz.'+a.instrument+'.tm00.bufr_d') )
        startTimes.append( dt - timedelta(hours=3) )
        endTimes.append( dt + timedelta(hours=3) )

    # 
    for i,fs in enumerate(filesToMake):
        print("Working on:{}".format( os.path.join(a.outpath,fs) ) )
        if( os.path.exists( os.path.join(a.outpath,fs) ) ):
            if( os.path.getsize( os.path.join(a.outpath,fs) ) != 0 ):
                print("skipping file! It exists already!")
                continue
 
        with open( os.path.join(a.outpath,fs), 'wb' ) as f:
            filesToDelete = []
            for b in allBufrFiles:
                bb = os.path.basename(b)
                isAfterStart = datetime.strptime(bb[pLen:pLen+12],'%Y%m%d%H%M') >= startTimes[i]
                isBeforeEnd = datetime.strptime(bb[pLen+12+bLen:pLen+12+bLen+12],'%Y%m%d%H%M') <= endTimes[i]   
                if isAfterStart and isBeforeEnd:
                    ff = open(b,'rb')
                    f.write( ff.read() )
                    ff.close()
                    filesToDelete.append( b )
            if(a.lucky):
                print("You have indicated you feel lucky. Good for you! I'm going to delete all the files I just concatenated... Still feel lucky?")
                for fdel in filesToDelete: os.remove(fdel)
                # since I just deleted these things, update the list.  
                allBufrFiles = glob.glob( os.path.join(a.inpath, a.prefix+"*.bufr") )
                allBufrFiles.sort()
        # if we just created an empty file, delete it. 
        if( os.path.getsize( os.path.join(a.outpath,fs) ) == 0 ): os.remove( os.path.join(a.outpath,fs) ) 

if __name__ == "__main__":
    cwd = os.getcwd()
    parser = argparse.ArgumentParser( description = "Cat together tiny BUFR files into 6hr chunks. --prefix needs to be exact  and --between doesn't need to be the exact text, just the lengths matter.")
    parser.add_argument('--start', help = 'start dtg YYYYMMDDhh', required = True, dest = 'start')
    parser.add_argument('--end', help = 'end dtg YYYYMMDDhh', required = True, dest = 'end')
    parser.add_argument('--prefix', help = "Everything in the input bufr file name \
                                            before the point where start time is tagged.",\
                        required = False, dest = 'prefix',default="NUCAPS-C0399_v2r0_j01_s")
    parser.add_argument('--instrument', help = "Instrument name to tag output file.", required = False, dest = 'instrument',default="cris")
    parser.add_argument('--input', help = "input path where tiny BUFR files are stored.", required = False, dest = 'inpath',default = cwd )
    parser.add_argument('--output', help = "output path where big BUFR files are stored.", required = False, dest = 'outpath',default = cwd )
    parser.add_argument('--between', help = "text between the start hour and the end year.", required = False, dest = 'between',default = '000_e' )
    parser.add_argument('--lucky', help="Do you feel lucky? Delete little files after doing cat.", dest='lucky', action='store_true' )

    a = parser.parse_args()
    go( a )
