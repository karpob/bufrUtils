#!/usr/bin/env python3
import glob, argparse, os, subprocess, sys


def go ( a ):
    """
    Main to this script: where most of the happens.
    Input --> a arguments for argparse
    """
    os.environ["AAPP_PREFIX"]='/discover/nobackup/projects/gmao/obsdev/bkarpowi/AAPP/AAPP_8.3/'
    os.environ["LD_LIBRARY_PATH"]='/discover/nobackup/bkarpowi/anaconda3/lib'
    aappBin = os.path.join(os.environ["AAPP_PREFIX"],'AAPP','bin')
    os.system('echo $PATH') 
    flist = glob.glob(os.path.join(a.inpath,'*.bufr'))
    for ifile in flist:
        p = subprocess.Popen(aappBin+'/eccodes_decodebufr_1c -i {} CRISFSR'.format(ifile), shell=True)
        p.wait() 

        ifile_l1c = os.path.join(a.inpath, os.path.basename(ifile).replace('bufr','l1c'))
        ofile_l1c = os.path.basename(ifile_l1c).replace('C2211','C431')
        ofile_l1c = os.path.join(a.outpath, ofile_l1c)
        p = subprocess.Popen(aappBin+'/cris_channels {} {}'.format(ifile_l1c,ofile_l1c), shell=True)
        p.wait()
        # remove now useless intermediate file.
        os.remove(ifile_l1c)
        # now take subset back to bufr.
        p = subprocess.Popen(aappBin+'/eccodes_encodebufr_1c -i {} CRISFSR'.format(ofile_l1c), shell=True)
        p.wait()
        # remove useless intermediate file
        os.remove(ofile_l1c) 
if __name__ == "__main__":
    cwd = os.getcwd()
    parser = argparse.ArgumentParser( description = ".")
    parser.add_argument('--input', help = "input path BUFR files are stored.", required = True, dest = 'inpath' )
    parser.add_argument('--output',help = " where you want to put the output files.", required = True, dest= 'outpath')
    a = parser.parse_args()
    go( a )
