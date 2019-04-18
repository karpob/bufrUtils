#!/usr/bin/env python3
import glob, argparse, os


def go ( a ):
    """
    Main to this script: where most of the happens.
    Input --> a arguments for argparse
    Output --> potentially lots of files, but fewer than came in.
    """


    # all source bufr files in the directory 
    allBufrFiles = glob.glob( os.path.join(a.inpath,"*.bufr_d") )
    allBufrFiles.sort()
    repoPath = os.path.dirname( os.path.realpath(__file__) ) 
    h = open( os.path.join(repoPath, 'CrIS_BUFR_Table.hdr') ,'rb')
    head = h.read()
    h.close()

    for f in allBufrFiles:
        print( "Working on: {}".format(f) )
        # open the file and read in data.
        fin = open(f, 'rb')
        din = fin.read()
        fin.close()
        # open it again, but overwrite. Write header, then data.
        fout = open(f,'wb')
        fout.write( head )
        fout.write( din ) 
        fout.close()


if __name__ == "__main__":
    cwd = os.getcwd()
    parser = argparse.ArgumentParser( description = "Tack on header to CrIS 399 files")
    parser.add_argument('--input', help = "input path where tiny BUFR files are stored.", required = False, dest = 'inpath',default = cwd )

    a = parser.parse_args()
    go( a )
