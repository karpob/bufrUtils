#!/usr/bin/env python3
import os,glob,argparse


def main(a):

    # store working directory so you can go back later.
    path2here = os.getcwd()

    # if outpath doesn't exits make it.
    if (not os.path.isdir(a.outpath) ):os.makedirs(a.outpath)

    # get nc files in ncpath directory
    ncfiles = glob.glob( os.path.join(a.ncpath,'*.nc') )

    # change the workpath to the specified workpath.
    os.chdir(a.workpath)

    #iterate over netcdf files, write input file for npr_main.
    for ncfile in ncfiles:
        with open('npr.filenames','w') as f:
            f.write( 'CRIS_CCR\n' )
            f.write( 'NC2BUFR\n' )
            f.write( ncfile+"\n" )
            f.write( os.path.join(a.outpath,os.path.split(ncfile)[-1].replace('.nc','.bufr'))+"\n" )
            f.write( os.path.join(a.tablepath,'CrIS_NSR_BUFR_Table')+"\n" )
            f.write( '29817\n' )
            f.write( '224\n' )
        os.system( a.mainpath ) # call the npr binary.
    # go back to the working directory you started in..
    os.chdir(path2here)
if __name__ =="__main__":
    parser = argparse.ArgumentParser( description = "Run NCEPs NC2BUFR (modified by Bryan Karpowicz to work with CrIS NSR).")
    parser.add_argument('--ncpath', help = 'input directory', required = True, dest = 'ncpath')
    parser.add_argument('--outpath', help = 'output directory', required = True, dest = 'outpath')
    parser.add_argument('--mainpath', help ='full path the npr_main', required = False, default='/discover/nobackup/projects/cris/bkarpowi/CrIS_CCR/code/main/main_npr')
    parser.add_argument('--workpath', help ='working directory path', required = False,dest = 'workpath', default=os.getcwd() ) 
    parser.add_argument('--tablepath', help ='working directory path', required = False,dest = 'tablepath', default=os.getcwd() ) 
    a = parser.parse_args()
    main(a)
    print('Done!')  
