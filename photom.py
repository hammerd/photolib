'''
ABOUT:
This module contains generic photometry routines such as applying PAM, source detection with daofind,
photometry with daophot, and other basic utilities that should work for many detectors (but tested primarily for WFC3)
and object types.

DEPENDS:
Python 2.5.4
stwcs 1.1.1

AUTHOR:
D. HAMMER 2013

HISTORY:
Nov. 2013: Original library.


FUTURE IMPROVEMENTS:


USE:
import photom.py

'''

__author__='D.M. HAMMER'
__version__= 0.1


import pyraf, os, glob, argparse, pdb, pyfits, pylab, fileinput, shutil, scipy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import interpolate
import numpy as np
from pyraf import iraf
from iraf import noao, digiphot, daophot
from astropy.table import Table
from astropy.io import ascii
from stwcs.wcsutil import HSTWCS


def make_PAMcorr_image(image, outfile='default'):
    pamdir = '/Users/hammer/Research/STScI/WFC3_TEAM/photometry/Monitoring/'

    # -- Parse output filename & save a copy to file
    if outfile == 'default': outfile = image.split('.fits')[0] + '_PAM.fits'
    shutil.copy(image,outfile)

    # -- Read in fits image and header (assume flt/flc - should handle both full- & sub-arrays)
    prihdr = pyfits.getheader(outfile)
    fdata = pyfits.getdata(outfile, ext=1)
    exptime = prihdr['EXPTIME']
    detector = prihdr['detector']

    # -- Cycle through each SCI extension
    hdulist = pyfits.open(image,mode='update')
    for ff in xrange(len(hdulist)):
        if hdulist[ff].name == 'SCI':

            # -- read in header and data info
            scihdr = hdulist[ff].header
            data = hdulist[ff].data
            chip = scihdr['CCDCHIP']
            sizaxis1 = scihdr['SIZAXIS1']
            sizaxis2 = scihdr['SIZAXIS2']
            x0 = np.abs(scihdr['LTV1'])
            y0 = np.abs(scihdr['LTV2'])
            x1 = x0 + sizaxis1
            y1 = y0 + sizaxis2

            # -- apply the PAM
            if detector == 'UVIS':
                if chip == 1:
                    pam=pyfits.getdata(pamdir+'UVIS1wfc3_map.fits')
                    hdulist[ff].data = data * pam[y0:y1,x0:x1]
                elif chip == 2:
                    pam=pyfits.getdata(pamdir+'UVIS2wfc3_map.fits')
                    hdulist[ff].data = data * pam[y0:y1,x0:x1]
                else: raise Exception('Chip case not handled.')
            elif detector == 'IR':
                pam=pyfits.getdata(pamdir+'ir_wfc3_map.fits')
                hdulist[ff].data = data * pam[y0:y1,x0:x1]
            else: raise Exception('Detector '+detector+' not covered in our case list.')
            hdulist.close()

    return outfile
