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
            sizaxis1 = scihdr['NAXIS1']
            sizaxis2 = scihdr['NAXIS2']
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


def make_counts_image(image, outfile='default'):
    '''FUNCTION TO CONVERT CNTS/SEC IMAGE TO COUNTS (IF NECESSARY)'''

    # -- parse output filename & save a copy to file (NOTE: if outfile == input image, data is overwritten).
    if (image != outfile):
        if outfile == 'default': outfile = image.split('.fits')[0] + '_cnts.fits'
        shutil.copy(image,outfile)
    else: print 'OVERWRITING DATA FOR IMAGE: '+image+'.'


    # -- determine if image is flt/flc, crclean, or drz/drc
    prihdr = pyfits.getheader(outfile,ext=0)
    pscale = prihdr.get('D001SCAL',default='NA')
    if pscale != 'NA': imtype = 'drz'
    elif len(image.split('crclean.fits')) > 1: imtype = 'crclean'
    else: imtype = 'flt'


    # -- initialize a few required parameters
    detector = prihdr['DETECTOR']
    exptime = prihdr['EXPTIME']


    # -- multiply by exposure time (only if image is already in cnts/sec)
    #      [notes] -- IR crcleans are actually in cnts, but "BUNIT" still says per second (can't trust bunit for now).
    #              -- We assume drz are cnts/sec, but this does not have to be true.
    if imtype == 'drz':
        # -- save background & pixel scale info
        if prihdr['extend'] == True: back = pyfits.getval(outfile,'MDRIZSKY',ext=1)  
        else: back = pyfits.getval(outfile,'MDRIZSKY',ext=0)
        pscale_nat = pyfits.getval(outfile,'D001ISCL',ext=0)
        pscale_img = pyfits.getval(outfile,'D001SCAL',ext=0)

        # -- assign the number of chips associated with this image
        if (prihdr['detector'] == 'IR'): nchips = 1.0                                          # IR
        elif (prihdr['subarray'] == True) and (len(prihdr['CCDAMP']) == 1): nchips = 1.0      # UVIS sub-array
        elif (prihdr['detector'] == 'UVIS') and (prihdr['subarray'] == False): nchips = 2.0   # UVIS full-frame
        else: raise exception('Image type is not defined.')

        # -- add background and correct for different pixel scale (original backgrd is measured in raw images)
        fdata = pyfits.getdata(outfile,ext=0)
        fdata_cnts = np.copy(fdata) * exptime + np.sum(back)/nchips * (pscale_img/pscale_nat)**2
        hdulist = pyfits.open(outfile,mode='update')
        hdulist[0].data = fdata_cnts
        hdulist.close()

    elif ((detector == 'IR') & (imtype == 'flt')):
        hdulist = pyfits.open(outfile,mode='update')
        for ff in xrange(len(hdulist)):
            if hdulist[ff].name == 'SCI': hdulist[ff].data = hdulist[ff].data * exptime
        hdulist.close()

    else: print 'IMAGE SHOULD ALREADY BE IN UNITS OF COUNTS. RETURNING...'

    return outfile
    

def run_daofind(image, outfile='default', dthreshold=3.0, fwhmpsf=2.5, backsigma=None,rdnoise=None):
    '''RUN DAOFIND ON INPUT IMAGE'''

    # -- parse output filename
    if outfile == 'default': outfile = image+'.coo'


    # -- extract header info 
    prihdr = pyfits.getheader(image)
    exptime = prihdr['exptime']
    instrum = prihdr['INSTRUME']
    detector = prihdr['DETECTOR']
    SUBARRAY = prihdr['SUBARRAY']
    ccdamp = prihdr['CCDAMP']


    # -- record filter name, native pixel scale, and no. of chips
    if instrum == 'WFC3':
        if detector == 'UVIS':
            pscale_nat = 0.03962
            if ((SUBARRAY == True) & (len(ccdamp) == 1)): nchips = 1.0
            elif SUBARRAY == False: nchips = 2.0
            else: raise Exception('Image type is not defined.')
        elif detector == 'IR':
            pscale_nat = 0.12825
            nchips = 1.0
        else: raise Exception('Detector '+detector+' not covered in our case list.')
    elif instrum == 'ACS':
        if detector == 'WFC':
            pscale_nat = 0.049
            if ((SUBARRAY == True) & (len(ccdamp) == 1)): nchips = 1.0
            elif SUBARRAY == False: nchips = 2.0
            else: raise Exception('Image type is not defined.')
        else: raise Exception('Detector '+detector+' not covered in our case list.')
    else: raise Exception('Instrument '+instrum+' not covered in our case list.')


    # -- record pixel scale of current image, image type, and number of flts
    sciext = []
    pscale_img = prihdr.get('D001SCAL',default='NA')
    if pscale_img == 'NA':
        imtype = 'flt'            # we dont distinguish between flt/crclean, i.e., assume pscales are equal
        pscale_img = pscale_nat
        num_flts = 1.0
        # -- record location of science extension
        hdulist = pyfits.open(image)
        for ext in xrange(len(hdulist)):
            if hdulist[ext].name == 'SCI': sciext.append(ext)
        hdulist.close()
        if len(sciext) != 1: raise Exception('We do not handle images with '+str(len(sciext))+' SCI extensions.')
    else:
        imtype ='drz'
        num_flts = prihdr['NDRIZIM']/nchips
        sciext.append(0)


    # -- estimate read noise
    if rdnoise == None:
        rdnoise = np.zeros(len(ccdamp))
        for namp in xrange(len(ccdamp)): rdnoise[namp] = prihdr['READNSE'+ccdamp[namp]]
    rdnoise_corr = np.sqrt(num_flts * (np.average(rdnoise) * pixscale_img/pixscale_nat)**2)


    # -- perform rough background noise calculation
    if backsigma == None:
        backstats=iraf.imstatistics(image+'['+str(sciext[0])+']', fields='stddev', lower = -100, upper = 100, nclip=5, \
                                    lsigma=3.0, usigma=3.0, cache='yes', format='no',Stdout=1)
        backsigma=float(backstats[0])


    # -- remove old daofind files/run daofind
    file_query = os.access(outfile, os.R_OK)        
    if file_query == True: os.remove(outfile)
    iraf.daofind.unlearn()
    iraf.daofind(image=image+'['+str(sciext[0])+']', interactive='no', verify='no',output=outfile, fwhmpsf=fwhmpsf, \
                 sigma=backsigma, readnoise=rdnoise_corr, itime=exptime, threshold=dthreshold, datamin=-10, datamax=100000)


    return outfile
    

def run_daophot(image, outfile='default', coordfile='NA', backmethod='mean', backval=None, backsigma=None,rdnoise=None,\
                apertures='1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,65,70', cbox=3.0, \
                annulus=17.0, dannulus=3.0, calgorithm='centroid', salgorithm='median', fwhmpsf=2.5, epadu=1.0):

    '''THIS PROCEDURE RUNS DAOPHOT ON INPUT IMAGE'''

    # Parse input parameters
    if outfile == 'default': outfile = image + '0.mag'
    if coordfile == 'NA': coordfile = image + '0.coo'

    # -- extract header info
    prihdr = pyfits.getheader(image)
    exptime = prihdr['exptime']
    instrum = prihdr['INSTRUME']
    detector = prihdr['DETECTOR']
    SUBARRAY = prihdr['SUBARRAY']
    ccdamp = prihdr['CCDAMP']

    if instrum == 'WFC3': filter = prihdr['FILTER']
    elif instrum == 'ACS':
        filter = prihdr['FILTER1']
        if filter[0] == 'C': filter == prihdr['FILTER2']
    else: raise Exception('Instrument '+instrum+' not covered in our case list.')

    # -- record native pixel scale and no. of chips
    if instrum == 'WFC3':
        if detector == 'UVIS':
            pscale_nat = 0.03962
            if ((SUBARRAY == True) & (len(ccdamp) == 1)): nchips = 1.0
            elif SUBARRAY == False: nchips = 2.0
            else: raise Exception('Image type is not defined.')
        elif detector == 'IR':
            pscale_nat = 0.12825
            nchips = 1.0
        else: raise Exception('WFC3 Detector '+detector+' not covered in our case list.')
    elif instrum == 'ACS':
        if detector == 'WFC':
            pscale_nat = 0.049
            if ((SUBARRAY == True) & (len(ccdamp) == 1)): nchips = 1.0
            elif SUBARRAY == False: nchips = 2.0
            else: raise Exception('Image type is not defined.')
        else: raise Exception('ACS Detector '+detector+' not covered in our case list.')
    else: raise Exception('Instrument '+instr+' not covered in our case list.')

    # -- record pixel scale of current image, image type, image axis lengths, and number of flts
    sciext = []
    pscale_img = prihdr.get('D001SCAL',default='NA')
    if pscale_img == 'NA':
        imtype = 'flt'            # we dont distinguish between flt/crclean, i.e., assume pscales are equal
        pscale_img = pscale_nat
        num_flts = 1.0
        naxis1 = pyfits.getval(image,'NAXIS1',ext=('SCI',1))
        naxis2 = pyfits.getval(image,'NAXIS2',ext=('SCI',1))
        # -- record location of science extension
        hdulist = pyfits.open(image)
        for ext in xrange(len(hdulist)):
            if hdulist[ext].name == 'SCI': sciext.append(ext)
        hdulist.close()
        if len(sciext) != 1: raise Exception('We do not handle images with '+str(len(sciext))+' SCI extensions.')
    else:
        imtype ='drz'
        num_flts = prihdr['NDRIZIM']/nchips
        naxis1 = pyfits.getval(image,'NAXIS1',ext=0)
        naxis2 = pyfits.getval(image,'NAXIS2',ext=0)
        sciext.append(0)

    # -- get zeropoints
    if instrum == 'WFC3': zeropt = get_wfc3_zeropoint(filter)
    elif instrum == 'ACS' and imtype == 'drz': zeropt = get_acs_zeropoint(prihdr)
    elif instrum == 'ACS' and imtype == 'flt': zeropt = get_acs_zeropoint(pyfits.getheader(image,ext=('SCI',1)))

    # -- estimate read noise
    if rdnoise == None:
        rdnoise = np.zeros(len(ccdamp))
        for namp in xrange(len(ccdamp)): rdnoise[namp] = prihdr['READNSE'+ccdamp[namp]]
    rdnoise_corr = np.sqrt(num_flts * (np.average(rdnoise) * pscale_img/pscale_nat)**2)

    # -- measure the background and noise
    if ((backval == None) | (backsigma == None)):
        # -- read in the x/y center of the source
        xc,yc = np.loadtxt(coordfile, unpack=True, usecols = (0,1))

        # -- create temporary image for bckgrd measurement that masks sources out to 80 pixels (assign a very low number)
        tmp_image = image+'.back.fits'
        shutil.copy(image, tmp_image)
        hdulist = pyfits.open(tmp_image, mode='update')
        maskim = hdulist[sciext[0]].data
        if detector == 'IR': maskrad = 30
        else: maskrad = 80
        maskim[circular_mask(maskim.shape, maskrad, x_offset=(xc-naxis1/2.0), y_offset=(yc-naxis2/2.0))] = -99999.0

        # -- Also mask out sources with zero effective exposure [WE ELIMINATE PIXELS WITHIN 20 OF IMAGE BORDER]
        maskim[:,0:20] = -99999.0
        maskim[:,-20:] = -99999.0
        maskim[0:20,:] = -99999.0
        maskim[-20:,:] = -99999.0

        # -- generate initial guess for lower/upper limits (use 10 sigma)
        fmaskim = np.ndarray.flatten(maskim)
        llim = -100
        ulim = 10000.0
        init_median,init_rms = meanclip(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],maxiter=7,return_median=1)
        llim = init_median - 10.0*init_rms
        ulim = init_median + 10.0*init_rms

        # -- measure background and rms
        if backmethod.lower() == 'mean':
            back,backrms=meanclip(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],maxiter=7)
        elif backmethod.lower() == 'median':
            back,backrms = meanclip(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],maxiter=7,return_median=1)
        elif backmethod.lower() == 'mode':
            backmean,backrms = meanclip(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],maxiter=7)
            nbins = np.ceil(80.0/(0.1*backrms))
            cc,bb,pp = pylab.hist(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],log=True,bins=nbins,range=(-40.0,40.0))
            back = bb[cc.argmax()] + (bb.max()-bb.min())/(2.0*(len(bb)-1))
        else: raise Exception('Background statistical method '+backmethod+' is not covered in our case list.')

        if backval == None: backval = back
        if backsigma == None: backsigma = backrms

        print '\n BACKGROUND =  '+str(backval)
        print ' BACKGROUND RMS =  '+str(backsigma)+' \n'

    # Case of no aperture size given (we select aperture sizes of UVIS=0.2/0.4", IR=0.27/0.4", ACS=0.25/0.5")
    if apertures == '':
        if instrum == 'WFC3' and detector == 'IR': apertures=str(0.27/pscale_img)+','+str(0.4/pscale_img)
        elif instrum == 'WFC3' and detector == 'UVIS': apertures=str(0.2/pscale_img)+','+str(0.4/pscale_img)
        elif instrum == 'ACS' and detector == 'WFC': apertures=str(0.25/pscale_img)+','+str(0.5/pscale_img)
        else: raise exception('Instrument/Detector '+instrum+'/'+detector+' not covered in case list.')

    # Remove old phot output files
    file_query = os.access(outfile, os.R_OK)      
    if file_query == True: os.remove(outfile)

    # Run phot
    iraf.phot.unlearn()         # reset daophot parameters to default values
    iraf.phot(image=image+'['+str(sciext[0])+']', interactive='no', verify='no', coords=coordfile, output=outfile, \
              fwhmpsf=fwhmpsf, sigma=backsigma, readnoise=rdnoise_corr, itime=exptime, calgorithm=calgorithm, \
              cbox=cbox, skyvalue=backval,apertures=apertures,zmag=zeropt,salgorithm='constant')
              #annulus=annulus, dannulus=dannulus


    return backval,backsigma    # return computed background stats for image



