#! /usr/bin/env python

'''
ABOUT:
This program performs point-source photometry using the IRAF tasks daofind and daoimage
on a single HST image or a list of images (handles flt,crclean,single_sci or drz).

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2012

HISTORY:
Sept. 2012: Original script (v0.1).
Oct. 2012: Added more robust handling of images with single chip (e.g., IR, calibration).

FUTURE IMPROVEMENTS:
PAM-correction for ACS?


USE:
python EE_stepping_program.py
'''

__author__='D.M. HAMMER'
__version__= 0.2


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
    #
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
        if backmethod.lowercase() == 'mean': back,backrms=meanclip(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],maxiter=7)
        elif backmethod.lowercase() == 'median': back,backrms = meanclip(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],maxiter=7,return_median=1)
        elif backmethod.lowercase() == 'mode':
            backmean,backrms = meanclip(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],maxiter=7)
            nbins = np.ceil(80.0/(0.1*backrms))
            cc,bb,pp = pylab.hist(fmaskim[(fmaskim > llim) & (fmaskim < ulim)],log=True,bins=nbins,range=(-40.0,40.0))
            back = bb[cc.argmax()] + (bb.max()-bb.min())/(2.0*(len(bb)-1))
        else: raise Exception('Background statistical method '+backmethod+' is not covered in our case list.')

        if backval == None: backval = back
        if backsigma == None: backsigma == backrms

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
	          cbox=cbox, skyvalue=backmean,apertures=apertures,zmag=zeropt,salgorithm='constant')
	          #annulus=annulus, dannulus=dannulus


	return backval,backsigma    # return computed background stats for image



def circular_mask(arr_shape, r, x_offset=0, y_offset=0):
    """
    Generate circular mask for 2D image.

    Parameters
    ----------
    arr_shape : tuple of int
        Shape of the array to use the mask.

    r : int
        Radius of the mask in pixels.

    x_offset, y_offset : int or float, optional
        Mask offset relative to image center.

    Returns
    -------
    Numpy indices of the mask, rounded to nearest
    integer.

    References
    ----------
    http://mail.scipy.org/pipermail/numpy-discussion/2011-January/054470.html

    """
    assert len(arr_shape) == 2, 'Image is not 2-D'

    ny, nx = arr_shape
    assert nx > 1 and ny > 1, 'Image is too small'

    assert isinstance(r, (int, long)) and r > 0, 'Radius must be int > 0'

    xcen = np.round(0.5 * nx - 0.5 + x_offset).astype('int')
    ycen = np.round(0.5 * ny - 0.5 + y_offset).astype('int')

    x1, x2 = xcen - r, xcen + r
    y1, y2 = ycen - r, ycen + r

    assert y1 >= 0 and y2 < ny and x1 >= 0 and x2 < nx, 'Mask falls outside image bounds'

    y, x = np.ogrid[-r:r, -r:r]
    i = np.where(x**2 + y**2 <= r**2)

    a = np.zeros(arr_shape).astype('bool')
    a[y1:y2, x1:x2][i] = True

    return np.where(a)



def replace_filevalue(file, orgval, newval):

    ''' REPLACE UNWANTED VALUES IN EXTERNAL FILE '''

    for line in fileinput.input(file, inplace = 1):
        print line.replace(str(orgval), str(newval)),
    fileinput.close()



def get_wfc3_zeropoint(filter):
    # array of WFC3 AB zeropoints (not all here - add as needed)
    zp = {'F225W':24.0403, 'F275W':24.1305, 'F336W':24.6682, 'F390W':25.3562, 'F438W':24.8206, 'F475W':25.6799, \
          'F547M':24.7467, 'F555W':25.7906, 'F606W':26.0691, 'F656N':20.4827, 'F658N':21.0287, 'F814W':25.0985, \
          'F850LP':23.8338,'F105W':26.2687, 'F125W':26.2303, 'F140W':26.4524, 'F160W':25.9463}

    if zp.has_key(filter.upper()): return zp[filter.upper()]
    else: raise Exception('Zeropoint is not specified for this filter: '+filter)


def get_acs_zeropoint(hdr):
    PHOTFLAM = hdr['PHOTFLAM']
    PHOTPLAM = hdr['PHOTPLAM']
    zeropt = -2.5*np.log10(PHOTFLAM) - 5.0*np.log10(PHOTPLAM) - 2.408  #AB
    return zeropt


def get_modelpsf_uvis(wave_eval, rad_eval=[-9999.0]):
	''' RETRIEVE MODEL PSF FROM G. HARTIG 2009 ISR (wave_eval is a scalar -- only 1 wavelength permitted)'''	

	# Checks
	if len(np.array([wave_eval])) != 1: raise Exception('PSF may be evaluated at only 1 wavelength.')

	# read UVIS EE vs radius from Hartig
	extfile = 'Hartig_EE_model.dat'
	alldata = np.loadtxt(extfile)
	wave = alldata[0,1:]
	aper_rad = alldata[1:,0]

	# Match input wavelength to table wavelength (MUST MATCH)
	gdw = np.where(wave_eval == wave)[0]
	if len(gdw) == 1: data = alldata[1:,(gdw[0]+1)]
	else: raise Exception('Unique wavelength not found in external table.')

	# set radius positions to be evaluated & perform boundary checks
        nrads = data.size
        rmin=aper_rad[0]
        rmax=aper_rad[nrads-1]
	if rad_eval[0] < -999.0: rad_eval = aper_rad       	#default is to use George's radius values
	if len(np.array(rad_eval)[rad_eval > rmax]) > 0 or len(np.array(rad_eval)[rad_eval < rmin]) > 0: raise Exception('Requested aperture radius is outside table boundaries.')

	# Evaluate table values at requested radii
	interp = scipy.interpolate.InterpolatedUnivariateSpline(aper_rad, data)   # establish class for spline fitted model
	if len(rad_eval) == 1: flux_eval = np.array([interp(rad_eval)])
	else: flux_eval = interp(rad_eval)

	# Check for wave input parameters outside exisiting boundaries or rules (OLD METHOD FOR ISR TABLE)
	#ny,nx=data.shape
        #xmin=wave[0]
        #xmax=wave[nx-1]
	#xeval = [wave_eval for x in xrange(len(rad_eval))]
	#if (wave_eval > xmax) or (wave_eval < xmin): raise Exception('Requested wavelength is outside table boundaries.')
	# Evaluate table at requested wave/radius by interpolating (APPLIES TO TABLE FROM ISR--NOT NEW TABLE)
	#interp = scipy.interpolate.RectBivariateSpline(aper_rad, wave, data)	# establish class for spline fitted model
	#flux_eval = interp.ev(rad_eval, xeval)

	return [rad_eval, flux_eval]



def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.02, verbose=1, return_array=0, return_median=0):
   """
   Computes an iteratively sigma-clipped mean on a
   data set. Clipping is done about median, but mean
   is returned by default (use return_median to return
   the clipped median value).

   .. note:: MYMEANCLIP routine from ACS library.

   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.
       * 08/01/2013 Added option to return the array indices of non-clipped pixels. DMH
		    Added option to return median of the clipped array. DMH

   Examples
   --------
   >>> mean, sigma = meanclip(indata)

   Parameters
   ----------
   indata: array_like
       Input data.

   clipsig: float
       Number of sigma at which to clip.

   maxiter: int
       Ceiling on number of clipping iterations.

   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.

   verbose: {0, 1}
       Print messages to screen?

   return_array: {0, 1}
        Return the final array indices that were used to compute statistics.

   Returns
   -------
   mean: float
       N-sigma clipped mean.

   sigma: float
       Standard deviation of remaining pixels.

   """
   # Flatten array
   skpix = indata.reshape( indata.size, )

   # initialize array to store indices of pixels used to compute stats
   arrind = np.arange(0,skpix.size)

   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0

   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = np.median(skpix)
       sig = np.std(skpix)
       wsm = np.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]
           arrind = arrind[wsm]

       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop

   if return_median:
        val = np.median(skpix)
        val_type = 'median'
   else:
        val  = np.mean( skpix )
        val_type = 'mean'
   sigma = np.std( skpix )

   if verbose:
       if return_median:
           prf = 'MEDIANCLIP:'
           print '%s %.1f-sigma clipped median' % (prf, clipsig)
           print '%s Median computed in %i iterations' % (prf, iter)
           print '%s Median = %.6f, sigma = %.6f' % (prf, val, sigma)

       else:
           prf = 'MEANCLIP:'
           print '%s %.1f-sigma clipped mean' % (prf, clipsig)
           print '%s Mean computed in %i iterations' % (prf, iter)
           print '%s Mean = %.6f, sigma = %.6f' % (prf, val, sigma)
	
   if return_array: return np.copy(arrind)
   else: return val, sigma


def get_pivot_wavel(hdr):
	'''RETURN THE PIVOT WAVELENGTH (THUS FAR, ONLY WORKS FOR WFC3-F275W,F438W,F606W, & F814W) '''

	pivot = {}
	pivot['WFC3'] = {'F275W':0.2704, 'F438W':0.4325, 'F606W':0.5887, 'F814W':0.8024}

	#--get instrument/filter information
	instrum = hdr['INSTRUME']

	if hdr['INSTRUME'] == 'WFC3': filter = hdr['FILTER']
	else: #assuming ACS
        	filter = fheader['FILTER1']
                if filter[0] == 'C': filter == hdr['FILTER2']

	#--return pivot wavelength (if available)
	if pivot.has_key(instrum) and pivot[instrum].has_key(filter): return pivot[instrum][filter]
	else:	return -9999.0
		#raise exception('Instrument and/or filter not currently supported (manually add to pivot list).'
		#return -9999.0


if __name__=='__main__':
        # Parse input parameters
        parser = argparse.ArgumentParser(description='Measure the EE from Abhi Stepping Program.')
        parser.add_argument('-im', '--images',default='*fl?.fits', type=str, help='Input fits file(s). \
                              Default is all FLT science images in current working directory.')
	parser.add_argument('-back', '--background', default='mean', type=str, choices=['mean','median','mode'], \
	                     help='Method of background subtraction for photometry (default=mean).')
        options = parser.parse_args()
        backmeth=options.background

	# Initialize filename and aperture correction variables
	file_list = glob.glob(options.images)
	cnts_name = []
	for file, ff in zip(file_list, xrange(len(file_list))):
		tmp=file.split('.fits')
		cnts_name.append(tmp[0] + '_cnts.fits')

        find_name = [cnts_name[x]+'0.coo' for x in xrange(len(cnts_name))]
        phot_name = [cnts_name[x]+'0.mag' for x in xrange(len(cnts_name))]


	# Initialize a dictionary (structure) to hold date for each image
	data = {}
	find_sharp = []
	find_round = []
	find_ground = []
	# Generate source catalogs: src detection & photometry within same image (i.e., not color catalogs)
        for file, ff in zip(file_list,xrange(len(file_list))):
		make_PAMcorr_image(file, outfile=cnts_name[ff])					
		run_daofind(cnts_name[ff], outfile=find_name[ff], dthreshold=800.0)
                replace_filevalue(find_name[ff], 'INDEF',0.0)   # all were associated with "sharp" parameter -- setting it to value that keeps the source within accepted boundary
		xx,yy,mm,sharp,round,ground,id = np.loadtxt(find_name[ff],unpack=True)

		# --Run through a decision tree to select the correct object if more than deteced by find
		if xx.size == 0:
                    print 'NO OBJECTS WERE DETECTED.'
		elif xx.size > 1:
		   # If more than 1 object detected - first remove any objects near borders as these are not expected
		    xxx = xx[(xx > 150.) & (xx < 410.) & (yy > 150.) & (yy < 365.)]
		    yyy = yy[(xx > 150.) & (xx < 410.) & (yy > 150.) & (yy < 365.)]
		    mmm = mm[(xx > 150.) & (xx < 410.) & (yy > 150.) & (yy < 365.)]
                    sxx = sharp[(xx > 150.) & (xx < 410.) & (yy > 150.) & (yy < 365.)]
                    rxx = round[(xx > 150.) & (xx < 410.) & (yy > 150.) & (yy < 365.)]
                    gxx = ground[(xx > 150.) & (xx < 410.) & (yy > 150.) & (yy < 365.)]
		    if xxx.size == 0:
                        print 'NO OBJECTS MATCHED IMAGE LOCATION CRITERIA.'
			if cnts_name[ff] == 'iaaug1jlq_flt.clean_cnts.fits':
			    xxx = np.array([xx[1]])
			    yyy = np.array([yy[1]])
			    sxx = np.array([sharp[1]])
			    rxx = np.array([round[1]])
			    gxx = np.array([ground[1]])
			    np.savetxt(find_name[ff],zip(xxx,yyy),fmt='%0.5f')
                        else: pdb.set_trace()
		    elif xxx.size == 1: np.savetxt(find_name[ff],zip(xxx,yyy),fmt='%0.5f')
		    else:
			# select the objects by find geometric parameters
			gd = np.where((sxx < 0.82) & (np.abs(gxx) < 0.4))[0]
			if gd.size == 1:
			    np.savetxt(find_name[ff],zip(xxx[gd],yyy[gd]),fmt='%0.5f')
			    sxx = sxx[gd]
			    rxx = rxx[gd]
			    gxx = gxx[gd]
			elif gd.size > 1:
			    # select the object with the brightest magnitude
			    #dist = np.sqrt((xxx[gd]-255.0)**2 + (yyy[gd]-255.0)**2)
			    mag = mmm[gd]
                            xxx = xxx[gd][mag == np.min(mag)]
                            yyy = yyy[gd][mag == np.min(mag)]
			    sxx = sxx[gd][mag == np.min(mag)]
                            rxx = rxx[gd][mag == np.min(mag)]
                            gxx = gxx[gd][mag == np.min(mag)]
                            np.savetxt(find_name[ff],zip(xxx,yyy),fmt='%0.5f')
			else:
			    print 'NO OBJECTS MATCHED SHARPNESS/ROUNDNESS CRITERIA.'
                            pdb.set_trace()
		else:
		    sxx = np.array([sharp])
		    rxx = np.array([round])
		    gxx = np.array([ground])


		find_sharp.append(sxx[0])
		find_round.append(rxx[0])
		find_ground.append(gxx[0])

		back, backrms = run_daophot(cnts_name[ff], coordfile=find_name[ff], outfile=phot_name[ff], calgorithm='gauss', backmethod=backmeth, cbox=10.)
		replace_filevalue(phot_name[ff], 'INDEF',-9999.0)
		iraf.txdump(phot_name[ff],'xcenter,ycenter,flux, mag', 'yes',Stdout=phot_name[ff]+'.trimmed')
		xc,yc,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26, \
		m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26=np.loadtxt(phot_name[ff]+'.trimmed',unpack=True)

		# record various properties of observation (chip,amp,mjd, etc.)
		dd=pyfits.open(file)
		hdr0=dd[0].header
		hdr1=dd[1].header
		im = dd[1].data
		dd.close()
		amps = hdr0['CCDAMP']
		biaslev = hdr0['BIASLEV'+amps]
		filter = hdr0['FILTER']
		LTV1 = hdr1['LTV1']
		LTV2 = hdr1['LTV2']
		chip = hdr1['CCDCHIP']
		shut = hdr0['SHUTRPOS']
		expt = hdr0['EXPTIME']
		mjd_avg = (hdr0['EXPEND'] - hdr0['EXPSTART'])/2. + hdr0['EXPSTART']
		mjd_deltat = (hdr0['EXPEND'] - hdr0['EXPSTART'])*24.0*60.0		# time between observation starts in minutes
		sizaxis1 = hdr1['NAXIS1']
		sizaxis2 = hdr1['NAXIS2']


 		xcp = xc - LTV1
		ycp = yc - LTV2


                # bit-wise or (add unique bits) DQ array across 3-pixel radius aperture
                dq = pyfits.getdata(file,ext=3)
                subdq = np.int32(dq[circular_mask(dq.shape,3, x_offset=(xc-sizaxis1/2.0),  y_offset=(yc-sizaxis2/2.0))])
                bitor = np.int32(0)
                for subind in subdq: bitor |= subind


		data[ff] = {'filename':file, 'amp':amps,'shutter':shut,'mjd_avg':mjd_avg, 'mjd_deltat': mjd_deltat, \
				'chip': chip, 'axis1':sizaxis1, 'axis2':sizaxis2, 'xc':xc, 'yc':yc,'xcp':xcp, 'ycp':ycp, \
				'background': back, 'background_rms':backrms, 'exptime': expt, 'biaslevel': biaslev, 'dqflag':bitor, \
				'flux':[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26], \
				'mag':[m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26]}


		# save image of source with measured barycenter on top
		if ff == 0:
		    fig = pylab.figure()
                    fig.subplots_adjust(wspace=0.4)
		pylab.clf()


		# plot #1 - object position
		sz=50.0
		x0=np.round(xc)-sz/2.
		x1=np.round(xc)+sz/2.
		y0=np.round(yc)-sz/2.
                y1=np.round(yc)+sz/2.
		ax1 = pylab.subplot(2,2,1)
		ax1.imshow(np.log10(im[y0:y1,x0:x1]),interpolation='nearest')
		ax1.autoscale(axis='both',enable=False)
		ax1.scatter([xc-x0-1.0], [yc-y0-1.0], marker='x', s=200., color='w')
		pylab.title('X = '+str(xc)+'  Y = '+str(yc),fontsize='small')


                # plot #2 - dq array at object position
                sz=50.0
                x0=np.round(xc)-sz/2.
                x1=np.round(xc)+sz/2.
                y0=np.round(yc)-sz/2.
                y1=np.round(yc)+sz/2.
                ax2 = pylab.subplot(2,2,3)
                ax2.imshow(dq[y0:y1,x0:x1],interpolation='nearest')
                ax2.autoscale(axis='both',enable=False)
                ax2.scatter([xc-x0-1.0], [yc-y0-1.0], marker='x', s=200., color='w',alpha=0.1)
                pylab.title('DQ ARRAY',fontsize='small')


		# plot #3 - background histogram
		tmp_image=glob.glob('*back.fits')[0]
		backim = pyfits.getdata(tmp_image)
		#--measure back statistics (mean and mode via IRAF)
		initback = iraf.imstatistics(tmp_image+'[0]', fields='mode,stddev', lower = -100, upper = 10000, nclip=7,lsigma=3.0, usigma=3.0, cache='yes', format='no',Stdout=1)
                llim = float(initback[0].split('  ')[0]) - 10.0*float(initback[0].split('  ')[1])
                ulim = float(initback[0].split('  ')[0]) + 10.0*float(initback[0].split('  ')[1])	
		backstats=iraf.imstatistics(tmp_image+'[0]', fields='mean,mode,stddev', lower = llim, upper = ulim, nclip=7,lsigma=3.0, usigma=3.0, cache='yes', format='no',Stdout=1)
		backmean=float(backstats[0].split('  ')[0])
                backmode=float(backstats[0].split('  ')[1])
		backrms=float(backstats[0].split('  ')[2])
		fbackim= np.ndarray.flatten(backim)
		gd=np.where((fbackim > llim) & (fbackim < ulim))[0]
		backmedian=meanclip(fbackim[gd],maxiter=7,return_median=1)[0]

                ax3 = pylab.subplot(1,2,2)
		nbins = np.ceil(80.0/(0.1*backrms))
                pylab.hist(fbackim[gd],log=True,bins=nbins,range=(-40.0,40.0))
		pylab.ylim(0.5,600000)
		pylab.xlim(-20,20)
		pylab.plot([backmode,backmode],[0.5,600000],ls='-',color='red',label='mode')
		pylab.plot([backmedian,backmedian],[0.5,600000],ls='--',color='aqua',label='median')
		pylab.plot([backmean,backmean],[0.5,600000],ls=':',color='black',label='mean')
		pylab.legend(loc=2, handletextpad=0.0,borderpad=0.0,frameon=False,handlelength=1.)
		pylab.title('Histogram of Background Pixels')
		pylab.xlabel('Background [e-]')
		pylab.ylabel('Number of Pixels')
		pylab.annotate('chip '+str(data[ff]['chip']),[0.77,0.95],xycoords='axes fraction')
		pylab.savefig(file.split('.fits')[0]+'_srcloc.pdf')


                tmp = glob.glob('*cnts.fit*')
                for tt in tmp: os.remove(tt)

	# SAVE ASCII CATALOG FOR SOURCES
	fnarr = [data[ff]['filename'] for ff in xrange(len(data))]
	amparr = [data[ff]['amp'] for ff in xrange(len(data))]
        shutarr = [data[ff]['shutter'] for ff in xrange(len(data))]
        mjdarr = [data[ff]['mjd_avg'] for ff in xrange(len(data))]
        mjddeltarr = [data[ff]['mjd_deltat'] for ff in xrange(len(data))]
        chiparr = [data[ff]['chip'] for ff in xrange(len(data))]
        axis1arr = [data[ff]['axis1'] for ff in xrange(len(data))]
        axis2arr = [data[ff]['axis2'] for ff in xrange(len(data))]
        xcarr = [data[ff]['xc'] for ff in xrange(len(data))]
        ycarr = [data[ff]['yc'] for ff in xrange(len(data))]
        xcparr = [data[ff]['xcp'] for ff in xrange(len(data))]
        ycparr = [data[ff]['ycp'] for ff in xrange(len(data))]
        backarr = [data[ff]['background'] for ff in xrange(len(data))]
        backrmsarr = [data[ff]['background_rms'] for ff in xrange(len(data))]
        exptimearr = [data[ff]['exptime'] for ff in xrange(len(data))]
	biaslevelarr = [data[ff]['biaslevel'] for ff in xrange(len(data))]
	dqflagarr = [data[ff]['dqflag'] for ff in xrange(len(data))]
	f1 = [data[ff]['flux'][0] for ff in xrange(len(data))]
        f2 = [data[ff]['flux'][1] for ff in xrange(len(data))]
        f3 = [data[ff]['flux'][2] for ff in xrange(len(data))]
        f4 = [data[ff]['flux'][3] for ff in xrange(len(data))]
        f5 = [data[ff]['flux'][4] for ff in xrange(len(data))]
        f6 = [data[ff]['flux'][5] for ff in xrange(len(data))]
        f7 = [data[ff]['flux'][6] for ff in xrange(len(data))]
        f8 = [data[ff]['flux'][7] for ff in xrange(len(data))]
        f9 = [data[ff]['flux'][8] for ff in xrange(len(data))]
        f10 = [data[ff]['flux'][9] for ff in xrange(len(data))]
        f12 = [data[ff]['flux'][10] for ff in xrange(len(data))]
        f14 = [data[ff]['flux'][11] for ff in xrange(len(data))]
        f16 = [data[ff]['flux'][12] for ff in xrange(len(data))]
        f18 = [data[ff]['flux'][13] for ff in xrange(len(data))]
        f20 = [data[ff]['flux'][14] for ff in xrange(len(data))]
        f24 = [data[ff]['flux'][15] for ff in xrange(len(data))]
        f28 = [data[ff]['flux'][16] for ff in xrange(len(data))]
        f32 = [data[ff]['flux'][17] for ff in xrange(len(data))]
        f36 = [data[ff]['flux'][18] for ff in xrange(len(data))]
        f40 = [data[ff]['flux'][19] for ff in xrange(len(data))]
        f45 = [data[ff]['flux'][20] for ff in xrange(len(data))]
        f50 = [data[ff]['flux'][21] for ff in xrange(len(data))]
        f55 = [data[ff]['flux'][22] for ff in xrange(len(data))]
        f60 = [data[ff]['flux'][23] for ff in xrange(len(data))]
        f65 = [data[ff]['flux'][24] for ff in xrange(len(data))]
        f70 = [data[ff]['flux'][25] for ff in xrange(len(data))]

        m1 = [data[ff]['mag'][0] for ff in xrange(len(data))]
        m2 = [data[ff]['mag'][1] for ff in xrange(len(data))]
        m3 = [data[ff]['mag'][2] for ff in xrange(len(data))]
        m4 = [data[ff]['mag'][3] for ff in xrange(len(data))]
        m5 = [data[ff]['mag'][4] for ff in xrange(len(data))]
        m6 = [data[ff]['mag'][5] for ff in xrange(len(data))]
        m7 = [data[ff]['mag'][6] for ff in xrange(len(data))]
        m8 = [data[ff]['mag'][7] for ff in xrange(len(data))]
        m9 = [data[ff]['mag'][8] for ff in xrange(len(data))]
        m10 = [data[ff]['mag'][9] for ff in xrange(len(data))]
        m12 = [data[ff]['mag'][10] for ff in xrange(len(data))]
        m14 = [data[ff]['mag'][11] for ff in xrange(len(data))]
        m16 = [data[ff]['mag'][12] for ff in xrange(len(data))]        
        m18 = [data[ff]['mag'][13] for ff in xrange(len(data))]
        m20 = [data[ff]['mag'][14] for ff in xrange(len(data))]
        m24 = [data[ff]['mag'][15] for ff in xrange(len(data))]
        m28 = [data[ff]['mag'][16] for ff in xrange(len(data))]
        m32 = [data[ff]['mag'][17] for ff in xrange(len(data))]
        m36 = [data[ff]['mag'][18] for ff in xrange(len(data))]
        m40 = [data[ff]['mag'][19] for ff in xrange(len(data))]
        m45 = [data[ff]['mag'][20] for ff in xrange(len(data))]
        m50 = [data[ff]['mag'][21] for ff in xrange(len(data))]
        m55 = [data[ff]['mag'][22] for ff in xrange(len(data))]
        m60 = [data[ff]['mag'][23] for ff in xrange(len(data))]
        m65 = [data[ff]['mag'][24] for ff in xrange(len(data))]
        m70 = [data[ff]['mag'][25] for ff in xrange(len(data))]

	tt = {'#filename':fnarr, 'amp':amparr, 'shutter':shutarr, 'mjd_avg':mjdarr, 'mjd_deltat':mjddeltarr, 'chip':chiparr, 'axis1':axis1arr, 'axis2':axis2arr, \
	      'xc':xcarr, 'yc':ycarr, 'xcp':xcparr, 'ycp':ycparr, 'background':backarr, 'background_rms':backrmsarr, 'exptime':exptimearr, 'biaslevel': biaslevelarr, 'dqflag': dqflagarr, \
	      'f1':f1, 'f2':f2, 'f3':f3,'f4':f4,'f5':f5,'f6':f6,'f7':f7,'f8':f8,'f9':f9,'f10':f10,'f12':f12,'f14':f14,'f16':f16,'f18':f18,'f20':f20,\
	      'f24':f24,'f28':f28,'f32':f32,'f36':f36,'f40':f40,'f45':f45,'f50':f50,'f55':f55,'f60':f60,'f65':f65,'f70':f70, \
              'm1':m1, 'm2':m2, 'm3':m3,'m4':m4,'m5':m5,'m6':m6,'m7':m7,'m8':m8,'m9':m9,'m10':m10,'m12':m12,'m14':m14,'m16':m16,'m18':m18,'m20':m20,\
              'm24':m24,'m28':m28,'m32':m32,'m36':m36,'m40':m40,'m45':m45,'m50':m50,'m55':m55,'m60':m60,'m65':m65,'m70':m70}


	ascii.write(tt, filter+'_photcat.dat', names=['#filename','amp','shutter','mjd_avg','mjd_deltat','chip','axis1','axis2','xc','yc','xcp','ycp','background','background_rms','exptime', 'biaslevel', 'dqflag', \
                                                'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f12','f14','f16',\
                                                'f18','f20','f24','f28','f32','f36','f40','f45','f50','f55','f60','f65','f70',\
                                                'm1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m12','m14','m16',\
                                                'm18','m20','m24','m28','m32','m36','m40','m45','m50','m55','m60','m65','m70'], \
						formats={'#filename':'%s','amp':'%s','shutter':'%s','mjd_avg':'%9.4f', 'mjd_deltat':'%6.4f', 'chip':'%i', 'axis1':'%i', 'axis2':'%i',\
						'xc':'%8.3f', 'yc':'%8.3f', 'xcp':'%8.3f', 'ycp':'%8.3f', 'background':'%0.5f','background_rms':'%0.5f', 'exptime':'%0.2f', 'biaslevel': '%0.4f', 'dqflag': '%i', \
						'f1':'%0.2f', 'f2':'%0.2f','f3':'%0.2f','f4':'%0.2f','f5':'%0.2f','f6':'%0.2f','f7':'%0.2f','f8':'%0.2f','f9':'%0.2f','f10':'%0.2f','f12':'%0.2f','f14':'%0.2f','f16':'%0.2f',\
						'f18':'%0.2f', 'f20':'%0.2f','f24':'%0.2f','f28':'%0.2f','f32':'%0.2f','f36':'%0.2f','f40':'%0.2f','f45':'%0.2f','f50':'%0.2f','f55':'%0.2f','f60':'%0.2f','f65':'%0.2f','f70':'%0.2f',\
						'm1':'%0.2f', 'm2':'%0.2f','m3':'%0.2f','m4':'%0.2f','m5':'%0.2f','m6':'%0.2f','m7':'%0.2f','m8':'%0.2f','m9':'%0.2f','m10':'%0.2f','m12':'%0.2f','m14':'%0.2f','m16':'%0.2f',\
                                                'm18':'%0.2f','m20':'%0.2f','m24':'%0.2f','m28':'%0.2f','m32':'%0.2f','m36':'%0.2f','m40':'%0.2f','m45':'%0.2f','m50':'%0.2f','m55':'%0.2f','m60':'%0.2f','m65':'%0.2f','m70':'%0.2f'})


	# Make diagnostic plots to determine which find parameters do a good job of identifying the standard star
	fnarr = np.array(fnarr)
	mjd_avg=np.array(mjdarr)
	find_sharp = np.array(find_sharp)
        find_round = np.array(find_round)
        find_ground = np.array(find_ground)

        #badf = ['ibbs25fgq_flt.fits', 'ibbsr1qqq_flt.fits', 'ibuc51aeq_flt.fits', 'ibbsy4n1q_flt.fits']	# F814W bad
	#badfind = np.where((fnarr == badf[0]) | (fnarr == badf[1]) | (fnarr == badf[2]) | (fnarr == badf[3]))[0]
        pylab.figure()
	pylab.scatter(mjd_avg,find_sharp)
	#pylab.scatter(mjd_avg[badfind],find_sharp[badfind],color='r',marker='o',s=100)
        pylab.figure()
        pylab.scatter(mjd_avg,find_round)
        #pylab.scatter(mjd_avg[badfind],find_round[badfind],color='r',marker='o',s=100)
        pylab.figure()
        pylab.scatter(mjd_avg,find_ground)
        #pylab.scatter(mjd_avg[badfind],find_ground[badfind],color='r',marker='o',s=100)

        pdb.set_trace()
