from __future__ import print_function
#DAO_SEXPARAMS = '/Users/david/Dropbox/research/photpipe/config/PS1/GPC1v3/idldaophot.default.sex'
#DAO_NNW = '/Users/david/Dropbox/research/photpipe/Cfiles/src/sex/config/default.nnw'
#HOSTPHOT_SEXOUTCOLS = '/Users/david/Dropbox/research/photpipe/config/PS1/GPC1v3/hostphot.sex.param'
#DAO_KERNEL = '/Users/david/Dropbox/research/photpipe/Cfiles/src/sex/config/gauss_4.0_7x7.conv'
DAO_SEXPARAMS = 'sexparams/idldaophot.default.sex'
DAO_NNW = 'sexparams/default.nnw'
HOSTPHOT_SEXOUTCOLS = 'sexparams/hostphot.sex.param'
DAO_KERNEL = 'sexparams/gauss_4.0_7x7.conv'
import os
import numpy as np
import sys
from subprocess import check_call
from txtobj import txtobj
from astropy.io import fits
from astropy import wcs

def getHostMatchfromSex(imfile,eventra,eventdec,sextable=None,maxRparam=6,returnShapePars=False,ext=0,zpt=None):
    print('getting host nearest position %s, %s'%(eventra,eventdec))

    if not sextable:
        sextable = runsex(imfile,zpt=zpt)
    hdr = fits.getheader(imfile,ext=ext)
    ImWCS = wcs.WCS(hdr)

    xylist = ImWCS.wcs_world2pix(eventra,eventdec,0)
    snx,sny = xylist

    C_xx = np.cos(sextable.THETA_IMAGE*np.pi/180.)**2./sextable.A_IMAGE**2. + \
           np.sin(sextable.THETA_IMAGE*np.pi/180.)**2./sextable.B_IMAGE**2.
    C_yy = np.sin(sextable.THETA_IMAGE*np.pi/180.)**2./sextable.A_IMAGE**2. + \
           np.cos(sextable.THETA_IMAGE*np.pi/180.)**2./sextable.B_IMAGE**2.
    C_xy = 2*np.cos(sextable.THETA_IMAGE*np.pi/180.)*np.sin(sextable.THETA_IMAGE*np.pi/180.)*(1./sextable.A_IMAGE**2. + 1/sextable.B_IMAGE**2.)
    x_r = snx - sextable.X_IMAGE
    y_r = sny - sextable.Y_IMAGE

    Rpar = np.sqrt(C_xx*x_r**2. + C_yy*y_r**2. + np.abs(C_xy*x_r*y_r))
    if np.min(Rpar) < maxRparam:
        iRpar = np.where(Rpar == np.min(Rpar))[0]
        hostra,hostdec = sextable.X_WORLD[iRpar],sextable.Y_WORLD[iRpar]
        if returnShapePars:
            return(hostra,hostdec,np.min(Rpar),sextable.A_IMAGE[iRpar],sextable.B_IMAGE[iRpar],sextable.THETA_IMAGE[iRpar],sextable.MAG_AUTO[iRpar])
        else:
            return(hostra,hostdec)
    else:
        if returnShapePars:
            return(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
        else:
            return(np.nan,np.nan)


def runsex(imfile, wtfile=None, maskfile=None, zpt=None, segmapname=None):

    imroot, imext = os.path.splitext(imfile)
    tmpcatname = "%s.sexcat.txt" % (imroot)
    if os.path.exists(tmpcatname):
        os.system('rm %s' % tmpcatname)

    sexcommand = ["sex "]
    sexcommand.append("%s " % imfile)
    sexcommand.append("-c %s " % DAO_SEXPARAMS)
    sexcommand.append("-STARNNW_NAME %s " % DAO_NNW)
    sexcommand.append("-PARAMETERS_NAME %s " % HOSTPHOT_SEXOUTCOLS)
    sexcommand.append("-FILTER_NAME %s " % DAO_KERNEL)
    sexcommand.append("-DETECT_THRESH 2.0 ")
    sexcommand.append("-GAIN 1.0 ")
    sexcommand.append("-DEBLEND_MINCONT 0.5 ")
    sexcommand.append("-CATALOG_NAME %s " % tmpcatname)
    if segmapname:
        sexcommand.append("-CHECKIMAGE_TYPE SEGMENTATION ")
        sexcommand.append("-CHECKIMAGE_NAME %s " % segmapname)
    if zpt:
        sexcommand.append("-MAG_ZEROPOINT %s " % zpt)
    sexcommand = ''.join(sexcommand)

    if wtfile:
        if maskfile:
            fname, fext = os.path.splitext(wtfile)
            noiseimfilename_tmp = '%s_tmp%s' % (fname, fext)
            mask = fits.getdata(maskfile)
            hdu = fits.open(wtfile)
            hdu[1].data[np.where(mask > 0)] = 1e16
            hdu.writeto(noiseimfilename_tmp, overwrite=True)
            sexcommand += ' -WEIGHT_TYPE MAP_VAR -WEIGHT_IMAGE %s' % noiseimfilename_tmp
        else:
            sexcommand += ' -WEIGHT_TYPE MAP_VAR -WEIGHT_IMAGE %s' % wtfile

            print('running SExtractor:')

    print(sexcommand)

    try:
        check_call(sexcommand, shell=True)
    except Exception as e:
        print('error: sextractor invocation failed', file=sys.stderr)
        print('command was:', sexcommand, file=sys.stderr)
        #os.system('rm %s %s'%(tmpcatname,noiseimfilename_tmp))
        raise

    if maskfile and wtfile:
        print('removing temporary noise image %s' % noiseimfilename_tmp)
        os.system('rm %s' % noiseimfilename_tmp)

    # Read in the catalog
    sextable = txtobj(tmpcatname, sexheader=True)

    ## DC - commenting out below for a test...
    # os.system('rm %s' % tmpcatname)

    return(sextable)
