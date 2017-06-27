#!/usr/bin/env python2
#########################################################################################
#                                                                                       #							
#   Dr. Andry Rajoelimanana                                                             #	
#                                                                                       #
#   SALT/RSS longslit reduction                                                         #
#                                                                                       #		
#########################################################################################

import yaml
from astropy.io import fits 
import os, sys
from pyraf import iraf
import numpy as np
iraf.noao(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.longslit(_doprint=0)
iraf.apextract(dispaxis='1')

def header(images, param):
    hdrs = fits.getheader(str(images))    
    if param in hdrs.keys():
        return hdrs[param]
    else:
        return False
        
# Add headers        
def putheader(images, param, value):    
    iraf.hedit(images, fields=param, value=value, add='yes', addonly='yes',
               verify='no',show='yes')        

# Get shell output        
def getsh(command):
        child = os.popen(command)
        data = child.read()
        err = child.close()
        if err:
                raise RuntimeError, '%s failed with error code %d' \
                        % (command, err)
                print '\nTHIS IS VERY BAD!!!'
                sys.exit(1)
        return data

# Make calibaration file using RSS setups 
def makesetups(images, maketype):     
       LAMPID = header(images,'LAMPID').replace(' ','')
       GRATING = header(images,'GRATING')
       GRANGLE = "%.2f" % header(images,'GR-ANGLE')
       ARANGLE = "%.2f" % header(images,'AR-ANGLE')
       BIN = header(images,'CCDSUM').replace(' ','x')     
       if (maketype == 'std'):
           sts = "Std_"+GRATING+"_"+ARANGLE+"_"+GRANGLE+"_"+BIN
       else:
           sts = LAMPID+"_"+GRATING+"_"+ARANGLE+"_"+GRANGLE+"_"+BIN
       return sts


# Remove existing file
def rmexist(inlist):
    for infile in inlist:
        if os.path.isfile(infile):
            os.remove(infile)

# Merge configuration file
def config_merge(custom, default):
    if isinstance(custom, dict) and isinstance(default, dict):
        for k,v in default.iteritems():
            if k not in custom:
                custom[k] = v
            else:
                custom[k] = config_merge(custom[k], v)
    return custom

# Load configuration file
def loadparam(config, sections):
    for section_name in sections:
        params = config[section_name].items()
        for param_id in params: 
            eval(section_name).setParam(param_id[0], param_id[1])    

# Load configuration file and merge it with the default one            
def load_config(config_file, config_default='/home/rajoelimananaa/software/saltpy/rssconf.yml'):
    config = yaml.safe_load(open(config_default))
    custom_config = {}
    if os.path.exists(config_file):
        custom_config = yaml.safe_load(open(config_file))
        
    config = config_merge(custom_config, config)

    sections = ['iraf.identify',
            'iraf.reidentify',
            'iraf.fitcoords',
            'iraf.transform',
            'iraf.fit1d',
            'iraf.apall',
            ]

    for section_name in sections:
        params = config[section_name].items()
        for param_id in params:
            eval(section_name).setParam(param_id[0], param_id[1])
    return config    


# Prepare images for reduction
def prepare_image(images, do_error=False):

    hdrs = fits.getheader(images)
    if 'OBSERVAT' not in hdrs.keys(): putheader(images, 'OBSERVAT', 'SAAO')
    if 'DISPAXIS' not in hdrs.keys(): putheader(images, 'DISPAXIS', '1')
    if 'AIRMASS' not in hdrs.keys(): rssairmass(images)
    
    if do_error:
        error_file = os.path.splitext(images)[0]+'_var.fits'
        hdrerrs =  fits.getheader(error_file)
        if 'OBSERVAT' not in hdrerrs.keys(): putheader(error_file, 'OBSERVAT', 'SAAO')
        if 'DISPAXIS' not in hdrerrs.keys(): putheader(error_file, 'DISPAXIS', 1)
        if 'AIRMASS' not in hdrerrs.keys(): rssairmass(error_file)


# Calculate airmass
def rssairmass(images):

    airms = header(images, 'TELALT')
    airmass = "%.5f" % (1./np.cos(np.radians(90.0-airms)))
    putheader(images, 'AIRMASS', airmass)


# Get approximate wavelength range
def rsswave(images):

    grating = header(images,'GRATING')[0:6]
    grang = float(header(images,'GR-ANGLE'))
    arang = float(header(images,'CAMANG'))
    cbin = int(header(images,'CCDSUM').split(' ')[0])
    cols = int(header(images,'NAXIS1'))
    if (grating == 'PG0300'):
        grat=300
        gamma0 = 0.0
    elif grating=='PG0900':
        grat=903.2
        gamma0 = -0.265
    elif grating=='PG1300':
        grat=1299.76
        gamma0 = -0.265
    elif grating=='PG1800':
        grat=1801.89
        gamma0 = -0.265
    elif grating=='PG2300':
        grat=2302.60
        gamma0 = -0.265
    elif grating=='PG3000':
        grat=3000.55
        gamma0 = -0.265
    T2Con=-5.00	
    T3Con=-1.00
    FCampoly = np.asarray([-0.0023,0.0365,-0.2100,0.5061,-0.1861,328.697])
    Grat0 = 1.407798403
    ArtErr = -4.2E-05
    Home0 = -0.063025809
    alpha_r = np.radians(grang+Grat0)
    beta0_r = np.radians(arang*(1+ArtErr)+Home0)-alpha_r
    gam0_r = np.radians(gamma0)
    lam0 = 1e7*np.cos(gam0_r)*(np.sin(alpha_r) + np.sin(beta0_r))/grat
    ww = lam0/1000. - 4.
    fcam = np.polyval(FCampoly,ww)
    lmm = grat
    disp = (1e7*np.cos(gam0_r)*np.cos(beta0_r)/lmm)/(fcam/.015)
    dfcam = 3.162*disp*np.polyval([FCampoly[x]*(5-x) for x in range(5)],ww)
    T2 = -0.25*(1e7*np.cos(gam0_r)*np.sin(beta0_r)/lmm)/(fcam/47.43)**2 + T2Con*disp*dfcam
    T3 = (-1./24.)*3162.*disp/(fcam/47.43)**2 + T3Con*disp
    T0 = lam0 + T2 
    T1 = 3162.*disp + 3*T3
    X = (np.array(range(cols))+1-cols/2)*cbin/3162.
    lam_X = T0+T1*X+T2*(2*X**2-1)+T3*(4*X**3-3*X)
    iraf.hedit(images, "CRVAL1", min(lam_X), add='yes', addonly='yes', verify='no' )
    iraf.hedit(images, "CRPIX1", 1., add='yes', addonly='yes', verify='no' )
    iraf.hedit(images, "CDELT1", (max(lam_X)-min(lam_X))/cols, add='yes', addonly='yes', verify='no' )
    iraf.hedit(images, "CD1_1", (max(lam_X)-min(lam_X))/cols, add='yes', addonly='yes', verify='no' )
    
    return
 
 
