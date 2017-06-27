#! /usr/bin/env python									                             
#########################################################################################
#                                                                                       #							
#   Dr. Andry Rajoelimanana                                                             #	
#                                                                                       #
#   SALT/RSS longslit reduction                                                         #
#                                                                                       #		
#########################################################################################

import os
from pyraf import iraf
from astropy.io import fits

import rsstools as rt
import numpy as np

iraf.noao(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.longslit(_doprint=0)
iraf.apextract(dispaxis='1')
iraf.stsdas(_doprint=0)
iraf.analysis(_doprint=0)
iraf.fitting(_doprint=0)

# Wavelength Calibration of the ARC spectra       
def rssidentify(arclist, config):

    saltreddir = config['rss']['rssdir'] 
    ref_spec= os.path.splitext(arclist)[0]
    reference= rt.makesetups(arclist,"arc")
    referfits = reference+'.fits'  
    referid_db = saltreddir+"/database/id"+str(reference)
    referfits_db = saltreddir+'/'+referfits   
    coordlist= reference.split('_')[0].lower()+'.dat'       
    coordfile = config['rss']['linelist_dir']+str(coordlist)

    rt.loadparam(config, ['iraf.identify', 'iraf.reidentify', 'iraf.fitcoords'])
      
    if reference == 'NONE_NONE_NONE_NONE_NONE':
        print "***** rssidentify error"
        print "***** We can not get your spectral setup *****"
        return


    if os.path.isfile(coordfile):
        rt.getsh('cp %s ./' %  (coordfile))
    else :
        print "***** rssidentify error *** The coordinate list for this arc IS NOT in our directory *****"
        return            
    
    if os.path.isfile(referfits_db) and os.path.isfile(referid_db):
        rt.getsh('cp %s ./'  % (referfits_db))
        if not os.path.exists('database'): os.system('mkdir database')
        rt.getsh('cp %s database/'  % (referid_db))

    else:
        print referfits_db, referid_db
        print "This setup is not in our database (setup = %s)" % (reference)
        rt.getsh('cp %s ./%s'  % (arclist, referfits))
        rt.rsswave(referfits)
        iraf.identify(reference, coordlist=coordlist, mode='h')
        iraf.reidentify(reference, reference, coordlist=coordlist, mode='h')
        iraf.fitcoord(reference, mode='h')

        os.system('cp '+referfits+' '+saltreddir+'/'+referfits)
        os.system('cp database/id'+reference+' '+saltreddir+'/database/id'+reference)
        os.system('cp database/fc'+reference+' '+saltreddir+'/database/fc'+reference)


    print "'***** Reducing of ARCs *****'"

    arc_id = "database/id"+str(ref_spec)
    arc_fc = "database/fc"+str(ref_spec)

    if not os.path.isfile(arc_id) or not os.path.isfile(arc_fc):
        iraf.reidentify(reference, ref_spec, coordlist=coordlist, mode='h')
        iraf.reidentify(ref_spec, ref_spec, coordlist=coordlist, mode='h')        
        iraf.fitcoord(ref_spec, mode='h')
       
# Wavelength calibration, background substraction, and flux calibration of the target spectra        
def rssextractone(images, arcl, config, logfile, do_error=False):

    print "'***** Reducing of objects ***** '"
    infiles = os.path.splitext(images)[0]
    ref_spec = os.path.splitext(arcl)[0]
    if os.path.isfile('b'+images): 
        hduorig= fits.open('b'+images)
        os.remove(images)
    else:
        hduorig= fits.open(images)
        os.system('mv '+images+'  b'+images)

    hduimg = fits.PrimaryHDU(data=hduorig[0].data,header=hduorig[0].header)
    hduimg.writeto(images)
    hduvarimg = fits.PrimaryHDU(data=hduorig[1].data,header=hduorig[0].header)
    rt.rmexist([infiles+'_var.fits'])
    hduvarimg.writeto(infiles+'_var.fits')  
    rt.prepare_image(images, do_error=do_error)
    rssidentify(arcl, config)
    nxpix = rt.header(images, 'NAXIS1')
        
    # output wavelength calibrated, background substracted file
    wavefits = infiles+'w.fits'
    backfits = infiles+'ws.fits'
    tmpfit = infiles+'sky.fits'
    error_file = infiles+'_var.fits'
    wavefits_err = infiles+'w_var.fits'
    backfits_err = infiles+'ws_var.fits'

    rt.loadparam(config, ['iraf.transform', 'iraf.fit1d'])
    # check if the error file exist       
    if do_error:
        iraf.errorpars.errtype='uniform'
        iraf.errorpars.resample='no'
        iraf.samplepars.axis=1
        iraf.samplepars.setParam('naverage', config['iraf.fit1d']['naverage'])
        iraf.samplepars.setParam('low_reject', config['iraf.fit1d']['low_reject'])
        iraf.samplepars.setParam('high_reject', config['iraf.fit1d']['high_reject'])
        iraf.samplepars.setParam('niterate', config['iraf.fit1d']['niterate'])
        iraf.samplepars.setParam('grow', config['iraf.fit1d']['grow'])
               
    rt.rmexist([wavefits, backfits, wavefits_err, backfits_err])
    # Wavelength calibration	   
    iraf.transform(infiles, wavefits, fitnames=ref_spec, logfiles = logfile)
    iraf.fit1d(wavefits, tmpfit, "fit", axis=2)                    
    iraf.imarith(wavefits, "-", tmpfit, backfits)

    # Calibrated error file
    if do_error:
        iraf.transform(error_file, wavefits_err, fitnames=ref_spec, logfiles = logfile)
        iraf.gfit1d(wavefits, 'tmptab.tab', function=config['iraf.fit1d']['function'], order=config['iraf.fit1d']['order'], xmin=1, xmax=nxpix, interactive='no')		       
        iraf.tabpar('tmptab.tab',"rms", 1, Stdout=logfile)
        rmsval = float(iraf.tabpar.value)
        tmpval=rmsval*rmsval
        print 'RMS of the fit (for VAR plane propagation) = %s ' %rmsval
        iraf.imexpr("a+b",backfits_err, a=wavefits_err, b=tmpval)

    print "**** Done with image %s ****" % (images)  


# Flux calibration using input sensitivity curve    
def rssfluxdefine(images,calflux,config, do_error=False):

    print "We will use setups file: %s" %(calflux)
     
    obj_nofits = os.path.splitext(images)[0]+'ws'
    inobj = obj_nofits+'.fits'
    inobj_err = obj_nofits+'_var.fits'
    extinct = config['rss']['extdir']
    calfits = obj_nofits+'f.fits'   
    rt.rmexist([calfits])
    
    if not rt.header(inobj_err,'AIRMASS'): rt.rssairmass(inobj_err)

    iraf.longslit.calibrate(input=inobj, output=calfits, extinct='yes', flux='yes',
                            fnu='no',sensitivity=calflux, extinction=extinct)
    rt.putheader(calfits, 'STDNAME', calflux)
 
    if do_error: 
        calfits_err = obj_nofits+'f_var.fits'
        obj_err = obj_nofits+'f_err.fits'
        rt.rmexist([calfits_err, obj_err])
        # Convert variance to error
        iraf.imexpr("sqrt(a)", obj_err, a = inobj_err)
        iraf.longslit.calibrate(input=obj_err, output=calfits_err, extinct='yes',
                                flux='yes',observa='saao',fnu='no', sensitivity=calflux,
                                extinction=extinct)
             
# Estimate sensitivity curve using spectrophotometric standard
def rssstdred(stdlist,stdarclist, sens, config, logfile):
    
    instdlist = np.loadtxt(stdlist, unpack=True, ndmin=1, dtype=str)
    instdarclist = np.loadtxt(stdarclist, unpack=True, ndmin=1, dtype=str)
    stdarcfile=instdarclist[0]
    saltreddir = config['rss']['rssdir'] 
    caldir = config['rss']['caldir']
    extinct= config['rss']['extdir']


    inmake=str(instdlist[0])
    stdsetup = rt.makesetups(stdarcfile,'std')
    
    rssdate = inmake.split('P')[1][0:8]
    std_out = str(stdsetup)+"_"+rssdate
    sens_out = std_out+'.fits'
    std_fits= str(std_out)+".fits"
    print "Here is the standard setup file : %s" %(std_fits)
    rt.rmexist([std_out,std_fits])
           
    for instdl in instdlist: 

        rssextractone(instdl,stdarcfile, config, logfile, do_error=False)        
        obj_nofits = os.path.splitext(instdl)[0]+'ws'
        inobj = obj_nofits+'.fits'
        star_name= (rt.header(instdl,'OBJECT')).lower()
        star_database = caldir+str(star_name)+".dat"
        onedfits = obj_nofits+"1d.fits"
        onedfits1 = obj_nofits+".0001.fits"
        rt.rmexist([onedfits,onedfits1])

        if not os.path.isfile(star_database):
            print ""
            print "*****  Cannot find file ",star_database, "*****"
            print "*****  in the directory caldir *****"
            print ""
            return
        rt.loadparam(config, ['iraf.apall'])
        iraf.apall(inobj, b_niterate=5, b_sample='-100:-40,40:100', b_high_reject=1.5, b_low_reject=2)

        os.system('mv %s %s' %  (onedfits1,onedfits))

        if not rt.header(onedfits,'AIRMASS'): rt.rssairmass(onedfits)
        print "**** start running standard %s" % (std_out)

        iraf.standard(input=onedfits, output=std_out, caldir=caldir, interact='yes',
                      star_name=star_name, extinct=extinct)
        
    # Sensitivity files    
    iraf.sensfunc(std_out, sens_out, interactive= 'yes', extinct=extinct,
                  newextinction='extinct.dat')
    rt.getsh('cp %s %s' %  (sens_out,str(saltreddir))) 
    rt.putheader(sens_out, 'FITSNAME', instdlist)

    for instdl in instdlist:    
        obj_nofits = os.path.splitext(instdl)[0]+'ws'
        inobj = obj_nofits+'.fits'
        in_cal = obj_nofits+'1d.fits'
        out_cal = obj_nofits+'1dp.fits'
        rt.rmexist([out_cal])
        iraf.calibrate(in_cal, output=out_cal, sensitivity=std_fits, extinction=extinct)

    return std_fits
            
        
def read_fromfile(infile):
    listout=[]
    with open(infile) as FileObj:
        listout = FileObj.read().splitlines()
    return listout     

       
def reduce(objlist,objarclist,stdlist='', stdarclist='', logfile='saltlongslit.log', do_error=True):

    inobjlist = np.loadtxt(objlist, unpack=True, ndmin=1, dtype=str)
    inobjarclist = np.loadtxt(objarclist, unpack=True, ndmin=1, dtype=str)
    arcfile=inobjarclist[0]

    # get parameter files to dict    
    config = rt.load_config('rssconfig.yml')

    if not os.path.exists('database'): os.system('mkdir database')
    if not os.path.isfile('login.cl'): os.system('mkiraf')

    # check Flux calibration file
    
    if (stdlist != '' and os.path.splitext(stdlist)[0] == stdlist):
        stdlist = rssstdred(stdlist, stdarclist,'yes', config, logfile)
        print stdlist
        if (os.path.splitext(stdlist)[0] != stdlist and os.path.isfile(stdlist)): 
             print 'Good to go'
        else:
             return  
    for inobj in inobjlist:
             
        print "Starting with object %s and %s" %(inobj,arcfile)
        rssextractone(inobj, arcfile, config, logfile, do_error=do_error)
        if (os.path.splitext(stdlist)[0] != stdlist and os.path.isfile(stdlist)):
            print "***** Standard file found =  %s *****" %(stdlist)
            rssfluxdefine(inobj, stdlist, config, do_error=do_error)
            
            
            

 
