#! /usr/bin/env python	
'''								                             
#########################################################################################
#                                                                                       #							
#   Dr. Andry Rajoelimanana                                                             #	
#                                                                                       #
#   SALT/RSS longslit reduction                                                         #
#                                                                                       #		
#########################################################################################
'''

import os
from pyraf import iraf
from astropy.io import fits
from shutil import copy2
import rsstools as rt
import numpy as np

iraf.noao(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.longslit(_doprint=0)
iraf.apextract(dispaxis='1')
iraf.stsdas(_doprint=0)
iraf.analysis(_doprint=0)
iraf.fitting(_doprint=0)


class arc(object):
    def __init__(self, name):
        self.name = name 
        self.noext = os.path.splitext(self.name)[0]
        self.idname = rt.makesetups(self.name,"arc")
        self.idfits = self.idname+'.fits'
        self.linelist = self.idname.split('_')[0].lower()+'.dat'  
    def __str__(self):
        return self.name+' :  idname --> '+self.idname
            

class target(object):
    def __init__(self, name):
        rt.rmexist(['f'+name])
        copy2(name,'f'+name)
        self.name = 'f'+name 
        self.noext = os.path.splitext(self.name)[0]
        self.var = self.insert_key('_var.fits')
        self.wave = self.insert_key('w.fits')
        self.bkg = self.insert_key('ws.fits')
        self.flux2d = self.insert_key('wsf.fits')
        self.sky = self.insert_key('sky.fits')
        self.wave_var = self.insert_key('w_var.fits') 
        self.bkg_var = self.insert_key('ws_var.fits')
        self.bkg_err = self.insert_key('ws_err.fits')        
        self.flux2d_var = self.insert_key('wsf_var.fits')   
        self.flux2d_err = self.insert_key('wsf_err.fits')     
    def insert_key(self, key):
        return self.noext+key


class standard(target):
    def __init__(self,name):
        target.__init__(self,name)
        self.idname = rt.makesetups(self.name,'std')
        self.date = self.name.split('P')[1][0:8]
        self.sens = self.idname+'.fits'
        self.sensout = self.idname+'_'+self.date
        self.object= (rt.header(self.name,'OBJECT')).lower()
        self.oned0 = self.insert_key('ws.0001.fits')
        self.oned = self.insert_key('ws1d.fits')    
    def __str__(self):
        return self.name+' : sensitivity function --> '+self.sens        


class rsspipe(object):
    def __init__(self, objfile, objarcfile, stdfile=None, stdarcfile=None,\
                 rssconfigfile=None, logfile=None, do_error=True, \
                 rssdatadir = None):
        self.objectfile = objfile
        self.arcfile = objarcfile
        # read object and arc list 
        objlist = self.readfile(self.objectfile)
        self.objlist = []
        for obj in objlist:
            self.objlist.append(target(obj))
        self.arclist = self.readfile(self.arcfile)[0]
        self.arc = arc(self.arclist)
        # load configuration
        self.rssconfigfile = rssconfigfile if rssconfigfile is not None else 'rssconfig.yml'
        self.rssconfig = rt.load_config(self.rssconfigfile)
        # read standard spectrophotometric
        if stdfile is not None:
            self.stdfile = stdfile
            stdlist = self.readfile(self.stdfile)
            self.stdlist = []
            for std in stdlist:
                self.stdlist.append(standard(std))
            if stdarcfile is not None:
                self.stdarcfile = stdarcfile
                stdarclist = self.readfile(stdarcfile)
                self.stdarclist = standard(stdarclist[0])                
            else:
                self.stdarclist = self.arc 
        # loading config path file
        self.logfile = logfile if logfile is not None else 'saltlongslit.log'
        self.set_configfile(rssdatadir)
        self.do_error = do_error
        
    def readfile(self, filename):
        with open(filename) as f:
            lines = f.read().splitlines()
        return lines
        
    def set_configfile(self, rssdatadir):
        if rssdatadir is not None:
            rssdatadir = rssdatadir
        else: 
            rssdatadir = '/Users/ando/bin/salt_red/saltdatabase'
        while not os.path.exists(rssdatadir):
            print("RSS data directory "+str(rssdatadir)+" IS NOT valid")
            rssdatadir = raw_input('Enter a correct path for RSS data directory:  ')
        self.rssdatadir = rssdatadir
        self.caldir = self.rssdatadir+'/caldir/'
        self.extinction = self.caldir+'/suth_extinct_burgh.dat'
        self.linelist_dir = self.rssdatadir+'/linelists/' 
               
    def rssidentify(self, arc):
        if not os.path.exists('database'): os.system('mkdir database')
        rt.loadparam(self.rssconfig, ['iraf.identify', 'iraf.reidentify', 'iraf.fitcoords'])
        coordfile = self.linelist_dir+arc.linelist
        while not os.path.exists(coordfile):
            print("***** rssidentify error *** The coordinate "+str(coordfile)+" IS NOT in our directory *****")
            coordfile = raw_input('Enter a correct path for coordinate files:  ')
        arc.coordfile = coordfile
        copy2(coordfile, './')
        # check if the setup has been already identified before
        if os.path.isfile(self.rssdatadir+"/database/id"+arc.idname) and \
                        os.path.isfile(self.rssdatadir+"/"+arc.idfits):
            copy2(self.rssdatadir+"/database/id"+arc.idname, 'database/')
            copy2(self.rssdatadir+"/"+arc.idfits,'./')
        else:
            copy2(arc.name, arc.idfits)
            rt.rsswave(arc.idfits)
            iraf.identify(arc.idname, coordlist=coordfile, mode='h')
            iraf.reidentify(arc.idname, arc.idname, coordlist=coordfile, mode='h')  
            iraf.fitcoord(arc.idname, mode='h')
            copy2(arc.idfits, self.rssdatadir+'/'+arc.idfits)
            copy2('database/id'+arc.idname, self.rssdatadir+'/database/id'+arc.idname)
            copy2('database/fc'+arc.idname, self.rssdatadir+'/database/fc'+arc.idname)
        print("***** Reducing of ARCs *****")
        if not os.path.isfile("database/id"+arc.noext) or not os.path.isfile("database/fc"+arc.noext):
            iraf.reidentify(arc.idname, arc.noext, coordlist=coordfile, mode='h')
            iraf.reidentify(arc.noext, arc.noext, coordlist=coordfile, mode='h')        
            iraf.fitcoord(arc.noext, mode='h')
            
    def reduce2d(self, obj, arc, do_error=None):
        images = obj.name  
        do_error = do_error if do_error is not None else self.do_error
        hduorig= fits.open(images)
        os.remove(images)
        hduimg = fits.PrimaryHDU(data=hduorig[0].data,header=hduorig[0].header)
        hduimg.writeto(images)
        hduvarimg = fits.PrimaryHDU(data=hduorig[1].data,header=hduorig[0].header)
        rt.rmexist([obj.var])
        hduvarimg.writeto(obj.var)    
        rt.prepare_image(images, do_error=do_error)
        self.rssidentify(arc)
        nxpix = rt.header(images, 'NAXIS1')
        rt.loadparam(self.rssconfig, ['iraf.transform', 'iraf.fit1d'])
        config= self.rssconfig
        if self.do_error:
            iraf.errorpars.errtype='uniform'
            iraf.errorpars.resample='no'
            iraf.samplepars.axis=1
            iraf.samplepars.setParam('naverage', config['iraf.fit1d']['naverage'])
            iraf.samplepars.setParam('low_reject', config['iraf.fit1d']['low_reject'])
            iraf.samplepars.setParam('high_reject', config['iraf.fit1d']['high_reject'])
            iraf.samplepars.setParam('niterate', config['iraf.fit1d']['niterate'])
            iraf.samplepars.setParam('grow', config['iraf.fit1d']['grow'])                   
        rt.rmexist([obj.wave, obj.bkg, obj.wave_var, obj.bkg_var])
        # Wavelength calibration	   
        iraf.transform(obj.noext, obj.wave, fitnames = arc.noext, logfiles = self.logfile)
        iraf.fit1d(obj.wave, obj.sky, "fit", axis=2)                    
        iraf.imarith(obj.wave, "-", obj.sky, obj.bkg)
        # Calibrated error file
        if do_error:
            iraf.transform(obj.var, obj.wave_var, fitnames = arc.noext, logfiles = self.logfile)
            iraf.gfit1d(obj.wave, 'tmptab.tab', function=config['iraf.fit1d']['function'], order=config['iraf.fit1d']['order'], xmin=1, xmax=nxpix, interactive='no')		       
            iraf.tabpar('tmptab.tab',"rms", 1, Stdout=self.logfile)
            rmsval = float(iraf.tabpar.value)
            tmpval=rmsval*rmsval
            print('RMS of the fit (for VAR plane propagation) = %s ' %rmsval)
            iraf.imexpr("a+b",obj.bkg_var, a=obj.wave_var, b=tmpval)
            print("**** Done with image %s ****" % (images)) 
            
    # Estimate sensitivity curve using spectrophotometric standard
    def rssstdred(self, obj, arc):
        rt.rmexist([obj[0].sensout,obj[0].sensout+'.fits'])
        for std in obj:
            self.reduce2d(std, arc, do_error=False)
            rt.loadparam(self.rssconfig, ['iraf.apall'])
            std_data = self.caldir+'/'+std.object+".dat"
            rt.rmexist([std.oned,std.oned0])
            while not os.path.isfile(std_data):
                print("Cannot find file ", std_data)
                std_data = raw_input('Enter a star database:  ')
            iraf.apall(std.bkg, b_niterate=5, b_sample='-100:-40,40:100', b_high_reject=1.5, b_low_reject=2)
            os.system('mv %s %s' %  (std.oned0,std.oned))
            if not rt.header(std.oned,'AIRMASS'): rt.rssairmass(std.oned)
            iraf.standard(input=std.oned, output=std.sensout, caldir=self.caldir, interact='yes',
                      star_name=std.object, extinct=self.extinction)
        # Sensitivity files    
        iraf.sensfunc(std.sensout, std.sensout+'.fits', interactive= 'yes', extinct=self.extinction,
                      newextinction='extinct.dat')
        #rt.putheader(std.sensout+'.fits', 'FITSNAME', ','.join([std.name for std in obj]))
        rt.getsh('cp %s %s' %  (std.sensout+'.fits',self.rssdatadir)) 
        return std.sensout+'.fits'
        
    def rssflux(self, obj, sens, do_error=None):
        do_error = do_error if do_error is not None else self.do_error    
        rt.rmexist([obj.flux2d])
        if not rt.header(obj.bkg_var,'AIRMASS'): rt.rssairmass(obj.bkg_var)
        iraf.longslit.calibrate(input=obj.bkg, output=obj.flux2d, extinct='yes',
                                flux='yes', fnu='no',sensitivity=sens, 
                                extinction=self.extinction)
        rt.putheader(obj.flux2d, 'STDNAME', sens)
        if do_error:
            rt.rmexist([obj.flux2d_var, obj.flux2d_err]) 
            # Convert variance to error
            iraf.imexpr("sqrt(a)", obj.bkg_err, a = obj.bkg_var)
            iraf.longslit.calibrate(input=obj.bkg_err, output=obj.flux2d_err, extinct='yes',
                                    flux='yes',observa='saao',fnu='no', sensitivity=sens,
                                    extinction=self.extinction)
                                    
#import sys
#sys.path.insert(0,'/Users/ando/andry/research/Make_software/saltrss')
#test = rsspipe('obj1', 'arc1', stdfile='std1')
#bbb = test.rssstdred(test.stdlist, test.stdarclist)
#for inobj in test.objlist:
#    test.reduce2d(inobj, test.arc)
#    test.rssflux(inobj, bbb)
