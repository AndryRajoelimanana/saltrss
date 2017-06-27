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
import pyfits
import numpy as np
import string
import rsstools as rt


def extract(objlist, apallconf = 'apall_config.yml'):
    config = rt.apall_config(apallconf,config_default='/Users/ando/andry/software/pythonpath/rsspipe/apall_conf.yml')
    inobjlist = np.loadtxt(objlist,dtype="str",unpack=True, ndmin=1)
    for inobj in range(len(inobjlist)):
        print "Starting with object %s " %(inobjlist[inobj])
        objj = os.path.splitext(inobjlist[inobj])[0]
        if os.path.isfile(objj+'wsf.fits'):
            inobjfile = objj+'wsf.fits'
        else:
            inobjfile = objj+'ws.fits'            
        salt_apall(inobjfile)
    salt_scombine(objlist)

def salt_apall(obj):
    outname= string.split(obj,sep='.')[0]+"1d.fits"
    iraf.apall(input=obj,output=' ')
    os.system('mv '+string.split(obj,sep='.')[0]+'.0001.fits  '+outname)

    err_obj=string.split(obj,sep='.')[0]+'_var.fits'
    err_obj1d=string.split(obj,sep='.')[0]+'1d_var.fits'

    if os.path.isfile(err_obj):
         iraf.apall.references=obj        
         iraf.apall.interactive='no'
         iraf.apall.find='no'
         iraf.apall.recenter='no'
         iraf.apall.resize='no'
         iraf.apall.edit='no'
         iraf.apall.trace='no'
         iraf.apall.fittrace='no'
         iraf.apall.extract='yes'
         iraf.apall.extras='no'
         iraf.apall.review='no'
         iraf.apall(input=err_obj,output=err_obj1d)
         os.system('mv '+string.split(obj,sep='.')[0]+'1d_var.0001.fits  '+err_obj1d)


def salt_scombine(objlist):

    obj = np.loadtxt(objlist,dtype="str",unpack=True, ndmin=1)
    objlist1d=objlist+'1d'
    objj = os.path.splitext(obj[0])[0]
    if os.path.isfile(objj+'wsf.fits'):
        os.system('sed s/.fits/wsf1d.fits/ '+objlist+' > '+objlist1d)
        objnofits=string.split(obj[0],sep='.')[0]+'wsf'
        os.system('sed s/.fits/wsf1d_var.fits/ '+objlist+' > templist')        
    else:
        os.system('sed s/.fits/ws1d.fits/ '+objlist+' > '+objlist1d)           
        objnofits=string.split(obj[0],sep='.')[0]+'ws'
        os.system('sed s/.fits/ws1d_var.fits/ '+objlist+' > templist')        

    low_ap= float(rt.getsh("grep -i 'low\t' database/ap"+objnofits+" |  awk '{print $3}'"))
    high_ap= float(rt.getsh("grep -i 'high\t' database/ap"+objnofits+" |  awk '{print $3}'"))

    apradius = abs(low_ap)+abs(high_ap)
    readnoise = 2.45 * np.sqrt(apradius)
    print "ggg = %s" %readnoise

    iraf.scombine.combine='average'
    iraf.scombine.rdnoise=readnoise
    iraf.scombine.reject='ccdclip'
    iraf.scombine.scale='median'

    outnam=''.join(rt.header(obj[0],'OBJECT').split()).lower()
    outname=outnam+'_'+objlist
    if os.path.isfile(outname+'.fits'):
       	 os.system('rm '+outname+'.fits')

    inname = str('@'+objlist1d)

    iraf.scombine(input=inname,output = outname,first='yes')

    tempfile='temp.fits'
    if os.path.isfile(tempfile):
       	 os.system('rm '+tempfile)
 
    iraf.imsum('@templist', tempfile)
    nspec=len(obj)
    outname_var=outname+'_var'
    if os.path.isfile(outname_var+'.fits'):
       	 os.system('rm '+outname_var+'.fits')

    nspec=nspec*nspec

    iraf.imarith(tempfile,'/',nspec,outname_var)

    os.system('rm *temp.fits')
    os.system('rm templist')
    iraf.dispcor(outname_var,output='',table=outname)
    outnametxt=outname+'.txt'
    outname_vartxt=outname_var+'.txt'
    iraf.wspectext(outname,outnametxt,header='no')
    iraf.wspectext(outname_var,outname_vartxt,header='no')

    rssdate = rt.header(obj[0],'DATE-OBS')
    rssdate = string.replace(rssdate,'-','')
    pastefinal=outname+'_final_'+rssdate+'.txt'
    os.system('paste '+outnametxt+' '+outname_vartxt+'  > '+pastefinal)


