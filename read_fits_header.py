import sys
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import astropy.utils as utils
from astropy.nddata import Cutout2D
from astropy.table import Table
import numpy as np

detector = []
detector += ['mos1']
detector += ['mos2']

for det in range(0,len(detector)):
    filename = 'analysis/%s-fov-img.fits'%(detector[det])
    hdu_list = fits.open(filename)
    print (filename)
    hdu_header = hdu_list[0].header
    
    print ('[deg] RA of target = %s'%(hdu_header['RA_OBJ']))
    print ('[deg] DEC of target = %s'%(hdu_header['DEC_OBJ']))
    
    output_file = open('analysis/%s_target_ra.txt'%(detector[det]), "w")
    output_file.write('%s\n'%(hdu_header['RA_OBJ']))
    output_file.close()
    output_file = open('analysis/%s_target_ra_offset.txt'%(detector[det]), "w")
    output_file.write('%s\n'%(float(hdu_header['RA_OBJ'])+0.2))
    output_file.close()
    output_file = open('analysis/%s_target_dec.txt'%(detector[det]), "w")
    output_file.write('%s\n'%(hdu_header['DEC_OBJ']))
    output_file.close()
    output_file = open('analysis/%s_target_dec_offset.txt'%(detector[det]), "w")
    output_file.write('%s\n'%(float(hdu_header['DEC_OBJ'])+0.2))
    output_file.close()
