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
    
    for idx_x in range(0,2):
        delta_x = 0.01*float(idx_x)
        output_file = open('analysis/sky_det_ref/%s_target_ra_%s.txt'%(detector[det],idx_x), "w")
        output_file.write('%s\n'%(float(hdu_header['RA_OBJ'])-0.+delta_x))
        output_file.close()
        output_file = open('analysis/sky_det_ref/%s_target_dec_%s.txt'%(detector[det],idx_x), "w")
        output_file.write('%s\n'%(float(hdu_header['DEC_OBJ'])-0.+delta_x))
        output_file.close()
