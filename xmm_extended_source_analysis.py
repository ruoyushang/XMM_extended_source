
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import astropy.utils as utils
from astropy.nddata import Cutout2D
from astropy.table import Table
import numpy as np
from numpy.linalg import inv
from operator import itemgetter, attrgetter
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import time

fig, ax = plt.subplots()
figsize_x = 8.
figsize_y = 6.
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

on_sample = sys.argv[1]
on_obsID = sys.argv[2]

#diagnostic_plots = True
diagnostic_plots = False

#detector = 'mos1S001'
detector = 'mos2S002'

energy_array = [2000,3000,5000,8000,12000]
#energy_array = [2000,12000]
ana_ccd_bins = [1,2,3,4,5,6,7]
#ana_ccd_bins = [0]

ch_low = 200
ch_high = 12000
ch_scale = 200
t_low = 0
t_high = 1
t_scale = 0.05
detx_low = -20000
detx_high = 20000
detx_scale = 500
detr_low = 0
detr_high = 20000
detr_scale = 500

src_ra = 0
src_dec = 0
if 'Cas_A' in on_sample:
    src_ra = 350.85
    src_dec = 58.84

def get_conversion_mtx():
    radec_filename = '../%s/%s/analysis/pointing_coord.txt'%(on_sample,on_obsID)
    detxy_filename = '../%s/%s/analysis/det_coord_%s.txt'%(on_sample,on_obsID,detector)
    ref_ra = []
    ref_dec = []
    radec_file = open(radec_filename)
    for line in radec_file:
        sky_ra = line.split()[1]
        sky_dec = line.split()[2]
        ref_ra += [float(sky_ra)]
        ref_dec += [float(sky_dec)]
    ref_detx = []
    ref_dety = []
    detxy_file = open(detxy_filename)
    for line in detxy_file:
        det_x = line.split()[1]
        det_y = line.split()[2]
        ref_detx += [float(det_x)]
        ref_dety += [float(det_y)]
    delta_ra = []
    delta_dec = []
    delta_detx = []
    delta_dety = []
    for point in range(1,len(ref_ra)):
        delta_ra += [ref_ra[point]-ref_ra[0]]
        delta_dec += [ref_dec[point]-ref_dec[0]]
        delta_detx += [ref_detx[point]-ref_detx[0]]
        delta_dety += [ref_dety[point]-ref_dety[0]]
    mtx_delta_sky = np.matrix([[delta_ra[0],delta_dec[0]],[delta_ra[1],delta_dec[1]]])
    mtx_delta_det = np.matrix([[delta_detx[0],delta_dety[0]],[delta_detx[1],delta_dety[1]]])
    att_ra = ref_ra[0]
    att_dec = ref_dec[0]
    att_detx = ref_detx[0]
    att_dety = ref_dety[0]
    inv_mtx_delta_det = inv(mtx_delta_det)
    inv_mtx_delta_sky = inv(mtx_delta_sky)
    mtx_conv_sky_2_det = np.matmul(mtx_delta_det,inv_mtx_delta_sky)
    mtx_conv_det_2_sky = np.matmul(mtx_delta_sky,inv_mtx_delta_det)
    att_detxy = np.array([att_detx,att_dety])
    delta_att_radec = np.matmul(mtx_conv_det_2_sky,att_detxy)
    att_radec = [delta_att_radec[0,0]+ref_ra[0], delta_att_radec[0,1]+ref_dec[0]]

    return att_radec, mtx_conv_sky_2_det, mtx_conv_det_2_sky

global_att_radec, global_mtx_conv_sky_2_det, global_mtx_conv_det_2_sky = get_conversion_mtx()

print ('global_att_radec:')
print (global_att_radec)

pix_size = (1./60.)/1200. # degree
if not 'extragalactic' in on_sample:
    detr_low = pow(pow(src_ra-global_att_radec[0],2)+pow(src_dec-global_att_radec[1],2),0.5)-20000*pix_size
    detr_low = max(0.,detr_low)
    #detr_high = pow(pow(src_ra-global_att_radec[0],2)+pow(src_dec-global_att_radec[1],2),0.5)+20000*pix_size
    detr_high = detr_low+60000*pix_size
    detr_scale = detr_scale*pix_size

fwc_2_sci_ratio = 0.4

class MyArray2D:

    def __init__(self):
        start_x = -19474.5
        start_y = -19474.5
        pixel_scale = 500
        image_size = 39000
        nrows = int(image_size/pixel_scale)
        ncols = int(image_size/pixel_scale)
        array_shape = (nrows,ncols)

        self.xaxis = np.zeros(array_shape[0])
        self.yaxis = np.zeros(array_shape[1])
        self.zaxis = np.zeros(array_shape)
        self.zerr = np.zeros(array_shape)
        for idx in range(0,len(self.xaxis)):
            self.xaxis[idx] = start_x + idx*pixel_scale
        for idx in range(0,len(self.yaxis)):
            self.yaxis[idx] = start_y + idx*pixel_scale

    def reset(self):
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                self.zaxis[idx_x,idx_y] = 0.
                self.zerr[idx_x,idx_y] = 0.
    def scale(self, factor):
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                self.zaxis[idx_x,idx_y] = self.zaxis[idx_x,idx_y]*factor
                self.zerr[idx_x,idx_y] = self.zerr[idx_x,idx_y]*factor
    def add(self, add_array, factor=1.):
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                self.zaxis[idx_x,idx_y] = self.zaxis[idx_x,idx_y]+add_array.zaxis[idx_x,idx_y]*factor
                self.zerr[idx_x,idx_y] = pow(pow(self.zerr[idx_x,idx_y],2)+pow(add_array.zerr[idx_x,idx_y]*factor,2),0.5)
    def fill(self, value_x, value_y):
        key_idx_x = 0
        key_idx_y = 0
        for idx_x in range(0,len(self.xaxis)-1):
            if self.xaxis[idx_x]<=value_x and self.xaxis[idx_x+1]>value_x:
                key_idx_x = idx_x
        for idx_y in range(0,len(self.yaxis)-1):
            if self.yaxis[idx_y]<=value_y and self.yaxis[idx_y+1]>value_y:
                key_idx_y = idx_y
        self.zaxis[key_idx_x,key_idx_y] += 1.
    def get_bin_content(self, value_x, value_y):
        key_idx_x = 0
        key_idx_y = 0
        for idx_x in range(0,len(self.xaxis)-1):
            if self.xaxis[idx_x]<=value_x and self.xaxis[idx_x+1]>value_x:
                key_idx_x = idx_x
        for idx_y in range(0,len(self.yaxis)-1):
            if self.yaxis[idx_y]<=value_y and self.yaxis[idx_y+1]>value_y:
                key_idx_y = idx_y
        return self.zaxis[key_idx_x,key_idx_y]
    def get_pixel_fraction(self):
        pix_on_cnt = 0
        pix_off_cnt = 0
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                pix_off_cnt += 1
                if self.zaxis[idx_x,idx_y]!=0:
                    pix_on_cnt += 1
        return pix_on_cnt/pix_off_cnt
    def get_error(self):
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                self.zerr[idx_x,idx_y] = pow(max(self.zaxis[idx_x,idx_y],1.),0.5)
    def integral(self):
        total_integral = 0.
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                cnt = self.zaxis[idx_x,idx_y]
                total_integral += cnt
        return total_integral
    def calc_significance(self, src_array, bkg_array):
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                src = src_array.zaxis[idx_x,idx_y]
                bkg = bkg_array.zaxis[idx_x,idx_y]
                if not src>0.: continue
                if not bkg>0.: continue
                alpha = 1.
                z_score = pow(2,0.5)*pow(src*math.log((1+alpha)/alpha*(src/(src+bkg))) + bkg*math.log((1+alpha)*bkg/(src+bkg)),0.5)
                if (src-bkg)>0.:
                    self.zaxis[idx_x,idx_y] = z_score
                else:
                    self.zaxis[idx_x,idx_y] = -z_score

class MyArray1D:

    def __init__(self,bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale):
        nchannels = int((bin_end-bin_start)/pixel_scale)
        array_shape = (nchannels)
        self.yaxis = np.zeros(array_shape)
        self.xaxis = np.zeros(array_shape)
        self.yerr = np.zeros(array_shape)
        self.xerr = np.zeros(array_shape)
        for idx in range(0,len(self.xaxis)):
            self.xaxis[idx] = bin_start + idx*pixel_scale

    def reset(self):
        for entry in range(0,len(self.yaxis)):
            self.yaxis[entry] = 0.
            self.yerr[entry] = 0.
    def set_floor(self):
        for entry in range(0,len(self.yaxis)):
            self.yaxis[entry] = max(0.,self.yaxis[entry])
    def scale(self, factor):
        for entry in range(0,len(self.yaxis)):
            self.yaxis[entry] = self.yaxis[entry]*factor
            self.yerr[entry] = self.yerr[entry]*factor
    def add(self, add_array, factor=1.):
        for entry in range(0,len(self.yaxis)):
            self.yaxis[entry] = self.yaxis[entry]+add_array.yaxis[entry]*factor
            self.yerr[entry] = pow(pow(self.yerr[entry],2)+pow(add_array.yerr[entry]*factor,2),0.5)
    def divide(self, target_array, factor=1.):
        for entry in range(0,len(self.yaxis)):
            if target_array.yaxis[entry]!=0.:
                self.yaxis[entry] = self.yaxis[entry]/(target_array.yaxis[entry]*factor)
                self.yerr[entry] = self.yerr[entry]/(target_array.yaxis[entry]*factor)
            else:
                self.yaxis[entry] = 0.
                self.yerr[entry] = 0.
    def multiply(self, target_array, factor=1.):
        for entry in range(0,len(self.yaxis)):
                self.yaxis[entry] = self.yaxis[entry]*(target_array.yaxis[entry]*factor)
                self.yerr[entry] = self.yerr[entry]*(target_array.yaxis[entry]*factor)
    def integral(self,integral_range=[]):
        integral_cnt = 0.
        if integral_range==[]:
            integral_range=[self.xaxis[0],self.xaxis[len(self.yaxis)-1]]
        for entry in range(0,len(self.xaxis)):
            if self.xaxis[entry]<integral_range[0]: continue
            if self.xaxis[entry]>integral_range[1]: continue
            integral_cnt += self.yaxis[entry]
        return integral_cnt
    def fill(self, value):
        for entry in range(0,len(self.xaxis)-1):
            if self.xaxis[entry]<=value and self.xaxis[entry+1]>value:
                self.yaxis[entry] += 1.
    def get_pixel_fraction(self):
        pix_on_cnt = 0
        pix_off_cnt = 0
        for entry in range(0,len(self.yaxis)):
            pix_off_cnt += 1
            if self.yaxis[entry]!=0:
                pix_on_cnt += 1
        return pix_on_cnt/pix_off_cnt
    def get_error(self):
        for entry in range(0,len(self.yaxis)):
            self.yerr[entry] = pow(max(self.yaxis[entry],1.),0.5)
            self.xerr[entry] = 0.5*(self.xaxis[1]-self.xaxis[0])
    def get_bin_content(self, value_x):
        key_idx_x = 0
        for idx_x in range(0,len(self.xaxis)-1):
            if self.xaxis[idx_x]<=value_x and self.xaxis[idx_x+1]>value_x:
                key_idx_x = idx_x
        return self.yaxis[key_idx_x]
    def normalize(self):
        reference = self.yaxis[0]
        for entry in range(0,len(self.yaxis)):
            self.yaxis[entry] = self.yaxis[entry]/reference
            self.yerr[entry] = self.yerr[entry]/reference
    def calc_significance(self, src_array, bkg_array):
        for idx_x in range(0,len(self.xaxis)):
            src = src_array.yaxis[idx_x]
            bkg = bkg_array.yaxis[idx_x]
            if not src>0.: continue
            if not bkg>0.: continue
            alpha = 1.
            z_score = pow(2,0.5)*pow(src*math.log((1+alpha)/alpha*(src/(src+bkg))) + bkg*math.log((1+alpha)*bkg/(src+bkg)),0.5)
            if (src-bkg)>0.:
                self.yaxis[idx_x] = z_score
            else:
                self.yaxis[idx_x] = -z_score

def read_event_file(filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[200,12000],ccd_id=0):

    # how to read events:
    # https://docs.astropy.org/en/stable/generated/examples/io/fits-tables.html#accessing-data-stored-as-a-table-in-a-multi-extension-fits-mef-file

    hdu_list = fits.open(filename)
    #print (filename)
    #print (hdu_list.info())

    events = Table.read(filename, hdu=1)
    #print('total events:')
    #print(len(events))
    #print('events.columns:')
    #print(events.columns)

    evt_count = 0.
    lightcurve_array = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
    pattern_array = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
    spectrum_array = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    detx_array = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
    image_array = MyArray2D()
    detr_array = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)

    time_start = events[0]['TIME']
    time_end = events[len(events)-1]['TIME']
    obs_duration = time_end-time_start

    for evt in range(0,len(events)):

        if (evt%10)!=0: continue

        evt_time = events[evt]['TIME']
        evt_pattern = events[evt]['PATTERN']
        evt_pi = events[evt]['PI']

        evt_detx = events[evt]['DETX']
        evt_dety = events[evt]['DETY']
        evt_detr = pow(pow(evt_detx,2)+pow(evt_dety,2),0.5)
        if not 'extragalactic' in on_sample:
            evt_detxy = np.array([evt_detx,evt_dety])
            evt_xy = np.matmul(global_mtx_conv_det_2_sky,evt_detxy)
            evt_x = evt_xy[0,0]+global_att_radec[0]
            evt_y = evt_xy[0,1]+global_att_radec[1]
            evt_detr = pow(pow(evt_x-src_ra,2)+pow(evt_y-src_dec,2),0.5)

        if evt_pi<energy_range[0]: continue
        if evt_pi>=energy_range[1]: continue

        # if 'sp-free' in evt_filter:
        #if evt_pattern==0:continue
        #if evt_pattern==1:continue
        #if evt_pattern==3:continue

        #if evt_pattern<=0:continue
        #if evt_pattern>12:continue
        if evt_pattern>4:continue

        #if 'Al' in evt_filter:
        #    if evt_pi<1200: continue
        #    if evt_pi>1900: continue
        #if evt_pi>1200 and evt_pi<1900: continue
        #if evt_pi<1900: continue

        if 'source' in evt_filter:
            #if abs(evt_detx)>1000: continue
            #if abs(evt_dety)>1000: continue
            if abs(evt_detx)>6000: continue
            if abs(evt_dety)>6000: continue
        if 'ring' in evt_filter:
            #if abs(evt_detx)<1500 and abs(evt_dety)<1500: continue
            #if abs(evt_detx)>6000: continue
            #if abs(evt_dety)>6000: continue
            if abs(evt_detx)<6000 and abs(evt_dety)<6000: continue
            if abs(evt_detx)>12000: continue
            if abs(evt_dety)>12000: continue

        if not mask_lc==None:
            zscore = mask_lc.get_bin_content((evt_time-time_start)/(time_end-time_start))
            if 'sp-veto' in evt_filter:
                if zscore>0: continue
            if 'sp-select' in evt_filter:
                if zscore==0: continue
        if not mask_map==None:
            zscore = mask_map.get_bin_content(evt_detx,evt_dety)
            if zscore>=2.0: continue
            if zscore<=-1.0: continue

        evt_ccd_id = 0
        if evt_detx<7000 and evt_detx>=-7000 and evt_dety<7000 and evt_dety>=-7000:
            evt_ccd_id = 1
        if evt_detx>=7000 and evt_dety<7000 and evt_dety>=-7000:
            evt_ccd_id = 3
        if evt_detx<-7000 and evt_dety<7000 and evt_dety>=-7000:
            evt_ccd_id = 6
        if evt_detx>=0 and evt_dety>=7000:
            evt_ccd_id = 4
        if evt_detx<0 and evt_dety>=7000:
            evt_ccd_id = 5
        if evt_detx>=0 and evt_dety<-7000:
            evt_ccd_id = 2
        if evt_detx<0 and evt_dety<-7000:
            evt_ccd_id = 7

        if not 'cor-evt' in filename and not ccd_id==0:
            if evt_ccd_id!=ccd_id: continue

        evt_count += 1.
        lightcurve_array.fill((evt_time-time_start)/(time_end-time_start))
        pattern_array.fill(evt_pattern)
        spectrum_array.fill(evt_pi)
        detx_array.fill(evt_detx)
        image_array.fill(evt_detx,evt_dety)
        detr_array.fill(evt_detr)

    return [obs_duration, evt_count, lightcurve_array, pattern_array, spectrum_array, detx_array, image_array, detr_array]

def fit_pattern(data_pattern,xray_pattern,spf_pattern,qpb_pattern):

    sp_free_data = 0.
    sp_free_xray = 0.
    sp_free_qpb = 0.
    for idx in range(0,len(data_pattern.xaxis)):
        if idx!=2: continue
        if idx>4: continue
        sp_free_data += data_pattern.yaxis[idx]
        sp_free_xray += xray_pattern.yaxis[idx]
        sp_free_qpb  += qpb_pattern.yaxis[idx]

    xray_scale = 0.
    if sp_free_xray>0.:
        xray_scale = (sp_free_data-sp_free_qpb)/sp_free_xray
    xray_scale = max(0.,xray_scale)

    sp_pattern_data = 0.
    sp_pattern_xray = 0.
    sp_pattern_qpb = 0.
    sp_pattern_spf = 0.
    for idx in range(0,len(data_pattern.xaxis)):
        if idx==2: continue
        if idx>4: continue
        sp_pattern_data += data_pattern.yaxis[idx]
        sp_pattern_xray += xray_pattern.yaxis[idx]
        sp_pattern_qpb  += qpb_pattern.yaxis[idx]
        sp_pattern_spf  += spf_pattern.yaxis[idx]
    sp_pattern_xray = sp_pattern_xray*xray_scale
    
    sp_scale = 0.
    if sp_pattern_spf>0.:
        sp_scale = (sp_pattern_data-sp_pattern_xray-sp_pattern_qpb)/sp_pattern_spf
    sp_scale = max(0.,sp_scale)


    # Define the function to minimize (residuals)
    def residuals(A):
        return (data_pattern.yaxis - A[0]*xray_pattern.yaxis - A[1]*spf_pattern.yaxis - qpb_pattern.yaxis) / data_pattern.yerr  # Include the weighted residuals
        #return (data_pattern.yaxis - A[0]*xray_pattern.yaxis - A[1]*spf_pattern.yaxis - qpb_pattern.yaxis)  # Include the weighted residuals

    data_pattern.get_error()
    # Perform the least squares fitting
    result = least_squares(residuals, x0=[xray_scale,sp_scale], bounds=([0.5*xray_scale,0.5*sp_scale], [2.*(xray_scale+1e-4),2.*(sp_scale+1e-4)]))  # x0 is the initial guess for A

    xray_scale = result.x[0]
    sp_scale = result.x[1]


    return sp_scale, xray_scale
    
def make_timecut_mask(lightcurve_data,lightcurve_bkgd):

    lightcurve_sci_fov_mask = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
    #lightcurve_sci_fov_mask.calc_significance(lightcurve_data,lightcurve_bkgd)
    lightcurve_sci_fov_mask.add(lightcurve_data)
    lightcurve_sci_fov_mask.divide(lightcurve_bkgd)
    for idx in range(0,len(lightcurve_sci_fov_mask.xaxis)):
        if lightcurve_sci_fov_mask.yaxis[idx]>1.5:
            lightcurve_sci_fov_mask.yaxis[idx] = 1.
        else:
            lightcurve_sci_fov_mask.yaxis[idx] = 0.
    
    return lightcurve_sci_fov_mask

start_time = time.time()

map_color = 'coolwarm'


SP_flare_mask = True
#point_source_mask = True
#SP_flare_mask = False
point_source_mask = False

off_sample = 'extragalactic'
off_obsID = 'ID0827241101' # Percentage of flaring time: 15.2%
#off_sample = on_sample
#off_obsID = on_obsID

#xray_sample = 'extragalactic'
#xray_obsID = 'ID0690900101'
#xray_obsID = 'ID0827251001'
xray_sample = 'Cas_A'
xray_obsID = 'ID0412180101'

plot_tag = detector+'_'+on_obsID

xray_sci_fov_evt_filename = '../%s/%s/analysis/%s-fov-evt.fits'%(xray_sample,xray_obsID,detector)
xray_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-fov-evt.fits'%(xray_sample,xray_obsID,detector)
xray_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(xray_sample,xray_obsID,detector)
xray_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(xray_sample,xray_obsID,detector)

off_sci_fov_evt_filename = '../%s/%s/analysis/%s-fov-evt.fits'%(off_sample,off_obsID,detector)
off_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-fov-evt.fits'%(off_sample,off_obsID,detector)
off_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(off_sample,off_obsID,detector)
off_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(off_sample,off_obsID,detector)

on_sci_fov_evt_filename = '../%s/%s/analysis/%s-fov-evt.fits'%(on_sample,on_obsID,detector)
on_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-fov-evt.fits'%(on_sample,on_obsID,detector)
on_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(on_sample,on_obsID,detector)
on_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(on_sample,on_obsID,detector)


def analyze_a_ccd_chip(energy_range=[200,12000],ccd_id=0):

    print ('prepare xray sample time cuts')
    
    output_all_fov = read_event_file(xray_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    output_all_cor = read_event_file(xray_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)

    xray_duration_all_fov = output_all_fov[0]
    xray_evt_count_all_fov = output_all_fov[1]
    xray_lightcurve_all_fov = output_all_fov[2]
    xray_pattern_all_fov = output_all_fov[3]
    xray_spectrum_all_fov = output_all_fov[4]
    xray_detx_all_fov = output_all_fov[5]
    xray_image_all_fov = output_all_fov[6]
    
    xray_duration_all_cor = output_all_cor[0]
    xray_evt_count_all_cor = output_all_cor[1]
    xray_lightcurve_all_cor = output_all_cor[2]
    xray_pattern_all_cor = output_all_cor[3]
    xray_spectrum_all_cor = output_all_cor[4]
    xray_detx_all_cor = output_all_cor[5]
    xray_image_all_cor = output_all_cor[6]

    area_pix_frac_all_fov = xray_image_all_fov.get_pixel_fraction()
    area_pix_frac_all_cor = xray_image_all_cor.get_pixel_fraction()
    
    time_pix_frac_all_fov = xray_lightcurve_all_fov.get_pixel_fraction()
    time_expo_all_fov = xray_duration_all_fov*time_pix_frac_all_fov
    
    time_pix_frac_all_cor = xray_lightcurve_all_cor.get_pixel_fraction()
    time_expo_all_cor = xray_duration_all_cor*time_pix_frac_all_cor
    
    xray_lightcurve_all_cor.scale(area_pix_frac_all_fov/area_pix_frac_all_cor)

    mask_lc = make_timecut_mask(xray_lightcurve_all_fov,xray_lightcurve_all_cor) 

    time_pix_frac_mask = mask_lc.get_pixel_fraction()
    if time_pix_frac_mask>0.8:
        mask_lc = None

    print ('apply x-ray sample space and time masks')
    
    output_sci_source = read_event_file(xray_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto source',energy_range=energy_range,ccd_id=0)
    output_sci_ring   = read_event_file(xray_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto ring',energy_range=energy_range,ccd_id=0)
    
    xray_duration_sci_source = output_sci_source[0]
    xray_evt_count_sci_source = output_sci_source[1]
    xray_lightcurve_sci_source = output_sci_source[2]
    xray_pattern_sci_source = output_sci_source[3]
    xray_spectrum_sci_source = output_sci_source[4]
    xray_detx_sci_source = output_sci_source[5]
    xray_image_sci_source = output_sci_source[6]
    
    xray_duration_sci_ring = output_sci_ring[0]
    xray_evt_count_sci_ring = output_sci_ring[1]
    xray_lightcurve_sci_ring = output_sci_ring[2]
    xray_pattern_sci_ring = output_sci_ring[3]
    xray_spectrum_sci_ring = output_sci_ring[4]
    xray_detx_sci_ring = output_sci_ring[5]
    xray_image_sci_ring = output_sci_ring[6]

    area_pix_frac_sci_source = xray_image_sci_source.get_pixel_fraction()
    area_pix_frac_sci_ring = xray_image_sci_ring.get_pixel_fraction()
    xray_pattern_sci_ring.scale(area_pix_frac_sci_source/area_pix_frac_sci_ring)

    xray_pattern_template = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
    xray_pattern_template.add(xray_pattern_sci_source)
    xray_pattern_template.add(xray_pattern_sci_ring,factor=-1.)

    print ('prepare off sample time cuts')
    
    output_all_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    output_all_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    
    off_duration_all_fov = output_all_fov[0]
    off_evt_count_all_fov = output_all_fov[1]
    off_lightcurve_all_fov = output_all_fov[2]
    off_pattern_all_fov = output_all_fov[3]
    off_spectrum_all_fov = output_all_fov[4]
    off_detx_all_fov = output_all_fov[5]
    off_image_all_fov = output_all_fov[6]
    
    off_duration_all_cor = output_all_cor[0]
    off_evt_count_all_cor = output_all_cor[1]
    off_lightcurve_all_cor = output_all_cor[2]
    off_pattern_all_cor = output_all_cor[3]
    off_spectrum_all_cor = output_all_cor[4]
    off_detx_all_cor = output_all_cor[5]
    off_image_all_cor = output_all_cor[6]
    
    off_duration_fwc_fov = output_fwc_fov[0]
    off_evt_count_fwc_fov = output_fwc_fov[1]
    off_lightcurve_fwc_fov = output_fwc_fov[2]
    off_pattern_fwc_fov = output_fwc_fov[3]
    off_spectrum_fwc_fov = output_fwc_fov[4]
    off_detx_fwc_fov = output_fwc_fov[5]
    off_image_fwc_fov = output_fwc_fov[6]
    
    off_duration_fwc_cor = output_fwc_cor[0]
    off_evt_count_fwc_cor = output_fwc_cor[1]
    off_lightcurve_fwc_cor = output_fwc_cor[2]
    off_pattern_fwc_cor = output_fwc_cor[3]
    off_spectrum_fwc_cor = output_fwc_cor[4]
    off_detx_fwc_cor = output_fwc_cor[5]
    off_image_fwc_cor = output_fwc_cor[6]
    
    area_pix_frac_all_fov = off_image_all_fov.get_pixel_fraction()
    area_pix_frac_all_cor = off_image_all_cor.get_pixel_fraction()
    area_pix_frac_fwc_fov = off_image_fwc_fov.get_pixel_fraction()
    area_pix_frac_fwc_cor = off_image_fwc_cor.get_pixel_fraction()
    
    time_pix_frac_all_fov = off_lightcurve_all_fov.get_pixel_fraction()
    time_expo_all_fov = off_duration_all_fov*time_pix_frac_all_fov
    time_pix_frac_fwc_fov = off_lightcurve_fwc_fov.get_pixel_fraction()
    time_expo_fwc_fov = off_duration_fwc_fov*time_pix_frac_fwc_fov
    
    time_pix_frac_all_cor = off_lightcurve_all_cor.get_pixel_fraction()
    time_expo_all_cor = off_duration_all_cor*time_pix_frac_all_cor
    time_pix_frac_fwc_cor = off_lightcurve_fwc_cor.get_pixel_fraction()
    time_expo_fwc_cor = off_duration_fwc_cor*time_pix_frac_fwc_cor
    
    off_lightcurve_all_cor.scale(area_pix_frac_all_fov/area_pix_frac_all_cor)
    off_lightcurve_fwc_fov.scale(fwc_2_sci_ratio*time_expo_all_fov/time_expo_fwc_fov)

    mask_lc = make_timecut_mask(off_lightcurve_all_fov,off_lightcurve_all_cor) 

    print ('apply off sample space and time masks')
    
    output_sci_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    output_sci_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    
    off_duration_sci_fov = output_sci_fov[0]
    off_evt_count_sci_fov = output_sci_fov[1]
    off_lightcurve_sci_fov = output_sci_fov[2]
    off_pattern_sci_fov = output_sci_fov[3]
    off_spectrum_sci_fov = output_sci_fov[4]
    off_detx_sci_fov = output_sci_fov[5]
    off_image_sci_fov = output_sci_fov[6]
    
    off_duration_sci_cor = output_sci_cor[0]
    off_evt_count_sci_cor = output_sci_cor[1]
    off_lightcurve_sci_cor = output_sci_cor[2]
    off_pattern_sci_cor = output_sci_cor[3]
    off_spectrum_sci_cor = output_sci_cor[4]
    off_detx_sci_cor = output_sci_cor[5]
    off_image_sci_cor = output_sci_cor[6]
    
    output_spf_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_range,ccd_id=ccd_id)
    output_spf_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_range,ccd_id=ccd_id)
    
    off_duration_spf_fov = output_spf_fov[0]
    off_evt_count_spf_fov = output_spf_fov[1]
    off_lightcurve_spf_fov = output_spf_fov[2]
    off_pattern_spf_fov = output_spf_fov[3]
    off_spectrum_spf_fov = output_spf_fov[4]
    off_detx_spf_fov = output_spf_fov[5]
    off_image_spf_fov = output_spf_fov[6]
    
    off_duration_spf_cor = output_spf_cor[0]
    off_evt_count_spf_cor = output_spf_cor[1]
    off_lightcurve_spf_cor = output_spf_cor[2]
    off_pattern_spf_cor = output_spf_cor[3]
    off_spectrum_spf_cor = output_spf_cor[4]
    off_detx_spf_cor = output_spf_cor[5]
    off_image_spf_cor = output_spf_cor[6]
    
    output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    
    off_duration_fwc_fov = output_fwc_fov[0]
    off_evt_count_fwc_fov = output_fwc_fov[1]
    off_lightcurve_fwc_fov = output_fwc_fov[2]
    off_pattern_fwc_fov = output_fwc_fov[3]
    off_spectrum_fwc_fov = output_fwc_fov[4]
    off_detx_fwc_fov = output_fwc_fov[5]
    off_image_fwc_fov = output_fwc_fov[6]
    
    off_duration_fwc_cor = output_fwc_cor[0]
    off_evt_count_fwc_cor = output_fwc_cor[1]
    off_lightcurve_fwc_cor = output_fwc_cor[2]
    off_pattern_fwc_cor = output_fwc_cor[3]
    off_spectrum_fwc_cor = output_fwc_cor[4]
    off_detx_fwc_cor = output_fwc_cor[5]
    off_image_fwc_cor = output_fwc_cor[6]

    area_pix_frac_sci_fov = off_image_sci_fov.get_pixel_fraction()
    area_pix_frac_sci_cor = off_image_sci_cor.get_pixel_fraction()
    area_pix_frac_spf_fov = off_image_spf_fov.get_pixel_fraction()
    area_pix_frac_spf_cor = off_image_spf_cor.get_pixel_fraction()
    area_pix_frac_fwc_fov = off_image_fwc_fov.get_pixel_fraction()
    area_pix_frac_fwc_cor = off_image_fwc_cor.get_pixel_fraction()
    
    time_pix_frac_sci_fov = off_lightcurve_sci_fov.get_pixel_fraction()
    time_expo_sci_fov = off_duration_sci_fov*time_pix_frac_sci_fov
    time_pix_frac_sci_cor = off_lightcurve_sci_cor.get_pixel_fraction()
    time_expo_sci_cor = off_duration_sci_cor*time_pix_frac_sci_cor
    
    time_pix_frac_spf_fov = off_lightcurve_spf_fov.get_pixel_fraction()
    time_expo_spf_fov = off_duration_spf_fov*time_pix_frac_spf_fov
    time_pix_frac_spf_cor = off_lightcurve_spf_cor.get_pixel_fraction()
    time_expo_spf_cor = off_duration_spf_cor*time_pix_frac_spf_cor
    
    time_pix_frac_fwc_fov = off_lightcurve_fwc_fov.get_pixel_fraction()
    time_expo_fwc_fov = off_duration_fwc_fov*time_pix_frac_fwc_fov
    time_pix_frac_fwc_cor = off_lightcurve_fwc_cor.get_pixel_fraction()
    time_expo_fwc_cor = off_duration_fwc_cor*time_pix_frac_fwc_cor
    
    off_lightcurve_sci_cor.scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    off_lightcurve_spf_fov.scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
    off_lightcurve_spf_cor.scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
    off_lightcurve_fwc_fov.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    off_lightcurve_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
    off_spectrum_sci_cor.scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    off_spectrum_spf_fov.scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
    off_spectrum_spf_cor.scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
    off_spectrum_fwc_fov.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    off_spectrum_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
    off_pattern_sci_cor.scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    off_pattern_spf_fov.scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
    off_pattern_spf_cor.scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
    off_pattern_fwc_fov.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    off_pattern_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)

    if diagnostic_plots:
        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_lightcurve_all_fov.xaxis,off_lightcurve_all_fov.yaxis,yerr=off_lightcurve_all_fov.yerr,color='k',label='All FoV')
        axbig.errorbar(off_lightcurve_sci_fov.xaxis,off_lightcurve_sci_fov.yaxis,yerr=off_lightcurve_sci_fov.yerr,color='red',label='Sci FoV')
        axbig.errorbar(off_lightcurve_all_cor.xaxis,off_lightcurve_all_cor.yaxis,yerr=off_lightcurve_all_cor.yerr,color='green',label='All Cor')
        axbig.set_yscale('log')
        axbig.set_xlabel('Time')
        axbig.legend(loc='best')
        fig.savefig("../output_plots/lightcurve_off_timecut_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_spectrum_sci_fov.xaxis,off_spectrum_sci_fov.yaxis,yerr=off_spectrum_sci_fov.yerr,label='Sci FoV')
        axbig.errorbar(off_spectrum_spf_fov.xaxis,off_spectrum_spf_fov.yaxis,yerr=off_spectrum_spf_fov.yerr,label='SPF FoV')
        axbig.errorbar(off_spectrum_fwc_fov.xaxis,off_spectrum_fwc_fov.yaxis,yerr=off_spectrum_fwc_fov.yerr,label='FWC FoV')
        axbig.set_xlabel('Channel')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("../output_plots/off_spectrum_fov_raw_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_pattern_sci_fov.xaxis,off_pattern_sci_fov.yaxis,yerr=off_pattern_sci_fov.yerr,label='Sci FoV')
        axbig.errorbar(off_pattern_spf_fov.xaxis,off_pattern_spf_fov.yaxis,yerr=off_pattern_spf_fov.yerr,label='SPF FoV')
        axbig.errorbar(off_pattern_fwc_fov.xaxis,off_pattern_fwc_fov.yaxis,yerr=off_pattern_fwc_fov.yerr,label='FWC FoV')
        axbig.set_xlabel('Channel')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("../output_plots/off_pattern_fov_raw_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

    sp_pattern_template = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
    sp_spectrum_template = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    
    sp_pattern_template.add(off_pattern_spf_fov)
    sp_pattern_template.add(off_pattern_sci_fov,factor=-1.)
    sp_spectrum_template.add(off_spectrum_spf_fov)
    sp_spectrum_template.add(off_spectrum_sci_fov,factor=-1.)

    print ('prepare on sample time cuts')
    
    output_all_fov = read_event_file(on_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    output_all_cor = read_event_file(on_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',ccd_id=0)
    
    on_duration_all_fov = output_all_fov[0]
    on_evt_count_all_fov = output_all_fov[1]
    on_lightcurve_all_fov = output_all_fov[2]
    on_pattern_all_fov = output_all_fov[3]
    on_spectrum_all_fov = output_all_fov[4]
    on_detx_all_fov = output_all_fov[5]
    on_image_all_fov = output_all_fov[6]
    
    on_duration_all_cor = output_all_cor[0]
    on_evt_count_all_cor = output_all_cor[1]
    on_lightcurve_all_cor = output_all_cor[2]
    on_pattern_all_cor = output_all_cor[3]
    on_spectrum_all_cor = output_all_cor[4]
    on_detx_all_cor = output_all_cor[5]
    on_image_all_cor = output_all_cor[6]
    
    on_duration_fwc_fov = output_fwc_fov[0]
    on_evt_count_fwc_fov = output_fwc_fov[1]
    on_lightcurve_fwc_fov = output_fwc_fov[2]
    on_pattern_fwc_fov = output_fwc_fov[3]
    on_spectrum_fwc_fov = output_fwc_fov[4]
    on_detx_fwc_fov = output_fwc_fov[5]
    on_image_fwc_fov = output_fwc_fov[6]
    
    on_duration_fwc_cor = output_fwc_cor[0]
    on_evt_count_fwc_cor = output_fwc_cor[1]
    on_lightcurve_fwc_cor = output_fwc_cor[2]
    on_pattern_fwc_cor = output_fwc_cor[3]
    on_spectrum_fwc_cor = output_fwc_cor[4]
    on_detx_fwc_cor = output_fwc_cor[5]
    on_image_fwc_cor = output_fwc_cor[6]

    area_pix_frac_all_fov = on_image_all_fov.get_pixel_fraction()
    area_pix_frac_all_cor = on_image_all_cor.get_pixel_fraction()
    area_pix_frac_fwc_fov = on_image_fwc_fov.get_pixel_fraction()
    area_pix_frac_fwc_cor = on_image_fwc_cor.get_pixel_fraction()
    
    time_pix_frac_all_fov = on_lightcurve_all_fov.get_pixel_fraction()
    time_expo_all_fov = on_duration_all_fov*time_pix_frac_all_fov
    time_pix_frac_fwc_fov = on_lightcurve_fwc_fov.get_pixel_fraction()
    time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
    
    time_pix_frac_all_cor = on_lightcurve_all_cor.get_pixel_fraction()
    time_expo_all_cor = on_duration_all_cor*time_pix_frac_all_cor
    time_pix_frac_fwc_cor = on_lightcurve_fwc_cor.get_pixel_fraction()
    time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor
    
    on_lightcurve_all_cor.scale(area_pix_frac_all_fov/area_pix_frac_all_cor)
    on_lightcurve_fwc_fov.scale(fwc_2_sci_ratio*time_expo_all_fov/time_expo_fwc_fov)

    mask_lc = make_timecut_mask(on_lightcurve_all_fov,on_lightcurve_all_cor) 
    
    time_pix_frac_mask = mask_lc.get_pixel_fraction()
    if time_pix_frac_mask>0.8:
        mask_lc = None
    if not SP_flare_mask:
        mask_lc = None

    print ('apply space and time masks')
    
    output_sci_fov = read_event_file(on_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    output_sci_cor = read_event_file(on_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    
    on_duration_sci_fov = output_sci_fov[0]
    on_evt_count_sci_fov = output_sci_fov[1]
    on_lightcurve_sci_fov = output_sci_fov[2]
    on_pattern_sci_fov = output_sci_fov[3]
    on_spectrum_sci_fov = output_sci_fov[4]
    on_detx_sci_fov = output_sci_fov[5]
    on_image_sci_fov = output_sci_fov[6]
    on_detr_sci_fov = output_sci_fov[7]
    
    on_duration_sci_cor = output_sci_cor[0]
    on_evt_count_sci_cor = output_sci_cor[1]
    on_lightcurve_sci_cor = output_sci_cor[2]
    on_pattern_sci_cor = output_sci_cor[3]
    on_spectrum_sci_cor = output_sci_cor[4]
    on_detx_sci_cor = output_sci_cor[5]
    on_image_sci_cor = output_sci_cor[6]
    
    output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    
    on_duration_fwc_fov = output_fwc_fov[0]
    on_evt_count_fwc_fov = output_fwc_fov[1]
    on_lightcurve_fwc_fov = output_fwc_fov[2]
    on_pattern_fwc_fov = output_fwc_fov[3]
    on_spectrum_fwc_fov = output_fwc_fov[4]
    on_detx_fwc_fov = output_fwc_fov[5]
    on_image_fwc_fov = output_fwc_fov[6]
    on_detr_fwc_fov = output_fwc_fov[7]
    
    on_duration_fwc_cor = output_fwc_cor[0]
    on_evt_count_fwc_cor = output_fwc_cor[1]
    on_lightcurve_fwc_cor = output_fwc_cor[2]
    on_pattern_fwc_cor = output_fwc_cor[3]
    on_spectrum_fwc_cor = output_fwc_cor[4]
    on_detx_fwc_cor = output_fwc_cor[5]
    on_image_fwc_cor = output_fwc_cor[6]
    
    area_pix_frac_sci_fov = on_image_sci_fov.get_pixel_fraction()
    area_pix_frac_sci_cor = on_image_sci_cor.get_pixel_fraction()
    area_pix_frac_fwc_fov = on_image_fwc_fov.get_pixel_fraction()
    area_pix_frac_fwc_cor = on_image_fwc_cor.get_pixel_fraction()
    
    time_pix_frac_sci_fov = on_lightcurve_sci_fov.get_pixel_fraction()
    time_expo_sci_fov = on_duration_sci_fov*time_pix_frac_sci_fov
    time_pix_frac_sci_cor = on_lightcurve_sci_cor.get_pixel_fraction()
    time_expo_sci_cor = on_duration_sci_cor*time_pix_frac_sci_cor
    
    time_pix_frac_fwc_fov = on_lightcurve_fwc_fov.get_pixel_fraction()
    time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
    time_pix_frac_fwc_cor = on_lightcurve_fwc_cor.get_pixel_fraction()
    time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor
    
    on_lightcurve_sci_cor.scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    on_lightcurve_fwc_fov.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    on_lightcurve_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
    on_spectrum_sci_cor.scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    on_spectrum_fwc_fov.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    on_spectrum_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
    on_pattern_sci_cor.scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    on_pattern_fwc_fov.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    on_pattern_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)

    if diagnostic_plots:
        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_lightcurve_all_fov.xaxis,on_lightcurve_all_fov.yaxis,yerr=on_lightcurve_all_fov.yerr,color='k',label='All FoV')
        axbig.errorbar(on_lightcurve_sci_fov.xaxis,on_lightcurve_sci_fov.yaxis,yerr=on_lightcurve_sci_fov.yerr,color='red',label='Sci FoV')
        axbig.errorbar(on_lightcurve_all_cor.xaxis,on_lightcurve_all_cor.yaxis,yerr=on_lightcurve_all_cor.yerr,color='green',label='All Cor')
        axbig.set_yscale('log')
        axbig.set_xlabel('Time')
        axbig.legend(loc='best')
        fig.savefig("../output_plots/lightcurve_on_timecut_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_spectrum_sci_fov.xaxis,on_spectrum_sci_fov.yaxis,yerr=on_spectrum_sci_fov.yerr,label='Sci FoV')
        axbig.errorbar(on_spectrum_fwc_fov.xaxis,on_spectrum_fwc_fov.yaxis,yerr=on_spectrum_fwc_fov.yerr,label='FWC FoV')
        axbig.set_xlabel('Channel')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("../output_plots/on_spectrum_fov_raw_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_spectrum_sci_cor.xaxis,on_spectrum_sci_cor.yaxis,yerr=on_spectrum_sci_cor.yerr,label='Sci Cor')
        axbig.errorbar(on_spectrum_fwc_cor.xaxis,on_spectrum_fwc_cor.yaxis,yerr=on_spectrum_fwc_cor.yerr,label='FWC Cor')
        axbig.set_xlabel('Channel')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("../output_plots/on_spectrum_cor_raw_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_pattern_sci_fov.xaxis,on_pattern_sci_fov.yaxis,yerr=on_pattern_sci_fov.yerr,label='Sci FoV')
        axbig.errorbar(on_pattern_fwc_fov.xaxis,on_pattern_fwc_fov.yaxis,yerr=on_pattern_fwc_fov.yerr,label='FWC FoV')
        axbig.set_xlabel('Channel')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("../output_plots/on_pattern_fov_raw_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_pattern_sci_cor.xaxis,on_pattern_sci_cor.yaxis,yerr=on_pattern_sci_cor.yerr,label='Sci Cor')
        axbig.errorbar(on_pattern_fwc_cor.xaxis,on_pattern_fwc_cor.yaxis,yerr=on_pattern_fwc_cor.yerr,label='FWC Cor')
        axbig.set_xlabel('Channel')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("../output_plots/on_pattern_cor_raw_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

    print ('predict SP background')
    
    qpb_pattern_template = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
    qpb_pattern_template.add(on_pattern_sci_cor)
    qpb_spectrum_template = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    qpb_spectrum_template.add(on_spectrum_sci_cor)

    spf_scaling, xray_scaling = fit_pattern(on_pattern_sci_fov,xray_pattern_template,sp_pattern_template,qpb_pattern_template)
    print ('xray_scaling = %s'%(xray_scaling))
    print ('spf_scaling = %s'%(spf_scaling))
        
    xray_pattern_template.scale(xray_scaling)
    sp_pattern_template.scale(spf_scaling)
        
    sp_spectrum_template.scale(spf_scaling)
    on_spectrum_sci_fov_bkg = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    on_spectrum_sci_fov_bkg.add(sp_spectrum_template)
    on_spectrum_sci_fov_bkg.add(qpb_spectrum_template)

    if diagnostic_plots:
        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_pattern_sci_fov.xaxis,on_pattern_sci_fov.yaxis,yerr=on_pattern_sci_fov.yerr,color='k',label='Data')
        axbig.errorbar(xray_pattern_template.xaxis,xray_pattern_template.yaxis,yerr=xray_pattern_template.yerr,color='g',label='X-ray')
        axbig.errorbar(qpb_pattern_template.xaxis,qpb_pattern_template.yaxis,yerr=qpb_pattern_template.yerr,color='b',label='QPB')
        axbig.bar(sp_pattern_template.xaxis,sp_pattern_template.yaxis,color='orange',label='SP Flare')
        axbig.set_yscale('log')
        axbig.set_xlabel('PATTERN')
        axbig.legend(loc='best')
        fig.savefig("../output_plots/pattern_on_fov_scaled_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_spectrum_sci_fov.xaxis,on_spectrum_sci_fov.yaxis,yerr=on_spectrum_sci_fov.yerr,color='k',label='Data')
        axbig.plot(on_spectrum_sci_fov_bkg.xaxis,on_spectrum_sci_fov_bkg.yaxis,color='red',label='Bkg')
        axbig.plot(sp_spectrum_template.xaxis,sp_spectrum_template.yaxis,label='SP')
        axbig.plot(qpb_spectrum_template.xaxis,qpb_spectrum_template.yaxis,label='QPB')
        axbig.set_yscale('log')
        axbig.set_xlabel('Channel')
        axbig.legend(loc='best')
        fig.savefig("../output_plots/spectrum_on_fit_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        xray_pattern_template.normalize()
        qpb_pattern_template.normalize()
        sp_pattern_template.normalize()
        
        fig.clf()
        axbig = fig.add_subplot()
        axbig.plot(xray_pattern_template.xaxis,xray_pattern_template.yaxis,label='X-ray')
        axbig.plot(qpb_pattern_template.xaxis,qpb_pattern_template.yaxis,label='QPB')
        axbig.plot(sp_pattern_template.xaxis,sp_pattern_template.yaxis,label='SP')
        axbig.set_yscale('log')
        axbig.set_xlabel('Pattern')
        axbig.legend(loc='best')
        fig.savefig("../output_plots/pattern_template_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()


    qpb_detx_template = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
    qpb_detx_template.add(on_detx_fwc_fov)
    sp_detx_template = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
    sp_detx_template.add(on_detx_fwc_fov)
    qpb_detr_template = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)
    qpb_detr_template.add(on_detr_fwc_fov)
    sp_detr_template = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)
    sp_detr_template.add(on_detr_fwc_fov)
    qpb_image_template = MyArray2D()
    qpb_image_template.add(on_image_fwc_fov)
    sp_image_template = MyArray2D()
    sp_image_template.add(on_image_fwc_fov)

    qpb_template_scaling = qpb_spectrum_template.integral()/qpb_detx_template.integral()
    sp_template_scaling = sp_spectrum_template.integral()/sp_detx_template.integral()

    qpb_detx_template.scale(qpb_template_scaling)
    sp_detx_template.scale(sp_template_scaling)
    on_detx_sci_fov_bkg = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
    on_detx_sci_fov_bkg.add(sp_detx_template)
    on_detx_sci_fov_bkg.add(qpb_detx_template)
    
    qpb_detr_template.scale(qpb_template_scaling)
    sp_detr_template.scale(sp_template_scaling)
    on_detr_sci_fov_bkg = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)
    on_detr_sci_fov_bkg.add(sp_detr_template)
    on_detr_sci_fov_bkg.add(qpb_detr_template)
    
    qpb_image_template.scale(qpb_template_scaling)
    sp_image_template.scale(sp_template_scaling)
    on_image_sci_fov_bkg = MyArray2D()
    on_image_sci_fov_bkg.add(sp_image_template)
    on_image_sci_fov_bkg.add(qpb_image_template)

    if diagnostic_plots:
        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_detx_sci_fov.xaxis,on_detx_sci_fov.yaxis,yerr=on_detx_sci_fov.yerr,color='k',label='Data')
        axbig.plot(on_detx_sci_fov_bkg.xaxis,on_detx_sci_fov_bkg.yaxis,color='red',label='Bkg')
        axbig.plot(sp_detx_template.xaxis,sp_detx_template.yaxis,label='SP')
        axbig.plot(qpb_detx_template.xaxis,qpb_detx_template.yaxis,label='QPB')
        axbig.set_yscale('log')
        axbig.set_ylim(bottom=1)
        axbig.set_xlabel('DETX')
        axbig.legend(loc='best')
        fig.savefig("../output_plots/detx_on_fit_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()
        
        axbig = fig.add_subplot()
        axbig.errorbar(on_detr_sci_fov.xaxis,on_detr_sci_fov.yaxis,yerr=on_detr_sci_fov.yerr,color='k',label='Data')
        axbig.plot(on_detr_sci_fov_bkg.xaxis,on_detr_sci_fov_bkg.yaxis,color='red',label='Bkg')
        axbig.plot(sp_detr_template.xaxis,sp_detr_template.yaxis,label='SP')
        axbig.plot(qpb_detr_template.xaxis,qpb_detr_template.yaxis,label='QPB')
        axbig.set_yscale('log')
        axbig.set_ylim(bottom=1)
        axbig.set_xlabel('DETR')
        axbig.legend(loc='best')
        fig.savefig("../output_plots/detr_on_fit_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()
        
        
        fig.clf()
        axbig = fig.add_subplot()
        label_x = 'DETY'
        label_y = 'DETX'
        axbig.set_xlabel(label_x)
        axbig.set_ylabel(label_y)
        xmin = on_image_sci_fov.xaxis.min()
        xmax = on_image_sci_fov.xaxis.max()
        ymin = on_image_sci_fov.yaxis.min()
        ymax = on_image_sci_fov.yaxis.max()
        axbig.imshow(on_image_sci_fov.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
        fig.savefig("../output_plots/skymap_on_sci_fov_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()
        
        fig.clf()
        axbig = fig.add_subplot()
        label_x = 'DETY'
        label_y = 'DETX'
        axbig.set_xlabel(label_x)
        axbig.set_ylabel(label_y)
        xmin = on_image_sci_fov_bkg.xaxis.min()
        xmax = on_image_sci_fov_bkg.xaxis.max()
        ymin = on_image_sci_fov_bkg.yaxis.min()
        ymax = on_image_sci_fov_bkg.yaxis.max()
        axbig.imshow(on_image_sci_fov_bkg.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
        fig.savefig("../output_plots/skymap_on_sci_fov_bkg_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()
        
        on_image_sci_fov_excess = MyArray2D()
        on_image_sci_fov_excess.add(on_image_sci_fov)
        on_image_sci_fov_excess.add(on_image_sci_fov_bkg,factor=-1.)
        
        fig.clf()
        axbig = fig.add_subplot()
        label_x = 'DETY'
        label_y = 'DETX'
        axbig.set_xlabel(label_x)
        axbig.set_ylabel(label_y)
        xmin = on_image_sci_fov_excess.xaxis.min()
        xmax = on_image_sci_fov_excess.xaxis.max()
        ymin = on_image_sci_fov_excess.yaxis.min()
        ymax = on_image_sci_fov_excess.yaxis.max()
        axbig.imshow(on_image_sci_fov_excess.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
        fig.savefig("../output_plots/skymap_on_sci_fov_excess_ch%s_ccd%s_%s.png"%(energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

    output_spectrum = [on_spectrum_sci_fov,on_spectrum_sci_fov_bkg,sp_spectrum_template,qpb_spectrum_template]
    output_detx = [on_detx_sci_fov,on_detx_sci_fov_bkg,sp_detx_template,qpb_detx_template]
    output_detr = [on_detr_sci_fov,on_detr_sci_fov_bkg,sp_detr_template,qpb_detr_template]

    return output_spectrum, output_detx, output_detr


on_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
bkg_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
qpb_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
sp_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)

on_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
bkg_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
qpb_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
sp_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)

on_detr_sum = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)
bkg_detr_sum = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)
qpb_detr_sum = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)
sp_detr_sum = MyArray1D(bin_start=detr_low,bin_end=detr_high,pixel_scale=detr_scale)

for eng in range(0,len(energy_array)-1):
    for ccd in range(0,len(ana_ccd_bins)):
        ana_output_spectrum, ana_output_detx, ana_output_detr = analyze_a_ccd_chip(energy_range=[energy_array[eng],energy_array[eng+1]],ccd_id=ana_ccd_bins[ccd])
        on_spectrum_sum.add(ana_output_spectrum[0])
        bkg_spectrum_sum.add(ana_output_spectrum[1])
        sp_spectrum_sum.add(ana_output_spectrum[2])
        qpb_spectrum_sum.add(ana_output_spectrum[3])
        on_detx_sum.add(ana_output_detx[0])
        bkg_detx_sum.add(ana_output_detx[1])
        sp_detx_sum.add(ana_output_detx[2])
        qpb_detx_sum.add(ana_output_detx[3])
        on_detr_sum.add(ana_output_detr[0])
        bkg_detr_sum.add(ana_output_detr[1])
        sp_detr_sum.add(ana_output_detr[2])
        qpb_detr_sum.add(ana_output_detr[3])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_spectrum_sum.xaxis,on_spectrum_sum.yaxis,yerr=on_spectrum_sum.yerr,color='k',label='Data')
axbig.plot(bkg_spectrum_sum.xaxis,bkg_spectrum_sum.yaxis,color='red',label='Bkg')
axbig.plot(sp_spectrum_sum.xaxis,sp_spectrum_sum.yaxis,label='SP')
axbig.plot(qpb_spectrum_sum.xaxis,qpb_spectrum_sum.yaxis,label='QPB')
axbig.set_yscale('log')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
fig.savefig("../output_plots/spectrum_on_fit_sum_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_detx_sum.xaxis,on_detx_sum.yaxis,yerr=on_detx_sum.yerr,color='k',label='Data')
axbig.plot(bkg_detx_sum.xaxis,bkg_detx_sum.yaxis,color='red',label='Bkg')
axbig.plot(sp_detx_sum.xaxis,sp_detx_sum.yaxis,label='SP')
axbig.plot(qpb_detx_sum.xaxis,qpb_detx_sum.yaxis,label='QPB')
axbig.set_yscale('log')
axbig.set_ylim(bottom=1)
axbig.set_xlabel('DETX')
axbig.legend(loc='best')
fig.savefig("../output_plots/detx_on_fit_sum_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_detr_sum.xaxis,on_detr_sum.yaxis,yerr=on_detr_sum.yerr,color='k',label='Data')
axbig.plot(bkg_detr_sum.xaxis,bkg_detr_sum.yaxis,color='red',label='Bkg')
axbig.plot(sp_detr_sum.xaxis,sp_detr_sum.yaxis,label='SP')
axbig.plot(qpb_detr_sum.xaxis,qpb_detr_sum.yaxis,label='QPB')
axbig.set_yscale('log')
axbig.set_ylim(bottom=1)
axbig.set_xlabel('R')
axbig.legend(loc='best')
fig.savefig("../output_plots/detr_on_fit_sum_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()



























end_time = time.time()
elapsed_time = end_time - start_time
print('elapsed_time = %s'%(elapsed_time))

exit()


