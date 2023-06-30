
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
from operator import itemgetter, attrgetter
from scipy.optimize import curve_fit
import time

fig, ax = plt.subplots()
figsize_x = 8.
figsize_y = 6.
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

energy_lower = 250
energy_upper = 12000

energy_bins = [250,12000]

ch_low = 200
ch_high = 12000
ch_scale = 200
t_low = 0
t_high = 1
t_scale = 0.05
detx_low = -20000
detx_high = 20000
detx_scale = 500

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

def read_event_file(filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[200,12000]):

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

    evt_count = []
    lightcurve_array = []
    pattern_array = []
    spectrum_array = []
    detx_array = []
    image_array = []

    for ebin in range(0,len(energy_bins)):
        evt_count += [0]
        lightcurve_array += [MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)]
        pattern_array += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
        spectrum_array += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
        detx_array += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
        image_array += [MyArray2D()]

    time_start = events[0]['TIME']
    time_end = events[len(events)-1]['TIME']
    obs_duration = time_end-time_start
    for evt in range(0,len(events)):
        evt_time = events[evt]['TIME']
        evt_pattern = events[evt]['PATTERN']
        evt_pi = events[evt]['PI']
        evt_detx = events[evt]['DETX']
        evt_dety = events[evt]['DETY']

        #if abs(evt_detx)>20000: continue
        #if abs(evt_dety)>20000: continue
        #if abs(evt_dety)>7000: continue
        #if evt_detx<7000: continue
        #if evt_detx>-7000: continue

        if evt_pi<energy_range[0]: continue
        if evt_pi>energy_range[1]: continue

        # if 'sp-free' in evt_filter:
        #if evt_pattern==0:continue
        #if evt_pattern==1:continue
        #if evt_pattern==3:continue

        #if evt_pattern<=0:continue
        #if evt_pattern>12:continue

        #if 'Al' in evt_filter:
        #    if evt_pi<1200: continue
        #    if evt_pi>1900: continue
        #if evt_pi>1200 and evt_pi<1900: continue
        #if evt_pi<1900: continue

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

        evt_ebin = 0
        for ebin in range(0,len(energy_bins)-1):
            if evt_pi>=energy_bins[ebin] and evt_pi<energy_bins[ebin+1]:
                evt_ebin = ebin+1

        evt_count[0] += 1.
        lightcurve_array[0].fill((evt_time-time_start)/(time_end-time_start))
        pattern_array[0].fill(evt_pattern)
        spectrum_array[0].fill(evt_pi)
        detx_array[0].fill(evt_detx)
        image_array[0].fill(evt_detx,evt_dety)

        evt_count[evt_ebin] += 1.
        lightcurve_array[evt_ebin].fill((evt_time-time_start)/(time_end-time_start))
        pattern_array[evt_ebin].fill(evt_pattern)
        spectrum_array[evt_ebin].fill(evt_pi)
        detx_array[evt_ebin].fill(evt_detx)
        image_array[evt_ebin].fill(evt_detx,evt_dety)

    return [obs_duration, evt_count, lightcurve_array, pattern_array, spectrum_array, detx_array, image_array]

def fit_pattern(data_pattern,xray_pattern,spf_pattern,qpb_pattern):

    ## Define the function to minimize (residuals)
    #def residuals(A):
    #    return (data_pattern.yaxis - A[0]*xray_pattern.yaxis - A[1]*spf_pattern.yaxis - qpb_pattern.yaxis) / data_pattern.error  # Include the weighted residuals

    #data_pattern.get_error()
    ## Perform the least squares fitting
    #result = least_squares(residuals, x0=[0.1,0.1], bounds=([0.,0.], [1.0,1.0]))  # x0 is the initial guess for A

    #xray_pattern.scale(result.x[0])
    #spf_pattern.scale(result.x[1])


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

#detector = 'mos1S001'
detector = 'mos2S002'


SP_flare_mask = True
#point_source_mask = True
#SP_flare_mask = False
point_source_mask = False

on_sample = sys.argv[1]
on_obsID = sys.argv[2]
#on_sample = 'extragalactic'
#on_obsID = 'ID0505460501' # Percentage of flaring time: 60%
#on_sample = 'extragalactic'
#on_obsID = 'ID0820560101' # Percentage of flaring time: 48.8%
#on_sample = 'extragalactic'
#on_obsID = 'ID0803160301' # Percentage of flaring time: 28.4%
#on_sample = 'extragalactic'
#on_obsID = 'ID0827241101' # Percentage of flaring time: 15.2%
#on_sample = 'extragalactic'
#on_obsID = 'ID0827200401' # Percentage of flaring time: 11.9%
#on_sample = 'extragalactic'
#on_obsID = 'ID0827200501' # Percentage of flaring time: 6%

#on_sample = 'extragalactic'
#on_obsID = 'ID0827251101' # Percentage of flaring time: 0.0%

#on_sample = 'RX_J1241.5_p3250'
#on_obsID = 'ID0056020901' # Percentage of flaring time: %
#on_sample = 'RX_J0256.5_p0006'
#on_obsID = 'ID0056020301' # Percentage of flaring time: 60.0%
#on_sample = 'RX_J1713.7_m3941'
#on_obsID = 'ID0093670301'
#on_sample = 'Cas_A'
#on_obsID = 'ID0412180101'
#on_sample = '3HWC_J1928_p178'
#on_obsID = 'ID0902120101'

off_sample = 'extragalactic'
off_obsID = 'ID0827241101' # Percentage of flaring time: 15.2%
#off_sample = on_sample
#off_obsID = on_obsID

plot_tag = detector+'_'+on_obsID

off_sci_fov_evt_filename = '../%s/%s/analysis/%s-fov-evt.fits'%(off_sample,off_obsID,detector)
off_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-fov-evt.fits'%(off_sample,off_obsID,detector)
off_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(off_sample,off_obsID,detector)
off_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(off_sample,off_obsID,detector)

on_sci_fov_evt_filename = '../%s/%s/analysis/%s-fov-evt.fits'%(on_sample,on_obsID,detector)
on_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-fov-evt.fits'%(on_sample,on_obsID,detector)
on_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(on_sample,on_obsID,detector)
on_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(on_sample,on_obsID,detector)

print ('prepare off sample time cuts')

output_all_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_all_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])

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

area_pix_frac_all_fov = off_image_all_fov[0].get_pixel_fraction()
area_pix_frac_all_cor = off_image_all_cor[0].get_pixel_fraction()
area_pix_frac_fwc_fov = off_image_fwc_fov[0].get_pixel_fraction()
area_pix_frac_fwc_cor = off_image_fwc_cor[0].get_pixel_fraction()

time_pix_frac_all_fov = off_lightcurve_all_fov[0].get_pixel_fraction()
time_expo_all_fov = off_duration_all_fov*time_pix_frac_all_fov
time_pix_frac_fwc_fov = off_lightcurve_fwc_fov[0].get_pixel_fraction()
time_expo_fwc_fov = off_duration_fwc_fov*time_pix_frac_fwc_fov

time_pix_frac_all_cor = off_lightcurve_all_cor[0].get_pixel_fraction()
time_expo_all_cor = off_duration_all_cor*time_pix_frac_all_cor
time_pix_frac_fwc_cor = off_lightcurve_fwc_cor[0].get_pixel_fraction()
time_expo_fwc_cor = off_duration_fwc_cor*time_pix_frac_fwc_cor

for ebin in range(0,len(energy_bins)):
    off_lightcurve_all_cor[ebin].scale(area_pix_frac_all_fov/area_pix_frac_all_cor)
    off_lightcurve_fwc_fov[ebin].scale(fwc_2_sci_ratio*time_expo_all_fov/time_expo_fwc_fov)

#mask_lc = make_timecut_mask(off_lightcurve_all_cor[0],off_lightcurve_fwc_fov[0]) 
mask_lc = make_timecut_mask(off_lightcurve_all_fov[0],off_lightcurve_all_cor[0]) 


print ('apply space and time masks')

output_sci_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=[energy_lower,energy_upper])
output_sci_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=[energy_lower,energy_upper])

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

output_spf_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=[energy_lower,energy_upper])
output_spf_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=[energy_lower,energy_upper])

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

output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])

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

area_pix_frac_sci_fov = off_image_sci_fov[0].get_pixel_fraction()
area_pix_frac_sci_cor = off_image_sci_cor[0].get_pixel_fraction()
area_pix_frac_spf_fov = off_image_spf_fov[0].get_pixel_fraction()
area_pix_frac_spf_cor = off_image_spf_cor[0].get_pixel_fraction()
area_pix_frac_fwc_fov = off_image_fwc_fov[0].get_pixel_fraction()
area_pix_frac_fwc_cor = off_image_fwc_cor[0].get_pixel_fraction()


time_pix_frac_sci_fov = off_lightcurve_sci_fov[0].get_pixel_fraction()
time_expo_sci_fov = off_duration_sci_fov*time_pix_frac_sci_fov
time_pix_frac_sci_cor = off_lightcurve_sci_cor[0].get_pixel_fraction()
time_expo_sci_cor = off_duration_sci_cor*time_pix_frac_sci_cor

time_pix_frac_spf_fov = off_lightcurve_spf_fov[0].get_pixel_fraction()
time_expo_spf_fov = off_duration_spf_fov*time_pix_frac_spf_fov
time_pix_frac_spf_cor = off_lightcurve_spf_cor[0].get_pixel_fraction()
time_expo_spf_cor = off_duration_spf_cor*time_pix_frac_spf_cor

time_pix_frac_fwc_fov = off_lightcurve_fwc_fov[0].get_pixel_fraction()
time_expo_fwc_fov = off_duration_fwc_fov*time_pix_frac_fwc_fov
time_pix_frac_fwc_cor = off_lightcurve_fwc_cor[0].get_pixel_fraction()
time_expo_fwc_cor = off_duration_fwc_cor*time_pix_frac_fwc_cor

print ('area_pix_frac_spf_fov')
print (area_pix_frac_spf_fov)
print ('time_expo_spf_fov')
print (time_expo_spf_fov)

for ebin in range(0,len(energy_bins)):

    off_lightcurve_sci_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    off_lightcurve_spf_fov[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
    off_lightcurve_spf_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
    off_lightcurve_fwc_fov[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    off_lightcurve_fwc_cor[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)

    off_spectrum_sci_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    off_spectrum_spf_fov[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
    off_spectrum_spf_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
    off_spectrum_fwc_fov[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    off_spectrum_fwc_cor[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
    off_pattern_sci_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    off_pattern_spf_fov[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
    off_pattern_spf_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
    off_pattern_fwc_fov[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    off_pattern_fwc_cor[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)


fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_lightcurve_all_fov[0].xaxis,off_lightcurve_all_fov[0].yaxis,yerr=off_lightcurve_all_fov[0].yerr,color='k',label='All FoV')
axbig.errorbar(off_lightcurve_sci_fov[0].xaxis,off_lightcurve_sci_fov[0].yaxis,yerr=off_lightcurve_sci_fov[0].yerr,color='red',label='Sci FoV')
axbig.errorbar(off_lightcurve_all_cor[0].xaxis,off_lightcurve_all_cor[0].yaxis,yerr=off_lightcurve_all_cor[0].yerr,color='green',label='All Cor')
axbig.errorbar(off_lightcurve_fwc_fov[0].xaxis,off_lightcurve_fwc_fov[0].yaxis,yerr=off_lightcurve_fwc_fov[0].yerr,color='blue',label='FWC FoV')
axbig.set_yscale('log')
axbig.set_xlabel('Time')
axbig.legend(loc='best')
fig.savefig("../output_plots/lightcurve_off_timecut_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_spectrum_sci_fov[0].xaxis,off_spectrum_sci_fov[0].yaxis,yerr=off_spectrum_sci_fov[0].yerr,label='Sci FoV')
axbig.errorbar(off_spectrum_spf_fov[0].xaxis,off_spectrum_spf_fov[0].yaxis,yerr=off_spectrum_spf_fov[0].yerr,label='SPF FoV')
axbig.errorbar(off_spectrum_fwc_fov[0].xaxis,off_spectrum_fwc_fov[0].yaxis,yerr=off_spectrum_fwc_fov[0].yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_spectrum_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_spectrum_sci_cor[0].xaxis,off_spectrum_sci_cor[0].yaxis,yerr=off_spectrum_sci_cor[0].yerr,label='Sci Cor')
axbig.errorbar(off_spectrum_spf_cor[0].xaxis,off_spectrum_spf_cor[0].yaxis,yerr=off_spectrum_spf_cor[0].yerr,label='SPF Cor')
axbig.errorbar(off_spectrum_fwc_cor[0].xaxis,off_spectrum_fwc_cor[0].yaxis,yerr=off_spectrum_fwc_cor[0].yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_spectrum_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_pattern_sci_fov[0].xaxis,off_pattern_sci_fov[0].yaxis,yerr=off_pattern_sci_fov[0].yerr,label='Sci FoV')
axbig.errorbar(off_pattern_spf_fov[0].xaxis,off_pattern_spf_fov[0].yaxis,yerr=off_pattern_spf_fov[0].yerr,label='SPF FoV')
axbig.errorbar(off_pattern_fwc_fov[0].xaxis,off_pattern_fwc_fov[0].yaxis,yerr=off_pattern_fwc_fov[0].yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_pattern_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_pattern_sci_cor[0].xaxis,off_pattern_sci_cor[0].yaxis,yerr=off_pattern_sci_cor[0].yerr,label='Sci Cor')
axbig.errorbar(off_pattern_spf_cor[0].xaxis,off_pattern_spf_cor[0].yaxis,yerr=off_pattern_spf_cor[0].yerr,label='SPF Cor')
axbig.errorbar(off_pattern_fwc_cor[0].xaxis,off_pattern_fwc_cor[0].yaxis,yerr=off_pattern_fwc_cor[0].yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_pattern_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()


sp_pattern_template = []
sp_spectrum_template = []
sp_detx_template = []
for ebin in range(0,len(energy_bins)):

    sp_pattern_template += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
    sp_pattern_template[ebin].add(off_pattern_spf_fov[ebin])
    sp_pattern_template[ebin].add(off_pattern_sci_fov[ebin],factor=-1.)
    
    sp_spectrum_template += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
    sp_spectrum_template[ebin].add(off_spectrum_spf_fov[ebin])
    sp_spectrum_template[ebin].add(off_spectrum_sci_fov[ebin],factor=-1.)
    
    sp_detx_template += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
    sp_detx_template[ebin].add(off_detx_spf_fov[ebin])
    sp_detx_template[ebin].add(off_detx_sci_fov[ebin],factor=-1.)




print ('prepare on sample time cuts')

output_all_fov = read_event_file(on_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_all_cor = read_event_file(on_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])

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

area_pix_frac_all_fov = on_image_all_fov[0].get_pixel_fraction()
area_pix_frac_all_cor = on_image_all_cor[0].get_pixel_fraction()
area_pix_frac_fwc_fov = on_image_fwc_fov[0].get_pixel_fraction()
area_pix_frac_fwc_cor = on_image_fwc_cor[0].get_pixel_fraction()

time_pix_frac_all_fov = on_lightcurve_all_fov[0].get_pixel_fraction()
time_expo_all_fov = on_duration_all_fov*time_pix_frac_all_fov
time_pix_frac_fwc_fov = on_lightcurve_fwc_fov[0].get_pixel_fraction()
time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov

time_pix_frac_all_cor = on_lightcurve_all_cor[0].get_pixel_fraction()
time_expo_all_cor = on_duration_all_cor*time_pix_frac_all_cor
time_pix_frac_fwc_cor = on_lightcurve_fwc_cor[0].get_pixel_fraction()
time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor


for ebin in range(0,len(energy_bins)):
    on_lightcurve_all_cor[ebin].scale(area_pix_frac_all_fov/area_pix_frac_all_cor)
    on_lightcurve_fwc_fov[ebin].scale(fwc_2_sci_ratio*time_expo_all_fov/time_expo_fwc_fov)

#mask_lc = make_timecut_mask(on_lightcurve_all_cor[0],on_lightcurve_fwc_fov[0]) 
mask_lc = make_timecut_mask(on_lightcurve_all_fov[0],on_lightcurve_all_cor[0]) 

time_pix_frac_mask = mask_lc.get_pixel_fraction()
if time_pix_frac_mask>0.8:
    mask_lc = None
if not SP_flare_mask:
    mask_lc = None


print ('apply space and time masks')

output_sci_fov = read_event_file(on_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=[energy_lower,energy_upper])
output_sci_cor = read_event_file(on_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=[energy_lower,energy_upper])

on_duration_sci_fov = output_sci_fov[0]
on_evt_count_sci_fov = output_sci_fov[1]
on_lightcurve_sci_fov = output_sci_fov[2]
on_pattern_sci_fov = output_sci_fov[3]
on_spectrum_sci_fov = output_sci_fov[4]
on_detx_sci_fov = output_sci_fov[5]
on_image_sci_fov = output_sci_fov[6]

on_duration_sci_cor = output_sci_cor[0]
on_evt_count_sci_cor = output_sci_cor[1]
on_lightcurve_sci_cor = output_sci_cor[2]
on_pattern_sci_cor = output_sci_cor[3]
on_spectrum_sci_cor = output_sci_cor[4]
on_detx_sci_cor = output_sci_cor[5]
on_image_sci_cor = output_sci_cor[6]

output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])
output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=[energy_lower,energy_upper])

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

area_pix_frac_sci_fov = on_image_sci_fov[0].get_pixel_fraction()
area_pix_frac_sci_cor = on_image_sci_cor[0].get_pixel_fraction()
print ('area_pix_frac_sci_fov = %s'%(area_pix_frac_sci_fov))
print ('area_pix_frac_sci_cor = %s'%(area_pix_frac_sci_cor))
area_pix_frac_fwc_fov = on_image_fwc_fov[0].get_pixel_fraction()
area_pix_frac_fwc_cor = on_image_fwc_cor[0].get_pixel_fraction()
print ('area_pix_frac_fwc_fov = %s'%(area_pix_frac_fwc_fov))
print ('area_pix_frac_fwc_cor = %s'%(area_pix_frac_fwc_cor))

time_pix_frac_sci_fov = on_lightcurve_sci_fov[0].get_pixel_fraction()
time_expo_sci_fov = on_duration_sci_fov*time_pix_frac_sci_fov
time_pix_frac_sci_cor = on_lightcurve_sci_cor[0].get_pixel_fraction()
time_expo_sci_cor = on_duration_sci_cor*time_pix_frac_sci_cor

time_pix_frac_fwc_fov = on_lightcurve_fwc_fov[0].get_pixel_fraction()
time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
time_pix_frac_fwc_cor = on_lightcurve_fwc_cor[0].get_pixel_fraction()
time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor


for ebin in range(0,len(energy_bins)):

    on_lightcurve_sci_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    on_lightcurve_fwc_fov[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    on_lightcurve_fwc_cor[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)

    on_spectrum_sci_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    on_spectrum_fwc_fov[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    on_spectrum_fwc_cor[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
    on_pattern_sci_cor[ebin].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
    on_pattern_fwc_fov[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
    on_pattern_fwc_cor[ebin].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)


fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_lightcurve_all_fov[0].xaxis,on_lightcurve_all_fov[0].yaxis,yerr=on_lightcurve_all_fov[0].yerr,color='k',label='All FoV')
axbig.errorbar(on_lightcurve_sci_fov[0].xaxis,on_lightcurve_sci_fov[0].yaxis,yerr=on_lightcurve_sci_fov[0].yerr,color='red',label='Sci FoV')
axbig.errorbar(on_lightcurve_all_cor[0].xaxis,on_lightcurve_all_cor[0].yaxis,yerr=on_lightcurve_all_cor[0].yerr,color='green',label='All Cor')
axbig.errorbar(on_lightcurve_fwc_fov[0].xaxis,on_lightcurve_fwc_fov[0].yaxis,yerr=on_lightcurve_fwc_fov[0].yerr,color='blue',label='FWC FoV')
axbig.set_yscale('log')
axbig.set_xlabel('Time')
axbig.legend(loc='best')
fig.savefig("../output_plots/lightcurve_on_timecut_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_spectrum_sci_fov[0].xaxis,on_spectrum_sci_fov[0].yaxis,yerr=on_spectrum_sci_fov[0].yerr,label='Sci FoV')
axbig.errorbar(on_spectrum_fwc_fov[0].xaxis,on_spectrum_fwc_fov[0].yaxis,yerr=on_spectrum_fwc_fov[0].yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_spectrum_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_spectrum_sci_cor[0].xaxis,on_spectrum_sci_cor[0].yaxis,yerr=on_spectrum_sci_cor[0].yerr,label='Sci Cor')
axbig.errorbar(on_spectrum_fwc_cor[0].xaxis,on_spectrum_fwc_cor[0].yaxis,yerr=on_spectrum_fwc_cor[0].yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_spectrum_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_pattern_sci_fov[0].xaxis,on_pattern_sci_fov[0].yaxis,yerr=on_pattern_sci_fov[0].yerr,label='Sci FoV')
axbig.errorbar(on_pattern_fwc_fov[0].xaxis,on_pattern_fwc_fov[0].yaxis,yerr=on_pattern_fwc_fov[0].yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_pattern_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_pattern_sci_cor[0].xaxis,on_pattern_sci_cor[0].yaxis,yerr=on_pattern_sci_cor[0].yerr,label='Sci Cor')
axbig.errorbar(on_pattern_fwc_cor[0].xaxis,on_pattern_fwc_cor[0].yaxis,yerr=on_pattern_fwc_cor[0].yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_pattern_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()



print ('predict SP background')

xray_pattern_template = []
qpb_pattern_template = []
qpb_spectrum_template = []
qpb_detx_template = []
for ebin in range(0,len(energy_bins)):
    xray_pattern_template += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
    xray_pattern_template[ebin].add(on_pattern_sci_cor[ebin])
    qpb_pattern_template += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
    qpb_pattern_template[ebin].add(on_pattern_sci_cor[ebin])
    qpb_spectrum_template += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
    qpb_spectrum_template[ebin].add(on_spectrum_sci_cor[ebin])
    qpb_detx_template += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
    qpb_detx_template[ebin].add(on_detx_fwc_fov[ebin])


on_spectrum_sci_fov_bkg = [] 
on_detx_sci_fov_bkg = []
for ebin in range(0,len(energy_bins)):

    spf_scaling, xray_scaling = fit_pattern(on_pattern_sci_fov[ebin],xray_pattern_template[ebin],sp_pattern_template[ebin],qpb_pattern_template[ebin])
    print ('xray_scaling = %s'%(xray_scaling))
    print ('spf_scaling = %s'%(spf_scaling))

    xray_pattern_template[ebin].scale(xray_scaling)
    sp_pattern_template[ebin].scale(spf_scaling)

    sp_spectrum_template[ebin].scale(spf_scaling)
    on_spectrum_sci_fov_bkg += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
    on_spectrum_sci_fov_bkg[ebin].add(sp_spectrum_template[ebin])
    on_spectrum_sci_fov_bkg[ebin].add(qpb_spectrum_template[ebin])

    qpb_detx_template[ebin].scale(qpb_spectrum_template[ebin].integral()/qpb_detx_template[ebin].integral())
    sp_detx_template[ebin].scale(sp_spectrum_template[ebin].integral()/sp_detx_template[ebin].integral())
    on_detx_sci_fov_bkg += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
    on_detx_sci_fov_bkg[ebin].add(sp_detx_template[ebin])
    on_detx_sci_fov_bkg[ebin].add(qpb_detx_template[ebin])

sp_spectrum_template[0].reset()
qpb_spectrum_template[0].reset()
sp_pattern_template[0].reset()
qpb_pattern_template[0].reset()
sp_detx_template[0].reset()
qpb_detx_template[0].reset()
on_spectrum_sci_fov_bkg[0].reset()
on_detx_sci_fov_bkg[0].reset()
for ebin in range(1,len(energy_bins)):
    sp_spectrum_template[0].add(sp_spectrum_template[ebin])
    qpb_spectrum_template[0].add(qpb_spectrum_template[ebin])
    sp_pattern_template[0].add(sp_pattern_template[ebin])
    qpb_pattern_template[0].add(qpb_pattern_template[ebin])
    sp_detx_template[0].add(sp_detx_template[ebin])
    qpb_detx_template[0].add(qpb_detx_template[ebin])
    on_spectrum_sci_fov_bkg[0].add(on_spectrum_sci_fov_bkg[ebin])
    on_detx_sci_fov_bkg[0].add(on_detx_sci_fov_bkg[ebin])

for ebin in range(0,len(energy_bins)):

    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(on_pattern_sci_fov[ebin].xaxis,on_pattern_sci_fov[ebin].yaxis,yerr=on_pattern_sci_fov[ebin].yerr,color='k',label='Data')
    axbig.errorbar(xray_pattern_template[ebin].xaxis,xray_pattern_template[ebin].yaxis,yerr=xray_pattern_template[ebin].yerr,label='X-ray')
    axbig.errorbar(qpb_pattern_template[ebin].xaxis,qpb_pattern_template[ebin].yaxis,yerr=qpb_pattern_template[ebin].yerr,label='QPB')
    axbig.bar(sp_pattern_template[ebin].xaxis,sp_pattern_template[ebin].yaxis,label='SP Flare')
    axbig.set_yscale('log')
    axbig.set_xlabel('PATTERN')
    axbig.legend(loc='best')
    fig.savefig("../output_plots/pattern_on_fov_scaled_%s_E%s.png"%(plot_tag,ebin),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(on_detx_sci_fov[ebin].xaxis,on_detx_sci_fov[ebin].yaxis,yerr=on_detx_sci_fov[ebin].yerr,color='k',label='Data')
    axbig.plot(on_detx_sci_fov_bkg[ebin].xaxis,on_detx_sci_fov_bkg[ebin].yaxis,color='red',label='Bkg')
    axbig.plot(sp_detx_template[ebin].xaxis,sp_detx_template[ebin].yaxis,label='SP')
    axbig.plot(qpb_detx_template[ebin].xaxis,qpb_detx_template[ebin].yaxis,label='QPB')
    axbig.set_yscale('log')
    axbig.set_xlabel('Channel')
    axbig.legend(loc='best')
    fig.savefig("../output_plots/detx_on_fit_%s_E%s.png"%(plot_tag,ebin),bbox_inches='tight')
    axbig.remove()



fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_spectrum_sci_fov[0].xaxis,on_spectrum_sci_fov[0].yaxis,yerr=on_spectrum_sci_fov[0].yerr,color='k',label='Data')
axbig.plot(on_spectrum_sci_fov_bkg[0].xaxis,on_spectrum_sci_fov_bkg[0].yaxis,color='red',label='Bkg')
axbig.plot(sp_spectrum_template[0].xaxis,sp_spectrum_template[0].yaxis,label='SP')
axbig.plot(qpb_spectrum_template[0].xaxis,qpb_spectrum_template[0].yaxis,label='QPB')
axbig.set_yscale('log')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
fig.savefig("../output_plots/spectrum_on_fit_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()



















end_time = time.time()
elapsed_time = end_time - start_time
print('elapsed_time = %s'%(elapsed_time))

exit()


