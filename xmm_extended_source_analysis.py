
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
from scipy.optimize import least_squares
import time

fig, ax = plt.subplots()
figsize_x = 8.
figsize_y = 6.
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

energy_array = [2000,3000,5000,8000,12000]

all_ccd_bins = [1,2,3,4,5,6,7]
ana_ccd_bins = [1,2,3,4,5,6,7]
xray_ccd_bins = [1,2,3,4,5,6,7]

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

def read_event_file(filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=[200,12000],ccd=all_ccd_bins):

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

    for ebin in range(0,len(energy_range)-1):
        evt_count_dE = []
        lightcurve_array_dE = []
        pattern_array_dE = []
        spectrum_array_dE = []
        detx_array_dE = []
        image_array_dE = []
        for ccdid in range(0,len(ccd)):
            evt_count_dE += [0]
            lightcurve_array_dE += [MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)]
            pattern_array_dE += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
            spectrum_array_dE += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
            detx_array_dE += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
            image_array_dE += [MyArray2D()]
        evt_count += [evt_count_dE]
        lightcurve_array += [lightcurve_array_dE]
        pattern_array += [pattern_array_dE]
        spectrum_array += [spectrum_array_dE]
        detx_array += [detx_array_dE]
        image_array += [image_array_dE]

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
        if evt_pi>energy_range[len(energy_range)-1]: continue

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

        evt_ebin = 0
        for ebin in range(0,len(energy_range)-1):
            if evt_pi>=energy_range[ebin] and evt_pi<energy_range[ebin+1]:
                evt_ebin = ebin

        ccd_id = 0
        if evt_detx<7000 and evt_detx>=-7000 and evt_dety<7000 and evt_dety>=-7000:
            ccd_id = 1
        if evt_detx>=7000 and evt_dety<7000 and evt_dety>=-7000:
            ccd_id = 3
        if evt_detx<-7000 and evt_dety<7000 and evt_dety>=-7000:
            ccd_id = 6
        if evt_detx>=0 and evt_dety>=7000:
            ccd_id = 4
        if evt_detx<0 and evt_dety>=7000:
            ccd_id = 5
        if evt_detx>=0 and evt_dety<-7000:
            ccd_id = 2
        if evt_detx<0 and evt_dety<-7000:
            ccd_id = 7

        evt_ccdid = -1
        for ccd_bin in range(0,len(ccd)):
            if ccd_id==ccd[ccd_bin]:
                evt_ccdid = ccd_bin
        if not 'cor-evt' in filename:
            if evt_ccdid==-1: continue

        if not 'cor-evt' in filename:
            evt_count[evt_ebin][evt_ccdid] += 1.
            lightcurve_array[evt_ebin][evt_ccdid].fill((evt_time-time_start)/(time_end-time_start))
            pattern_array[evt_ebin][evt_ccdid].fill(evt_pattern)
            spectrum_array[evt_ebin][evt_ccdid].fill(evt_pi)
            detx_array[evt_ebin][evt_ccdid].fill(evt_detx)
            image_array[evt_ebin][evt_ccdid].fill(evt_detx,evt_dety)
        else:
            for ccd_bin in range(0,len(ccd)):
                evt_count[evt_ebin][ccd_bin] += 1.
                lightcurve_array[evt_ebin][ccd_bin].fill((evt_time-time_start)/(time_end-time_start))
                pattern_array[evt_ebin][ccd_bin].fill(evt_pattern)
                spectrum_array[evt_ebin][ccd_bin].fill(evt_pi)
                detx_array[evt_ebin][ccd_bin].fill(evt_detx)
                image_array[evt_ebin][ccd_bin].fill(evt_detx,evt_dety)

    return [obs_duration, evt_count, lightcurve_array, pattern_array, spectrum_array, detx_array, image_array]

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



print ('prepare xray sample time cuts')

output_all_fov = read_event_file(xray_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)
output_all_cor = read_event_file(xray_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)

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

for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
    
        area_pix_frac_all_fov = xray_image_all_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_all_cor = xray_image_all_cor[ebin][ccdid].get_pixel_fraction()
        
        time_pix_frac_all_fov = xray_lightcurve_all_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_all_fov = xray_duration_all_fov*time_pix_frac_all_fov
        
        time_pix_frac_all_cor = xray_lightcurve_all_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_all_cor = xray_duration_all_cor*time_pix_frac_all_cor
    
        xray_lightcurve_all_cor[ebin][ccdid].scale(area_pix_frac_all_fov/area_pix_frac_all_cor)

xray_lightcurve_all_fov_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
xray_lightcurve_all_cor_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        xray_lightcurve_all_fov_sum.add(xray_lightcurve_all_fov[ebin][ccdid])
        xray_lightcurve_all_cor_sum.add(xray_lightcurve_all_cor[ebin][ccdid])

mask_lc = make_timecut_mask(xray_lightcurve_all_fov_sum,xray_lightcurve_all_cor_sum) 

time_pix_frac_mask = mask_lc.get_pixel_fraction()
if time_pix_frac_mask>0.8:
    mask_lc = None

print ('apply x-ray sample space and time masks')

output_sci_source = read_event_file(xray_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto source',energy_range=energy_array,ccd=xray_ccd_bins)
output_sci_ring = read_event_file(xray_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto ring',energy_range=energy_array,ccd=xray_ccd_bins)

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

for ebin in range(0,len(energy_array)-1):
    area_pix_frac_sci_source = 0.
    area_pix_frac_sci_ring = 0.
    for ccdid in range(0,len(xray_ccd_bins)):
        area_pix_frac_sci_source += xray_image_sci_source[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_sci_ring += xray_image_sci_ring[ebin][ccdid].get_pixel_fraction()
    for ccdid in range(0,len(xray_ccd_bins)):
        xray_pattern_sci_ring[ebin][ccdid].scale(area_pix_frac_sci_source/area_pix_frac_sci_ring)

xray_pattern_template = []
for ebin in range(0,len(energy_array)-1):
    xray_pattern_template_dE = []
    for ccdid in range(0,len(ana_ccd_bins)):
        xray_pattern_template_dE += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
    xray_pattern_template += [xray_pattern_template_dE]
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(ana_ccd_bins)):
        xray_pattern_template[0][0].add(xray_pattern_sci_source[ebin][ccdid])
        xray_pattern_template[0][0].add(xray_pattern_sci_ring[ebin][ccdid],factor=-1.)

for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(ana_ccd_bins)):
        if ebin==0 and ccdid==0: continue
        xray_pattern_template[ebin][ccdid].add(xray_pattern_template[0][0])


print ('prepare off sample time cuts')

output_all_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)
output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)
output_all_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)
output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)

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

for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
    
        area_pix_frac_all_fov = off_image_all_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_all_cor = off_image_all_cor[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_fov = off_image_fwc_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_cor = off_image_fwc_cor[ebin][ccdid].get_pixel_fraction()
        
        time_pix_frac_all_fov = off_lightcurve_all_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_all_fov = off_duration_all_fov*time_pix_frac_all_fov
        time_pix_frac_fwc_fov = off_lightcurve_fwc_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_fov = off_duration_fwc_fov*time_pix_frac_fwc_fov
        
        time_pix_frac_all_cor = off_lightcurve_all_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_all_cor = off_duration_all_cor*time_pix_frac_all_cor
        time_pix_frac_fwc_cor = off_lightcurve_fwc_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_cor = off_duration_fwc_cor*time_pix_frac_fwc_cor
    
        off_lightcurve_all_cor[ebin][ccdid].scale(area_pix_frac_all_fov/area_pix_frac_all_cor)
        off_lightcurve_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*time_expo_all_fov/time_expo_fwc_fov)

off_lightcurve_all_fov_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
off_lightcurve_all_cor_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        off_lightcurve_all_fov_sum.add(off_lightcurve_all_fov[ebin][ccdid])
        off_lightcurve_all_cor_sum.add(off_lightcurve_all_cor[ebin][ccdid])

mask_lc = make_timecut_mask(off_lightcurve_all_fov_sum,off_lightcurve_all_cor_sum) 


print ('apply off sample space and time masks')

output_sci_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_array)
output_sci_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_array)

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

output_spf_fov = read_event_file(off_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_array)
output_spf_cor = read_event_file(off_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_array)

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

output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_array)
output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_array)

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

off_image_sci_fov_sum = MyArray2D()
off_image_sci_cor_sum = MyArray2D()
off_image_spf_fov_sum = MyArray2D()
off_image_spf_cor_sum = MyArray2D()
off_image_fwc_fov_sum = MyArray2D()
off_image_fwc_cor_sum = MyArray2D()
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        off_image_sci_fov_sum.add(off_image_sci_fov[ebin][ccdid])
        off_image_sci_cor_sum.add(off_image_sci_cor[ebin][ccdid])
        off_image_spf_fov_sum.add(off_image_spf_fov[ebin][ccdid])
        off_image_spf_cor_sum.add(off_image_spf_cor[ebin][ccdid])
        off_image_fwc_fov_sum.add(off_image_fwc_fov[ebin][ccdid])
        off_image_fwc_cor_sum.add(off_image_fwc_cor[ebin][ccdid])

area_pix_frac_sci_fov = off_image_sci_fov_sum.get_pixel_fraction()
area_pix_frac_sci_cor = off_image_sci_cor_sum.get_pixel_fraction()
area_pix_frac_spf_fov = off_image_spf_fov_sum.get_pixel_fraction()
area_pix_frac_spf_cor = off_image_spf_cor_sum.get_pixel_fraction()
area_pix_frac_fwc_fov = off_image_fwc_fov_sum.get_pixel_fraction()
area_pix_frac_fwc_cor = off_image_fwc_cor_sum.get_pixel_fraction()


for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
    
        time_pix_frac_sci_fov = off_lightcurve_sci_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_sci_fov = off_duration_sci_fov*time_pix_frac_sci_fov
        time_pix_frac_sci_cor = off_lightcurve_sci_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_sci_cor = off_duration_sci_cor*time_pix_frac_sci_cor
        
        time_pix_frac_spf_fov = off_lightcurve_spf_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_spf_fov = off_duration_spf_fov*time_pix_frac_spf_fov
        time_pix_frac_spf_cor = off_lightcurve_spf_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_spf_cor = off_duration_spf_cor*time_pix_frac_spf_cor
        
        time_pix_frac_fwc_fov = off_lightcurve_fwc_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_fov = off_duration_fwc_fov*time_pix_frac_fwc_fov
        time_pix_frac_fwc_cor = off_lightcurve_fwc_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_cor = off_duration_fwc_cor*time_pix_frac_fwc_cor
    
        off_lightcurve_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        off_lightcurve_spf_fov[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
        off_lightcurve_spf_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
        off_lightcurve_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        off_lightcurve_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
        off_spectrum_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        off_spectrum_spf_fov[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
        off_spectrum_spf_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
        off_spectrum_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        off_spectrum_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
        
        off_pattern_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        off_pattern_spf_fov[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_spf_fov*time_expo_sci_fov/time_expo_spf_fov)
        off_pattern_spf_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_spf_cor*time_expo_sci_fov/time_expo_spf_cor)
        off_pattern_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        off_pattern_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)

off_lightcurve_sci_fov_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
off_lightcurve_sci_cor_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        off_lightcurve_sci_fov_sum.add(off_lightcurve_sci_fov[ebin][ccdid])
        off_lightcurve_sci_cor_sum.add(off_lightcurve_sci_cor[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_lightcurve_all_fov_sum.xaxis,off_lightcurve_all_fov_sum.yaxis,yerr=off_lightcurve_all_fov_sum.yerr,color='k',label='All FoV')
axbig.errorbar(off_lightcurve_sci_fov_sum.xaxis,off_lightcurve_sci_fov_sum.yaxis,yerr=off_lightcurve_sci_fov_sum.yerr,color='red',label='Sci FoV')
axbig.errorbar(off_lightcurve_all_cor_sum.xaxis,off_lightcurve_all_cor_sum.yaxis,yerr=off_lightcurve_all_cor_sum.yerr,color='green',label='All Cor')
axbig.set_yscale('log')
axbig.set_xlabel('Time')
axbig.legend(loc='best')
fig.savefig("../output_plots/lightcurve_off_timecut_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

off_spectrum_sci_fov_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
off_spectrum_spf_fov_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
off_spectrum_fwc_fov_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        off_spectrum_sci_fov_sum.add(off_spectrum_sci_fov[ebin][ccdid])
        off_spectrum_spf_fov_sum.add(off_spectrum_spf_fov[ebin][ccdid])
        off_spectrum_fwc_fov_sum.add(off_spectrum_fwc_fov[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_spectrum_sci_fov_sum.xaxis,off_spectrum_sci_fov_sum.yaxis,yerr=off_spectrum_sci_fov_sum.yerr,label='Sci FoV')
axbig.errorbar(off_spectrum_spf_fov_sum.xaxis,off_spectrum_spf_fov_sum.yaxis,yerr=off_spectrum_spf_fov_sum.yerr,label='SPF FoV')
axbig.errorbar(off_spectrum_fwc_fov_sum.xaxis,off_spectrum_fwc_fov_sum.yaxis,yerr=off_spectrum_fwc_fov_sum.yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_spectrum_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

off_spectrum_sci_cor_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
off_spectrum_spf_cor_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
off_spectrum_fwc_cor_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        off_spectrum_sci_cor_sum.add(off_spectrum_sci_cor[ebin][ccdid])
        off_spectrum_spf_cor_sum.add(off_spectrum_spf_cor[ebin][ccdid])
        off_spectrum_fwc_cor_sum.add(off_spectrum_fwc_cor[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_spectrum_sci_cor_sum.xaxis,off_spectrum_sci_cor_sum.yaxis,yerr=off_spectrum_sci_cor_sum.yerr,label='Sci Cor')
axbig.errorbar(off_spectrum_spf_cor_sum.xaxis,off_spectrum_spf_cor_sum.yaxis,yerr=off_spectrum_spf_cor_sum.yerr,label='SPF Cor')
axbig.errorbar(off_spectrum_fwc_cor_sum.xaxis,off_spectrum_fwc_cor_sum.yaxis,yerr=off_spectrum_fwc_cor_sum.yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_spectrum_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

off_pattern_sci_fov_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
off_pattern_spf_fov_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
off_pattern_fwc_fov_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        off_pattern_sci_fov_sum.add(off_pattern_sci_fov[ebin][ccdid])
        off_pattern_spf_fov_sum.add(off_pattern_spf_fov[ebin][ccdid])
        off_pattern_fwc_fov_sum.add(off_pattern_fwc_fov[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_pattern_sci_fov_sum.xaxis,off_pattern_sci_fov_sum.yaxis,yerr=off_pattern_sci_fov_sum.yerr,label='Sci FoV')
axbig.errorbar(off_pattern_spf_fov_sum.xaxis,off_pattern_spf_fov_sum.yaxis,yerr=off_pattern_spf_fov_sum.yerr,label='SPF FoV')
axbig.errorbar(off_pattern_fwc_fov_sum.xaxis,off_pattern_fwc_fov_sum.yaxis,yerr=off_pattern_fwc_fov_sum.yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_pattern_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

off_pattern_sci_cor_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
off_pattern_spf_cor_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
off_pattern_fwc_cor_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        off_pattern_sci_cor_sum.add(off_pattern_sci_cor[ebin][ccdid])
        off_pattern_spf_cor_sum.add(off_pattern_spf_cor[ebin][ccdid])
        off_pattern_fwc_cor_sum.add(off_pattern_fwc_cor[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(off_pattern_sci_cor_sum.xaxis,off_pattern_sci_cor_sum.yaxis,yerr=off_pattern_sci_cor_sum.yerr,label='Sci Cor')
axbig.errorbar(off_pattern_spf_cor_sum.xaxis,off_pattern_spf_cor_sum.yaxis,yerr=off_pattern_spf_cor_sum.yerr,label='SPF Cor')
axbig.errorbar(off_pattern_fwc_cor_sum.xaxis,off_pattern_fwc_cor_sum.yaxis,yerr=off_pattern_fwc_cor_sum.yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/off_pattern_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()


sp_pattern_template = []
sp_spectrum_template = []
for ebin in range(0,len(energy_array)-1):
    sp_pattern_template_dE = []
    sp_spectrum_template_dE = []
    for ccdid in range(0,len(all_ccd_bins)):
        sp_pattern_template_dE += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
        sp_spectrum_template_dE += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
    sp_pattern_template += [sp_pattern_template_dE]
    sp_spectrum_template += [sp_spectrum_template_dE]
    
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        sp_pattern_template[ebin][ccdid].add(off_pattern_spf_fov_sum)
        sp_pattern_template[ebin][ccdid].add(off_pattern_sci_fov_sum,factor=-1.)
        sp_spectrum_template[ebin][ccdid].add(off_spectrum_spf_fov_sum)
        sp_spectrum_template[ebin][ccdid].add(off_spectrum_sci_fov_sum,factor=-1.)


print ('prepare on sample time cuts')

output_all_fov = read_event_file(on_sci_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)
output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)
output_all_cor = read_event_file(on_sci_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)
output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,mask_lc=None,mask_map=None,evt_filter='',energy_range=energy_array)

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


for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
    
        area_pix_frac_all_fov = on_image_all_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_all_cor = on_image_all_cor[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_fov = on_image_fwc_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_cor = on_image_fwc_cor[ebin][ccdid].get_pixel_fraction()
        
        time_pix_frac_all_fov = on_lightcurve_all_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_all_fov = on_duration_all_fov*time_pix_frac_all_fov
        time_pix_frac_fwc_fov = on_lightcurve_fwc_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
        
        time_pix_frac_all_cor = on_lightcurve_all_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_all_cor = on_duration_all_cor*time_pix_frac_all_cor
        time_pix_frac_fwc_cor = on_lightcurve_fwc_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor
    
        on_lightcurve_all_cor[ebin][ccdid].scale(area_pix_frac_all_fov/area_pix_frac_all_cor)
        on_lightcurve_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*time_expo_all_fov/time_expo_fwc_fov)

on_lightcurve_all_fov_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
on_lightcurve_all_cor_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        on_lightcurve_all_fov_sum.add(on_lightcurve_all_fov[ebin][ccdid])
        on_lightcurve_all_cor_sum.add(on_lightcurve_all_cor[ebin][ccdid])

mask_lc = make_timecut_mask(on_lightcurve_all_fov_sum,on_lightcurve_all_cor_sum) 

time_pix_frac_mask = mask_lc.get_pixel_fraction()
if time_pix_frac_mask>0.8:
    mask_lc = None
if not SP_flare_mask:
    mask_lc = None


print ('apply space and time masks')

output_sci_fov = read_event_file(on_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_array,ccd=ana_ccd_bins)
output_sci_cor = read_event_file(on_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_array,ccd=ana_ccd_bins)

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

output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_array,ccd=ana_ccd_bins)
output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_array,ccd=ana_ccd_bins)

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


for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(ana_ccd_bins)):
    
        area_pix_frac_sci_fov = on_image_sci_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_sci_cor = on_image_sci_cor[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_fov = on_image_fwc_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_cor = on_image_fwc_cor[ebin][ccdid].get_pixel_fraction()
        
        time_pix_frac_sci_fov = on_lightcurve_sci_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_sci_fov = on_duration_sci_fov*time_pix_frac_sci_fov
        time_pix_frac_sci_cor = on_lightcurve_sci_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_sci_cor = on_duration_sci_cor*time_pix_frac_sci_cor
        
        time_pix_frac_fwc_fov = on_lightcurve_fwc_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
        time_pix_frac_fwc_cor = on_lightcurve_fwc_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor
    
        on_lightcurve_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        on_lightcurve_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        on_lightcurve_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
        on_spectrum_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        on_spectrum_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        on_spectrum_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
        
        on_pattern_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        on_pattern_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        on_pattern_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)

on_lightcurve_sci_fov_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
on_lightcurve_sci_cor_sum = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        on_lightcurve_sci_fov_sum.add(on_lightcurve_sci_fov[ebin][ccdid])
        on_lightcurve_sci_cor_sum.add(on_lightcurve_sci_cor[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_lightcurve_all_fov_sum.xaxis,on_lightcurve_all_fov_sum.yaxis,yerr=on_lightcurve_all_fov_sum.yerr,color='k',label='All FoV')
axbig.errorbar(on_lightcurve_sci_fov_sum.xaxis,on_lightcurve_sci_fov_sum.yaxis,yerr=on_lightcurve_sci_fov_sum.yerr,color='red',label='Sci FoV')
axbig.errorbar(on_lightcurve_all_cor_sum.xaxis,on_lightcurve_all_cor_sum.yaxis,yerr=on_lightcurve_all_cor_sum.yerr,color='green',label='All Cor')
axbig.set_yscale('log')
axbig.set_xlabel('Time')
axbig.legend(loc='best')
fig.savefig("../output_plots/lightcurve_on_timecut_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

on_spectrum_sci_fov_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
on_spectrum_fwc_fov_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        on_spectrum_sci_fov_sum.add(on_spectrum_sci_fov[ebin][ccdid])
        on_spectrum_fwc_fov_sum.add(on_spectrum_fwc_fov[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_spectrum_sci_fov_sum.xaxis,on_spectrum_sci_fov_sum.yaxis,yerr=on_spectrum_sci_fov_sum.yerr,label='Sci FoV')
axbig.errorbar(on_spectrum_fwc_fov_sum.xaxis,on_spectrum_fwc_fov_sum.yaxis,yerr=on_spectrum_fwc_fov_sum.yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_spectrum_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

on_spectrum_sci_cor_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
on_spectrum_fwc_cor_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        on_spectrum_sci_cor_sum.add(on_spectrum_sci_cor[ebin][ccdid])
        on_spectrum_fwc_cor_sum.add(on_spectrum_fwc_cor[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_spectrum_sci_cor_sum.xaxis,on_spectrum_sci_cor_sum.yaxis,yerr=on_spectrum_sci_cor_sum.yerr,label='Sci Cor')
axbig.errorbar(on_spectrum_fwc_cor_sum.xaxis,on_spectrum_fwc_cor_sum.yaxis,yerr=on_spectrum_fwc_cor_sum.yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_spectrum_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

on_pattern_sci_fov_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
on_pattern_fwc_fov_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        on_pattern_sci_fov_sum.add(on_pattern_sci_fov[ebin][ccdid])
        on_pattern_fwc_fov_sum.add(on_pattern_fwc_fov[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_pattern_sci_fov_sum.xaxis,on_pattern_sci_fov_sum.yaxis,yerr=on_pattern_sci_fov_sum.yerr,label='Sci FoV')
axbig.errorbar(on_pattern_fwc_fov_sum.xaxis,on_pattern_fwc_fov_sum.yaxis,yerr=on_pattern_fwc_fov_sum.yerr,label='FWC FoV')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_pattern_fov_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

on_pattern_sci_cor_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
on_pattern_fwc_cor_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(all_ccd_bins)):
        on_pattern_sci_cor_sum.add(on_pattern_sci_cor[ebin][ccdid])
        on_pattern_fwc_cor_sum.add(on_pattern_fwc_cor[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_pattern_sci_cor_sum.xaxis,on_pattern_sci_cor_sum.yaxis,yerr=on_pattern_sci_cor_sum.yerr,label='Sci Cor')
axbig.errorbar(on_pattern_fwc_cor_sum.xaxis,on_pattern_fwc_cor_sum.yaxis,yerr=on_pattern_fwc_cor_sum.yerr,label='FWC Cor')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
axbig.set_yscale('log')
fig.savefig("../output_plots/on_pattern_cor_raw_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()



print ('predict SP background')

qpb_pattern_template = []
qpb_spectrum_template = []
for ebin in range(0,len(energy_array)-1):
    qpb_pattern_template_dE = []
    qpb_spectrum_template_dE = []
    for ccdid in range(0,len(ana_ccd_bins)):
        qpb_pattern_template_dE += [MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)]
        qpb_pattern_template_dE[ccdid].add(on_pattern_sci_cor[ebin][ccdid])
        qpb_spectrum_template_dE += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
        qpb_spectrum_template_dE[ccdid].add(on_spectrum_sci_cor[ebin][ccdid])
    qpb_pattern_template += [qpb_pattern_template_dE]
    qpb_spectrum_template += [qpb_spectrum_template_dE]

on_spectrum_sci_fov_bkg = [] 
spf_scaling_array = []
xray_scaling_array = []
for ebin in range(0,len(energy_array)-1):
    on_spectrum_sci_fov_bkg_dE = [] 
    spf_scaling_array_dE = []
    xray_scaling_array_dE = []
    for ccdid in range(0,len(ana_ccd_bins)):
    
        spf_scaling, xray_scaling = fit_pattern(on_pattern_sci_fov[ebin][ccdid],xray_pattern_template[ebin][ccdid],sp_pattern_template[ebin][ccdid],qpb_pattern_template[ebin][ccdid])
        print ('xray_scaling = %s'%(xray_scaling))
        print ('spf_scaling = %s'%(spf_scaling))
        spf_scaling_array_dE += [spf_scaling]
        xray_scaling_array_dE += [xray_scaling]
    
        xray_pattern_template[ebin][ccdid].scale(xray_scaling)
        sp_pattern_template[ebin][ccdid].scale(spf_scaling)
    
        sp_spectrum_template[ebin][ccdid].scale(spf_scaling)
        on_spectrum_sci_fov_bkg_dE += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
        on_spectrum_sci_fov_bkg_dE[ccdid].add(sp_spectrum_template[ebin][ccdid])
        on_spectrum_sci_fov_bkg_dE[ccdid].add(qpb_spectrum_template[ebin][ccdid])
    on_spectrum_sci_fov_bkg += [on_spectrum_sci_fov_bkg_dE]
    spf_scaling_array += [spf_scaling_array_dE]
    xray_scaling_array += [xray_scaling_array_dE]


sp_spectrum_template_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
qpb_spectrum_template_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
sp_pattern_template_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
qpb_pattern_template_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
xray_pattern_template_sum = MyArray1D(bin_start=0,bin_end=25,pixel_scale=1)
on_spectrum_sci_fov_bkg_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(ana_ccd_bins)):
        sp_spectrum_template_sum.add(sp_spectrum_template[ebin][ccdid])
        qpb_spectrum_template_sum.add(qpb_spectrum_template[ebin][ccdid])
        sp_pattern_template_sum.add(sp_pattern_template[ebin][ccdid])
        qpb_pattern_template_sum.add(qpb_pattern_template[ebin][ccdid])
        xray_pattern_template_sum.add(xray_pattern_template[ebin][ccdid])
        on_spectrum_sci_fov_bkg_sum.add(on_spectrum_sci_fov_bkg[ebin][ccdid])

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_pattern_sci_fov_sum.xaxis,on_pattern_sci_fov_sum.yaxis,yerr=on_pattern_sci_fov_sum.yerr,color='k',label='Data')
axbig.errorbar(xray_pattern_template_sum.xaxis,xray_pattern_template_sum.yaxis,yerr=xray_pattern_template_sum.yerr,color='g',label='X-ray')
axbig.errorbar(qpb_pattern_template_sum.xaxis,qpb_pattern_template_sum.yaxis,yerr=qpb_pattern_template_sum.yerr,color='b',label='QPB')
axbig.bar(sp_pattern_template_sum.xaxis,sp_pattern_template_sum.yaxis,color='orange',label='SP Flare')
axbig.set_yscale('log')
axbig.set_xlabel('PATTERN')
axbig.legend(loc='best')
fig.savefig("../output_plots/pattern_on_fov_scaled_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()


fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_spectrum_sci_fov_sum.xaxis,on_spectrum_sci_fov_sum.yaxis,yerr=on_spectrum_sci_fov_sum.yerr,color='k',label='Data')
axbig.plot(on_spectrum_sci_fov_bkg_sum.xaxis,on_spectrum_sci_fov_bkg_sum.yaxis,color='red',label='Bkg')
axbig.plot(sp_spectrum_template_sum.xaxis,sp_spectrum_template_sum.yaxis,label='SP')
axbig.plot(qpb_spectrum_template_sum.xaxis,qpb_spectrum_template_sum.yaxis,label='QPB')
axbig.set_yscale('log')
axbig.set_xlabel('Channel')
axbig.legend(loc='best')
fig.savefig("../output_plots/spectrum_on_fit_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()


xray_pattern_template_sum.normalize()
qpb_pattern_template_sum.normalize()
sp_pattern_template_sum.normalize()

fig.clf()
axbig = fig.add_subplot()
axbig.plot(xray_pattern_template_sum.xaxis,xray_pattern_template_sum.yaxis,label='X-ray')
axbig.plot(qpb_pattern_template_sum.xaxis,qpb_pattern_template_sum.yaxis,label='QPB')
axbig.plot(sp_pattern_template_sum.xaxis,sp_pattern_template_sum.yaxis,label='SP')
axbig.set_yscale('log')
axbig.set_xlabel('Pattern')
axbig.legend(loc='best')
fig.savefig("../output_plots/pattern_template_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()



print ('get science data')

output_sci_fov = read_event_file(on_sci_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_array,ccd=ana_ccd_bins)
output_sci_cor = read_event_file(on_sci_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_array,ccd=ana_ccd_bins)

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

output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_array,ccd=ana_ccd_bins)
output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_array,ccd=ana_ccd_bins)

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

for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(ana_ccd_bins)):
    
        area_pix_frac_sci_fov = on_image_sci_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_sci_cor = on_image_sci_cor[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_fov = on_image_fwc_fov[ebin][ccdid].get_pixel_fraction()
        area_pix_frac_fwc_cor = on_image_fwc_cor[ebin][ccdid].get_pixel_fraction()
        
        time_pix_frac_sci_fov = on_lightcurve_sci_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_sci_fov = on_duration_sci_fov*time_pix_frac_sci_fov
        time_pix_frac_sci_cor = on_lightcurve_sci_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_sci_cor = on_duration_sci_cor*time_pix_frac_sci_cor
        
        time_pix_frac_fwc_fov = on_lightcurve_fwc_fov[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
        time_pix_frac_fwc_cor = on_lightcurve_fwc_cor[ebin][ccdid].get_pixel_fraction()
        time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor
    
        on_lightcurve_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        on_lightcurve_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        on_lightcurve_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
    
        on_spectrum_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        on_spectrum_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        on_spectrum_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)
        
        on_pattern_sci_cor[ebin][ccdid].scale(area_pix_frac_sci_fov/area_pix_frac_sci_cor)
        on_pattern_fwc_fov[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_fov*time_expo_sci_fov/time_expo_fwc_fov)
        on_pattern_fwc_cor[ebin][ccdid].scale(fwc_2_sci_ratio*area_pix_frac_sci_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_fwc_cor)


qpb_detx_template = []
sp_detx_template = []
qpb_image_template = []
sp_image_template = []
for ebin in range(0,len(energy_array)-1):
    qpb_detx_template_dE = []
    sp_detx_template_dE = []
    qpb_image_template_dE = []
    sp_image_template_dE = []
    for ccdid in range(0,len(ana_ccd_bins)):
        qpb_detx_template_dE += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
        qpb_detx_template_dE[ccdid].add(on_detx_fwc_fov[ebin][ccdid])
        sp_detx_template_dE += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
        sp_detx_template_dE[ccdid].add(on_detx_fwc_fov[ebin][ccdid])
        qpb_image_template_dE += [MyArray2D()]
        qpb_image_template_dE[ccdid].add(on_image_fwc_fov[ebin][ccdid])
        sp_image_template_dE += [MyArray2D()]
        sp_image_template_dE[ccdid].add(on_image_fwc_fov[ebin][ccdid])
    qpb_detx_template += [qpb_detx_template_dE]
    sp_detx_template += [sp_detx_template_dE]
    qpb_image_template += [qpb_image_template_dE]
    sp_image_template += [sp_image_template_dE]

on_detx_sci_fov_bkg = []
on_image_sci_fov_bkg = []
for ebin in range(0,len(energy_array)-1):
    on_detx_sci_fov_bkg_dE = []
    on_image_sci_fov_bkg_dE = []
    for ccdid in range(0,len(ana_ccd_bins)):
    
        spf_scaling = spf_scaling_array[ebin][ccdid]
        xray_scaling = xray_scaling_array[ebin][ccdid]
    
        qpb_detx_template[ebin][ccdid].scale(qpb_spectrum_template[ebin][ccdid].integral()/qpb_detx_template[ebin][ccdid].integral())
        sp_detx_template[ebin][ccdid].scale(sp_spectrum_template[ebin][ccdid].integral()/sp_detx_template[ebin][ccdid].integral())
        on_detx_sci_fov_bkg_dE += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
        on_detx_sci_fov_bkg_dE[ccdid].add(sp_detx_template[ebin][ccdid])
        on_detx_sci_fov_bkg_dE[ccdid].add(qpb_detx_template[ebin][ccdid])
    
        qpb_image_template[ebin][ccdid].scale(qpb_spectrum_template[ebin][ccdid].integral()/qpb_image_template[ebin][ccdid].integral())
        sp_image_template[ebin][ccdid].scale(sp_spectrum_template[ebin][ccdid].integral()/sp_image_template[ebin][ccdid].integral())
        on_image_sci_fov_bkg_dE += [MyArray2D()]
        on_image_sci_fov_bkg_dE[ccdid].add(sp_image_template[ebin][ccdid])
        on_image_sci_fov_bkg_dE[ccdid].add(qpb_image_template[ebin][ccdid])

    on_detx_sci_fov_bkg += [on_detx_sci_fov_bkg_dE]
    on_image_sci_fov_bkg += [on_image_sci_fov_bkg_dE]

sp_detx_template_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
qpb_detx_template_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
on_detx_sci_fov_bkg_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
on_detx_sci_fov_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
sp_image_template_sum = MyArray2D()
qpb_image_template_sum = MyArray2D()
on_image_sci_fov_bkg_sum = MyArray2D()
on_image_sci_fov_sum = MyArray2D()
for ebin in range(0,len(energy_array)-1):
    for ccdid in range(0,len(ana_ccd_bins)):
        sp_detx_template_sum.add(sp_detx_template[ebin][ccdid])
        qpb_detx_template_sum.add(qpb_detx_template[ebin][ccdid])
        on_detx_sci_fov_bkg_sum.add(on_detx_sci_fov_bkg[ebin][ccdid])
        on_detx_sci_fov_sum.add(on_detx_sci_fov[ebin][ccdid])
        sp_image_template_sum.add(sp_image_template[ebin][ccdid])
        qpb_image_template_sum.add(qpb_image_template[ebin][ccdid])
        on_image_sci_fov_bkg_sum.add(on_image_sci_fov_bkg[ebin][ccdid])
        on_image_sci_fov_sum.add(on_image_sci_fov[ebin][ccdid])


fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(on_detx_sci_fov_sum.xaxis,on_detx_sci_fov_sum.yaxis,yerr=on_detx_sci_fov_sum.yerr,color='k',label='Data')
axbig.plot(on_detx_sci_fov_bkg_sum.xaxis,on_detx_sci_fov_bkg_sum.yaxis,color='red',label='Bkg')
axbig.plot(sp_detx_template_sum.xaxis,sp_detx_template_sum.yaxis,label='SP')
axbig.plot(qpb_detx_template_sum.xaxis,qpb_detx_template_sum.yaxis,label='QPB')
axbig.set_yscale('log')
axbig.set_ylim(bottom=1)
axbig.set_xlabel('DETX')
axbig.legend(loc='best')
fig.savefig("../output_plots/detx_on_fit_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()


fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = on_image_sci_fov_sum.xaxis.min()
xmax = on_image_sci_fov_sum.xaxis.max()
ymin = on_image_sci_fov_sum.yaxis.min()
ymax = on_image_sci_fov_sum.yaxis.max()
axbig.imshow(on_image_sci_fov_sum.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("../output_plots/skymap_on_sci_fov_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = on_image_sci_fov_bkg_sum.xaxis.min()
xmax = on_image_sci_fov_bkg_sum.xaxis.max()
ymin = on_image_sci_fov_bkg_sum.yaxis.min()
ymax = on_image_sci_fov_bkg_sum.yaxis.max()
axbig.imshow(on_image_sci_fov_bkg_sum.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("../output_plots/skymap_on_sci_fov_bkg_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

on_image_sci_fov_excess = MyArray2D()
on_image_sci_fov_excess.add(on_image_sci_fov_sum)
on_image_sci_fov_excess.add(on_image_sci_fov_bkg_sum,factor=-1.)

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
fig.savefig("../output_plots/skymap_on_sci_fov_excess_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()


















end_time = time.time()
elapsed_time = end_time - start_time
print('elapsed_time = %s'%(elapsed_time))

exit()


