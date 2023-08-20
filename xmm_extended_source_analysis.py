
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import wcs
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
from scipy.optimize import minimize
import time

import common_functions

MyArray1D = common_functions.MyArray1D
MyArray2D = common_functions.MyArray2D
pattern_low = common_functions.pattern_low
pattern_high = common_functions.pattern_high
pattern_scale = common_functions.pattern_scale
ch_low = common_functions.ch_low
ch_high = common_functions.ch_high
ch_scale = common_functions.ch_scale
t_low = common_functions.t_low
t_high = common_functions.t_high
t_scale = common_functions.t_scale
detx_low = common_functions.detx_low
detx_high = common_functions.detx_high
detx_scale = common_functions.detx_scale
detr_low = common_functions.detr_low
detr_high = common_functions.detr_high
detr_scale = common_functions.detr_scale
sky_ra_low = common_functions.sky_ra_low
sky_ra_high = common_functions.sky_ra_high
sky_dec_low = common_functions.sky_dec_low
sky_dec_high = common_functions.sky_dec_high
sky_scale = common_functions.sky_scale


fig, ax = plt.subplots()
figsize_x = 8.
figsize_y = 6.
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

on_sample = sys.argv[1]
on_obsID = sys.argv[2]

diagnostic_plots = True
#diagnostic_plots = False
fast_test = True
#fast_test = False

#detector = 'mos1'
#detector = 'mos2'
detector = sys.argv[3]

on_filter = sys.argv[4]
#on_filter = 'fov'
#on_filter = 'reg1'
#on_filter = 'reg2'
#on_filter = 'reg3'
#on_filter = 'reg4'

#energy_array = [200,1000,2000,3000,5000,8000,12000]
energy_array = [1000,2000,4000,6000,8000,10000,12000]
#energy_array = [2000,12000]
#ana_ccd_bins = [1,2,3,4,5,6,7]
ana_ccd_bins = [0]

output_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/'+on_obsID

source_radius = 1000
ring_inner_radius = 1500
ring_outer_radius = 6000


src_ra = 0
src_dec = 0
if 'Cas_A' in on_sample:
    #src_ra = 350.8500000
    #src_dec = 58.8150000
    src_ra = 350.95
    src_dec = 58.71

src_detx = 0
src_dety = 0
if 'ID0400210101' in on_obsID:
    src_detx = 13000
    src_dety = 13000

def SaveImage(image_input,filename):

    new_image_data = image_input.zaxis

    my_hdu = fits.PrimaryHDU(new_image_data)
    my_hdul = fits.HDUList([my_hdu])
    my_hdul.writeto('%s/%s.fits'%(output_dir,filename), overwrite=True)
    
def SaveSpectrum(spec_input,filename):

    count = fits.Column(name='Count', array=spec_input.yaxis, format='D')
    energy = fits.Column(name='Energy', array=spec_input.xaxis, format='D')
    my_table = fits.BinTableHDU.from_columns([count,energy],name='SPECTRUM')
    my_table.writeto('%s/%s.fits'%(output_dir,filename), overwrite=True)

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

#global_att_radec, global_mtx_conv_sky_2_det, global_mtx_conv_det_2_sky = get_conversion_mtx()
#
#print ('global_att_radec:')
#print (global_att_radec)
#
#pix_size = (1./60.)/1200. # degree
#if not 'extragalactic' in on_sample:
#    detr_low = pow(pow(src_ra-global_att_radec[0],2)+pow(src_dec-global_att_radec[1],2),0.5)-20000*pix_size
#    detr_low = max(0.,detr_low)
#    #detr_high = pow(pow(src_ra-global_att_radec[0],2)+pow(src_dec-global_att_radec[1],2),0.5)+20000*pix_size
#    detr_high = detr_low+60000*pix_size
#    detr_scale = detr_scale*pix_size
#    sky_ra_low = global_att_radec[0]-0.25
#    sky_ra_high = global_att_radec[0]+0.25
#    sky_dec_low = global_att_radec[1]-0.25
#    sky_dec_high = global_att_radec[1]+0.25

fwc_2_sci_ratio = 0.3


def read_event_file(filename,rmf_name,mask_lc=None,mask_map=None,write_events=False,evt_filter='',energy_range=[200,12000],ccd_id=0):

    # how to read events:
    # https://docs.astropy.org/en/stable/generated/examples/io/fits-tables.html#accessing-data-stored-as-a-table-in-a-multi-extension-fits-mef-file

    #rmf_table = Table.read(rmf_name, hdu=1)

    print ('read file %s'%(filename))
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
    pattern_array = []
    spectrum_array = []
    detx_array = []
    image_array = []
    for ch in range(0,len(energy_range)):
        pattern_array += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
        spectrum_array += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
        detx_array += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
        image_array += [MyArray2D()]

    time_start = events[0]['TIME']
    time_end = events[len(events)-1]['TIME']
    obs_duration = time_end-time_start

    evt_detx_list = []
    evt_dety_list = []
    evt_pi_list = []

    for evt in range(0,len(events)):

        if fast_test:
            if (evt%10)!=0: continue

        evt_time = events[evt]['TIME']
        evt_pattern = events[evt]['PATTERN']

        #evt_pha = events[evt]['PHA']
        #if evt_pha>=2400: continue
        #energy_lo = rmf_table[evt_pha]['ENERG_LO']
        #energy_hi = rmf_table[evt_pha]['ENERG_HI']
        #evt_energy = 0.5*(energy_lo+energy_hi)*1000.

        evt_pi = events[evt]['PI']
        #if evt_pi<energy_range[0]: continue
        #if evt_pi>energy_range[len(energy_range)-1]: continue

        channel = 0
        for ch in range(1,len(energy_range)):
            if evt_pi>=energy_range[ch-1] and evt_pi<energy_range[ch]:
                channel = ch

        evt_detx = events[evt]['DETX']
        evt_dety = events[evt]['DETY']
        #evt_x = 0.
        #evt_y = 0.
        #if not 'extragalactic' in on_sample:
        #    evt_detxy = np.array([evt_detx,evt_dety])
        #    evt_xy = np.matmul(global_mtx_conv_det_2_sky,evt_detxy)
        #    evt_x = evt_xy[0,0]+global_att_radec[0]
        #    evt_y = evt_xy[0,1]+global_att_radec[1]

        if evt_pi<energy_range[0]: continue
        if evt_pi>=energy_range[len(energy_range)-1]: continue

        if evt_pattern>4:continue
        #if evt_pattern==0:continue

        if 'source' in evt_filter:
            if abs(evt_detx)>source_radius: continue
            if abs(evt_dety)>source_radius: continue
            #if abs(evt_detx)>6000: continue
            #if abs(evt_dety)>6000: continue
        elif 'ring' in evt_filter:
            if abs(evt_detx)<ring_inner_radius and abs(evt_dety)<ring_inner_radius: continue
            if abs(evt_detx)>ring_outer_radius: continue
            if abs(evt_dety)>ring_outer_radius: continue
            #if abs(evt_detx)<6000 and abs(evt_dety)<6000: continue
            #if abs(evt_detx)>12000: continue
            #if abs(evt_dety)>12000: continue
        elif 'halo' in evt_filter:
            #distance = pow(pow(evt_x-src_ra,2)+pow(evt_y-src_dec,2),0.5)
            #if distance<0.06: continue
            distance = pow(pow(evt_detx-src_detx,2)+pow(evt_dety-src_dety,2),0.5)
            if distance<6000: continue
            if distance>12000: continue

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

        evt_detx_list += [evt_detx]
        evt_dety_list += [evt_dety]
        evt_pi_list += [evt_pi]

        evt_count += 1.
        lightcurve_array.fill((evt_time-time_start)/(time_end-time_start))
        pattern_array[0].fill(evt_pattern)
        spectrum_array[0].fill(evt_pi)
        detx_array[0].fill(evt_detx)
        image_array[0].fill(evt_detx,evt_dety)
        pattern_array[channel].fill(evt_pattern)
        spectrum_array[channel].fill(evt_pi)
        detx_array[channel].fill(evt_detx)
        image_array[channel].fill(evt_detx,evt_dety)

    if write_events:
        col_detx = fits.Column(name='DETX', array=evt_detx_list, format='D')
        col_dety = fits.Column(name='DETY', array=evt_dety_list, format='D')
        col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
        my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_pi],name='EVENT')
        my_table.writeto('%s/sci_events_%s_ccd%s.fits'%(output_dir,detector,ccd_id), overwrite=True)

    return [obs_duration, evt_count, lightcurve_array, pattern_array, spectrum_array, detx_array, image_array]

def func_powerlaw(x,A,B):

    return A*pow(x/1000.,B)

def fit_pattern(energy_range,data_pattern,xray_pattern,spf_pattern,qpb_pattern):

    xray_qpb_ratio = []
    spfree_data_qpb_ratio = []
    data_qpb_ratio = []
    for ch in range(1,len(energy_range)):
        sp_free_data = 0.
        sp_free_qpb = 0.
        for idx in range(0,len(data_pattern[ch].xaxis)):
            #if idx!=2 and idx!=4: continue
            if idx!=2: continue
            sp_free_data += data_pattern[ch].yaxis[idx]
            sp_free_qpb  += qpb_pattern[ch].yaxis[idx]
        xray_qpb_ratio += [max(sp_free_data/sp_free_qpb-1.,0.)]
        spfree_data_qpb_ratio += [sp_free_data/sp_free_qpb]
        data = 0.
        qpb = 0.
        for idx in range(0,len(data_pattern[ch].xaxis)):
            if idx==2: continue
            if idx==4: continue
            if idx>4: continue
            data += data_pattern[ch].yaxis[idx]
            qpb  += qpb_pattern[ch].yaxis[idx]
        data_qpb_ratio += [data/qpb]

    sp_scale_list = []
    xray_scale_list = []
    for ch in range(1,len(energy_range)):
        xray_scale_list += [xray_qpb_ratio[ch-1]]
        data_cnt = data_pattern[ch].integral()
        qpb_cnt = qpb_pattern[ch].integral()
        xray_cnt = xray_qpb_ratio[ch-1]*xray_pattern[ch].integral()
        raw_spf_cnt = spf_pattern[ch].integral()
        sp_scale = max((data_cnt-qpb_cnt-xray_cnt)/raw_spf_cnt,0.)
        sp_scale_list += [sp_scale]

    min_data_qpb_diff = 1e10
    min_data_qpb_diff_ch = 0
    for ch in range(1,len(energy_range)):
        data_qpb_diff = abs(spfree_data_qpb_ratio[ch-1]-1.)
        if min_data_qpb_diff>data_qpb_diff:
            min_data_qpb_diff = data_qpb_diff
            min_data_qpb_diff_ch = ch

    qpb_cnt_list = []
    spf_cnt_list = []
    xray_cnt_list = []
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Initial estimate')
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for ch in range(1,len(energy_range)):
        data_cnt = data_pattern[ch].integral()
        qpb_cnt = qpb_pattern[ch].integral()
        spf_cnt = sp_scale_list[ch-1]*spf_pattern[ch].integral()
        xray_cnt = xray_scale_list[ch-1]*xray_pattern[ch].integral()
        model_cnt = qpb_cnt + spf_cnt + xray_cnt
        qpb_cnt_list += [qpb_cnt]
        spf_cnt_list += [spf_cnt]
        xray_cnt_list += [xray_cnt]
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print ('ch = %s'%(ch))
        print ('sp_scale_list[ch-1] = %s'%(sp_scale_list[ch-1]))
        print ('data_cnt = %s'%(data_cnt))
        print ('model_cnt = %s'%(model_cnt))
        print ('xray_cnt = %s'%(xray_cnt))
        print ('spf_cnt = %s'%(spf_cnt))
        print ('qpb_cnt = %s'%(qpb_cnt))

    #return xray_cnt_list, spf_cnt_list, qpb_cnt_list

    # Define the function to minimize (residuals)
    def residuals(A):
        chi2 = 0.
        for ch in range(1,len(energy_range)):
            channel_weight = 1./abs(spfree_data_qpb_ratio[ch-1]-1.)
            data_cnt_sum = 0.
            qpb_model_sum = 0.
            spf_model_sum = 0.
            xray_raw_sum = 0.
            #spf_scale = A[0]
            spf_scale = A[0]*np.exp(energy_range[ch]/1000.*(-1.*A[1]))
            penalty_spf = (spf_scale-sp_scale_list[ch-1])*spf_pattern[ch].integral()*channel_weight
            for pattern in range(0,5):
                data_cnt = data_pattern[ch].yaxis[pattern]
                qpb_model = qpb_pattern[ch].yaxis[pattern]
                spf_model = spf_scale*spf_pattern[ch].yaxis[pattern]
                xray_raw = xray_pattern[ch].yaxis[pattern]
                data_cnt_sum += data_cnt
                qpb_model_sum += qpb_model
                spf_model_sum += spf_model
                xray_raw_sum += xray_raw
            xray_scale_sum = (data_cnt_sum-qpb_model_sum-spf_model_sum)/xray_raw_sum
            data_cnt_sum = 0.
            qpb_model_sum = 0.
            spf_model_sum = 0.
            xray_model_sum = 0.
            for pattern in range(0,5):
                data_cnt = data_pattern[ch].yaxis[pattern]
                qpb_model = qpb_pattern[ch].yaxis[pattern]
                spf_model = spf_scale*spf_pattern[ch].yaxis[pattern]
                xray_model = xray_scale_sum*xray_pattern[ch].yaxis[pattern]
                data_cnt_sum += data_cnt
                qpb_model_sum += qpb_model
                spf_model_sum += spf_model
                xray_model_sum += xray_model
                if data_cnt==0: continue
                pattern_weight = data_pattern[ch].integral()/data_cnt
                weight = pattern_weight
                chi = (data_cnt - xray_model - spf_model - qpb_model)*weight
                chi2 += chi*chi
            penalty_xray = 0.
            if xray_model_sum<0.:
                penalty_xray = 100.*abs(xray_model_sum)
            penalty_spf_neg = 0.
            if spf_model_sum<0.:
                penalty_spf_neg = 100.*abs(spf_model_sum)
            chi2 += penalty_spf*penalty_spf + penalty_xray*penalty_xray + penalty_spf_neg*penalty_spf_neg
        return chi2

    print ('min_data_qpb_diff_ch = %s'%(min_data_qpb_diff_ch))
    print ('param_init = %s'%(sp_scale_list[min_data_qpb_diff_ch-1]))
    qpb_scale = 1.
    param_init = [sp_scale_list[min_data_qpb_diff_ch-1],0.]
    param_bound_upper = [10.,2]
    param_bound_lower = [0.,-2]
    param_init = np.array(param_init)
    param_bound_upper = np.array(param_bound_upper)
    param_bound_lower = np.array(param_bound_lower)
    #result = least_squares(residuals, x0=param_init, bounds=(param_bound_lower,param_bound_upper))  # x0 is the initial guess for A
    result = minimize(residuals, x0=param_init, method='Nelder-Mead')  # x0 is the initial guess for A

    qpb_scale_norm = 1.

    qpb_cnt_list = []
    spf_cnt_list = []
    xray_cnt_list = []
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Optimized estimate')
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for ch in range(1,len(energy_range)):
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        sp_scale_norm = result.x[0]*np.exp(energy_range[ch]/1000.*(-1.*result.x[1]))
        print ('ch = %s'%(ch))
        print ('sp_scale_list[ch-1] = %s'%(sp_scale_list[ch-1]))
        print ('sp_scale_fit = %s'%(sp_scale_norm))
        spfree_data_qpb_diff = abs(spfree_data_qpb_ratio[ch-1]-1.)
        data_qpb_diff = abs(data_qpb_ratio[ch-1]-1.)
        print ('spfree_data_qpb_diff = %s'%(spfree_data_qpb_diff))
        print ('data_qpb_diff = %s'%(data_qpb_diff))
        #if spfree_data_qpb_diff<1.0 and data_qpb_diff>1.0:
        #    print ('use initial estimate')
        #    sp_scale_norm = sp_scale_list[ch-1] 
        data_cnt_sum = 0.
        qpb_model_sum = 0.
        spf_model_sum = 0.
        xray_raw_sum = 0.
        for pattern in range(0,5):
            data_cnt = data_pattern[ch].yaxis[pattern]
            data_cnt_sum += data_cnt
            qpb_model = qpb_pattern[ch].yaxis[pattern]
            qpb_model_sum += qpb_model
            spf_model = sp_scale_norm*spf_pattern[ch].yaxis[pattern]
            spf_model_sum += spf_model
            xray_raw = xray_pattern[ch].yaxis[pattern]
            xray_raw_sum += xray_raw
        xray_scale_norm = (data_cnt_sum-qpb_model_sum-spf_model_sum)/xray_raw_sum
        data_cnt = data_pattern[ch].integral()
        qpb_cnt = qpb_scale_norm*qpb_pattern[ch].integral()
        spf_cnt = sp_scale_norm*spf_pattern[ch].integral()
        xray_cnt = xray_scale_norm*xray_pattern[ch].integral()
        model_cnt = qpb_cnt + spf_cnt + xray_cnt
        qpb_cnt_list += [qpb_cnt]
        spf_cnt_list += [spf_cnt]
        xray_cnt_list += [xray_cnt]
        print ('data_cnt = %s'%(data_cnt))
        print ('model_cnt = %s'%(model_cnt))
        print ('xray_cnt = %s'%(xray_cnt))
        print ('spf_cnt = %s'%(spf_cnt))
        print ('qpb_cnt = %s'%(qpb_cnt))

    return xray_cnt_list, spf_cnt_list, qpb_cnt_list

    

    
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

xray_sample = 'extragalactic'
xray_obsID = 'ID0690900101'
#xray_obsID = 'ID0827251001'
#xray_sample = 'Cas_A'
#xray_obsID = 'ID0412180101'

plot_tag = detector+'_'+on_filter+'_'+on_obsID

xray_sci_fov_evt_filename = '../%s/%s/analysis/%s-fov-evt.fits'%(xray_sample,xray_obsID,detector)
xray_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-fov-evt.fits'%(xray_sample,xray_obsID,detector)
xray_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(xray_sample,xray_obsID,detector)
xray_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(xray_sample,xray_obsID,detector)
xray_rmf_filename = '../%s/%s/analysis/%s-src-rmf.fits'%(xray_sample,xray_obsID,detector)

off_sci_fov_evt_filename = '../%s/%s/analysis/%s-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
off_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
off_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(off_sample,off_obsID,detector)
off_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(off_sample,off_obsID,detector)
off_rmf_filename = '../%s/%s/analysis/%s-src-rmf.fits'%(off_sample,off_obsID,detector)

on_sci_fov_evt_filename = '../%s/%s/analysis/%s-%s-evt.fits'%(on_sample,on_obsID,detector,on_filter)
on_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-%s-evt.fits'%(on_sample,on_obsID,detector,on_filter)
on_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(on_sample,on_obsID,detector)
on_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(on_sample,on_obsID,detector)
on_rmf_filename = '../%s/%s/analysis/%s-src-rmf.fits'%(on_sample,on_obsID,detector)


def analyze_a_ccd_chip(energy_range=[200,12000],ccd_id=0):

    #print ('prepare xray sample time cuts')
    #
    #output_all_fov = read_event_file(xray_sci_fov_evt_filename,xray_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    #output_all_cor = read_event_file(xray_sci_cor_evt_filename,xray_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    #output_fwc_fov = read_event_file(xray_fwc_fov_evt_filename,xray_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    #output_fwc_cor = read_event_file(xray_fwc_cor_evt_filename,xray_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')

    #xray_duration_all_fov = output_all_fov[0]
    #xray_evt_count_all_fov = output_all_fov[1]
    #xray_lightcurve_all_fov = output_all_fov[2]
    #xray_pattern_all_fov = output_all_fov[3]
    #xray_spectrum_all_fov = output_all_fov[4]
    #xray_detx_all_fov = output_all_fov[5]
    #xray_image_all_fov = output_all_fov[6]
    #
    #xray_duration_all_cor = output_all_cor[0]
    #xray_evt_count_all_cor = output_all_cor[1]
    #xray_lightcurve_all_cor = output_all_cor[2]
    #xray_pattern_all_cor = output_all_cor[3]
    #xray_spectrum_all_cor = output_all_cor[4]
    #xray_detx_all_cor = output_all_cor[5]
    #xray_image_all_cor = output_all_cor[6]

    #xray_duration_fwc_fov = output_fwc_fov[0]
    #xray_evt_count_fwc_fov = output_fwc_fov[1]
    #xray_lightcurve_fwc_fov = output_fwc_fov[2]
    #xray_pattern_fwc_fov = output_fwc_fov[3]
    #xray_spectrum_fwc_fov = output_fwc_fov[4]
    #xray_detx_fwc_fov = output_fwc_fov[5]
    #xray_image_fwc_fov = output_fwc_fov[6]
    #
    #xray_duration_fwc_cor = output_fwc_cor[0]
    #xray_evt_count_fwc_cor = output_fwc_cor[1]
    #xray_lightcurve_fwc_cor = output_fwc_cor[2]
    #xray_pattern_fwc_cor = output_fwc_cor[3]
    #xray_spectrum_fwc_cor = output_fwc_cor[4]
    #xray_detx_fwc_cor = output_fwc_cor[5]
    #xray_image_fwc_cor = output_fwc_cor[6]

    #area_pix_frac_fwc_fov = xray_image_fwc_fov[0].integral()
    #area_pix_frac_fwc_cor = xray_image_fwc_cor[0].integral()
    #
    #time_pix_frac_all_fov = xray_lightcurve_all_fov.get_pixel_fraction()
    #time_expo_all_fov = xray_duration_all_fov*time_pix_frac_all_fov
    #
    #time_pix_frac_all_cor = xray_lightcurve_all_cor.get_pixel_fraction()
    #time_expo_all_cor = xray_duration_all_cor*time_pix_frac_all_cor
    #
    #xray_lightcurve_all_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)

    #mask_lc = make_timecut_mask(xray_lightcurve_all_fov,xray_lightcurve_all_cor) 

    #time_pix_frac_mask = mask_lc.get_pixel_fraction()
    #if time_pix_frac_mask>0.8:
    #    mask_lc = None

    #print ('apply x-ray sample space and time masks')
    #
    #output_sci_source = read_event_file(xray_sci_fov_evt_filename,xray_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto source',energy_range=energy_range,ccd_id=0)
    #output_sci_ring   = read_event_file(xray_sci_fov_evt_filename,xray_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto ring',energy_range=energy_range,ccd_id=0)
    #
    #xray_duration_sci_source = output_sci_source[0]
    #xray_evt_count_sci_source = output_sci_source[1]
    #xray_lightcurve_sci_source = output_sci_source[2]
    #xray_pattern_sci_source = output_sci_source[3]
    #xray_spectrum_sci_source = output_sci_source[4]
    #xray_detx_sci_source = output_sci_source[5]
    #xray_image_sci_source = output_sci_source[6]
    #
    #xray_duration_sci_ring = output_sci_ring[0]
    #xray_evt_count_sci_ring = output_sci_ring[1]
    #xray_lightcurve_sci_ring = output_sci_ring[2]
    #xray_pattern_sci_ring = output_sci_ring[3]
    #xray_spectrum_sci_ring = output_sci_ring[4]
    #xray_detx_sci_ring = output_sci_ring[5]
    #xray_image_sci_ring = output_sci_ring[6]

    ##area_pix_frac_sci_source = xray_image_sci_source[0].get_pixel_fraction()
    ##area_pix_frac_sci_ring = xray_image_sci_ring[0].get_pixel_fraction()
    #area_pix_frac_sci_source = source_radius*source_radius
    #area_pix_frac_sci_ring = ring_outer_radius*ring_outer_radius - ring_inner_radius*ring_inner_radius
    #for ch in range(0,len(energy_range)):
    #    xray_pattern_sci_ring[ch].scale(area_pix_frac_sci_source/area_pix_frac_sci_ring)

    #xray_pattern_fov_template = []
    #for ch in range(0,len(energy_range)):
    #    xray_pattern_fov_template += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
    #for ch in range(0,len(energy_range)):
    #    xray_pattern_fov_template[ch].add(xray_pattern_sci_source[0])
    #    xray_pattern_fov_template[ch].add(xray_pattern_sci_ring[0],factor=-1.)

    print ('prepare off sample time cuts')
    
    output_all_fov = read_event_file(off_sci_fov_evt_filename,off_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,off_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_all_cor = read_event_file(off_sci_cor_evt_filename,off_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,off_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    
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
    
    area_pix_frac_fwc_fov = off_image_fwc_fov[0].integral()
    area_pix_frac_fwc_cor = off_image_fwc_cor[0].integral()
    
    time_pix_frac_all_fov = off_lightcurve_all_fov.get_pixel_fraction()
    time_expo_all_fov = off_duration_all_fov*time_pix_frac_all_fov
    time_pix_frac_fwc_fov = off_lightcurve_fwc_fov.get_pixel_fraction()
    time_expo_fwc_fov = off_duration_fwc_fov*time_pix_frac_fwc_fov
    
    time_pix_frac_all_cor = off_lightcurve_all_cor.get_pixel_fraction()
    time_expo_all_cor = off_duration_all_cor*time_pix_frac_all_cor
    time_pix_frac_fwc_cor = off_lightcurve_fwc_cor.get_pixel_fraction()
    time_expo_fwc_cor = off_duration_fwc_cor*time_pix_frac_fwc_cor
    
    fwc_2_sci_ratio = off_image_all_cor[0].integral()/off_image_fwc_cor[0].integral()
    print ('Before timecut, OFF data fwc_2_sci_ratio = %s'%(fwc_2_sci_ratio))

    off_lightcurve_all_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    off_lightcurve_fwc_fov.scale(fwc_2_sci_ratio)

    mask_lc = make_timecut_mask(off_lightcurve_all_fov,off_lightcurve_all_cor) 

    print ('apply off sample space and time masks')
    
    output_sci_fov = read_event_file(off_sci_fov_evt_filename,off_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range)
    output_sci_cor = read_event_file(off_sci_cor_evt_filename,off_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range)
    
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
    
    output_spf_fov = read_event_file(off_sci_fov_evt_filename,off_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_range)
    output_spf_cor = read_event_file(off_sci_cor_evt_filename,off_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_range)
    
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
    
    output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,off_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range)
    output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,off_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range)
    
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

    area_pix_frac_fwc_fov = off_image_fwc_fov[0].integral()
    area_pix_frac_fwc_cor = off_image_fwc_cor[0].integral()
    
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
    
    fwc_2_sci_ratio = off_image_sci_cor[0].integral()/off_image_fwc_cor[0].integral()
    print ('After timecut, OFF data fwc_2_sci_ratio = %s'%(fwc_2_sci_ratio))

    off_lightcurve_sci_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    off_lightcurve_spf_fov.scale(time_expo_sci_fov/time_expo_spf_fov)
    off_lightcurve_spf_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_spf_cor)
    off_lightcurve_fwc_fov.scale(fwc_2_sci_ratio)
    off_lightcurve_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    
    for ch in range(0,len(energy_range)):

        #fwc_2_sci_ratio = off_image_sci_cor[ch].integral()/off_image_fwc_cor[ch].integral()
        #print ('After timecut, ch = %s, OFF data fwc_2_sci_ratio = %s'%(ch,fwc_2_sci_ratio))

        off_spectrum_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        off_spectrum_spf_fov[ch].scale(time_expo_sci_fov/time_expo_spf_fov)
        off_spectrum_spf_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_spf_cor)
        off_spectrum_fwc_fov[ch].scale(fwc_2_sci_ratio)
        off_spectrum_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    
        off_pattern_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        off_pattern_spf_fov[ch].scale(time_expo_sci_fov/time_expo_spf_fov)
        off_pattern_spf_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_spf_cor)
        off_pattern_fwc_fov[ch].scale(fwc_2_sci_ratio)
        off_pattern_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)

    if diagnostic_plots:
        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_lightcurve_all_fov.xaxis,off_lightcurve_all_fov.yaxis,yerr=off_lightcurve_all_fov.yerr,color='k',label='All FoV')
        axbig.errorbar(off_lightcurve_sci_fov.xaxis,off_lightcurve_sci_fov.yaxis,yerr=off_lightcurve_sci_fov.yerr,color='red',label='Sci FoV')
        axbig.errorbar(off_lightcurve_all_cor.xaxis,off_lightcurve_all_cor.yaxis,yerr=off_lightcurve_all_cor.yerr,color='green',label='All Cor')
        axbig.errorbar(off_lightcurve_fwc_cor.xaxis,off_lightcurve_fwc_cor.yaxis,yerr=off_lightcurve_fwc_cor.yerr,color='blue',label='FWC Cor')
        axbig.set_yscale('log')
        axbig.set_xlabel('Time')
        axbig.legend(loc='best')
        fig.savefig("%s/lightcurve_off_timecut_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_spectrum_sci_fov[0].xaxis,off_spectrum_sci_fov[0].yaxis,yerr=off_spectrum_sci_fov[0].yerr,label='Sci FoV')
        axbig.errorbar(off_spectrum_spf_fov[0].xaxis,off_spectrum_spf_fov[0].yaxis,yerr=off_spectrum_spf_fov[0].yerr,label='SPF FoV')
        axbig.errorbar(off_spectrum_fwc_fov[0].xaxis,off_spectrum_fwc_fov[0].yaxis,yerr=off_spectrum_fwc_fov[0].yerr,label='FWC FoV')
        axbig.set_xlabel('Energy [eV]')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("%s/off_spectrum_fov_raw_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_spectrum_sci_cor[0].xaxis,off_spectrum_sci_cor[0].yaxis,yerr=off_spectrum_sci_cor[0].yerr,label='Sci Cor')
        axbig.errorbar(off_spectrum_spf_cor[0].xaxis,off_spectrum_spf_cor[0].yaxis,yerr=off_spectrum_spf_cor[0].yerr,label='SPF Cor')
        axbig.errorbar(off_spectrum_fwc_cor[0].xaxis,off_spectrum_fwc_cor[0].yaxis,yerr=off_spectrum_fwc_cor[0].yerr,label='FWC Cor')
        axbig.set_xlabel('Energy [eV]')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("%s/off_spectrum_cor_raw_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_pattern_sci_fov[0].xaxis,off_pattern_sci_fov[0].yaxis,yerr=off_pattern_sci_fov[0].yerr,label='Sci FoV')
        axbig.errorbar(off_pattern_spf_fov[0].xaxis,off_pattern_spf_fov[0].yaxis,yerr=off_pattern_spf_fov[0].yerr,label='SPF FoV')
        axbig.errorbar(off_pattern_fwc_fov[0].xaxis,off_pattern_fwc_fov[0].yaxis,yerr=off_pattern_fwc_fov[0].yerr,label='FWC FoV')
        axbig.set_xlabel('Pattern')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("%s/off_pattern_fov_raw_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

    sp_pattern_fov_template = []
    sp_spectrum_fov_template = []
    for ch in range(0,len(energy_range)):
        sp_pattern_fov_template += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
        sp_spectrum_fov_template += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
    
    for ch in range(0,len(energy_range)):
        sp_spectrum_fov_template[ch].add(off_spectrum_spf_fov[ch])
        sp_spectrum_fov_template[ch].add(off_spectrum_sci_fov[ch],factor=-1.)
        ch_norm = sp_spectrum_fov_template[ch].integral()/off_pattern_spf_fov[0].integral()
        sp_pattern_fov_template[ch].add(off_pattern_spf_fov[0],factor=ch_norm)
        sp_pattern_fov_template[ch].add(off_pattern_sci_fov[0],factor=-1.*ch_norm)

    print ('prepare on sample time cuts')
    
    output_all_fov = read_event_file(on_sci_fov_evt_filename,on_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,on_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_all_cor = read_event_file(on_sci_cor_evt_filename,on_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,on_rmf_filename,mask_lc=None,mask_map=None,evt_filter='')
    
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

    area_pix_frac_fwc_fov = on_image_fwc_fov[0].integral()
    area_pix_frac_fwc_cor = on_image_fwc_cor[0].integral()
    
    time_pix_frac_all_fov = on_lightcurve_all_fov.get_pixel_fraction()
    time_expo_all_fov = on_duration_all_fov*time_pix_frac_all_fov
    time_pix_frac_fwc_fov = on_lightcurve_fwc_fov.get_pixel_fraction()
    time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
    
    time_pix_frac_all_cor = on_lightcurve_all_cor.get_pixel_fraction()
    time_expo_all_cor = on_duration_all_cor*time_pix_frac_all_cor
    time_pix_frac_fwc_cor = on_lightcurve_fwc_cor.get_pixel_fraction()
    time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor
    
    fwc_2_sci_ratio = on_image_all_cor[0].integral()/on_image_fwc_cor[0].integral()
    print ('Before timecut, ON data fwc_2_sci_ratio = %s'%(fwc_2_sci_ratio))

    on_lightcurve_all_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    on_lightcurve_fwc_fov.scale(fwc_2_sci_ratio)

    mask_lc = make_timecut_mask(on_lightcurve_all_fov,on_lightcurve_all_cor) 
    
    time_pix_frac_mask = mask_lc.get_pixel_fraction()
    if time_pix_frac_mask>0.8:
        mask_lc = None
    if not SP_flare_mask:
        mask_lc = None

    print ('apply space and time masks')
    
    output_sci_fov = read_event_file(on_sci_fov_evt_filename,on_rmf_filename,mask_lc=mask_lc,mask_map=None,write_events=True,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    output_sci_cor = read_event_file(on_sci_cor_evt_filename,on_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    
    on_duration_sci_fov = output_sci_fov[0]
    on_evt_count_sci_fov = output_sci_fov[1]
    on_lightcurve_sci_fov = output_sci_fov[2]
    on_pattern_sci_fov = output_sci_fov[3]
    on_spectrum_sci_fov = output_sci_fov[4]
    on_detx_sci_fov = output_sci_fov[5]
    on_image_sci_fov = output_sci_fov[6]
    #on_image_sky_sci_fov = output_sci_fov[7]
    
    on_duration_sci_cor = output_sci_cor[0]
    on_evt_count_sci_cor = output_sci_cor[1]
    on_lightcurve_sci_cor = output_sci_cor[2]
    on_pattern_sci_cor = output_sci_cor[3]
    on_spectrum_sci_cor = output_sci_cor[4]
    on_detx_sci_cor = output_sci_cor[5]
    on_image_sci_cor = output_sci_cor[6]
    
    output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,on_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,on_rmf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    
    on_duration_fwc_fov = output_fwc_fov[0]
    on_evt_count_fwc_fov = output_fwc_fov[1]
    on_lightcurve_fwc_fov = output_fwc_fov[2]
    on_pattern_fwc_fov = output_fwc_fov[3]
    on_spectrum_fwc_fov = output_fwc_fov[4]
    on_detx_fwc_fov = output_fwc_fov[5]
    on_image_fwc_fov = output_fwc_fov[6]
    #on_image_sky_fwc_fov = output_fwc_fov[7]
    
    on_duration_fwc_cor = output_fwc_cor[0]
    on_evt_count_fwc_cor = output_fwc_cor[1]
    on_lightcurve_fwc_cor = output_fwc_cor[2]
    on_pattern_fwc_cor = output_fwc_cor[3]
    on_spectrum_fwc_cor = output_fwc_cor[4]
    on_detx_fwc_cor = output_fwc_cor[5]
    on_image_fwc_cor = output_fwc_cor[6]
    
    area_pix_frac_fwc_fov = on_image_fwc_fov[0].integral()
    area_pix_frac_fwc_cor = on_image_fwc_cor[0].integral()
    
    time_pix_frac_sci_fov = on_lightcurve_sci_fov.get_pixel_fraction()
    time_expo_sci_fov = on_duration_sci_fov*time_pix_frac_sci_fov
    time_pix_frac_sci_cor = on_lightcurve_sci_cor.get_pixel_fraction()
    time_expo_sci_cor = on_duration_sci_cor*time_pix_frac_sci_cor
    
    time_pix_frac_fwc_fov = on_lightcurve_fwc_fov.get_pixel_fraction()
    time_expo_fwc_fov = on_duration_fwc_fov*time_pix_frac_fwc_fov
    time_pix_frac_fwc_cor = on_lightcurve_fwc_cor.get_pixel_fraction()
    time_expo_fwc_cor = on_duration_fwc_cor*time_pix_frac_fwc_cor

    fwc_2_sci_ratio = on_image_sci_cor[0].integral()/on_image_fwc_cor[0].integral()
    print ('After timecut, ON data fwc_2_sci_ratio = %s'%(fwc_2_sci_ratio))
    
    on_lightcurve_sci_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    on_lightcurve_fwc_fov.scale(fwc_2_sci_ratio)
    on_lightcurve_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    
    for ch in range(0,len(energy_range)):

        #fwc_2_sci_ratio = on_image_sci_cor[ch].integral()/on_image_fwc_cor[ch].integral()
        #print ('After timecut, ch = %s, ON data fwc_2_sci_ratio = %s'%(ch,fwc_2_sci_ratio))

        on_spectrum_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_spectrum_fwc_fov[ch].scale(fwc_2_sci_ratio)
        on_spectrum_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_pattern_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_pattern_fwc_fov[ch].scale(fwc_2_sci_ratio)
        on_pattern_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_detx_fwc_fov[ch].scale(fwc_2_sci_ratio)
        on_image_fwc_fov[ch].scale(fwc_2_sci_ratio)

    if diagnostic_plots:

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_lightcurve_fwc_fov.xaxis,on_lightcurve_fwc_fov.yaxis,yerr=on_lightcurve_fwc_fov.yerr,color='k',label='FWC FoV')
        axbig.errorbar(on_lightcurve_fwc_cor.xaxis,on_lightcurve_fwc_cor.yaxis,yerr=on_lightcurve_fwc_cor.yerr,color='red',label='FWC Cor')
        axbig.set_yscale('log')
        axbig.set_xlabel('Time')
        axbig.legend(loc='best')
        fig.savefig("%s/lightcurve_on_fwc_timecut_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_lightcurve_sci_fov.xaxis,on_lightcurve_sci_fov.yaxis,yerr=on_lightcurve_sci_fov.yerr,color='k',label='Sci FoV')
        axbig.errorbar(on_lightcurve_sci_cor.xaxis,on_lightcurve_sci_cor.yaxis,yerr=on_lightcurve_sci_cor.yerr,color='red',label='Sci Cor')
        axbig.set_yscale('log')
        axbig.set_xlabel('Time')
        axbig.legend(loc='best')
        fig.savefig("%s/lightcurve_on_sci_timecut_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_lightcurve_all_fov.xaxis,on_lightcurve_all_fov.yaxis,yerr=on_lightcurve_all_fov.yerr,color='k',label='All FoV')
        axbig.errorbar(on_lightcurve_sci_fov.xaxis,on_lightcurve_sci_fov.yaxis,yerr=on_lightcurve_sci_fov.yerr,color='red',label='Sci FoV')
        axbig.errorbar(on_lightcurve_all_cor.xaxis,on_lightcurve_all_cor.yaxis,yerr=on_lightcurve_all_cor.yerr,color='green',label='All Cor')
        axbig.set_yscale('log')
        axbig.set_xlabel('Time')
        axbig.legend(loc='best')
        fig.savefig("%s/lightcurve_on_all_timecut_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_spectrum_sci_fov[0].xaxis,on_spectrum_sci_fov[0].yaxis,yerr=on_spectrum_sci_fov[0].yerr,label='Sci FoV')
        axbig.errorbar(on_spectrum_fwc_fov[0].xaxis,on_spectrum_fwc_fov[0].yaxis,yerr=on_spectrum_fwc_fov[0].yerr,label='FWC FoV')
        axbig.set_xlabel('Energy [eV]')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("%s/on_spectrum_fov_raw_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_spectrum_sci_cor[0].xaxis,on_spectrum_sci_cor[0].yaxis,yerr=on_spectrum_sci_cor[0].yerr,label='Sci Cor')
        axbig.errorbar(on_spectrum_fwc_cor[0].xaxis,on_spectrum_fwc_cor[0].yaxis,yerr=on_spectrum_fwc_cor[0].yerr,label='FWC Cor')
        axbig.set_xlabel('Energy [eV]')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("%s/on_spectrum_cor_raw_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_pattern_sci_fov[0].xaxis,on_pattern_sci_fov[0].yaxis,yerr=on_pattern_sci_fov[0].yerr,label='Sci FoV')
        axbig.errorbar(on_pattern_fwc_fov[0].xaxis,on_pattern_fwc_fov[0].yaxis,yerr=on_pattern_fwc_fov[0].yerr,label='FWC FoV')
        axbig.set_xlabel('Pattern')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("%s/on_pattern_fov_raw_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(on_pattern_sci_cor[0].xaxis,on_pattern_sci_cor[0].yaxis,yerr=on_pattern_sci_cor[0].yerr,label='Sci Cor')
        axbig.errorbar(on_pattern_fwc_cor[0].xaxis,on_pattern_fwc_cor[0].yaxis,yerr=on_pattern_fwc_cor[0].yerr,label='FWC Cor')
        axbig.set_xlabel('Pattern')
        axbig.legend(loc='best')
        axbig.set_yscale('log')
        fig.savefig("%s/on_pattern_cor_raw_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

    print ('predict SP background')
    
    qpb_pattern_fov_template = []
    qpb_spectrum_fov_template = []
    qpb_detx_fov_template = []
    qpb_image_fov_template = []
    for ch in range(0,len(energy_range)):
        qpb_pattern_fov_template += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
        qpb_spectrum_fov_template += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
        qpb_detx_fov_template += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
        qpb_image_fov_template += [MyArray2D()]
    for ch in range(0,len(energy_range)):
        qpb_spectrum_fov_template[ch].add(on_spectrum_fwc_fov[ch])
        qpb_detx_fov_template[ch].add(on_detx_fwc_fov[ch])
        qpb_image_fov_template[ch].add(on_image_fwc_fov[ch])
        ch_norm = on_spectrum_fwc_fov[ch].integral()/on_pattern_fwc_fov[ch].integral()
        qpb_pattern_fov_template[ch].add(on_pattern_fwc_fov[ch],factor=ch_norm)

    for ch in range(0,len(energy_range)):
        sp_rescale = qpb_pattern_fov_template[0].integral()/sp_pattern_fov_template[0].integral()
        sp_pattern_fov_template[ch].scale(sp_rescale)

    xray_pattern_fov_template = []
    for ch in range(0,len(energy_range)):
        xray_pattern_fov_template += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
    for ch in range(0,len(energy_range)):
        xray_pattern_fov_template[ch].add(qpb_pattern_fov_template[ch])

    xray_cnt_list, spf_cnt_list, qpb_cnt_list = fit_pattern(energy_range,on_pattern_sci_fov,xray_pattern_fov_template,sp_pattern_fov_template,qpb_pattern_fov_template)
        
    sp_detx_fov_template = []
    sp_image_fov_template = []
    for ch in range(0,len(energy_range)):
        sp_detx_fov_template += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
        sp_image_fov_template += [MyArray2D()]
    for ch in range(0,len(energy_range)):
        sp_pattern_norm = sp_pattern_fov_template[ch].integral()
        qpb_pattern_norm = qpb_pattern_fov_template[ch].integral()
        sp_detx_fov_template[ch].add(on_detx_fwc_fov[ch],sp_pattern_norm/qpb_pattern_norm)
        sp_image_fov_template[ch].add(on_image_fwc_fov[ch],sp_pattern_norm/qpb_pattern_norm)

    for ch in range(1,len(energy_range)):
        xray_scale = xray_cnt_list[ch-1]/xray_pattern_fov_template[ch].integral()
        xray_pattern_fov_template[ch].scale(xray_scale)
        qpb_scale = qpb_cnt_list[ch-1]/qpb_pattern_fov_template[ch].integral()
        qpb_pattern_fov_template[ch].scale(qpb_scale)
        sp_scale = spf_cnt_list[ch-1]/sp_pattern_fov_template[ch].integral()
        sp_pattern_fov_template[ch].scale(sp_scale)

        qpb_scale = qpb_cnt_list[ch-1]/qpb_detx_fov_template[ch].integral()
        qpb_detx_fov_template[ch].scale(qpb_scale)
        qpb_scale = qpb_cnt_list[ch-1]/qpb_image_fov_template[ch].integral()
        qpb_image_fov_template[ch].scale(qpb_scale)
        qpb_scale = qpb_cnt_list[ch-1]/qpb_spectrum_fov_template[ch].integral()
        qpb_spectrum_fov_template[ch].scale(qpb_scale)

        sp_scale = spf_cnt_list[ch-1]/sp_spectrum_fov_template[ch].integral()
        sp_spectrum_fov_template[ch].scale(sp_scale)
        sp_scale = spf_cnt_list[ch-1]/sp_detx_fov_template[ch].integral()
        sp_detx_fov_template[ch].scale(sp_scale)
        sp_scale = spf_cnt_list[ch-1]/sp_image_fov_template[ch].integral()
        sp_image_fov_template[ch].scale(sp_scale)


    on_spectrum_sci_fov_bkg_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    sp_spectrum_fov_template_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    qpb_spectrum_fov_template_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    for ch in range(1,len(energy_range)):
        sp_spectrum_fov_template_sum.add(sp_spectrum_fov_template[ch])
        qpb_spectrum_fov_template_sum.add(qpb_spectrum_fov_template[ch])
        on_spectrum_sci_fov_bkg_sum.add(sp_spectrum_fov_template[ch])
        on_spectrum_sci_fov_bkg_sum.add(qpb_spectrum_fov_template[ch])

    on_detx_sci_fov_bkg_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
    sp_detx_fov_template_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
    qpb_detx_fov_template_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
    for ch in range(1,len(energy_range)):
        sp_detx_fov_template_sum.add(sp_detx_fov_template[ch])
        qpb_detx_fov_template_sum.add(qpb_detx_fov_template[ch])
        on_detx_sci_fov_bkg_sum.add(sp_detx_fov_template[ch])
        on_detx_sci_fov_bkg_sum.add(qpb_detx_fov_template[ch])
    
    on_image_sci_fov_bkg_sum = MyArray2D()
    sp_image_fov_template_sum = MyArray2D()
    qpb_image_fov_template_sum = MyArray2D()
    for ch in range(1,len(energy_range)):
        sp_image_fov_template_sum.add(sp_image_fov_template[ch])
        qpb_image_fov_template_sum.add(qpb_image_fov_template[ch])
        on_image_sci_fov_bkg_sum.add(sp_image_fov_template[ch])
        on_image_sci_fov_bkg_sum.add(qpb_image_fov_template[ch])
    on_image_sci_fov_bkg_sum.flattening()
    sp_image_fov_template_sum.flattening()
    qpb_image_fov_template_sum.flattening()


    if diagnostic_plots:
        for ch in range(1,len(energy_range)):
            fig.clf()
            axbig = fig.add_subplot()
            axbig.errorbar(on_pattern_sci_fov[ch].xaxis,on_pattern_sci_fov[ch].yaxis,yerr=on_pattern_sci_fov[ch].yerr,color='k',label='Data')
            axbig.errorbar(xray_pattern_fov_template[ch].xaxis,xray_pattern_fov_template[ch].yaxis,yerr=xray_pattern_fov_template[ch].yerr,color='r',label='X-ray')
            axbig.errorbar(qpb_pattern_fov_template[ch].xaxis,qpb_pattern_fov_template[ch].yaxis,yerr=qpb_pattern_fov_template[ch].yerr,color='b',label='QPB')
            axbig.bar(sp_pattern_fov_template[ch].xaxis,sp_pattern_fov_template[ch].yaxis,color='orange',label='SP Flare')
            axbig.set_yscale('log')
            axbig.set_ylim(bottom=1)
            axbig.set_xlabel('PATTERN')
            axbig.legend(loc='best')
            fig.savefig("%s/pattern_on_fov_scaled_ch%s_ccd%s_%s.png"%(output_dir,ch,ccd_id,plot_tag),bbox_inches='tight')
            axbig.remove()

        for ch in range(1,len(energy_range)):
            bkg_spectrum_fov = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
            bkg_spectrum_fov.add(qpb_spectrum_fov_template[ch])
            bkg_spectrum_fov.add(sp_spectrum_fov_template[ch])
            fig.clf()
            axbig = fig.add_subplot()
            axbig.errorbar(on_spectrum_sci_fov[ch].xaxis,on_spectrum_sci_fov[ch].yaxis,yerr=on_spectrum_sci_fov[ch].yerr,color='k',label='Data')
            axbig.errorbar(qpb_spectrum_fov_template[ch].xaxis,qpb_spectrum_fov_template[ch].yaxis,yerr=qpb_spectrum_fov_template[ch].yerr,color='b',label='QPB')
            axbig.errorbar(sp_spectrum_fov_template[ch].xaxis,sp_spectrum_fov_template[ch].yaxis,yerr=sp_spectrum_fov_template[ch].yerr,color='orange',label='SP Flare')
            axbig.errorbar(bkg_spectrum_fov.xaxis,bkg_spectrum_fov.yaxis,yerr=bkg_spectrum_fov.yerr,color='red',label='BKG')
            axbig.set_yscale('log')
            axbig.set_ylim(bottom=1)
            axbig.set_xlabel('Energy [eV]')
            axbig.legend(loc='best')
            fig.savefig("%s/spectrum_on_fov_scaled_ch%s_ccd%s_%s.png"%(output_dir,ch,ccd_id,plot_tag),bbox_inches='tight')
            axbig.remove()

        for ch in range(1,len(energy_range)):
            xray_pattern_fov_template[ch].normalize()
            qpb_pattern_fov_template[ch].normalize()
            sp_pattern_fov_template[ch].normalize()
            
            fig.clf()
            axbig = fig.add_subplot()
            axbig.plot(xray_pattern_fov_template[ch].xaxis,xray_pattern_fov_template[ch].yaxis,label='X-ray')
            axbig.plot(qpb_pattern_fov_template[ch].xaxis,qpb_pattern_fov_template[ch].yaxis,label='QPB')
            axbig.plot(sp_pattern_fov_template[ch].xaxis,sp_pattern_fov_template[ch].yaxis,label='SP')
            axbig.set_yscale('log')
            axbig.set_xlabel('Pattern')
            axbig.legend(loc='best')
            fig.savefig("%s/pattern_fov_template_ch%s_ccd%s_%s.png"%(output_dir,ch,ccd_id,plot_tag),bbox_inches='tight')
            axbig.remove()


    output_spectrum = [on_spectrum_sci_fov[0],on_spectrum_sci_fov_bkg_sum,sp_spectrum_fov_template_sum,qpb_spectrum_fov_template_sum]
    output_detx = [on_detx_sci_fov[0],on_detx_sci_fov_bkg_sum,sp_detx_fov_template_sum,qpb_detx_fov_template_sum]
    output_image = [on_image_sci_fov[0],on_image_sci_fov_bkg_sum,sp_image_fov_template_sum,qpb_image_fov_template_sum]

    return output_spectrum, output_detx, output_image


sci_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
bkg_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
qpb_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spf_spectrum_sum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)

on_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
bkg_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
qpb_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
sp_detx_sum = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)

#excess_image_sky_sum = MyArray2D(coord='radec')
#on_image_sky_sum = MyArray2D(coord='radec')
#bkg_image_sky_sum = MyArray2D(coord='radec')
#qpb_image_sky_sum = MyArray2D(coord='radec')
#sp_image_sky_sum = MyArray2D(coord='radec')

excess_image_sum = MyArray2D()
on_image_sum = MyArray2D()
bkg_image_sum = MyArray2D()
qpb_image_sum = MyArray2D()
sp_image_sum = MyArray2D()

for ccd in range(0,len(ana_ccd_bins)):
    ana_output_spectrum, ana_output_detx, ana_output_image = analyze_a_ccd_chip(energy_range=energy_array,ccd_id=ana_ccd_bins[ccd])

    prob_bkg_spectrum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    prob_bkg_spectrum.add(ana_output_spectrum[1])
    prob_bkg_spectrum.normalize2()

    evt_detx_list = []
    evt_dety_list = []
    evt_pi_list = []
    for idx_x in range(0,len(ana_output_image[1].xaxis)):
        for idx_y in range(0,len(ana_output_image[1].yaxis)):
            evt_detx = ana_output_image[1].xaxis[idx_x]
            evt_dety = ana_output_image[1].yaxis[idx_y]
            sum_bkg = ana_output_image[1].zaxis[idx_x,idx_y]
            evt_pi = np.random.choice(a=prob_bkg_spectrum.xaxis, size=int(sum_bkg*10), p=prob_bkg_spectrum.yaxis)
            for evt in range(0,len(evt_pi)):
                evt_detx_list += [evt_detx]
                evt_dety_list += [evt_dety]
                evt_pi_list += [evt_pi[evt]]
    col_detx = fits.Column(name='DETX', array=evt_detx_list, format='D')
    col_dety = fits.Column(name='DETY', array=evt_dety_list, format='D')
    col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
    my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_pi],name='EVENT')
    my_table.writeto('%s/bkg_events_%s_ccd%s.fits'%(output_dir,detector,ana_ccd_bins[ccd]), overwrite=True)

    prob_spf_spectrum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    prob_spf_spectrum.add(ana_output_spectrum[2])
    prob_spf_spectrum.normalize2()

    evt_detx_list = []
    evt_dety_list = []
    evt_pi_list = []
    for idx_x in range(0,len(ana_output_image[2].xaxis)):
        for idx_y in range(0,len(ana_output_image[2].yaxis)):
            evt_detx = ana_output_image[2].xaxis[idx_x]
            evt_dety = ana_output_image[2].yaxis[idx_y]
            sum_bkg = ana_output_image[2].zaxis[idx_x,idx_y]
            evt_pi = np.random.choice(a=prob_spf_spectrum.xaxis, size=int(sum_bkg*10), p=prob_spf_spectrum.yaxis)
            for evt in range(0,len(evt_pi)):
                evt_detx_list += [evt_detx]
                evt_dety_list += [evt_dety]
                evt_pi_list += [evt_pi[evt]]
    col_detx = fits.Column(name='DETX', array=evt_detx_list, format='D')
    col_dety = fits.Column(name='DETY', array=evt_dety_list, format='D')
    col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
    my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_pi],name='EVENT')
    my_table.writeto('%s/spf_events_%s_ccd%s.fits'%(output_dir,detector,ana_ccd_bins[ccd]), overwrite=True)

    prob_qpb_spectrum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    prob_qpb_spectrum.add(ana_output_spectrum[3])
    prob_qpb_spectrum.normalize2()

    evt_detx_list = []
    evt_dety_list = []
    evt_pi_list = []
    for idx_x in range(0,len(ana_output_image[3].xaxis)):
        for idx_y in range(0,len(ana_output_image[3].yaxis)):
            evt_detx = ana_output_image[3].xaxis[idx_x]
            evt_dety = ana_output_image[3].yaxis[idx_y]
            sum_bkg = ana_output_image[3].zaxis[idx_x,idx_y]
            evt_pi = np.random.choice(a=prob_qpb_spectrum.xaxis, size=int(sum_bkg*10), p=prob_qpb_spectrum.yaxis)
            for evt in range(0,len(evt_pi)):
                evt_detx_list += [evt_detx]
                evt_dety_list += [evt_dety]
                evt_pi_list += [evt_pi[evt]]
    col_detx = fits.Column(name='DETX', array=evt_detx_list, format='D')
    col_dety = fits.Column(name='DETY', array=evt_dety_list, format='D')
    col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
    my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_pi],name='EVENT')
    my_table.writeto('%s/qpb_events_%s_ccd%s.fits'%(output_dir,detector,ana_ccd_bins[ccd]), overwrite=True)

    sci_spectrum_sum.add(ana_output_spectrum[0])
    bkg_spectrum_sum.add(ana_output_spectrum[1])
    spf_spectrum_sum.add(ana_output_spectrum[2])
    qpb_spectrum_sum.add(ana_output_spectrum[3])
    on_detx_sum.add(ana_output_detx[0])
    bkg_detx_sum.add(ana_output_detx[1])
    sp_detx_sum.add(ana_output_detx[2])
    qpb_detx_sum.add(ana_output_detx[3])
    #on_image_sky_sum.add(ana_output_image_sky[0])
    #bkg_image_sky_sum.add(ana_output_image_sky[1])
    #sp_image_sky_sum.add(ana_output_image_sky[2])
    #qpb_image_sky_sum.add(ana_output_image_sky[3])
    on_image_sum.add(ana_output_image[0])
    bkg_image_sum.add(ana_output_image[1])
    sp_image_sum.add(ana_output_image[2])
    qpb_image_sum.add(ana_output_image[3])

excess_image_sum.add(on_image_sum)
excess_image_sum.add(bkg_image_sum,-1.)
#excess_image_sky_sum.add(on_image_sky_sum)
#excess_image_sky_sum.add(bkg_image_sky_sum,-1.)

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(sci_spectrum_sum.xaxis,sci_spectrum_sum.yaxis,yerr=sci_spectrum_sum.yerr,color='k',label='Data')
axbig.plot(bkg_spectrum_sum.xaxis,bkg_spectrum_sum.yaxis,color='red',label='Bkg')
axbig.plot(spf_spectrum_sum.xaxis,spf_spectrum_sum.yaxis,label='SP')
axbig.plot(qpb_spectrum_sum.xaxis,qpb_spectrum_sum.yaxis,label='QPB')
axbig.set_yscale('log')
axbig.set_ylim(bottom=1)
axbig.set_xlabel('Energy [eV]')
axbig.legend(loc='best')
fig.savefig("%s/spectrum_on_fit_sum_%s.png"%(output_dir,plot_tag),bbox_inches='tight')
axbig.remove()

#SaveSpectrum(sci_spectrum_sum,'sci_spectrum_sum')
#SaveSpectrum(bkg_spectrum_sum,'bkg_spectrum_sum')
#SaveSpectrum(qpb_spectrum_sum,'qpb_spectrum_sum')
#SaveSpectrum(spf_spectrum_sum,'spf_spectrum_sum')

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
fig.savefig("%s/detx_on_fit_sum_%s.png"%(output_dir,plot_tag),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = on_image_sum.xaxis.min()
xmax = on_image_sum.xaxis.max()
ymin = on_image_sum.yaxis.min()
ymax = on_image_sum.yaxis.max()
axbig.imshow(excess_image_sum.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/image_det_sum_%s.png"%(output_dir,plot_tag),bbox_inches='tight')
axbig.remove()

SaveImage(excess_image_sum,'image_det_sum')



























end_time = time.time()
elapsed_time = end_time - start_time
print('elapsed_time = %s'%(elapsed_time))

exit()


