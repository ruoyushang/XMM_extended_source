
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
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
find_nearest_ref_det_idx = common_functions.find_nearest_ref_det_idx
LoadCoordinateMatrix = common_functions.LoadCoordinateMatrix


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

#energy_array = [2000,4000,6000,8000,10000,12000]
energy_array = [2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000]
#ana_ccd_bins = [1,2,3,4,5,6,7]
#ana_ccd_bins = [1]
ana_ccd_bins = [0]

output_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/'+on_obsID

source_radius = 1000
ring_inner_radius = 1500
ring_outer_radius = 6000

sample_scale = 10.

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

def find_point_sources(image_data,image_mask):

    avg_cnt = 0.
    pix_used = 0.
    for idx_x in range(0,len(image_mask.xaxis)):
        for idx_y in range(0,len(image_mask.yaxis)):
            mask = image_mask.zaxis[idx_x,idx_y]
            if mask==1: continue
            avg_cnt += image_data.zaxis[idx_x,idx_y]
            pix_used += 1.
    avg_cnt = avg_cnt/pix_used

    rms_cnt = 0.
    pix_used = 0.
    for idx_x in range(0,len(image_mask.xaxis)):
        for idx_y in range(0,len(image_mask.yaxis)):
            mask = image_mask.zaxis[idx_x,idx_y]
            if mask==1: continue
            rms_cnt += pow(image_data.zaxis[idx_x,idx_y]-avg_cnt,2)
            pix_used += 1.
    rms_cnt = pow(rms_cnt/pix_used,0.5)

    for idx_x in range(0,len(image_mask.xaxis)):
        for idx_y in range(0,len(image_mask.yaxis)):
            mask = image_mask.zaxis[idx_x,idx_y]
            if mask==1: continue
            significance = (image_data.zaxis[idx_x,idx_y]-avg_cnt)/rms_cnt
            if significance>5.:
                image_mask.zaxis[idx_x,idx_y] = 1

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*u.deg, b*u.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def ConvertRaDecToGalactic(ra, dec):
    my_sky = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg

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

def ConvertDet2Sky(target_det,origin_sky,origin_det,mtx_conv_det_2_sky):

    target_det = np.array(target_det)
    origin_sky = np.array(origin_sky)
    origin_det = np.array(origin_det)

    delta_det = target_det-origin_det
    delta_sky = np.matmul(mtx_conv_det_2_sky,delta_det)
    target_sky = delta_sky+origin_sky

    return target_sky[0,0], target_sky[0,1]

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

ref_sky = []
ref_det = []
mtx_det_2_sky = []
mtx_sky_2_det = []
for idx_ra in range(0,1):
    for idx_dec in range(0,1):
        ref_sky_local, ref_det_local, mtx_det_2_sky_local, mtx_sky_2_det_local = LoadCoordinateMatrix(idx_ra,idx_dec,on_sample,on_obsID,detector)
        ref_sky += [ref_sky_local]
        ref_det += [ref_det_local]
        mtx_det_2_sky += [mtx_det_2_sky_local]
        mtx_sky_2_det += [mtx_sky_2_det_local]
if on_obsID=='ID0412180101' or on_obsID=='ID0400210101':
    ref_sky_local, ref_det_local, mtx_det_2_sky_local, mtx_sky_2_det_local = LoadCoordinateMatrix(97,97,on_sample,on_obsID,detector)
    ref_sky += [ref_sky_local]
    ref_det += [ref_det_local]
    mtx_det_2_sky += [mtx_det_2_sky_local]
    mtx_sky_2_det += [mtx_sky_2_det_local]


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


def read_event_file(filename,arf_name,mask_lc=None,mask_map=None,write_events=False,evt_filter='',energy_range=[200,12000],ccd_id=0):

    # how to read events:
    # https://docs.astropy.org/en/stable/generated/examples/io/fits-tables.html#accessing-data-stored-as-a-table-in-a-multi-extension-fits-mef-file


    time_frac = 1.
    if not mask_lc==None:
        time_frac = 0.
        total_time = 0.
        select_time = 0.
        for t in range(0,len(mask_lc.xaxis)):
            total_time += 1.
            zscore = mask_lc.yaxis[t]
            if 'sp-veto' in evt_filter:
                if zscore>0: 
                    continue
                else:
                    select_time += 1.
            if 'sp-select' in evt_filter:
                if zscore==0: 
                    continue
                else:
                    select_time += 1.
        time_frac = select_time/total_time

    print ('read file %s'%(filename))
    print (f'{evt_filter}, time_frac = {time_frac}')
    hdu_list = fits.open(filename)
    #print (filename)
    #print (hdu_list.info())

    events = Table.read(filename, hdu=1)
    #print('total events:')
    #print(len(events))
    #print('events.columns:')
    #print(events.columns)

    evt_count = 0.
    lightcurve_array_all = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
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
    obs_duration = (time_end-time_start)*time_frac

    evt_detx_list = []
    evt_dety_list = []
    evt_ra_list = []
    evt_dec_list = []
    #evt_l_list = []
    #evt_b_list = []
    evt_pi_list = []

    for evt in range(0,len(events)):

        if fast_test:
            if (evt%10)!=0: continue

        evt_time = events[evt]['TIME']
        evt_pattern = events[evt]['PATTERN']

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

        if evt_pattern>4: continue
        if 'sp-free' in evt_filter:
            if evt_pattern!=2: continue

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

        lightcurve_array_all.fill((evt_time-time_start)/(time_end-time_start))
        if not mask_lc==None:
            zscore = mask_lc.get_bin_content((evt_time-time_start)/(time_end-time_start))
            if 'sp-veto' in evt_filter:
                if zscore>0: continue
            if 'sp-select' in evt_filter:
                if zscore==0: continue

        #ref_det_idx_x, ref_det_idx_y = find_nearest_ref_det_idx([evt_detx,evt_dety],ref_det)
        ref_det_idx = 0
        if on_obsID=='ID0412180101' or on_obsID=='ID0400210101':
            ref_det_idx = 1
        evt_ra, evt_dec = ConvertDet2Sky([evt_detx,evt_dety],ref_sky[ref_det_idx],ref_det[ref_det_idx],mtx_det_2_sky[ref_det_idx])
        #evt_l, evt_b = ConvertRaDecToGalactic(evt_ra,evt_dec)

        evt_detx_list += [evt_detx]
        evt_dety_list += [evt_dety]
        evt_ra_list += [evt_ra]
        evt_dec_list += [evt_dec]
        #evt_l_list += [evt_l]
        #evt_b_list += [evt_b]
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
        col_ra = fits.Column(name='RA', array=evt_ra_list, format='D')
        col_dec = fits.Column(name='DEC', array=evt_dec_list, format='D')
        #col_l = fits.Column(name='GalL', array=evt_l_list, format='D')
        #col_b = fits.Column(name='GalB', array=evt_b_list, format='D')
        col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
        #my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_l,col_b,col_pi],name='EVENT')
        my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_pi],name='EVENT')
        my_table.header['REF_RA'] = ref_sky[0][0]
        my_table.header['REF_DEC'] = ref_sky[0][1]
        my_table.header['obs_duration'] = obs_duration
        my_table.writeto('%s/sci_events_%s_ccd%s.fits'%(output_dir,detector,ccd_id), overwrite=True)
        #fits.writeto('%s/sci_events_%s_ccd%s.fits'%(output_dir,detector,ccd_id), my_table, overwrite=True)

        energy_list = []
        area_list = []
        arf_in_table = Table.read(arf_name, hdu=1)
        for e in range(0,len(energy_array)-1):
            energy_lo = energy_array[e]
            energy_hi = energy_array[e+1]
            entry_cnt = 0.
            avg_area = 0.
            for entry in range(0,len(arf_in_table)):
                energy = 0.5*(arf_in_table[entry]['ENERG_LO']+arf_in_table[entry]['ENERG_HI'])*1000.
                if energy<energy_lo: continue
                if energy>=energy_hi: continue
                area = arf_in_table[entry]['SPECRESP']
                avg_area += area
                entry_cnt += 1.
            if entry_cnt>0:
                avg_area = avg_area/entry_cnt
            energy_list += [energy_lo]
            area_list += [avg_area]
        col_energy = fits.Column(name='ENERGY', array=energy_list, format='D')
        col_area = fits.Column(name='AREA', array=area_list, format='D')
        arf_out_table = fits.BinTableHDU.from_columns([col_energy,col_area],name='ENTRY')
        arf_out_table.writeto('%s/sci_area_%s_ccd%s.fits'%(output_dir,detector,ccd_id), overwrite=True)

        time_list = []
        all_evt_list = []
        sci_evt_list = []
        for t in range(0,len(lightcurve_array_all.xaxis)):
            time_list += [lightcurve_array_all.xaxis[t]]
            all_evt_list += [lightcurve_array_all.yaxis[t]]
            sci_evt_list += [lightcurve_array.yaxis[t]]
        col_time = fits.Column(name='TIME', array=time_list, format='D')
        col_allevt = fits.Column(name='ALL_EVT', array=all_evt_list, format='D')
        col_scievt = fits.Column(name='SCI_EVT', array=sci_evt_list, format='D')
        lc_out_table = fits.BinTableHDU.from_columns([col_time,col_allevt,col_scievt],name='ENTRY')
        lc_out_table.writeto('%s/sci_lightcurve_%s_ccd%s.fits'%(output_dir,detector,ccd_id), overwrite=True)

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
        if sp_free_qpb==0.:
            xray_qpb_ratio += [0.]
            spfree_data_qpb_ratio += [0.]
        else:
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
        if qpb==0.:
            data_qpb_ratio += [0.]
        else:
            data_qpb_ratio += [data/qpb]

    sp_scale_list = []
    max_sp_scale_list = []
    xray_scale_list = []
    for ch in range(1,len(energy_range)):
        xray_scale_list += [xray_qpb_ratio[ch-1]]
        data_cnt = data_pattern[ch].integral()
        qpb_cnt = qpb_pattern[ch].integral()
        xray_cnt = xray_qpb_ratio[ch-1]*xray_pattern[ch].integral()
        raw_spf_cnt = spf_pattern[ch].integral()
        if raw_spf_cnt==0.:
            sp_scale = 0.
            max_sp_scale = 0.
        else:
            sp_scale = max((data_cnt-qpb_cnt-xray_cnt)/raw_spf_cnt,0.)
            max_sp_scale = max((data_cnt-qpb_cnt)/raw_spf_cnt,0.)
        sp_scale_list += [sp_scale]
        max_sp_scale_list += [max_sp_scale]

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
            data_cnt_sum = 0.
            qpb_model_sum = 0.
            spf_model_sum = 0.
            xray_raw_sum = 0.
            spf_scale = A[0]*sp_scale_list[min_data_qpb_diff_ch-1]
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
            chi_sum = 0.
            for pattern in range(0,5):
                if pattern!=2 and pattern!=4: continue
                data_cnt = data_pattern[ch].yaxis[pattern]
                qpb_model = qpb_pattern[ch].yaxis[pattern]
                spf_model = spf_scale*spf_pattern[ch].yaxis[pattern]
                xray_model = xray_scale_sum*xray_pattern[ch].yaxis[pattern]
                if data_cnt==0: continue
                pattern_weight = 1./pow(data_cnt,0.5)
                #pattern_weight = 1./data_cnt
                chi = (data_cnt - xray_model - spf_model - qpb_model)*pattern_weight
                chi_sum += chi

            #chi2 += chi_sum
            chi2 += chi_sum*chi_sum
            
            penalty_chi = 0.
            for pattern in range(0,5):
                xray_model = xray_scale_sum*xray_pattern[ch].yaxis[pattern]
                pattern_weight = 1./pow(data_cnt,0.5)
                chi = xray_model*pattern_weight
                penalty_chi += chi
            if penalty_chi<0.:
                chi2 += penalty_chi*penalty_chi

            penalty_chi = 0.
            for pattern in range(0,5):
                spf_model = spf_scale*spf_pattern[ch].yaxis[pattern]
                pattern_weight = 1./pow(data_cnt,0.5)
                chi = spf_model*pattern_weight
                penalty_chi += chi
            if penalty_chi<0.:
                chi2 += penalty_chi*penalty_chi

            #if energy_range[ch]>=11000:  # almost X-ray free because of low effective area
            #    penalty_chi = 0.
            #    for pattern in range(0,5):
            #        data_cnt = data_pattern[ch].yaxis[pattern]
            #        qpb_model = qpb_pattern[ch].yaxis[pattern]
            #        spf_model = spf_scale*spf_pattern[ch].yaxis[pattern]
            #        xray_model = xray_scale_sum*xray_pattern[ch].yaxis[pattern]
            #        if data_cnt==0: continue
            #        pattern_weight = 1./pow(data_cnt,0.5)
            #        chi = (data_cnt - spf_model - qpb_model)*pattern_weight
            #        penalty_chi += chi
            #    chi2 += penalty_chi*penalty_chi

        return chi2

    max_sp_scale = 0.
    for ch in range(1,len(energy_range)):
        max_sp_scale += max_sp_scale_list[ch-1]
    max_sp_scale = max_sp_scale/float(len(energy_range)-1)
    if sp_scale_list[min_data_qpb_diff_ch-1]>0.:
        max_sp_scale = max_sp_scale/sp_scale_list[min_data_qpb_diff_ch-1]
    else:
        max_sp_scale = 1.

    print ('min_data_qpb_diff_ch = %s'%(min_data_qpb_diff_ch))
    print ('sp_scale_list[min_data_qpb_diff_ch-1] = %s'%(sp_scale_list[min_data_qpb_diff_ch-1]))
    print ('max_sp_scale = %s'%(max_sp_scale))
    qpb_scale = 1.
    param_init = [min(1.,max_sp_scale)]
    param_bound_upper = [max_sp_scale]
    param_bound_lower = [0.]
    #param_init = [1.0,0.1,-1.]
    #param_bound_upper = [3.,3.,0.1]
    #param_bound_lower = [0.,0.,-0.1]
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
    print ('result.x[0] = %s'%(result.x[0]))
    #print ('result.x[1] = %s'%(result.x[1]))
    #print ('result.x[2] = %s'%(result.x[2]))
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for ch in range(1,len(energy_range)):
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        sp_scale_norm = result.x[0]*sp_scale_list[min_data_qpb_diff_ch-1]
        print ('ch = %s'%(ch))
        print ('sp_scale_list[ch-1] = %s'%(sp_scale_list[ch-1]))
        print ('sp_scale_fit = %s'%(sp_scale_norm))
        spfree_data_qpb_diff = abs(spfree_data_qpb_ratio[ch-1]-1.)
        data_qpb_diff = abs(data_qpb_ratio[ch-1]-1.)
        print ('spfree_data_qpb_diff = %s'%(spfree_data_qpb_diff))
        print ('data_qpb_diff = %s'%(data_qpb_diff))
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

    

    
def make_timecut_mask(lightcurve_sci,lightcurve_cor,lightcurve_qpb,lightcurve_sci_fov_mask):

    avg_rate_fov = 0.
    for idx in range(0,len(lightcurve_sci.xaxis)):
        if lightcurve_sci_fov_mask.yaxis[idx]==1: continue
        avg_rate_fov += (lightcurve_sci.yaxis[idx])
    avg_rate_fov = avg_rate_fov/float(len(lightcurve_sci.xaxis))

    avg_rate_cor = 0.
    for idx in range(0,len(lightcurve_cor.xaxis)):
        if lightcurve_sci_fov_mask.yaxis[idx]==1: continue
        avg_rate_cor += (lightcurve_cor.yaxis[idx])
    avg_rate_cor = avg_rate_cor/float(len(lightcurve_cor.xaxis))

    rms_rate_cor = 0.
    for idx in range(0,len(lightcurve_cor.xaxis)):
        if lightcurve_sci_fov_mask.yaxis[idx]==1: continue
        rms_rate_cor += pow((lightcurve_cor.yaxis[idx])-avg_rate_cor,2)
    rms_rate_cor = pow(rms_rate_cor/float(len(lightcurve_cor.xaxis)),0.5)
    rms_rate_cor = rms_rate_cor*avg_rate_fov/avg_rate_cor

    for idx in range(0,len(lightcurve_sci.xaxis)):
        sci_rate = (lightcurve_sci.yaxis[idx])
        if (sci_rate-avg_rate_fov)/rms_rate_cor>3.0:
            lightcurve_sci_fov_mask.yaxis[idx] = 1.

    return lightcurve_sci_fov_mask


start_time = time.time()

map_color = 'coolwarm'


SP_flare_mask = True
#point_source_mask = True
#SP_flare_mask = False
point_source_mask = False

plot_tag = detector+'_'+on_filter+'_'+on_obsID


def analyze_a_ccd_chip(energy_range=[200,12000],ccd_id=0):

    #off_sample = 'extragalactic'
    #off_obsID = 'ID0820560101' # Percentage of flaring time: 15.2%
    off_sample = on_sample
    off_obsID = on_obsID
    
    xray_sample = 'extragalactic'
    xray_obsID = 'ID0690900101'
    #xray_obsID = 'ID0827251001'
    #xray_sample = 'Cas_A'
    #xray_obsID = 'ID0412180101'

    xray_sci_fov_evt_filename = '../%s/%s/analysis/%s-fov-evt.fits'%(xray_sample,xray_obsID,detector)
    xray_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-fov-evt.fits'%(xray_sample,xray_obsID,detector)
    xray_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(xray_sample,xray_obsID,detector)
    xray_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(xray_sample,xray_obsID,detector)
    xray_arf_filename = '../%s/%s/analysis/%s-fov-arf.fits'%(xray_sample,xray_obsID,detector)
    
    off_sci_fov_evt_filename = '../%s/%s/analysis/%s-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
    off_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
    off_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(off_sample,off_obsID,detector)
    off_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(off_sample,off_obsID,detector)
    off_arf_filename = '../%s/%s/analysis/%s-fov-arf.fits'%(off_sample,off_obsID,detector)
    
    on_sci_fov_evt_filename = '../%s/%s/analysis/%s-%s-evt.fits'%(on_sample,on_obsID,detector,on_filter)
    on_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-%s-evt.fits'%(on_sample,on_obsID,detector,on_filter)
    on_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(on_sample,on_obsID,detector)
    on_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(on_sample,on_obsID,detector)
    on_arf_filename = '../%s/%s/analysis/%s-fov-arf.fits'%(on_sample,on_obsID,detector)

    #print ('prepare xray sample time cuts')
    #
    #output_all_fov = read_event_file(xray_sci_fov_evt_filename,xray_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    #output_all_cor = read_event_file(xray_sci_cor_evt_filename,xray_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    #output_fwc_fov = read_event_file(xray_fwc_fov_evt_filename,xray_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    #output_fwc_cor = read_event_file(xray_fwc_cor_evt_filename,xray_arf_filename,mask_lc=None,mask_map=None,evt_filter='')

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

    #mask_lc = make_timecut_mask(xray_lightcurve_all_cor,xray_lightcurve_fwc_cor) 

    #time_pix_frac_mask = mask_lc.get_pixel_fraction()
    #if time_pix_frac_mask>0.8:
    #    mask_lc = None

    #print ('apply x-ray sample space and time masks')
    #
    #output_sci_source = read_event_file(xray_sci_fov_evt_filename,xray_arf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto source',energy_range=energy_range,ccd_id=0)
    #output_sci_ring   = read_event_file(xray_sci_fov_evt_filename,xray_arf_filename,mask_lc=mask_lc,mask_map=None,evt_filter='sp-veto ring',energy_range=energy_range,ccd_id=0)
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

    print ('prepare on sample time cuts')
    
    output_all_nsp = read_event_file(on_sci_fov_evt_filename,on_arf_filename,mask_lc=None,mask_map=None,evt_filter='sp-free')
    output_fwc_nsp = read_event_file(on_fwc_fov_evt_filename,on_arf_filename,mask_lc=None,mask_map=None,evt_filter='sp-free')
    output_all_fov = read_event_file(on_sci_fov_evt_filename,on_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,on_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_all_cor = read_event_file(on_sci_cor_evt_filename,on_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,on_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    
    on_duration_all_nsp = output_all_nsp[0]
    on_evt_count_all_nsp = output_all_nsp[1]
    on_lightcurve_all_nsp = output_all_nsp[2]
    on_pattern_all_nsp = output_all_nsp[3]
    on_spectrum_all_nsp = output_all_nsp[4]
    on_detx_all_nsp = output_all_nsp[5]
    on_image_all_nsp = output_all_nsp[6]
    
    on_duration_fwc_nsp = output_fwc_nsp[0]
    on_evt_count_fwc_nsp = output_fwc_nsp[1]
    on_lightcurve_fwc_nsp = output_fwc_nsp[2]
    on_pattern_fwc_nsp = output_fwc_nsp[3]
    on_spectrum_fwc_nsp = output_fwc_nsp[4]
    on_detx_fwc_nsp = output_fwc_nsp[5]
    on_image_fwc_nsp = output_fwc_nsp[6]
    
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
    
    #fwc_2_sci_ratio = on_image_all_cor[0].integral()/on_image_fwc_cor[0].integral()
    fwc_2_sci_ratio = 0.45
    print ('Before timecut, ON data fwc_2_sci_ratio = %s'%(fwc_2_sci_ratio))

    on_lightcurve_all_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    on_lightcurve_fwc_fov.scale(fwc_2_sci_ratio)
    on_lightcurve_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    on_lightcurve_fwc_nsp.scale(fwc_2_sci_ratio)

    on_mask_lc = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
    make_timecut_mask(on_lightcurve_all_fov,on_lightcurve_all_cor,on_lightcurve_fwc_cor,on_mask_lc) 
    make_timecut_mask(on_lightcurve_all_fov,on_lightcurve_all_cor,on_lightcurve_fwc_cor,on_mask_lc) 

    output_sci_fov = read_event_file(on_sci_fov_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range)
    output_spf_fov = read_event_file(on_sci_fov_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_range)
    on_evt_count_sci_fov = output_sci_fov[1]
    on_evt_count_spf_fov = output_spf_fov[1]
    
    print ('on_evt_count_sci_fov = %s'%(on_evt_count_sci_fov))
    print ('on_evt_count_spf_fov = %s'%(on_evt_count_spf_fov))
    if on_evt_count_spf_fov<1.0*on_evt_count_sci_fov:
        print ('not enough SPF data, use template OFF data')
        off_sample = 'extragalactic'
        off_obsID = 'ID0820560101' # Percentage of flaring time: 15.2%
        off_sci_fov_evt_filename = '../%s/%s/analysis/%s-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
        off_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
        off_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(off_sample,off_obsID,detector)
        off_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(off_sample,off_obsID,detector)
        off_arf_filename = '../%s/%s/analysis/%s-fov-arf.fits'%(off_sample,off_obsID,detector)

    time_pix_frac_mask = on_mask_lc.get_pixel_fraction()
    print ('time_pix_frac_mask = %s'%(time_pix_frac_mask))
    if time_pix_frac_mask>0.8:
        on_mask_lc = None
        print ('not enough SPF data, use template OFF data')
        off_sample = 'extragalactic'
        off_obsID = 'ID0820560101' # Percentage of flaring time: 15.2%
        off_sci_fov_evt_filename = '../%s/%s/analysis/%s-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
        off_fwc_fov_evt_filename = '../%s/%s/analysis/%s-fwc-%s-evt.fits'%(off_sample,off_obsID,detector,on_filter)
        off_sci_cor_evt_filename = '../%s/%s/analysis/%s-cor-evt.fits'%(off_sample,off_obsID,detector)
        off_fwc_cor_evt_filename = '../%s/%s/analysis/%s-fwc-cor-evt.fits'%(off_sample,off_obsID,detector)
        off_arf_filename = '../%s/%s/analysis/%s-fov-arf.fits'%(off_sample,off_obsID,detector)
    if not SP_flare_mask:
        on_mask_lc = None

    print ('use off sample: %s'%(off_sci_fov_evt_filename))
    print ('prepare off sample time cuts')
    
    output_all_nsp = read_event_file(off_sci_fov_evt_filename,off_arf_filename,mask_lc=None,mask_map=None,evt_filter='sp-free')
    output_fwc_nsp = read_event_file(off_fwc_fov_evt_filename,off_arf_filename,mask_lc=None,mask_map=None,evt_filter='sp-free')
    output_all_fov = read_event_file(off_sci_fov_evt_filename,off_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,off_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_all_cor = read_event_file(off_sci_cor_evt_filename,off_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,off_arf_filename,mask_lc=None,mask_map=None,evt_filter='')
    
    off_duration_all_nsp = output_all_nsp[0]
    off_evt_count_all_nsp = output_all_nsp[1]
    off_lightcurve_all_nsp = output_all_nsp[2]
    off_pattern_all_nsp = output_all_nsp[3]
    off_spectrum_all_nsp = output_all_nsp[4]
    off_detx_all_nsp = output_all_nsp[5]
    off_image_all_nsp = output_all_nsp[6]
    
    off_duration_fwc_nsp = output_fwc_nsp[0]
    off_evt_count_fwc_nsp = output_fwc_nsp[1]
    off_lightcurve_fwc_nsp = output_fwc_nsp[2]
    off_pattern_fwc_nsp = output_fwc_nsp[3]
    off_spectrum_fwc_nsp = output_fwc_nsp[4]
    off_detx_fwc_nsp = output_fwc_nsp[5]
    off_image_fwc_nsp = output_fwc_nsp[6]
    
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
    
    #fwc_2_sci_ratio = off_image_all_cor[0].integral()/off_image_fwc_cor[0].integral()
    fwc_2_sci_ratio = 0.45
    print ('Before timecut, OFF data fwc_2_sci_ratio = %s'%(fwc_2_sci_ratio))

    off_lightcurve_all_cor.scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    off_lightcurve_fwc_fov.scale(fwc_2_sci_ratio)
    off_lightcurve_fwc_cor.scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    off_lightcurve_fwc_nsp.scale(fwc_2_sci_ratio)

    off_mask_lc = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
    make_timecut_mask(off_lightcurve_all_fov,off_lightcurve_all_cor,off_lightcurve_fwc_cor,off_mask_lc) 
    make_timecut_mask(off_lightcurve_all_fov,off_lightcurve_all_cor,off_lightcurve_fwc_cor,off_mask_lc) 

    if diagnostic_plots:
        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(off_lightcurve_all_fov.xaxis,off_lightcurve_all_fov.yaxis,yerr=off_lightcurve_all_fov.yerr,color='k',label='All FoV')
        axbig.errorbar(off_lightcurve_all_cor.xaxis,off_lightcurve_all_cor.yaxis,yerr=off_lightcurve_all_cor.yerr,color='green',label='All Cor')
        axbig.errorbar(off_lightcurve_fwc_cor.xaxis,off_lightcurve_fwc_cor.yaxis,yerr=off_lightcurve_fwc_cor.yerr,color='blue',label='FWC Cor')
        axbig.set_yscale('log')
        axbig.set_xlabel('Time')
        axbig.legend(loc='best')
        fig.savefig("%s/lightcurve_off_pre_timecut_ch%s_ccd%s_%s.png"%(output_dir,energy_range[0],ccd_id,plot_tag),bbox_inches='tight')
        axbig.remove()

    print ('apply off sample space and time masks')
    
    output_sci_fov = read_event_file(off_sci_fov_evt_filename,off_arf_filename,mask_lc=off_mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range)
    output_sci_cor = read_event_file(off_sci_cor_evt_filename,off_arf_filename,mask_lc=off_mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range)
    
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
    
    output_spf_fov = read_event_file(off_sci_fov_evt_filename,off_arf_filename,mask_lc=off_mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_range)
    output_spf_cor = read_event_file(off_sci_cor_evt_filename,off_arf_filename,mask_lc=off_mask_lc,mask_map=None,evt_filter='sp-select',energy_range=energy_range)
    
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
    
    output_fwc_fov = read_event_file(off_fwc_fov_evt_filename,off_arf_filename,mask_lc=off_mask_lc,mask_map=None,evt_filter='',energy_range=energy_range)
    output_fwc_cor = read_event_file(off_fwc_cor_evt_filename,off_arf_filename,mask_lc=off_mask_lc,mask_map=None,evt_filter='',energy_range=energy_range)
    
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

        fwc_2_sci_ratio = off_image_sci_cor[ch].integral()/off_image_fwc_cor[ch].integral()
        print ('After timecut, ch = %s, OFF data fwc_2_sci_ratio = %s'%(ch,fwc_2_sci_ratio))

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

        off_image_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        off_image_spf_fov[ch].scale(time_expo_sci_fov/time_expo_spf_fov)
        off_image_spf_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_spf_cor)
        off_image_fwc_fov[ch].scale(fwc_2_sci_ratio)
        off_image_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    
        off_detx_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        off_detx_spf_fov[ch].scale(time_expo_sci_fov/time_expo_spf_fov)
        off_detx_spf_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor*time_expo_sci_fov/time_expo_spf_cor)
        off_detx_fwc_fov[ch].scale(fwc_2_sci_ratio)
        off_detx_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
    
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
    sp_image_fov_template = []
    sp_detx_fov_template = []
    for ch in range(0,len(energy_range)):
        sp_pattern_fov_template += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
        sp_spectrum_fov_template += [MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)]
        sp_image_fov_template += [MyArray2D()]
        sp_detx_fov_template += [MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)]
    
    for ch in range(0,len(energy_range)):
        sp_spectrum_fov_template[ch].add(off_spectrum_spf_fov[ch])
        sp_spectrum_fov_template[ch].add(off_spectrum_sci_fov[ch],factor=-1.)
        sp_image_fov_template[ch].add(off_image_spf_fov[ch])
        sp_image_fov_template[ch].add(off_image_sci_fov[ch],factor=-1.)
        sp_detx_fov_template[ch].add(off_detx_spf_fov[ch])
        sp_detx_fov_template[ch].add(off_detx_sci_fov[ch],factor=-1.)
        ch_norm = sp_spectrum_fov_template[ch].integral()/off_pattern_spf_fov[0].integral()
        sp_pattern_fov_template[ch].add(off_pattern_spf_fov[0],factor=ch_norm)
        sp_pattern_fov_template[ch].add(off_pattern_sci_fov[0],factor=-1.*ch_norm)

    print ('apply on sample space and time masks')
    
    output_sci_fov = read_event_file(on_sci_fov_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=None,write_events=True,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    output_sci_cor = read_event_file(on_sci_cor_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=None,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    
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
    
    output_fwc_fov = read_event_file(on_fwc_fov_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    output_fwc_cor = read_event_file(on_fwc_cor_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=None,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)
    
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

    ### Remove strong X-ray hot spots ###
    image_det_mask = MyArray2D()
    find_point_sources(on_image_sci_fov[0],image_det_mask)
    find_point_sources(on_image_sci_fov[0],image_det_mask)

    output_sci_fov_mask = read_event_file(on_sci_fov_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=image_det_mask,evt_filter='sp-veto',energy_range=energy_range,ccd_id=ccd_id)
    output_fwc_fov_mask = read_event_file(on_fwc_fov_evt_filename,on_arf_filename,mask_lc=on_mask_lc,mask_map=image_det_mask,evt_filter='',energy_range=energy_range,ccd_id=ccd_id)

    on_duration_sci_fov_mask = output_sci_fov_mask[0]
    on_evt_count_sci_fov_mask = output_sci_fov_mask[1]
    on_lightcurve_sci_fov_mask = output_sci_fov_mask[2]
    on_pattern_sci_fov_mask = output_sci_fov_mask[3]
    on_spectrum_sci_fov_mask = output_sci_fov_mask[4]
    on_detx_sci_fov_mask = output_sci_fov_mask[5]
    on_image_sci_fov_mask = output_sci_fov_mask[6]
    
    on_duration_fwc_fov_mask = output_fwc_fov_mask[0]
    on_evt_count_fwc_fov_mask = output_fwc_fov_mask[1]
    on_lightcurve_fwc_fov_mask = output_fwc_fov_mask[2]
    on_pattern_fwc_fov_mask = output_fwc_fov_mask[3]
    on_spectrum_fwc_fov_mask = output_fwc_fov_mask[4]
    on_detx_fwc_fov_mask = output_fwc_fov_mask[5]
    on_image_fwc_fov_mask = output_fwc_fov_mask[6]
    ### Remove strong X-ray hot spots ###
    

    area_pix_frac_fwc_fov_mask = on_image_fwc_fov_mask[0].integral()
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

        fwc_2_sci_ratio = on_image_sci_cor[ch].integral()/on_image_fwc_cor[ch].integral()
        print ('After timecut, ch = %s, OFF data fwc_2_sci_ratio = %s'%(ch,fwc_2_sci_ratio))

        on_spectrum_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_spectrum_fwc_fov[ch].scale(fwc_2_sci_ratio)
        on_spectrum_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_pattern_sci_cor[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_pattern_fwc_fov[ch].scale(fwc_2_sci_ratio)
        on_pattern_fwc_cor[ch].scale(fwc_2_sci_ratio*area_pix_frac_fwc_fov/area_pix_frac_fwc_cor)
        on_detx_fwc_fov[ch].scale(fwc_2_sci_ratio)
        on_image_fwc_fov[ch].scale(fwc_2_sci_ratio)

        on_pattern_sci_fov_mask[ch].scale(area_pix_frac_fwc_fov/area_pix_frac_fwc_fov_mask)

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
        axbig.errorbar(on_lightcurve_fwc_cor.xaxis,on_lightcurve_fwc_cor.yaxis,yerr=on_lightcurve_fwc_cor.yerr,color='blue',label='FWC Cor')
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
        #ch_norm = on_spectrum_fwc_fov[ch].integral()/on_pattern_fwc_fov[ch].integral()
        ch_norm = 1.
        qpb_pattern_fov_template[ch].add(on_pattern_fwc_fov[ch],factor=ch_norm)

    for ch in range(0,len(energy_range)):
        sp_rescale = qpb_pattern_fov_template[0].integral()/sp_pattern_fov_template[0].integral()
        sp_pattern_fov_template[ch].scale(sp_rescale)

    xray_pattern_fov_template = []
    for ch in range(0,len(energy_range)):
        xray_pattern_fov_template += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
    for ch in range(0,len(energy_range)):
        xray_pattern_fov_template[ch].add(qpb_pattern_fov_template[ch])

    #xray_cnt_list, spf_cnt_list, qpb_cnt_list = fit_pattern(energy_range,on_pattern_sci_fov_mask,xray_pattern_fov_template,sp_pattern_fov_template,qpb_pattern_fov_template)
    xray_cnt_list, spf_cnt_list, qpb_cnt_list = fit_pattern(energy_range,on_pattern_sci_fov,xray_pattern_fov_template,sp_pattern_fov_template,qpb_pattern_fov_template)
        
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
    
    sp_image_fov_template_sum = MyArray2D()
    qpb_image_fov_template_sum = MyArray2D()
    for ch in range(1,len(energy_range)):
        sp_image_fov_template_sum.add(sp_image_fov_template[ch])
        qpb_image_fov_template_sum.add(qpb_image_fov_template[ch])
    #sp_image_fov_template_sum.flattening()
    qpb_image_fov_template_sum.flattening()

    on_image_sci_fov_bkg_sum = MyArray2D()
    for ch in range(1,len(energy_range)):
        on_image_sci_fov_bkg_sum.add(sp_image_fov_template[ch])
        on_image_sci_fov_bkg_sum.add(qpb_image_fov_template[ch])

    all_pattern_fov_template = []
    for ch in range(0,len(energy_range)):
        all_pattern_fov_template += [MyArray1D(bin_start=pattern_low,bin_end=pattern_high,pixel_scale=pattern_scale)]
        all_pattern_fov_template[ch].add(xray_pattern_fov_template[ch])
        all_pattern_fov_template[ch].add(qpb_pattern_fov_template[ch])
        all_pattern_fov_template[ch].add(sp_pattern_fov_template[ch])

    if diagnostic_plots:
        for ch in range(1,len(energy_range)):
            fig.clf()
            axbig = fig.add_subplot()
            axbig.errorbar(on_pattern_sci_fov[ch].xaxis,on_pattern_sci_fov[ch].yaxis,yerr=on_pattern_sci_fov[ch].yerr,color='k',label='Data')
            axbig.errorbar(xray_pattern_fov_template[ch].xaxis,xray_pattern_fov_template[ch].yaxis,yerr=xray_pattern_fov_template[ch].yerr,color='g',label='X-ray')
            axbig.errorbar(qpb_pattern_fov_template[ch].xaxis,qpb_pattern_fov_template[ch].yaxis,yerr=qpb_pattern_fov_template[ch].yerr,color='b',label='QPB')
            axbig.bar(sp_pattern_fov_template[ch].xaxis,sp_pattern_fov_template[ch].yaxis,color='orange',label='SP Flare')
            axbig.errorbar(all_pattern_fov_template[ch].xaxis,all_pattern_fov_template[ch].yaxis,yerr=all_pattern_fov_template[ch].yerr,color='r',label='All')
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

    #prob_bkg_spectrum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    #prob_bkg_spectrum.add(ana_output_spectrum[1])
    #prob_bkg_spectrum.normalize2()

    #evt_detx_list = []
    #evt_dety_list = []
    #evt_ra_list = []
    #evt_dec_list = []
    ##evt_l_list = []
    ##evt_b_list = []
    #evt_pi_list = []
    #for idx_x in range(0,len(ana_output_image[1].xaxis)):
    #    for idx_y in range(0,len(ana_output_image[1].yaxis)):
    #        pix_detx = ana_output_image[1].xaxis[idx_x]
    #        pix_dety = ana_output_image[1].yaxis[idx_y]
    #        sum_bkg = ana_output_image[1].zaxis[idx_x,idx_y]
    #        for evt in range(0,round(sample_scale*sum_bkg)):
    #            evt_detx = np.random.uniform(low=pix_detx, high=pix_detx+detx_scale)
    #            evt_dety = np.random.uniform(low=pix_dety, high=pix_dety+detx_scale)
    #            #ref_det_idx_x, ref_det_idx_y = find_nearest_ref_det_idx([evt_detx,evt_dety],ref_det)
    #            ref_det_idx = 0
    #            if on_obsID=='ID0412180101' or on_obsID=='ID0400210101':
    #                ref_det_idx = 1
    #            evt_ra, evt_dec = ConvertDet2Sky([evt_detx,evt_dety],ref_sky[ref_det_idx],ref_det[ref_det_idx],mtx_det_2_sky[ref_det_idx])
    #            #evt_l, evt_b = ConvertRaDecToGalactic(evt_ra,evt_dec)
    #            evt_pi = np.random.choice(a=prob_bkg_spectrum.xaxis, p=prob_bkg_spectrum.yaxis)
    #            evt_detx_list += [evt_detx]
    #            evt_dety_list += [evt_dety]
    #            evt_ra_list += [evt_ra]
    #            evt_dec_list += [evt_dec]
    #            #evt_l_list += [evt_l]
    #            #evt_b_list += [evt_b]
    #            evt_pi_list += [evt_pi]
    #col_detx = fits.Column(name='DETX', array=evt_detx_list, format='D')
    #col_dety = fits.Column(name='DETY', array=evt_dety_list, format='D')
    #col_ra = fits.Column(name='RA', array=evt_ra_list, format='D')
    #col_dec = fits.Column(name='DEC', array=evt_dec_list, format='D')
    ##col_l = fits.Column(name='GalL', array=evt_l_list, format='D')
    ##col_b = fits.Column(name='GalB', array=evt_b_list, format='D')
    #col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
    ##my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_l,col_b,col_pi],name='EVENT')
    #my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_pi],name='EVENT')
    #my_table.writeto('%s/bkg_events_%s_ccd%s.fits'%(output_dir,detector,ana_ccd_bins[ccd]), overwrite=True)

    prob_spf_spectrum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    prob_spf_spectrum.add(ana_output_spectrum[2])
    prob_spf_spectrum.normalize2()

    evt_detx_list = []
    evt_dety_list = []
    evt_ra_list = []
    evt_dec_list = []
    #evt_l_list = []
    #evt_b_list = []
    evt_pi_list = []
    for idx_x in range(0,len(ana_output_image[2].xaxis)):
        for idx_y in range(0,len(ana_output_image[2].yaxis)):
            pix_detx = ana_output_image[2].xaxis[idx_x]
            pix_dety = ana_output_image[2].yaxis[idx_y]
            sum_bkg = ana_output_image[2].zaxis[idx_x,idx_y]
            for evt in range(0,round(sample_scale*sum_bkg)):
                evt_detx = np.random.uniform(low=pix_detx, high=pix_detx+detx_scale)
                evt_dety = np.random.uniform(low=pix_dety, high=pix_dety+detx_scale)
                #ref_det_idx_x, ref_det_idx_y = find_nearest_ref_det_idx([evt_detx,evt_dety],ref_det)
                ref_det_idx = 0
                if on_obsID=='ID0412180101' or on_obsID=='ID0400210101':
                    ref_det_idx = 1
                evt_ra, evt_dec = ConvertDet2Sky([evt_detx,evt_dety],ref_sky[ref_det_idx],ref_det[ref_det_idx],mtx_det_2_sky[ref_det_idx])
                #evt_l, evt_b = ConvertRaDecToGalactic(evt_ra,evt_dec)
                evt_pi = np.random.choice(a=prob_spf_spectrum.xaxis, p=prob_spf_spectrum.yaxis)
                evt_detx_list += [evt_detx]
                evt_dety_list += [evt_dety]
                evt_ra_list += [evt_ra]
                evt_dec_list += [evt_dec]
                #evt_l_list += [evt_l]
                #evt_b_list += [evt_b]
                evt_pi_list += [evt_pi]
    col_detx = fits.Column(name='DETX', array=evt_detx_list, format='D')
    col_dety = fits.Column(name='DETY', array=evt_dety_list, format='D')
    col_ra = fits.Column(name='RA', array=evt_ra_list, format='D')
    col_dec = fits.Column(name='DEC', array=evt_dec_list, format='D')
    #col_l = fits.Column(name='GalL', array=evt_l_list, format='D')
    #col_b = fits.Column(name='GalB', array=evt_b_list, format='D')
    col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
    #my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_l,col_b,col_pi],name='EVENT')
    my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_pi],name='EVENT')
    my_table.writeto('%s/spf_events_%s_ccd%s.fits'%(output_dir,detector,ana_ccd_bins[ccd]), overwrite=True)

    prob_qpb_spectrum = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    prob_qpb_spectrum.add(ana_output_spectrum[3])
    prob_qpb_spectrum.normalize2()

    evt_detx_list = []
    evt_dety_list = []
    evt_ra_list = []
    evt_dec_list = []
    #evt_l_list = []
    #evt_b_list = []
    evt_pi_list = []
    for idx_x in range(0,len(ana_output_image[3].xaxis)):
        for idx_y in range(0,len(ana_output_image[3].yaxis)):
            pix_detx = ana_output_image[3].xaxis[idx_x]
            pix_dety = ana_output_image[3].yaxis[idx_y]
            sum_bkg = ana_output_image[3].zaxis[idx_x,idx_y]
            for evt in range(0,round(sample_scale*sum_bkg)):
                evt_detx = np.random.uniform(low=pix_detx, high=pix_detx+detx_scale)
                evt_dety = np.random.uniform(low=pix_dety, high=pix_dety+detx_scale)
                #ref_det_idx_x, ref_det_idx_y = find_nearest_ref_det_idx([evt_detx,evt_dety],ref_det)
                ref_det_idx = 0
                if on_obsID=='ID0412180101' or on_obsID=='ID0400210101':
                    ref_det_idx = 1
                evt_ra, evt_dec = ConvertDet2Sky([evt_detx,evt_dety],ref_sky[ref_det_idx],ref_det[ref_det_idx],mtx_det_2_sky[ref_det_idx])
                #evt_l, evt_b = ConvertRaDecToGalactic(evt_ra,evt_dec)
                evt_pi = np.random.choice(a=prob_qpb_spectrum.xaxis, p=prob_qpb_spectrum.yaxis)
                evt_detx_list += [evt_detx]
                evt_dety_list += [evt_dety]
                evt_ra_list += [evt_ra]
                evt_dec_list += [evt_dec]
                #evt_l_list += [evt_l]
                #evt_b_list += [evt_b]
                evt_pi_list += [evt_pi]
    col_detx = fits.Column(name='DETX', array=evt_detx_list, format='D')
    col_dety = fits.Column(name='DETY', array=evt_dety_list, format='D')
    col_ra = fits.Column(name='RA', array=evt_ra_list, format='D')
    col_dec = fits.Column(name='DEC', array=evt_dec_list, format='D')
    #col_l = fits.Column(name='GalL', array=evt_l_list, format='D')
    #col_b = fits.Column(name='GalB', array=evt_b_list, format='D')
    col_pi = fits.Column(name='PI', array=evt_pi_list, format='D')
    #my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_l,col_b,col_pi],name='EVENT')
    my_table = fits.BinTableHDU.from_columns([col_detx,col_dety,col_ra,col_dec,col_pi],name='EVENT')
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


