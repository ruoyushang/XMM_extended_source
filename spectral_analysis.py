import sys
from astropy.io import fits
from astropy.table import Table
from astropy import units as my_unit
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
import numpy as np
from numpy.linalg import inv
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import common_functions

#on_sample = 'extragalactic'
#on_obsID = 'ID0505460501'
#on_obsID = 'ID0690900101'
#on_obsID = 'ID0803160301'
#on_obsID = 'ID0820560101'
#on_obsID = 'ID0827200401'
#on_obsID = 'ID0827200501'
#on_obsID = 'ID0827241101'
#on_obsID = 'ID0827251001'
#on_obsID = 'ID0827251101'

on_sample = 'Cas_A'
on_obsID = 'ID0412180101'
#on_obsID = 'ID0400210101' # Cas A Northern lobe
#on_obsID = 'ID0782961401' # angular distance to Cas A: 34.7 arcmin

#on_sample = '3HWC_J1928_p178'
#on_obsID = 'ID0902120101'

#detector = 'mos1'
detector = 'mos2'

ana_ccd_bins = [0]
#ana_ccd_bins = [1,2,3,4,5,6,7]

#exclusion_inner = 0.
#exclusion_outer = 0.15
exclusion_inner = 0.15
exclusion_outer = 1e10

point_source_cut = False

energy_cut = 2000

energy_array = [2000,4000,6000,8000,10000,12000]

roi_ra = 350.85
roi_dec = 58.815

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
find_nearest_ref_sky_idx = common_functions.find_nearest_ref_sky_idx
LoadCoordinateMatrix = common_functions.LoadCoordinateMatrix

def ConvertDet2Sky(target_det,origin_sky,origin_det,mtx_conv_det_2_sky):

    target_det = np.array(target_det)
    origin_sky = np.array(origin_sky)
    origin_det = np.array(origin_det)

    delta_det = target_det-origin_det
    delta_sky = np.matmul(mtx_conv_det_2_sky,delta_det)
    target_sky = delta_sky+origin_sky

    return target_sky[0,0], target_sky[0,1]

def ConvertSky2Det(target_sky,origin_sky,origin_det,mtx_conv_sky_2_det):

    target_sky = np.array(target_sky)
    origin_sky = np.array(origin_sky)
    origin_det = np.array(origin_det)

    delta_sky = target_sky-origin_sky
    delta_det = np.matmul(mtx_conv_sky_2_det,delta_sky)
    target_det = delta_det+origin_det

    return target_det[0,0], target_det[0,1]


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

def DistanceToROI(ra,dec):

    ref_sky_idx = 0
    if on_obsID=='ID0412180101' or on_obsID=='ID0400210101':
        ref_sky_idx = 1

    evt_detx, evt_dety = ConvertSky2Det([ra,dec],ref_sky[ref_sky_idx],ref_det[ref_sky_idx],mtx_sky_2_det[ref_sky_idx])
    roi_detx, roi_dety = ConvertSky2Det([roi_ra,roi_dec],ref_sky[ref_sky_idx],ref_det[ref_sky_idx],mtx_sky_2_det[ref_sky_idx])

    return pow(pow(evt_detx-roi_detx,2)+pow(evt_dety-roi_dety,2),0.5)*0.05/(60.*60.)

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*my_unit.deg, b*my_unit.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def ConvertRaDecToGalactic(ra, dec):
    my_sky = SkyCoord(ra*my_unit.deg, dec*my_unit.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg

def ConvertRaDecMapToGalacticMap(radec_map, galactic_map):

    for idx_x in range(0,len(radec_map.xaxis)):
        for idx_y in range(0,len(radec_map.yaxis)):
            pix_ra = radec_map.xaxis[idx_x]
            pix_dec = radec_map.yaxis[idx_y]
            pix_content = radec_map.zaxis[idx_x,idx_y]
            pix_l, pix_b = ConvertRaDecToGalactic(pix_ra,pix_dec)
            galactic_map.fill(pix_l,pix_b,weight=pix_content)



sky_ra_center = 0.
sky_dec_center = 0.

input_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/'+on_obsID
output_dir = '/Users/rshang/xmm_analysis/output_plots/'

image_det_mask = MyArray2D(pixel_scale=500)
image_det_sci = MyArray2D(pixel_scale=500)
image_det_bkg = MyArray2D(pixel_scale=500)
for ccd in range(0,len(ana_ccd_bins)):

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    sci_hdu_list = fits.open(sci_filename)
    sci_table = Table.read(sci_filename, hdu=1)
    for entry in range(0,len(sci_table)):
        evt_pi = sci_table[entry]['PI']
        evt_detx = sci_table[entry]['DETX']
        evt_dety = sci_table[entry]['DETY']
        if evt_pi<energy_cut: continue
        image_det_sci.fill(evt_detx,evt_dety)

    sky_ra_center = sci_hdu_list[1].header['REF_RA']
    sky_dec_center = sci_hdu_list[1].header['REF_DEC']
    sky_l_center, sky_b_center = ConvertRaDecToGalactic(sky_ra_center, sky_dec_center)

    bkg_filename = '%s/bkg_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    bkg_hdu_list = fits.open(bkg_filename)
    bkg_table = Table.read(bkg_filename, hdu=1)
    for entry in range(0,len(bkg_table)):
        evt_pi = bkg_table[entry]['PI']
        evt_detx = bkg_table[entry]['DETX']
        evt_dety = bkg_table[entry]['DETY']
        if evt_pi<energy_cut: continue
        image_det_bkg.fill(evt_detx,evt_dety,weight=0.1)

sky_ra_low = sky_ra_center-0.28
sky_ra_high = sky_ra_center+0.28
sky_dec_low = sky_dec_center-0.28
sky_dec_high = sky_dec_center+0.28
sky_l_low = sky_l_center-0.28
sky_l_high = sky_l_center+0.28
sky_b_low = sky_b_center-0.28
sky_b_high = sky_b_center+0.28


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
            if significance>3.:
                image_mask.zaxis[idx_x,idx_y] = 1


if point_source_cut:
    find_point_sources(image_det_sci,image_det_mask)
    find_point_sources(image_det_sci,image_det_mask)


image_det_sci = MyArray2D()
image_det_bkg = MyArray2D()
image_icrs_sci = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=0.56,pixel_scale=500.*0.05/(60.*60.))
image_galactic_sci = MyArray2D(start_x=sky_l_low,start_y=sky_b_low,image_size=0.56,pixel_scale=500.*0.05/(60.*60.))
spectrum_sci = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_bkg = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_qpb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_spf = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
detx_sci = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_bkg = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_qpb = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_spf = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
for ccd in range(0,len(ana_ccd_bins)):

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    sci_hdu_list = fits.open(sci_filename)
    sci_table = Table.read(sci_filename, hdu=1)
    for entry in range(0,len(sci_table)):
        evt_pi = sci_table[entry]['PI']
        evt_detx = sci_table[entry]['DETX']
        evt_dety = sci_table[entry]['DETY']
        evt_ra = sci_table[entry]['RA']
        evt_dec = sci_table[entry]['DEC']
        distance_2_roi = DistanceToROI(evt_ra,evt_dec)
        if distance_2_roi<exclusion_inner or distance_2_roi>exclusion_outer: continue
        if evt_pi<energy_cut: continue
        mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if mask==1: continue
        image_det_sci.fill(evt_detx,evt_dety)
        image_icrs_sci.fill(evt_ra,evt_dec)
        spectrum_sci.fill(evt_pi)
        detx_sci.fill(evt_detx)

    bkg_filename = '%s/bkg_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    bkg_hdu_list = fits.open(bkg_filename)
    bkg_table = Table.read(bkg_filename, hdu=1)
    for entry in range(0,len(bkg_table)):
        evt_pi = bkg_table[entry]['PI']
        evt_detx = bkg_table[entry]['DETX']
        evt_dety = bkg_table[entry]['DETY']
        evt_ra = bkg_table[entry]['RA']
        evt_dec = bkg_table[entry]['DEC']
        distance_2_roi = DistanceToROI(evt_ra,evt_dec)
        if distance_2_roi<exclusion_inner or distance_2_roi>exclusion_outer: continue
        if evt_pi<energy_cut: continue
        mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if mask==1: continue
        image_det_bkg.fill(evt_detx,evt_dety,weight=0.1)
        spectrum_bkg.fill(evt_pi,weight=0.1)
        detx_bkg.fill(evt_detx,weight=0.1)

    qpb_filename = '%s/qpb_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    qpb_hdu_list = fits.open(qpb_filename)
    qpb_table = Table.read(qpb_filename, hdu=1)
    for entry in range(0,len(qpb_table)):
        evt_pi = qpb_table[entry]['PI']
        evt_detx = qpb_table[entry]['DETX']
        evt_dety = qpb_table[entry]['DETY']
        evt_ra = qpb_table[entry]['RA']
        evt_dec = qpb_table[entry]['DEC']
        distance_2_roi = DistanceToROI(evt_ra,evt_dec)
        if distance_2_roi<exclusion_inner or distance_2_roi>exclusion_outer: continue
        if evt_pi<energy_cut: continue
        mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if mask==1: continue
        spectrum_qpb.fill(evt_pi,weight=0.1)
        detx_qpb.fill(evt_detx,weight=0.1)

    spf_filename = '%s/spf_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    spf_hdu_list = fits.open(spf_filename)
    spf_table = Table.read(spf_filename, hdu=1)
    for entry in range(0,len(spf_table)):
        evt_pi = spf_table[entry]['PI']
        evt_detx = spf_table[entry]['DETX']
        evt_dety = spf_table[entry]['DETY']
        evt_ra = spf_table[entry]['RA']
        evt_dec = spf_table[entry]['DEC']
        distance_2_roi = DistanceToROI(evt_ra,evt_dec)
        if distance_2_roi<exclusion_inner or distance_2_roi>exclusion_outer: continue
        if evt_pi<energy_cut: continue
        mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if mask==1: continue
        spectrum_spf.fill(evt_pi,weight=0.1)
        detx_spf.fill(evt_detx,weight=0.1)

sci_bkg_ratio = []
for ch in range(0,len(energy_array)-1):
    sci_cnt = spectrum_sci.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    bkg_cnt = spectrum_bkg.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    if bkg_cnt>0.:
        sci_bkg_ratio += [sci_cnt/bkg_cnt]
print ('sci_bkg_ratio = %s'%(sci_bkg_ratio))

fig, ax = plt.subplots()
figsize_x = 8.
figsize_y = 6.
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

map_color = 'coolwarm'

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_det_mask.xaxis.min()
xmax = image_det_mask.xaxis.max()
ymin = image_det_mask.yaxis.min()
ymax = image_det_mask.yaxis.max()
axbig.imshow(image_det_mask.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_det_mask.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_det_sci.xaxis.min()
xmax = image_det_sci.xaxis.max()
ymin = image_det_sci.yaxis.min()
ymax = image_det_sci.yaxis.max()
axbig.imshow(image_det_sci.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_det_sci.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'RA'
label_y = 'DEC'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_icrs_sci.xaxis.min()
xmax = image_icrs_sci.xaxis.max()
ymin = image_icrs_sci.yaxis.min()
ymax = image_icrs_sci.yaxis.max()
axbig.imshow(image_icrs_sci.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_icrs_sci.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

ConvertRaDecMapToGalacticMap(image_icrs_sci,image_galactic_sci)
fig.clf()
axbig = fig.add_subplot()
label_x = 'Gal. l'
label_y = 'Gal. b'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_galactic_sci.xaxis.min()
xmax = image_galactic_sci.xaxis.max()
ymin = image_galactic_sci.yaxis.min()
ymax = image_galactic_sci.yaxis.max()
axbig.imshow(image_galactic_sci.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_galactic_sci.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_det_sci.xaxis.min()
xmax = image_det_sci.xaxis.max()
ymin = image_det_sci.yaxis.min()
ymax = image_det_sci.yaxis.max()
axbig.imshow(image_det_bkg.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_det_bkg.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(spectrum_sci.xaxis,spectrum_sci.yaxis,yerr=spectrum_sci.yerr,color='k',label='Data')
axbig.plot(spectrum_qpb.xaxis,spectrum_qpb.yaxis,color='blue',label='QPB')
axbig.plot(spectrum_spf.xaxis,spectrum_spf.yaxis,color='green',label='SPF')
axbig.plot(spectrum_bkg.xaxis,spectrum_bkg.yaxis,color='red',label='Bkg')
axbig.set_yscale('log')
#axbig.set_ylim(bottom=1)
axbig.set_xlabel('Energy [eV]')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_spectrum_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(detx_sci.xaxis,detx_sci.yaxis,yerr=detx_sci.yerr,color='k',label='Data')
axbig.plot(detx_qpb.xaxis,detx_qpb.yaxis,color='blue',label='QPB')
axbig.plot(detx_spf.xaxis,detx_spf.yaxis,color='green',label='SPF')
axbig.plot(detx_bkg.xaxis,detx_bkg.yaxis,color='red',label='Bkg')
axbig.set_yscale('log')
axbig.set_ylim(bottom=1)
axbig.set_xlabel('DETX')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_detx_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()



energy_threshold = 4000
energy_axis = []
res_count_axis = []
res_count_error = []
for entry in range(0,len(spectrum_sci.xaxis)):
    sci_entry_energy = spectrum_sci.xaxis[entry]
    sci_entry_count = spectrum_sci.yaxis[entry]
    bkg_entry_count = spectrum_bkg.yaxis[entry]
    if sci_entry_energy<energy_threshold: continue
    energy_axis += [sci_entry_energy]
    res_count_axis += [sci_entry_count-bkg_entry_count]
    res_count_error += [max(1.,pow(sci_entry_count,0.5))]
energy_axis = np.array(energy_axis)
res_count_axis = np.array(res_count_axis)
res_count_error = np.array(res_count_error)


def continuum_func(x,A,gamma):
    return A*np.exp(x/1000.*(-1.*gamma))

def line_emission_func(x,A,E0,sigma):
    return A*np.exp(-0.5*pow((x-E0)/sigma,2))

def emission_model(x,A_nontherm,gamma_nontherm,A_Fe_K_alpha,E_Fe_K_alpha,sigma):
    return continuum_func(x,A_nontherm,gamma_nontherm) + line_emission_func(x,A_Fe_K_alpha,E_Fe_K_alpha,sigma)

start = (1e3,2,200,6600,107)
limit_upper = (1e6,10,1e6,7000,117)
limit_lower = (0,0,0,6000,97)
popt, pcov = curve_fit(emission_model,energy_axis,res_count_axis,p0=start,sigma=res_count_error,absolute_sigma=True,bounds=(limit_lower,limit_upper))

model_fit = emission_model(energy_axis, *popt)

powerlaw_amp = popt[0]
powerlaw_idx = popt[1]
Fe_K_alpha_amp = popt[2]
Fe_K_alpha_E = popt[3]
E_sigma = popt[4]
print ('powerlaw_amp = %s'%(powerlaw_amp))
print ('powerlaw_idx = %s'%(powerlaw_idx))
print ('Fe_K_alpha_amp = %s'%(Fe_K_alpha_amp))
print ('Fe_K_alpha_E = %s'%(Fe_K_alpha_E))
print ('E_sigma = %s'%(E_sigma))

fig, ax = plt.subplots()
figsize_x = 8.
figsize_y = 6.
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(energy_axis,res_count_axis,yerr=res_count_error,color='k',label='Data')
axbig.plot(energy_axis,model_fit,color='red',label='Model')
axbig.set_yscale('log')
axbig.set_ylim(bottom=1)
axbig.set_xlabel('Energy [eV]')
axbig.legend(loc='best')
fig.savefig("../output_plots/fit_model.png",bbox_inches='tight')
axbig.remove()

print ('I am done.')
