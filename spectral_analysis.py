import sys
from astropy.io import fits
from astropy.table import Table
from astropy import units as my_unit
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
import numpy as np
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
#on_obsID = 'ID0412180101'
on_obsID = 'ID0400210101' # Cas A Northern lobe
#on_obsID = 'ID0782961401' # angular distance to Cas A: 34.7 arcmin

#on_sample = '3HWC_J1928_p178'
#on_obsID = 'ID0902120101'

#detector = 'mos1'
detector = 'mos2'

ana_ccd_bins = [0]
#ana_ccd_bins = [1,2,3,4,5,6,7]

exclusion_bins = 0
energy_cut = 2000
#ratio_cut = 0.75
ratio_cut = 1e10

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

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*my_unit.deg, b*my_unit.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def ConvertRaDecToGalactic(ra, dec):
    my_sky = SkyCoord(ra*my_unit.deg, dec*my_unit.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg

sky_ra_center = 0.
sky_dec_center = 0.

input_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/'+on_obsID
output_dir = '/Users/rshang/xmm_analysis/output_plots/'

image_det_mask = MyArray2D(pixel_scale=1000)
image_det_sci = MyArray2D(pixel_scale=1000)
image_det_bkg = MyArray2D(pixel_scale=1000)
for ccd in range(0,len(ana_ccd_bins)):

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    sci_hdu_list = fits.open(sci_filename)
    sci_table = Table.read(sci_filename, hdu=1)
    for entry in range(0,len(sci_table)):
        evt_pi = sci_table[entry]['PI']
        evt_detx = sci_table[entry]['DETX']
        evt_dety = sci_table[entry]['DETY']
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
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
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
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

image_det_mask.calc_ratio(image_det_sci,image_det_bkg)
ratio_dist = MyArray1D(bin_start=-1.,bin_end=3.,pixel_scale=0.1)
for idx_x in range(0,len(image_det_mask.xaxis)):
    for idx_y in range(0,len(image_det_mask.yaxis)):
        content = image_det_mask.zaxis[idx_x,idx_y]
        if content==0.: continue
        ratio_dist.fill(content)


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
        evt_l, evt_b = ConvertRaDecToGalactic(evt_ra,evt_dec)
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
        if evt_pi<energy_cut: continue
        zscore = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if zscore>ratio_cut: continue
        image_det_sci.fill(evt_detx,evt_dety)
        image_icrs_sci.fill(evt_ra,evt_dec)
        image_galactic_sci.fill(evt_l,evt_b)
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
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
        if evt_pi<energy_cut: continue
        zscore = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if zscore>ratio_cut: continue
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
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
        if evt_pi<energy_cut: continue
        zscore = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if zscore>ratio_cut: continue
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
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
        if evt_pi<energy_cut: continue
        zscore = image_det_mask.get_bin_content(evt_detx,evt_dety)
        if zscore>ratio_cut: continue
        spectrum_spf.fill(evt_pi,weight=0.1)
        detx_spf.fill(evt_detx,weight=0.1)

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
axbig.set_xlabel('DETX')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_detx_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.plot(ratio_dist.xaxis,ratio_dist.yaxis)
axbig.set_yscale('log')
axbig.set_xlabel('Ratio')
fig.savefig("%s/%s_%s_ratio_dist.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

#sci_filename = '../output_plots/sci_spectrum_sum.fits'
#bkg_filename = '../output_plots/bkg_spectrum_sum.fits'
#
#sci_hdu_list = fits.open(sci_filename)
#print (sci_filename)
#print (sci_hdu_list.info())
#
#
#sci_table = Table.read(sci_filename, hdu=1)
#bkg_table = Table.read(bkg_filename, hdu=1)
#print('sci_table.columns:')
#print(sci_table.columns)
#
#energy_threshold = 4000
#energy_axis = []
#res_count_axis = []
#res_count_error = []
#for entry in range(0,len(sci_table)):
#    sci_entry_energy = sci_table[entry]['Energy']
#    sci_entry_count = sci_table[entry]['Count']
#    bkg_entry_count = bkg_table[entry]['Count']
#    if sci_entry_energy<energy_threshold: continue
#    energy_axis += [sci_entry_energy]
#    res_count_axis += [sci_entry_count-bkg_entry_count]
#    res_count_error += [max(1.,pow(sci_entry_count,0.5))]
#    print('sci_entry_energy = %s'%(sci_entry_energy))
#    print('sci_entry_count = %s'%(sci_entry_count))
#    print('bkg_entry_count = %s'%(bkg_entry_count))
#energy_axis = np.array(energy_axis)
#res_count_axis = np.array(res_count_axis)
#res_count_error = np.array(res_count_error)
#
#
#def continuum_func(x,A,gamma):
#    #return A*pow(x/1000.,-1.*gamma)
#    return A*np.exp(x/1000.*(-1.*gamma))
#
#def line_emission_func(x,A,E0,sigma):
#    return A*np.exp(-0.5*pow((x-E0)/sigma,2))
#
#def emission_model(x,A_nontherm,gamma_nontherm,A_Fe_K_alpha,E_Fe_K_alpha,sigma):
#    return continuum_func(x,A_nontherm,gamma_nontherm) + line_emission_func(x,A_Fe_K_alpha,E_Fe_K_alpha,sigma)
#
#start = (1e3,2,200,6600,107)
#limit_upper = (1e6,10,1e6,7000,117)
#limit_lower = (0,0,0,6000,97)
#popt, pcov = curve_fit(emission_model,energy_axis,res_count_axis,p0=start,sigma=res_count_error,absolute_sigma=True,bounds=(limit_lower,limit_upper))
#
#model_fit = emission_model(energy_axis, *popt)
#
#powerlaw_amp = popt[0]
#powerlaw_idx = popt[1]
#Fe_K_alpha_amp = popt[2]
#Fe_K_alpha_E = popt[3]
#E_sigma = popt[4]
#print ('powerlaw_amp = %s'%(powerlaw_amp))
#print ('powerlaw_idx = %s'%(powerlaw_idx))
#print ('Fe_K_alpha_amp = %s'%(Fe_K_alpha_amp))
#print ('Fe_K_alpha_E = %s'%(Fe_K_alpha_E))
#print ('E_sigma = %s'%(E_sigma))
#
#fig, ax = plt.subplots()
#figsize_x = 8.
#figsize_y = 6.
#fig.set_figheight(figsize_y)
#fig.set_figwidth(figsize_x)
#
#fig.clf()
#axbig = fig.add_subplot()
#axbig.errorbar(energy_axis,res_count_axis,yerr=res_count_error,color='k',label='Data')
#axbig.plot(energy_axis,model_fit,color='red',label='Model')
#axbig.set_yscale('log')
#axbig.set_ylim(bottom=1)
#axbig.set_xlabel('Energy [eV]')
#axbig.legend(loc='best')
#fig.savefig("../output_plots/fit_model.png",bbox_inches='tight')
#axbig.remove()

print ('I am done.')
