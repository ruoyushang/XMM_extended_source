import sys
from astropy.io import fits
from astropy.table import Table
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import common_functions

#on_sample = 'extragalactic'
#on_obsID = 'ID0690900101'
on_sample = 'Cas_A'
on_obsID = 'ID0412180101'


image_sci = common_functions.MyArray2D()
image_bkg = common_functions.MyArray2D()
spectrum_sci = common_functions.MyArray1D(bin_start=common_functions.ch_low,bin_end=common_functions.ch_high,pixel_scale=common_functions.ch_scale)
spectrum_bkg = common_functions.MyArray1D(bin_start=common_functions.ch_low,bin_end=common_functions.ch_high,pixel_scale=common_functions.ch_scale)

ana_ccd_bins = [1,2,3,4,5,6,7]

exclusion_bins = 5000

output_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/'+on_obsID
for ccd in range(0,len(ana_ccd_bins)):

    sci_filename = '%s/sci_events_ccd%s.fits'%(output_dir,ana_ccd_bins[ccd])
    sci_hdu_list = fits.open(sci_filename)
    sci_table = Table.read(sci_filename, hdu=1)
    for entry in range(0,len(sci_table)):
        evt_pi = sci_table[entry]['PI']
        evt_detx = sci_table[entry]['DETX']
        evt_dety = sci_table[entry]['DETY']
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
        image_sci.fill(evt_detx,evt_dety)
        spectrum_sci.fill(evt_pi)

    bkg_filename = '%s/bkg_events_ccd%s.fits'%(output_dir,ana_ccd_bins[ccd])
    bkg_hdu_list = fits.open(bkg_filename)
    bkg_table = Table.read(bkg_filename, hdu=1)
    for entry in range(0,len(bkg_table)):
        evt_pi = bkg_table[entry]['PI']
        evt_detx = bkg_table[entry]['DETX']
        evt_dety = bkg_table[entry]['DETY']
        if abs(evt_detx)<exclusion_bins and abs(evt_dety)<exclusion_bins: continue
        image_bkg.fill(evt_detx,evt_dety,weight=0.1)
        spectrum_bkg.fill(evt_pi,weight=0.1)

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
xmin = image_sci.xaxis.min()
xmax = image_sci.xaxis.max()
ymin = image_sci.yaxis.min()
ymax = image_sci.yaxis.max()
axbig.imshow(image_sci.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/image_sci.png"%(output_dir),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_sci.xaxis.min()
xmax = image_sci.xaxis.max()
ymin = image_sci.yaxis.min()
ymax = image_sci.yaxis.max()
axbig.imshow(image_bkg.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/image_bkg.png"%(output_dir),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(spectrum_sci.xaxis,spectrum_sci.yaxis,yerr=spectrum_sci.yerr,color='k',label='Data')
axbig.plot(spectrum_bkg.xaxis,spectrum_bkg.yaxis,color='red',label='Bkg')
axbig.set_yscale('log')
#axbig.set_ylim(bottom=1)
axbig.set_xlabel('Energy [eV]')
axbig.legend(loc='best')
fig.savefig("%s/spectrum_model.png"%(output_dir),bbox_inches='tight')
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
