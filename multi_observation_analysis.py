
import sys
import shutil
from astropy.io import fits
from astropy.table import Table
from astropy import units as my_unit
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
import numpy as np
from numpy.linalg import inv
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from scipy import fftpack
import matplotlib.pyplot as plt
from matplotlib import colors
import pickle

import common_functions

on_sample = sys.argv[1]
on_obsID = sys.argv[2]
job = sys.argv[3]
ana_tag = sys.argv[4]


#measure_cxb = True
measure_cxb = False

#bkg_only = True
bkg_only = False

include_cxb = True
if measure_cxb:
    include_cxb = False

# all events
find_extended_src = False
select_mask_events = False

# select events in RoI
#find_extended_src = True
#select_mask_events = True

# background study
if bkg_only:
    find_extended_src = True
    select_mask_events = False
if measure_cxb:
    find_extended_src = True
    select_mask_events = False

show_log_spectrum = True

do_energy_flux = True # use this to enable acceptance correction
write_xspec_output = False
my_spectrum_unit = 'erg/cm2/s/sr'
#do_energy_flux = False # use this if you need xspec output
#write_xspec_output = True
#my_spectrum_unit = 'Photons /cm2/s/sr/keV'

ana_ccd_bins = [0]
#ana_ccd_bins = [1,2,3,4,5,6,7]

energy_cut_lower = 500
#energy_cut_lower = 1700
energy_cut_upper = 8000
#energy_cut_upper = 12000

on_exposure = 0.
total_spectral_volume = 0.

#energy_array = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000]
#delta_energy = energy_array[1]-energy_array[0]
#cxb_delta_energy = delta_energy/4.
#cxb_energy_array = []
#cxb_energy = energy_array[0]
#while cxb_energy<=energy_array[len(energy_array)-1]:
#    cxb_energy_array += [cxb_energy]
#    cxb_energy += cxb_delta_energy
#print (f'cxb_energy_array = {cxb_energy_array}')

sample_scale = 10.

map_rebin = 2
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
dety_low = common_functions.dety_low
dety_high = common_functions.dety_high
dety_scale = common_functions.dety_scale
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

fig, ax = plt.subplots()
figsize_x = 8.
figsize_y = 6.
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)
map_color = 'coolwarm'
color_list = ['salmon','yellowgreen','deepskyblue']

on_run_list = on_obsID.split('+')
len_run_list = float(len(on_run_list))

def ConvertSky2Det(target_sky,origin_det,origin_sky,mtx_conv_sky_2_det):

    print (f'target_sky = {target_sky}')
    print (f'origin_sky = {origin_sky}')

    target_sky = np.array(target_sky)
    origin_sky = np.array(origin_sky)
    origin_det = np.array(origin_det)

    delta_sky = target_sky-origin_sky
    delta_det = np.matmul(np.transpose(mtx_conv_sky_2_det),delta_sky)
    target_det = origin_det
    target_det[0] += 1.*delta_det[0,0]
    target_det[1] += 1.*delta_det[0,1]

    return target_det[0], target_det[1]

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*my_unit.deg, b*my_unit.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def ConvertRaDecToGalactic(ra, dec):
    my_sky = SkyCoord(ra*my_unit.deg, dec*my_unit.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg

def smooth_map(image_data,image_smooth,mode):

    bin_size = image_data.xaxis[1]-image_data.xaxis[0]
    print ('smoothing map...')
    kernel_radius = (detx_high-detx_low)/mode
    print (f'kernel_radius = {kernel_radius}')
    kernel_pix_size = int(kernel_radius/bin_size)
    print (f'kernel_pix_size = {kernel_pix_size}')
    gaus_norm = 2.*np.pi*kernel_radius*kernel_radius
    image_kernel = MyArray2D(pixel_scale=bin_size)
    central_bin_x = int(len(image_kernel.xaxis)/2)
    central_bin_y = int(len(image_kernel.yaxis)/2)
    for idx_x in range(0,len(image_kernel.xaxis)):
        for idx_y in range(0,len(image_kernel.yaxis)):
            pix_x = image_kernel.xaxis[idx_x]
            pix_y = image_kernel.yaxis[idx_y]
            distance = pow(pix_x*pix_x+pix_y*pix_y,0.5)
            pix_content = np.exp(-(distance*distance)/(2.*kernel_radius*kernel_radius))
            image_kernel.zaxis[idx_x,idx_y] = pix_content/gaus_norm

    for idx_x1 in range(0,len(image_data.xaxis)):
        for idx_y1 in range(0,len(image_data.yaxis)):
            image_smooth.zaxis[idx_x1,idx_y1] = 0.
            for idx_x2 in range(idx_x1-3*kernel_pix_size,idx_x1+3*kernel_pix_size):
                for idx_y2 in range(idx_y1-3*kernel_pix_size,idx_y1+3*kernel_pix_size):
                    if idx_x2<0: continue
                    if idx_y2<0: continue
                    if idx_x2>=len(image_data.xaxis): continue
                    if idx_y2>=len(image_data.yaxis): continue
                    old_content = image_data.zaxis[idx_x2,idx_y2]
                    scale = image_kernel.zaxis[central_bin_x+idx_x2-idx_x1,central_bin_y+idx_y2-idx_y1] 
                    image_smooth.zaxis[idx_x1,idx_y1] += old_content*scale

    print ('done smoothing map.')

def MakeRadialProjection(image_map_xry, profile_r):

    for idx_z in range(0,len(profile_r.xaxis)):
        profile_r.yaxis[idx_z] = 0.

    delta_pix_z = profile_r.xaxis[1]-profile_r.xaxis[0]
    for idx_z in range(0,len(profile_r.xaxis)):
        pix_z_lo = profile_r.xaxis[idx_z]
        pix_z_hi = profile_r.xaxis[idx_z]+delta_pix_z
        pix_integral_xry = 0.
        pix_count = 0.
        for idx_x in range(0,len(image_map_xry.xaxis)):
            delta_x = image_map_xry.xaxis[idx_x] - mask_ra
            for idx_y in range(0,len(image_map_xry.yaxis)):
                delta_y = image_map_xry.yaxis[idx_y] - mask_dec
                distance = pow(delta_x*delta_x+delta_y*delta_y,0.5)
                if distance<pix_z_lo: continue
                if distance>=pix_z_hi: continue
                pix_content_xry = image_map_xry.zaxis[idx_x,idx_y]
                pix_integral_xry += pix_content_xry
                pix_count += 1.
        if pix_count>0.:
            profile_r.yaxis[idx_z] = pix_integral_xry/pix_count

def calculate_pix_gravity(image_data,src_x1,src_y1,src_x2,src_y2):

    min_distance = image_data.xaxis[1]-image_data.xaxis[0]
    src_mass1 = image_data.get_bin_content(src_x1,src_y1)
    src_mass2 = image_data.get_bin_content(src_x2,src_y2)
    distance = pow(pow(src_x1-src_x2,2)+pow(src_y1-src_y2,2),0.5)
    distance = max(min_distance,distance)
    if src_mass1==0.:
        return 0.
    elif src_mass2==0.:
        return 0.
    else:
        gravity = (src_mass1*src_mass2)/distance
        return gravity

def image_cleaning_with_svd(image_data,image_svd_cleaned):

    k = 5
    U_signal, S_signal, V_signal = np.linalg.svd(image_data.zaxis,full_matrices=False)
    image_svd_cleaned.zaxis = U_signal[:, :k] @ np.diag(S_signal[:k]) @ V_signal[:k, :]

def fft_filter(image_mask,image_data,image_fft_filter,lofreq_cut,hifreq_cut):

    rows, cols = image_data.zaxis.shape

    image_data_mask = MyArray2D()
    image_data_mask.zaxis = image_data.zaxis
    for idx_x in range(0,len(image_data.xaxis)):
        for idx_y in range(0,len(image_data.yaxis)):
            mask = image_mask.zaxis[idx_x,idx_y]
            if mask==1: 
                image_data_mask.zaxis[idx_x,idx_y] = 0.

    # Compute 2D Fourier transform
    F_data = fftpack.fft2(image_data_mask.zaxis)

    # Shift the zero-frequency component to the center
    F_data_shifted_lo = fftpack.fftshift(F_data)
    F_data_shifted_hi = fftpack.fftshift(F_data)

    # Filter out high-frequency components
    crow, ccol = rows // 2, cols // 2  # Center coordinates
    n_lo = lofreq_cut  # Number of lowest-frequency components to remove
    n_hi = hifreq_cut  # Number of highest-frequency components to remove
    if n_hi==0:
        n_hi = crow
    F_data_shifted_lo[crow-n_lo:crow+n_lo+1, ccol-n_lo:ccol+n_lo+1] = 0.
    F_data_shifted_hi[crow-n_hi:crow+n_hi+1, ccol-n_hi:ccol+n_hi+1] = 0.

    # Shift the zero-frequency component back to the corner
    F_filtered = fftpack.ifftshift(F_data_shifted_lo-F_data_shifted_hi)

    # Compute the inverse Fourier transform to get the filtered image
    image_fft_filter.zaxis = np.real(fftpack.ifft2(F_filtered))

def find_point_sources_with_gravity(mode,image_data,image_mask,threshold=5.0):

    print ('calculating pixel gravity...')
    kernel_radius = (detx_high-detx_low)/mode
    bin_size = image_data.xaxis[1]-image_data.xaxis[0]
    kernel_pix_size = int(kernel_radius/bin_size)

    all_pix_gravity = []
    for idx_x1 in range(0,len(image_data.xaxis)):
        for idx_y1 in range(0,len(image_data.yaxis)):
            mask = image_mask.zaxis[idx_x1,idx_y1]
            if mask==1: continue
            content = image_data.zaxis[idx_x1,idx_y1]
            if content==0: continue
            sum_gravity = 0.
            for idx_x2 in range(idx_x1-3*kernel_pix_size,idx_x1+3*kernel_pix_size):
                for idx_y2 in range(idx_y1-3*kernel_pix_size,idx_y1+3*kernel_pix_size):
                    if idx_x2<0: continue
                    if idx_y2<0: continue
                    if idx_x2>=len(image_data.xaxis): continue
                    if idx_y2>=len(image_data.yaxis): continue
                    pix_x1 = image_data.xaxis[idx_x1]
                    pix_y1 = image_data.yaxis[idx_y1]
                    pix_x2 = image_data.xaxis[idx_x2]
                    pix_y2 = image_data.yaxis[idx_y2]
                    sum_gravity += calculate_pix_gravity(image_data,pix_x1,pix_y1,pix_x2,pix_y2)
            all_pix_gravity += [sum_gravity]


    mean_gravity = np.mean(all_pix_gravity)
    rms_gravity = 0.
    for entry in range(0,len(all_pix_gravity)):
        rms_gravity += pow(all_pix_gravity[entry]-mean_gravity,2)
    rms_gravity = pow(rms_gravity/len(all_pix_gravity),0.5)
    print (f'mean_gravity = {mean_gravity}')
    print (f'rms_gravity = {rms_gravity}')

    for idx_x1 in range(0,len(image_data.xaxis)):
        for idx_y1 in range(0,len(image_data.yaxis)):
            mask = image_mask.zaxis[idx_x1,idx_y1]
            if mask==1: continue
            content = image_data.zaxis[idx_x1,idx_y1]
            if content==0: continue
            sum_gravity = 0.
            for idx_x2 in range(idx_x1-3*kernel_pix_size,idx_x1+3*kernel_pix_size):
                for idx_y2 in range(idx_y1-3*kernel_pix_size,idx_y1+3*kernel_pix_size):
                    if idx_x2<0: continue
                    if idx_y2<0: continue
                    if idx_x2>=len(image_data.xaxis): continue
                    if idx_y2>=len(image_data.yaxis): continue
                    pix_x1 = image_data.xaxis[idx_x1]
                    pix_y1 = image_data.yaxis[idx_y1]
                    pix_x2 = image_data.xaxis[idx_x2]
                    pix_y2 = image_data.yaxis[idx_y2]
                    sum_gravity += calculate_pix_gravity(image_data,pix_x1,pix_y1,pix_x2,pix_y2)
            if sum_gravity>threshold*rms_gravity+mean_gravity:
                image_mask.zaxis[idx_x1,idx_y1] = 1

    max_gravity = max(all_pix_gravity)
    hist_gravity = MyArray1D(bin_start=mean_gravity-2.*rms_gravity,bin_end=mean_gravity+4.*rms_gravity,pixel_scale=6.*rms_gravity/100.)
    for entry in range(0,len(all_pix_gravity)):
        hist_gravity.fill(all_pix_gravity[entry])

    return hist_gravity

def draw_stacked_histogram(fig,hist_data,hist_bkg,bkg_colors,bkg_labels,xlabel,ylabel,plot_name,plot_name_ul,log_scale=False):

    n_bins = len(hist_data.xaxis)
    stack_bkg = []
    for h1 in range(0,len(hist_bkg)):
        stack_bkg += [MyArray1D(bin_start=hist_data.xaxis[0],bin_end=hist_data.xaxis[n_bins-1],pixel_scale=hist_data.xaxis[1]-hist_data.xaxis[0])]
        for h2 in range(h1,len(hist_bkg)):
            stack_bkg[h1].add(hist_bkg[h2])
    
    min_data = 1e10
    for x in range(0,len(hist_data.xaxis)):
        if hist_data.yaxis[x]<=0.: continue
        if hist_data.yaxis[x]<min_data:
            min_data = hist_data.yaxis[x]

    avg_qpb = 0.
    bin_count = 0.
    if log_scale:
        for x in range(0,len(hist_bkg[len(hist_bkg)-1].xaxis)):
            energy = hist_bkg[len(hist_bkg)-1].xaxis[x]
            if 'spectrum' in plot_name:
                if energy>2000:
                    avg_qpb += hist_bkg[len(hist_bkg)-1].yaxis[x]
                    bin_count += 1.
            else:
                avg_qpb += hist_bkg[len(hist_bkg)-1].yaxis[x]
                bin_count += 1.
        avg_qpb = avg_qpb/bin_count

    fig.clf()
    axbig = fig.add_subplot()
    for h in range(0,len(hist_bkg)):
        axbig.fill_between(stack_bkg[h].xaxis,0.,stack_bkg[h].yaxis,color=bkg_colors[h],label=bkg_labels[h])
    axbig.fill_between(stack_bkg[0].xaxis,stack_bkg[0].yaxis-stack_bkg[0].yerr,stack_bkg[0].yaxis+stack_bkg[0].yerr,color='k',alpha=0.2)
    axbig.errorbar(hist_data.xaxis,hist_data.yaxis,yerr=hist_data.yerr,color='k',label='Data')
    if log_scale:
        axbig.set_yscale('log')
        axbig.set_ylim(bottom=0.5*min_data)
    if 'spectrum' in plot_name:
        axbig.set_xscale('log')
    axbig.set_xlabel(xlabel)
    axbig.set_ylabel(ylabel)
    axbig.legend(loc='best')
    fig.savefig("%s"%(plot_name),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    axbig = fig.add_subplot()
    axbig.plot(stack_bkg[0].xaxis,3.*stack_bkg[0].yerr,color='k')
    if log_scale:
        axbig.set_yscale('log')
    if 'spectrum' in plot_name:
        axbig.set_xscale('log')
    axbig.set_xlabel(xlabel)
    axbig.set_ylabel(ylabel)
    fig.savefig("%s"%(plot_name_ul),bbox_inches='tight')
    axbig.remove()

def get_cxb_spectrum(cxb_output_dir,hist_cxb,obs_sky_l,obs_sky_b,detector,obsID=''):

    cxb_file_name = "%s/cxb_spectrum_%s.txt"%(cxb_output_dir,on_sample)
    print (f'read {cxb_file_name}')
    cxb_file = open(cxb_file_name)
    cxb_measurement = []
    cxb_measurement_weight = []
    use_this_data = True
    distance = 0.
    exposure = 0.
    qpb_frac = 0.
    for line in cxb_file:
        if '#' in line:
            line_split = line.split(';')
            line_obsid = line_split[1]
            line_detector = line_split[2]
            line_spf_frac = line_split[4].strip(' SPF')
            qpb_frac = 1./(float(line_spf_frac)+1.)
            line_exposure = line_split[5].strip('\n')
            exposure = float(line_exposure)
            line_coord = line_split[3]
            line_coord = line_coord.strip('( )')
            distance = pow(obs_sky_l-float(line_coord[0]),2)+pow(obs_sky_b-float(line_coord[1]),2)
            distance = pow(distance,0.5)
            if not obsID=='':
                if not obsID in line_obsid:
                    use_this_data = False
            else:
                if distance<0.3:
                    use_this_data = False
            #if not detector in line_detector:
            #    use_this_data = False
            continue
        if not use_this_data: 
            use_this_data = True
            continue
        line_strip = line.strip('\n')
        line_split = line_strip.split(' ')
        cxb_measurement_single = []
        for entry in range(0,len(line_split)):
            cxb_measurement_single += [float(line_split[entry])]
        cxb_measurement += [cxb_measurement_single]
        #cxb_measurement_weight += [1./distance]
        cxb_measurement_weight += [pow(exposure,0.5)*qpb_frac]
    cxb_file.close()

    if not obsID=='':
        print (f'Load data for obsID = {obsID}, detector = {detector}')
    print (f'{len(cxb_measurement)} CXB measurements are used.')

    hist_cxb_measurement = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    for ch in range(0,len(hist_cxb.xaxis)):
        n_samples = 0.
        avg_cxb = 0.
        total_weight = 0.
        for m in range(0,len(cxb_measurement)):
            n_samples += 1.
            total_weight += cxb_measurement_weight[m]
            avg_cxb += cxb_measurement[m][ch]*cxb_measurement_weight[m]
        avg_cxb = avg_cxb/total_weight
        hist_cxb_measurement.yaxis[ch] = avg_cxb
    for ch in range(0,len(hist_cxb.xaxis)):
        n_samples = 0.
        avg_cxb = hist_cxb_measurement.yaxis[ch]
        rms_cxb = 0.
        total_weight = 0.
        for m in range(0,len(cxb_measurement)):
            n_samples += 1.
            total_weight += cxb_measurement_weight[m]
            rms_cxb += pow(cxb_measurement[m][ch]-avg_cxb,2)*cxb_measurement_weight[m]
        rms_cxb = pow(rms_cxb/total_weight,0.5)
        hist_cxb_measurement.yerr[ch] = rms_cxb

    for x in range(0,len(hist_cxb.xaxis)):
        energy = hist_cxb.xaxis[x]
        if energy<energy_cut_lower: continue
        if energy>energy_cut_upper: continue
        hist_cxb.yaxis[x] = hist_cxb_measurement.get_bin_content(energy)
        hist_cxb.yerr[x] = hist_cxb_measurement.get_bin_error(energy)

def analyze_one_observation(obsID,detector):

    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print (f'analyze obsID {obsID} detector {detector}')


    input_dir = '/Users/rshang/xmm_analysis/output_extended_analysis/'+on_sample+'/'+ana_tag+'/ID'+obsID
    cxb_output_dir = '/Users/rshang/xmm_analysis/output_extended_analysis/'+on_sample+'/'+ana_tag
    output_dir = '/Users/rshang/xmm_analysis/output_plots/plot_%s/'%(on_sample)

    print (f'input_dir = {input_dir}')

    image_det_mask = MyArray2D(pixel_scale=map_rebin*detx_scale)
    image_det_init_sci = MyArray2D(pixel_scale=map_rebin*detx_scale)
    image_det_init_bkg = MyArray2D(pixel_scale=map_rebin*detx_scale)
    image_det_diff = MyArray2D(pixel_scale=map_rebin*detx_scale)
    image_det_diff_lofreq = MyArray2D(pixel_scale=map_rebin*detx_scale)

    for ccd in range(0,len(ana_ccd_bins)):
    
        sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        print (f'sci_filename = {sci_filename}')
        sci_hdu_list = fits.open(sci_filename)
        sci_table = Table.read(sci_filename, hdu=1)
        for entry in range(0,len(sci_table)):
            evt_pi = sci_table[entry]['PI']
            evt_detx = sci_table[entry]['DETX']
            evt_dety = sci_table[entry]['DETY']
            if evt_pi<energy_cut_lower: continue
            if evt_pi>energy_cut_upper: continue
            image_det_init_sci.fill(evt_detx,evt_dety)
    
        bkg_filename = '%s/qpb_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        bkg_hdu_list = fits.open(bkg_filename)
        bkg_table = Table.read(bkg_filename, hdu=1)
        for entry in range(0,len(bkg_table)):
            evt_pi = bkg_table[entry]['PI']
            evt_detx = bkg_table[entry]['DETX']
            evt_dety = bkg_table[entry]['DETY']
            if evt_pi<energy_cut_lower: continue
            if evt_pi>energy_cut_upper: continue
            image_det_init_bkg.fill(evt_detx,evt_dety,weight=1./sample_scale)
        bkg_filename = '%s/spf_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        bkg_hdu_list = fits.open(bkg_filename)
        bkg_table = Table.read(bkg_filename, hdu=1)
        for entry in range(0,len(bkg_table)):
            evt_pi = bkg_table[entry]['PI']
            evt_detx = bkg_table[entry]['DETX']
            evt_dety = bkg_table[entry]['DETY']
            if evt_pi<energy_cut_lower: continue
            if evt_pi>energy_cut_upper: continue
            image_det_init_bkg.fill(evt_detx,evt_dety,weight=1./sample_scale)

    fft_lofreq_mode = 100
    image_det_diff.zaxis = image_det_init_sci.zaxis - image_det_init_bkg.zaxis
    smooth_map(image_det_diff,image_det_diff_lofreq,fft_lofreq_mode)

    fig.clf()
    axbig = fig.add_subplot()
    label_x = 'DETY'
    label_y = 'DETX'
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    xmin = image_det_diff_lofreq.xaxis.min()
    xmax = image_det_diff_lofreq.xaxis.max()
    ymin = image_det_diff_lofreq.yaxis.min()
    ymax = image_det_diff_lofreq.yaxis.max()
    axbig.imshow(image_det_diff_lofreq.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
    fig.savefig("%s/image_det_diff_lofreq_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
    axbig.remove()

    if find_extended_src:
        hist_gravity = find_point_sources_with_gravity(fft_lofreq_mode,image_det_diff_lofreq,image_det_mask)
        hist_gravity = find_point_sources_with_gravity(fft_lofreq_mode,image_det_diff_lofreq,image_det_mask)
        hist_gravity = find_point_sources_with_gravity(fft_lofreq_mode,image_det_diff_lofreq,image_det_mask)
        fig.clf()
        axbig = fig.add_subplot()
        axbig.plot(hist_gravity.xaxis,hist_gravity.yaxis)
        axbig.set_yscale('log')
        axbig.set_xlabel('Gravity')
        axbig.legend(loc='best')
        fig.savefig("%s/gravity_lofreq_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
        axbig.remove()

    if abs(mask_detx)<20000 and abs(mask_dety)<20000: 
        for idx_x in range(0,len(image_det_mask.xaxis)):
            for idx_y in range(0,len(image_det_mask.yaxis)):
                pix_detx = image_det_mask.xaxis[idx_x]
                pix_dety = image_det_mask.yaxis[idx_y]
                distance = pow(pow(pix_detx-mask_detx,2)+pow(pix_dety-mask_dety,2),0.5)
                if distance<mask_inner_detr:
                    image_det_mask.zaxis[idx_x,idx_y] = 1.
                if distance>mask_outer_detr:
                    image_det_mask.zaxis[idx_x,idx_y] = 1.

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
    fig.savefig("%s/image_det_mask_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
    axbig.remove()


    image_det_all = MyArray2D(pixel_scale=detx_scale)
    image_det_sci = MyArray2D(pixel_scale=detx_scale)
    image_det_qpb = MyArray2D(pixel_scale=detx_scale)
    image_det_spf = MyArray2D(pixel_scale=detx_scale)
    image_det_cxb = MyArray2D(pixel_scale=detx_scale)
    image_det_xry = MyArray2D(pixel_scale=detx_scale)
    spectrum_sci = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    spectrum_qpb = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    spectrum_qpb_radial = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    spectrum_spf = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    spectrum_spf_radial = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    spectrum_cxb = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    spectrum_xry = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
    image_icrs_sci = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_xry = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_qpb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_spf = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_cxb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    radial_sci = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
    radial_qpb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
    radial_spf = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
    radial_cxb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))

    area_curve = []
    area_curve += [MyArray1D(bin_start=2000,bin_end=12000,pixel_scale=2000)]
    area_curve += [MyArray1D(bin_start=2000,bin_end=12000,pixel_scale=2000)]
    area_curve += [MyArray1D(bin_start=2000,bin_end=12000,pixel_scale=2000)]
    area_curve += [MyArray1D(bin_start=2000,bin_end=12000,pixel_scale=2000)]
    lightcurve_all = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)
    lightcurve_sci = MyArray1D(bin_start=t_low,bin_end=t_high,pixel_scale=t_scale)

    fov_size = 0.
    for idx_x in range(0,len(image_det_mask.xaxis)):
        for idx_y in range(0,len(image_det_mask.yaxis)):
            evt_cnt = image_det_init_sci.zaxis[idx_x,idx_y]
            mask = image_det_mask.zaxis[idx_x,idx_y]
            if evt_cnt==0.: continue
            if not select_mask_events:
                if mask==1: continue
            else:
                if mask!=1: continue
            fov_size += pow(map_rebin*detx_scale*0.05/(60.*60.),2)/3282.8 # steradian
    print (f'fov_size = {fov_size}')
    pix_size = pow(detx_scale*0.05/(60.*60.),2)/3282.8

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[0])
    sci_hdu_list = fits.open(sci_filename)
    sky_ra_center = sci_hdu_list[1].header['REF_RA']
    sky_dec_center = sci_hdu_list[1].header['REF_DEC']
    sky_l_center, sky_b_center = ConvertRaDecToGalactic(sky_ra_center, sky_dec_center)

    aeff_r = -1
    for r in range(0,len(area_curve)):
        arf_filename = '%s/sci_area_r%s_%s_ccd%s.fits'%(input_dir,r,detector,ana_ccd_bins[0])
        arf_hdu_list = fits.open(arf_filename)
        arf_table = Table.read(arf_filename, hdu=1)
        for entry in range(0,len(arf_table)):
            entry_energy = arf_table[entry]['ENERGY']
            entry_area = arf_table[entry]['AREA']
            print (f'energy = {entry_energy}, area = {entry_area}')
            area_curve[r].fill(entry_energy,weight=entry_area)
        print (f'area_curve[{r}] = {area_curve[r].yaxis[0]}')
        if area_curve[r].yaxis[0]>0. and aeff_r==-1:
            aeff_r = r

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[0])
    sci_hdu_list = fits.open(sci_filename)
    exposure = sci_hdu_list[1].header['EXPOSURE']
    keV_to_erg = 1.60218e-9
    spectral_volume = ch_scale/1000.*fov_size*exposure*area_curve[aeff_r].yaxis[0]*keV_to_erg
    spacial_volume = (energy_cut_upper-energy_cut_lower)/1000.*pix_size*exposure*area_curve[aeff_r].yaxis[0]*keV_to_erg
    print (f'area = {area_curve[aeff_r].yaxis[0]}')
    print (f'exposure = {exposure}')
    print (f'spectral_volume = {spectral_volume}')

    for ccd in range(0,len(ana_ccd_bins)):
    
        lc_filename = '%s/sci_lightcurve_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        lc_hdu_list = fits.open(lc_filename)
        lc_table = Table.read(lc_filename, hdu=1)
        for entry in range(0,len(lc_table)):
            entry_time = lc_table[entry]['TIME']
            entry_allevt = lc_table[entry]['ALL_EVT']
            entry_scievt = lc_table[entry]['SCI_EVT']
            lightcurve_all.fill(entry_time,weight=entry_allevt)
            lightcurve_sci.fill(entry_time,weight=entry_scievt)

        sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        sci_hdu_list = fits.open(sci_filename)
        sci_table = Table.read(sci_filename, hdu=1)
        for entry in range(0,len(sci_table)):
            evt_pi = sci_table[entry]['PI']
            evt_detx = sci_table[entry]['DETX']
            evt_dety = sci_table[entry]['DETY']
            evt_ra = sci_table[entry]['RA']
            evt_dec = sci_table[entry]['DEC']
            if evt_pi<energy_cut_lower: continue
            if evt_pi>energy_cut_upper: continue
            radius = int(pow(evt_detx*evt_detx+evt_dety*evt_dety,0.5)/5000.)
            convert_to_energy_flux = 1.
            unfolding_acceptance = 1.
            if do_energy_flux: 
                convert_to_energy_flux = pow(evt_pi/1000.,2)*keV_to_erg*keV_to_erg
                unfolding_acceptance = area_curve[aeff_r].yaxis[0]/area_curve[radius].get_bin_content(evt_pi)
            image_det_all.fill(evt_detx,evt_dety,weight=1./spacial_volume*convert_to_energy_flux)
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if not select_mask_events:
                if mask==1: continue
            else:
                if mask!=1: continue
            delta_ra = evt_ra - sky_ra_center
            delta_dec = evt_dec - sky_dec_center
            spectrum_sci.fill(evt_pi,weight=1./spectral_volume*convert_to_energy_flux*unfolding_acceptance)
            spectrum_comb_sci.fill(evt_pi,weight=1./spectral_volume*1./len_run_list*convert_to_energy_flux*unfolding_acceptance)
            if write_xspec_output:
                spectrum_comb_raw_sci.fill(evt_pi)
            image_icrs_sci.fill(evt_ra,evt_dec,weight=1./spacial_volume*convert_to_energy_flux*unfolding_acceptance)
            image_icrs_comb_sci.fill(evt_ra,evt_dec,weight=1./spacial_volume*1./len_run_list*convert_to_energy_flux*unfolding_acceptance)
            image_det_sci.fill(evt_detx,evt_dety,weight=1./spacial_volume*convert_to_energy_flux*unfolding_acceptance)

        qpb_filename = '%s/qpb_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        qpb_hdu_list = fits.open(qpb_filename)
        qpb_table = Table.read(qpb_filename, hdu=1)
        for entry in range(0,len(qpb_table)):
            evt_pi = qpb_table[entry]['PI']
            evt_detx = qpb_table[entry]['DETX']
            evt_dety = qpb_table[entry]['DETY']
            evt_ra = qpb_table[entry]['RA']
            evt_dec = qpb_table[entry]['DEC']
            if evt_pi<energy_cut_lower: continue
            if evt_pi>energy_cut_upper: continue
            radius = int(pow(evt_detx*evt_detx+evt_dety*evt_dety,0.5)/5000.)
            convert_to_energy_flux = 1.
            unfolding_acceptance = 1.
            if do_energy_flux: 
                convert_to_energy_flux = pow(evt_pi/1000.,2)*keV_to_erg*keV_to_erg
                unfolding_acceptance = area_curve[aeff_r].yaxis[0]/area_curve[radius].get_bin_content(evt_pi)
            radial_acceptance = area_curve[radius].get_bin_content(evt_pi)/area_curve[aeff_r].yaxis[0]
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if not select_mask_events:
                if mask==1: continue
            else:
                if mask!=1: continue
            delta_ra = evt_ra - sky_ra_center
            delta_dec = evt_dec - sky_dec_center
            spectrum_qpb.fill(evt_pi,weight=1./spectral_volume*1./sample_scale*convert_to_energy_flux*unfolding_acceptance)
            spectrum_comb_qpb.fill(evt_pi,weight=1./spectral_volume*1./sample_scale*1./len_run_list*convert_to_energy_flux*unfolding_acceptance)
            if write_xspec_output:
                spectrum_comb_raw_qpb.fill(evt_pi,weight=1./sample_scale)
            image_icrs_qpb.fill(evt_ra,evt_dec,weight=1./spacial_volume*1./sample_scale*convert_to_energy_flux*unfolding_acceptance)
            image_icrs_comb_qpb.fill(evt_ra,evt_dec,weight=1./spacial_volume*1./sample_scale*1./len_run_list*convert_to_energy_flux*unfolding_acceptance)
            image_icrs_cxb.fill(evt_ra,evt_dec,weight=1./spacial_volume*1./sample_scale*1./len_run_list*radial_acceptance*convert_to_energy_flux*unfolding_acceptance)
            spectrum_qpb_radial.fill(evt_pi,weight=1./spectral_volume*1./sample_scale*radial_acceptance*convert_to_energy_flux*unfolding_acceptance)
            image_det_qpb.fill(evt_detx,evt_dety,weight=1./spacial_volume*1./sample_scale*convert_to_energy_flux*unfolding_acceptance)
            image_det_cxb.fill(evt_detx,evt_dety,weight=1./spacial_volume*1./sample_scale*1./len_run_list*radial_acceptance*convert_to_energy_flux*unfolding_acceptance)

        spf_filename = '%s/spf_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        spf_hdu_list = fits.open(spf_filename)
        spf_table = Table.read(spf_filename, hdu=1)
        for entry in range(0,len(spf_table)):
            evt_pi = spf_table[entry]['PI']
            evt_detx = spf_table[entry]['DETX']
            evt_dety = spf_table[entry]['DETY']
            evt_ra = spf_table[entry]['RA']
            evt_dec = spf_table[entry]['DEC']
            if evt_pi<energy_cut_lower: continue
            if evt_pi>energy_cut_upper: continue
            radius = int(pow(evt_detx*evt_detx+evt_dety*evt_dety,0.5)/5000.)
            convert_to_energy_flux = 1.
            unfolding_acceptance = 1.
            if do_energy_flux: 
                convert_to_energy_flux = pow(evt_pi/1000.,2)*keV_to_erg*keV_to_erg
                unfolding_acceptance = area_curve[aeff_r].yaxis[0]/area_curve[radius].get_bin_content(evt_pi)
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if not select_mask_events:
                if mask==1: continue
            else:
                if mask!=1: continue
            delta_ra = evt_ra - sky_ra_center
            delta_dec = evt_dec - sky_dec_center
            spectrum_spf.fill(evt_pi,weight=1./spectral_volume*1./sample_scale*convert_to_energy_flux*unfolding_acceptance)
            spectrum_spf_radial.fill(evt_pi,weight=1./spectral_volume*1./sample_scale*convert_to_energy_flux*unfolding_acceptance)
            spectrum_comb_spf.fill(evt_pi,weight=1./spectral_volume*1./sample_scale*1./len_run_list*convert_to_energy_flux*unfolding_acceptance)
            if write_xspec_output:
                spectrum_comb_raw_spf.fill(evt_pi,weight=1./sample_scale)
            image_icrs_spf.fill(evt_ra,evt_dec,weight=1./spacial_volume*1./sample_scale*convert_to_energy_flux*unfolding_acceptance)
            image_det_spf.fill(evt_detx,evt_dety,weight=1./spacial_volume*1./sample_scale*convert_to_energy_flux*unfolding_acceptance)

    spf_cnt = spectrum_spf.integral()
    spf_radial_cnt = spectrum_spf_radial.integral()
    spf_radial_norm = 0.
    if spf_radial_cnt>0.:
        spf_radial_norm = spf_cnt/spf_radial_cnt
    image_det_spf.scale(spf_radial_norm)
    image_icrs_spf.scale(spf_radial_norm)
    image_icrs_comb_spf.add(image_icrs_spf,factor=1./len_run_list)

    qpb_cnt = spectrum_qpb.integral()
    spf_cnt = spectrum_spf.integral()
    spf_frac = spf_cnt/qpb_cnt
    print (f'spf_cnt/qpb_cnt = {spf_cnt/qpb_cnt}')

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[0])
    sci_hdu_list = fits.open(sci_filename)
    exposure = sci_hdu_list[1].header['EXPOSURE']

    cxb_spectrum_measurement_tmp = []
    for ch in range(0,len(spectrum_cxb.xaxis)):
        sci_cnt = spectrum_sci.yaxis[ch]
        spf_cnt = spectrum_spf.yaxis[ch]
        qpb_cnt = spectrum_qpb.yaxis[ch]
        bkg_cnt = spf_cnt+qpb_cnt
        cxb_spectrum_measurement_tmp += [(sci_cnt-bkg_cnt)]

    MakeRadialProjection(image_icrs_sci, radial_sci)
    MakeRadialProjection(image_icrs_qpb, radial_qpb)
    MakeRadialProjection(image_icrs_spf, radial_spf)

    if measure_cxb:
        edit_mode = 'w'
        if job!='0': edit_mode = 'a'

        cxb_file = open("%s/cxb_spectrum_%s.txt"%(cxb_output_dir,on_sample),edit_mode)
        cxb_file.write('# job %s; %s; %s; (%0.1f,%0.1f); %0.2f SPF; %0.1f \n'%(job,obsID,detector,sky_l_center,sky_b_center,spf_frac,exposure))
        for entry in range(0,len(cxb_spectrum_measurement_tmp)):
            cxb_file.write('%0.2e'%(cxb_spectrum_measurement_tmp[entry]))
            if entry!=len(cxb_spectrum_measurement_tmp)-1:
                cxb_file.write(' ')
            else:
                cxb_file.write('\n')
        cxb_file.close()

        get_cxb_spectrum(cxb_output_dir,spectrum_cxb,sky_l_center,sky_b_center,detector,obsID=obsID)

    elif include_cxb:
        #get_cxb_spectrum(cxb_output_dir,spectrum_cxb,sky_l_center,sky_b_center,detector,obsID=obsID)
        get_cxb_spectrum(cxb_output_dir,spectrum_cxb,sky_l_center,sky_b_center,detector)
        spectrum_comb_cxb.add(spectrum_cxb,factor=1./len_run_list)

    cxb_sum = spectrum_cxb.integral()
    qpb_sum = spectrum_qpb_radial.integral()
    cxb_2_qpb_ratio = cxb_sum/qpb_sum
    image_icrs_comb_cxb.add(image_icrs_cxb,factor=cxb_2_qpb_ratio)

    fig.clf()
    axbig = fig.add_subplot()
    axbig.plot(area_curve[0].xaxis,area_curve[0].yaxis,label='R0')
    axbig.plot(area_curve[1].xaxis,area_curve[1].yaxis,label='R1')
    axbig.plot(area_curve[2].xaxis,area_curve[2].yaxis,label='R2')
    axbig.plot(area_curve[3].xaxis,area_curve[3].yaxis,label='R3')
    axbig.set_xlabel('Energy [eV]')
    axbig.set_ylabel('Area [cm2]')
    axbig.legend(loc='best')
    fig.savefig("%s/area_curve_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    axbig = fig.add_subplot()
    axbig.plot(lightcurve_all.xaxis,lightcurve_all.yaxis,color='k')
    axbig.plot(lightcurve_sci.xaxis,lightcurve_sci.yaxis,color='r')
    axbig.set_yscale('log')
    axbig.set_xlabel('Time')
    axbig.set_ylabel('Count')
    axbig.legend(loc='best')
    fig.savefig("%s/lightcurve_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
    axbig.remove()

    font = {'family': 'serif',
            'color':  'white',
            'weight': 'normal',
            'size': 10,
            }
    fig.clf()
    axbig = fig.add_subplot()
    label_x = 'DETY'
    label_y = 'DETX'
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    xmin = image_det_all.xaxis.min()
    xmax = image_det_all.xaxis.max()
    ymin = image_det_all.yaxis.min()
    ymax = image_det_all.yaxis.max()
    axbig.imshow(image_det_all.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax), norm=colors.LogNorm())
    axbig.contour(image_det_mask.zaxis[:,:], np.arange(0.0, 1.5, 1.0), linestyles='solid', colors='red', extent=(xmin,xmax,ymin,ymax))
    axbig.text(-17000, 17000, 'Galac coord. = (%0.1f, %0.1f)'%(sky_l_center,sky_b_center), fontdict=font)
    fig.savefig("%s/image_det_all_log_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    axbig = fig.add_subplot()
    label_x = 'DETY'
    label_y = 'DETX'
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    xmin = image_det_all.xaxis.min()
    xmax = image_det_all.xaxis.max()
    ymin = image_det_all.yaxis.min()
    ymax = image_det_all.yaxis.max()
    axbig.imshow(image_det_all.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
    axbig.contour(image_det_mask.zaxis[:,:], np.arange(0.0, 1.5, 1.0), linestyles='solid', colors='red', extent=(xmin,xmax,ymin,ymax))
    axbig.text(-17000, 17000, 'Galac coord. = (%0.1f, %0.1f)'%(sky_l_center,sky_b_center), fontdict=font)
    fig.savefig("%s/image_det_all_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
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
    axbig.contour(image_det_mask.zaxis[:,:], np.arange(0.0, 1.5, 1.0), linestyles='solid', colors='red', extent=(xmin,xmax,ymin,ymax))
    axbig.text(-17000, 17000, 'Galac coord. = (%0.1f, %0.1f)'%(sky_l_center,sky_b_center), fontdict=font)
    fig.savefig("%s/image_det_sci_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region),bbox_inches='tight')
    axbig.remove()

    
    plot_data = spectrum_sci
    plot_bkg = []
    plot_color = []
    plot_label = []
    plot_bkg += [spectrum_cxb]
    plot_color += [color_list[2]]
    plot_label += ['CXB']
    plot_bkg += [spectrum_spf]
    plot_color += [color_list[1]]
    plot_label += ['SPF']
    plot_bkg += [spectrum_qpb]
    plot_color += [color_list[0]]
    plot_label += ['QPB']
    save_name = "%s/spectrum_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region)
    save_name_ul = "%s/spectrum_ul_job%s_%s_%s_%s.png"%(output_dir,job,obsID,detector,region)
    draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'Energy [eV]',my_spectrum_unit,save_name,save_name_ul,show_log_spectrum)

    return spectral_volume
    

def get_observation_pointings(run_list):

    sky_ra_first = 0.
    sky_dec_first = 0.
    sky_ra_avg = 0.
    sky_dec_avg = 0.
    sky_ra_min = 1e10
    sky_dec_min = 1e10
    sky_ra_max = -1e10
    sky_dec_max = -1e10
    n_obs = 0.
    for run in run_list:
        obsID = run.split('_')[0]
        detector = run.split('_')[1]
        input_dir = '/Users/rshang/xmm_analysis/output_extended_analysis/'+on_sample+'/'+ana_tag+'/ID'+obsID
        print (f'input_dir = {input_dir}')
        sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[0])
        sci_hdu_list = fits.open(sci_filename)
        sky_ra = sci_hdu_list[1].header['REF_RA']
        sky_dec = sci_hdu_list[1].header['REF_DEC']
        print (f'sky_ra = {sky_ra}')
        print (f'sky_dec = {sky_dec}')
        if sky_ra_min>sky_ra:
            sky_ra_min = sky_ra
        if sky_dec_min>sky_dec:
            sky_dec_min = sky_dec
        if sky_ra_max<sky_ra:
            sky_ra_max = sky_ra
        if sky_dec_max<sky_dec:
            sky_dec_max = sky_dec
        if n_obs==0:
            sky_ra_first = sky_ra
            sky_dec_first = sky_dec
        n_obs += 1.
        sky_ra_avg += sky_ra
        sky_dec_avg += sky_dec
    sky_ra_avg = sky_ra_avg/n_obs
    sky_dec_avg = sky_dec_avg/n_obs

    print (f'sky_ra_avg = {sky_ra_avg}')
    print (f'sky_dec_avg = {sky_dec_avg}')
    print (f'sky_ra_min = {sky_ra_min}')
    print (f'sky_dec_min = {sky_dec_min}')
    print (f'sky_ra_max = {sky_ra_max}')
    print (f'sky_dec_max = {sky_dec_max}')

    return sky_ra_min, sky_dec_min, sky_ra_max, sky_dec_max, sky_ra_first, sky_dec_first

def main():

    global mask_detx
    global mask_dety

    if not use_det_coord:
        for run in on_run_list:
            obsID = run.split('_')[0]
            detector = run.split('_')[1]
        
            ref_sky = []
            ref_det = []
            mtx_det_2_sky = []
            mtx_sky_2_det = []
            for idx_ra in range(0,1):
                for idx_dec in range(0,1):
                    ref_sky_local, ref_det_local, mtx_det_2_sky_local, mtx_sky_2_det_local = LoadCoordinateMatrix(idx_ra,idx_dec,on_sample,'ID'+obsID)
                    ref_sky += [ref_sky_local]
                    ref_det += [ref_det_local]
                    mtx_det_2_sky += [mtx_det_2_sky_local]
                    mtx_sky_2_det += [mtx_sky_2_det_local]
            mask_detx, mask_dety = ConvertSky2Det([mask_ra,mask_dec],ref_det[0],[sky_ra_center,sky_dec_center],mtx_sky_2_det[0])
            print (f'mask_ra = {mask_ra}')
            print (f'mask_dec = {mask_dec}')
            print (f'mask_detx = {mask_detx}')
            print (f'mask_dety = {mask_dety}')

    total_spectral_volume = 0.
    for run in on_run_list:
        obsID = run.split('_')[0]
        detector = run.split('_')[1]
        spectral_volume = analyze_one_observation(obsID,detector)
        total_spectral_volume += spectral_volume
    
    output_dir = '/Users/rshang/xmm_analysis/output_plots/plot_%s/'%(on_sample)
    common_functions.DrawSkyMap(fig,map_color,image_icrs_comb_sci,"%s/image_icrs_job%s_%s_%s.png"%(output_dir,job,obsID,region))
    common_functions.DrawSkyMap(fig,map_color,image_icrs_comb_sci,"%s/image_icrs_log_job%s_%s_%s.png"%(output_dir,job,obsID,region),log_scale=True)
    
    image_icrs_comb_xry.add(image_icrs_comb_sci)
    image_icrs_comb_xry.add(image_icrs_comb_qpb,factor=-1.)
    image_icrs_comb_xry.add(image_icrs_comb_spf,factor=-1.)
    image_icrs_comb_xry.add(image_icrs_comb_cxb,factor=-1.)
    common_functions.DrawSkyMap(fig,map_color,image_icrs_comb_xry,"%s/image_icrs_xry_log_job%s_%s_%s.png"%(output_dir,job,obsID,region))
    
    plot_data = spectrum_comb_sci
    plot_bkg = []
    plot_color = []
    plot_label = []
    plot_bkg += [spectrum_comb_cxb]
    plot_color += [color_list[2]]
    plot_label += ['CXB']
    plot_bkg += [spectrum_comb_spf]
    plot_color += [color_list[1]]
    plot_label += ['SPF']
    plot_bkg += [spectrum_comb_qpb]
    plot_color += [color_list[0]]
    plot_label += ['QPB']
    save_name = "%s/spectrum_comb_job%s_%s_%s.png"%(output_dir,job,obsID,region)
    save_name_ul = "%s/spectrum_comb_ul_job%s_%s_%s.png"%(output_dir,job,obsID,region)
    draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'Energy [eV]',my_spectrum_unit,save_name,save_name_ul,show_log_spectrum)
    
    MakeRadialProjection(image_icrs_comb_sci, radial_comb_sci)
    MakeRadialProjection(image_icrs_comb_qpb, radial_comb_qpb)
    MakeRadialProjection(image_icrs_comb_spf, radial_comb_spf)
    MakeRadialProjection(image_icrs_comb_cxb, radial_comb_cxb)
    
    cxb_error_ratio = np.sum(spectrum_comb_cxb.yerr)/np.sum(spectrum_comb_cxb.yaxis)
    for entry in range(0,len(radial_comb_cxb.yaxis)):
        radial_comb_cxb.yerr[entry] = radial_comb_cxb.yaxis[entry]*cxb_error_ratio
    plot_data = radial_comb_sci
    plot_bkg = []
    plot_color = []
    plot_label = []
    plot_bkg += [radial_comb_cxb]
    plot_color += [color_list[2]]
    plot_label += ['CXB']
    plot_bkg += [radial_comb_spf]
    plot_color += [color_list[1]]
    plot_label += ['SPF']
    plot_bkg += [radial_comb_qpb]
    plot_color += [color_list[0]]
    plot_label += ['QPB']
    save_name = "%s/radial_comb_job%s_%s_%s.png"%(output_dir,job,obsID,region)
    save_name_ul = "%s/radial_comb_ul_job%s_%s_%s.png"%(output_dir,job,obsID,region)
    draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'Angular distance [deg]',my_spectrum_unit,save_name,save_name_ul,True)
    
    spectrum_comb_raw_xry.add(spectrum_comb_raw_sci)
    spectrum_comb_raw_xry.add(spectrum_comb_raw_qpb,factor=-1.)
    spectrum_comb_raw_xry.add(spectrum_comb_raw_spf,factor=-1.)
    spectrum_comb_raw_xry.add(spectrum_comb_raw_cxb,factor=-1.)
    spectrum_comb_raw_bkg.add(spectrum_comb_raw_qpb)
    spectrum_comb_raw_bkg.add(spectrum_comb_raw_spf)
    spectrum_comb_raw_bkg.add(spectrum_comb_raw_cxb)
        
    output_filename = f'{output_dir}/spectrum_{obsID}_{region}.pkl'
    with open(output_filename,"wb") as file:
        pickle.dump(spectrum_comb_raw_xry, file)

    spectrum_comb_bkg.add(spectrum_comb_cxb)
    spectrum_comb_bkg.add(spectrum_comb_spf)
    spectrum_comb_bkg.add(spectrum_comb_qpb)
    list_chl = []
    list_ul = []
    for ch in range(0,len(spectrum_comb_sci.xaxis)):
        data = spectrum_comb_sci.yaxis[ch]
        bkg = spectrum_comb_bkg.yaxis[ch]
        err = spectrum_comb_bkg.yerr[ch]
        list_chl += [ch]
        list_ul += [3.*err]
    col_channel = fits.Column(name='CHANNEL', array=list_chl, format='I')
    col_upperlimit = fits.Column(name='UL', array=list_ul, format='J', unit=my_spectrum_unit)
    # Create a FITS header template from an XMM output
    input_filename = '/Users/rshang/xmm_analysis/'+on_sample+'/ID'+obsID+'/analysis/mos2-fov-r0-pi.fits'
    hdu_list = fits.open(input_filename)
    my_spec_header = hdu_list[1].header
    my_spec_table = fits.BinTableHDU.from_columns([col_channel,col_upperlimit],name='SPECTRUM',header=my_spec_header)
    hdul = fits.HDUList([hdu_list[0], my_spec_table])
    print ('write output to %s/table_upper_limit_%s_%s.fits'%(output_dir,obsID,region))
    hdul.writeto('%s/table_upper_limit_%s_%s.fits'%(output_dir,obsID,region), overwrite=True)
    
    if write_xspec_output:
    
        print (f'total_spectral_volume = {total_spectral_volume}')
        for ch in range(0,len(spectrum_comb_raw_cxb.xaxis)):
            ch_energy = spectrum_comb_raw_cxb.xaxis[ch]
            ch_rebin = spectrum_comb_cxb.get_bin(ch_energy)
            spectrum_comb_raw_cxb.yaxis[ch] = spectrum_comb_cxb.yaxis[ch_rebin]*total_spectral_volume
            spectrum_comb_raw_cxb.yaxis[ch] = max(0.,spectrum_comb_raw_cxb.yaxis[ch])
        
        plot_data = spectrum_comb_raw_sci
        plot_bkg = []
        plot_color = []
        plot_label = []
        plot_bkg += [spectrum_comb_raw_cxb]
        plot_color += [color_list[2]]
        plot_label += ['CXB']
        plot_bkg += [spectrum_comb_raw_spf]
        plot_color += [color_list[1]]
        plot_label += ['SPF']
        plot_bkg += [spectrum_comb_raw_qpb]
        plot_color += [color_list[0]]
        plot_label += ['QPB']
        save_name = "%s/spectrum_comb_raw_job%s_%s_%s.png"%(output_dir,job,obsID,region)
        save_name_ul = "%s/spectrum_comb_raw_ul_job%s_%s_%s.png"%(output_dir,job,obsID,region)
        draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'Energy [eV]',my_spectrum_unit,save_name,save_name_ul,show_log_spectrum)
        
        on_exposure = 0.
        for run in on_run_list:
            obsID = run.split('_')[0]
            input_dir = '/Users/rshang/xmm_analysis/output_extended_analysis/'+on_sample+'/'+ana_tag+'/ID'+obsID
            sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,'mos1',ana_ccd_bins[0])
            sci_hdu_list = fits.open(sci_filename)
            on_exposure += sci_hdu_list[1].header['EXPOSURE']
            sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,'mos2',ana_ccd_bins[0])
            sci_hdu_list = fits.open(sci_filename)
            on_exposure += sci_hdu_list[1].header['EXPOSURE']
        print (f'Total exposure = {on_exposure}')
        
        input_filename = '/Users/rshang/xmm_analysis/'+on_sample+'/ID'+obsID+'/analysis/mos2-fov-r0-arf.fits'
        hdu_list = fits.open(input_filename)
        mytable = Table.read(input_filename, hdu=1)
        original_nbins = len(mytable)
        spectrum_rebin_raw_xry = MyArray1D(bin_start=0,bin_end=original_nbins,pixel_scale=1)
        spectrum_rebin_raw_sci = MyArray1D(bin_start=0,bin_end=original_nbins,pixel_scale=1)
        spectrum_rebin_raw_bkg = MyArray1D(bin_start=0,bin_end=original_nbins,pixel_scale=1)
        for entry in range(0,len(mytable)):
            spectrum_rebin_raw_xry.xaxis[entry] = 0.5*(mytable[entry]['ENERG_LO']+mytable[entry]['ENERG_HI'])*1000.
            spectrum_rebin_raw_sci.xaxis[entry] = 0.5*(mytable[entry]['ENERG_LO']+mytable[entry]['ENERG_HI'])*1000.
            spectrum_rebin_raw_bkg.xaxis[entry] = 0.5*(mytable[entry]['ENERG_LO']+mytable[entry]['ENERG_HI'])*1000.
    
        # Create a FITS header template from an XMM output
        input_filename = '/Users/rshang/xmm_analysis/'+on_sample+'/ID'+obsID+'/analysis/mos2-fov-r0-pi.fits'
        hdu_list = fits.open(input_filename)
        mytable = Table.read(input_filename, hdu=1)
        original_delta_energy = spectrum_comb_raw_xry.xaxis[1]-spectrum_comb_raw_xry.xaxis[0]
        rebin_delta_energy = spectrum_rebin_raw_xry.xaxis[1]-spectrum_rebin_raw_xry.xaxis[0]
        rebin_energy = spectrum_comb_raw_xry.xaxis[0]
        ch_start = spectrum_rebin_raw_xry.get_bin(rebin_energy)
        for ch in range(ch_start,len(spectrum_rebin_raw_xry.xaxis)):
            ch_energy = spectrum_rebin_raw_xry.xaxis[ch]
            ch_rebin = spectrum_comb_raw_xry.get_bin(ch_energy)
            spectrum_rebin_raw_sci.yaxis[ch] = spectrum_comb_raw_sci.yaxis[ch_rebin]*rebin_delta_energy/original_delta_energy
            spectrum_rebin_raw_bkg.yaxis[ch] = spectrum_comb_raw_bkg.yaxis[ch_rebin]*rebin_delta_energy/original_delta_energy
            spectrum_rebin_raw_xry.yaxis[ch] = spectrum_comb_raw_xry.yaxis[ch_rebin]*rebin_delta_energy/original_delta_energy
        
        fig.clf()
        axbig = fig.add_subplot()
        axbig.plot(spectrum_rebin_raw_xry.xaxis,spectrum_rebin_raw_xry.yaxis,color='k')
        axbig.set_xlabel('Energy')
        axbig.set_ylabel('Count')
        fig.savefig("%s/spectrum_rebin_%s_%s.png"%(output_dir,obsID,region),bbox_inches='tight')
        axbig.remove()
    
        list_chl = []
        list_cnt = []
        list_bkg = []
        for ch in range(0,len(spectrum_rebin_raw_xry.xaxis)):
            list_chl += [ch]
            list_cnt += [spectrum_rebin_raw_sci.yaxis[ch]]
            list_bkg += [spectrum_rebin_raw_bkg.yaxis[ch]]
        
        # Create a FITS table extension
        col_channel = fits.Column(name='CHANNEL', array=list_chl, format='I')
        col_count = fits.Column(name='COUNTS', array=list_cnt, format='J', unit='count')
        col_bkgnd = fits.Column(name='COUNTS', array=list_bkg, format='J', unit='count')
        my_spec_header = hdu_list[1].header
        my_spec_table = fits.BinTableHDU.from_columns([col_channel,col_count],name='SPECTRUM',header=my_spec_header)
        my_bkgd_table = fits.BinTableHDU.from_columns([col_channel,col_bkgnd],name='BACKFILE',header=my_spec_header)
        
        arf_filename = '/Users/rshang/xmm_analysis/'+on_sample+'/ID'+obsID+'/analysis/mos2-fov-r0-arf.fits'
        my_arf_hdu = fits.open(arf_filename)[1]
        my_arf_hdu.name = 'ANCRFILE'
        rmf_filename = '/Users/rshang/xmm_analysis/'+on_sample+'/ID'+obsID+'/analysis/mos2-fov-r0-rmf.fits'
        my_rmf_hdu = fits.open(rmf_filename)[1]
        my_rmf_hdu.name = 'RESPFILE'
    
        # Combine the primary and table extensions
        print ('write output to %s/table_spectrum_data_%s_%s.fits'%(output_dir,obsID,region))
        hdul = fits.HDUList([hdu_list[0], my_spec_table])
        hdul.writeto('%s/table_spectrum_data_%s_%s.fits'%(output_dir,obsID,region), overwrite=True)
            
        print ('write output to %s/table_spectrum_bkgd_%s_%s.fits'%(output_dir,obsID,region))
        hdul = fits.HDUList([hdu_list[0], my_bkgd_table])
        hdul.writeto('%s/table_spectrum_bkgd_%s_%s.fits'%(output_dir,obsID,region), overwrite=True)
    
        input_filename = '/Users/rshang/xmm_analysis/'+on_sample+'/ID'+obsID+'/analysis/mos2-fov-r0-arf.fits'
        output_filename = '%s/table_arf_%s_%s.fits'%(output_dir,obsID,region)
        shutil.copyfile(input_filename, output_filename)
    
        input_filename = '/Users/rshang/xmm_analysis/'+on_sample+'/ID'+obsID+'/analysis/mos2-fov-r0-rmf.fits'
        output_filename = '%s/table_rmf_%s_%s.fits'%(output_dir,obsID,region)
        shutil.copyfile(input_filename, output_filename)

#use_det_coord = True
use_det_coord = False
mask_detx = -350.
mask_dety = -350.

region = 'r0'
#mask_inner_radius = 0.0
#mask_outer_radius = 0.015
mask_inner_radius = 0.015
mask_outer_radius = 0.3

list_region = []
list_mask_inner_radius = []
list_mask_outer_radius = []
list_region += [region]
list_mask_inner_radius += [mask_inner_radius]
list_mask_outer_radius += [mask_outer_radius]

# Geminga
mask_ra = 98.4760
mask_dec = 17.7689

# Cas A
#mask_ra = 350.8075+0.05
#mask_dec = 58.8072+0.01
#list_region += ['r0']
#list_mask_inner_radius += [0.0]
#list_mask_outer_radius += [0.05]
#list_region += ['r1']
#list_mask_inner_radius += [0.05]
#list_mask_outer_radius += [0.1]
#list_region += ['r2']
#list_mask_inner_radius += [0.1]
#list_mask_outer_radius += [0.15]
#list_region += ['r3']
#list_mask_inner_radius += [0.15]
#list_mask_outer_radius += [0.2]
#list_region += ['r4']
#list_mask_inner_radius += [0.2]
#list_mask_outer_radius += [0.25]
#list_region += ['r5']
#list_mask_inner_radius += [0.25]
#list_mask_outer_radius += [0.45]

# 3HWC J1928
#mask_ra = 292.1499583
#mask_dec = 17.9011111

# 3HWC J1930
#mask_ra = 292.6258333
#mask_dec = 18.8708889

# Be/X-ray binary system EXO 2030+375
#mask_ra = 308.0637083
#mask_dec = 37.6375000
#list_region += ['r0']
#list_mask_inner_radius += [0.0]
#list_mask_outer_radius += [0.3]
#list_region += ['r0']
#list_mask_inner_radius += [0.0]
#list_mask_outer_radius += [0.05]
#list_region += ['r1']
#list_mask_inner_radius += [0.05]
#list_mask_outer_radius += [0.1]
#list_region += ['r2']
#list_mask_inner_radius += [0.1]
#list_mask_outer_radius += [0.15]
#list_region += ['r3']
#list_mask_inner_radius += [0.15]
#list_mask_outer_radius += [0.2]
#list_region += ['r4']
#list_mask_inner_radius += [0.2]
#list_mask_outer_radius += [0.25]

sky_ra_lower, sky_dec_lower, sky_ra_upper, sky_dec_upper, sky_ra_center, sky_dec_center = get_observation_pointings(on_run_list)
sky_l_center, sky_b_center = ConvertRaDecToGalactic(sky_ra_center, sky_dec_center)
sky_l_lower, sky_b_lower = ConvertRaDecToGalactic(sky_ra_lower, sky_dec_lower)
sky_l_upper, sky_b_upper = ConvertRaDecToGalactic(sky_ra_upper, sky_dec_upper)

sky_ra_low = sky_ra_lower-0.28
sky_ra_high = sky_ra_upper+0.28
sky_dec_low = sky_dec_lower-0.28
sky_dec_high = sky_dec_upper+0.28
sky_l_low = sky_l_lower-0.28
sky_l_high = sky_l_upper+0.28
sky_b_low = sky_b_lower-0.28
sky_b_high = sky_b_upper+0.28

image_icrs_size = max(sky_ra_high-sky_ra_low,sky_dec_high-sky_dec_low)
image_icrs_comb_sci = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
image_icrs_comb_xry = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
image_icrs_comb_qpb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
image_icrs_comb_spf = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
image_icrs_comb_cxb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
spectrum_comb_sci = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_bkg = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_qpb = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_spf = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_cxb = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_xry = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_raw_sci = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_raw_qpb = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_raw_spf = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_raw_cxb = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_raw_xry = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
spectrum_comb_raw_bkg = MyArray1D(bin_start=energy_cut_lower,bin_end=energy_cut_upper,pixel_scale=ch_scale)
radial_comb_sci = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_qpb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_spf = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_cxb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_sci = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))

for r in range(0,len(list_region)):
    region = list_region[r]
    mask_inner_radius = list_mask_inner_radius[r]
    mask_outer_radius = list_mask_outer_radius[r]
    mask_inner_detr = mask_inner_radius/(0.05/(60.*60.))
    mask_outer_detr = mask_outer_radius/(0.05/(60.*60.))

    image_icrs_comb_sci.reset()
    image_icrs_comb_xry.reset()
    image_icrs_comb_qpb.reset()
    image_icrs_comb_spf.reset()
    image_icrs_comb_cxb.reset()
    spectrum_comb_sci.reset()
    spectrum_comb_bkg.reset()
    spectrum_comb_qpb.reset()
    spectrum_comb_spf.reset()
    spectrum_comb_cxb.reset()
    spectrum_comb_xry.reset()
    spectrum_comb_raw_sci.reset()
    spectrum_comb_raw_qpb.reset()
    spectrum_comb_raw_spf.reset()
    spectrum_comb_raw_cxb.reset()
    spectrum_comb_raw_xry.reset()
    spectrum_comb_raw_bkg.reset()
    radial_comb_sci.reset()
    radial_comb_qpb.reset()
    radial_comb_spf.reset()
    radial_comb_cxb.reset()
    radial_comb_sci.reset()

    main()


