
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
from scipy import fftpack
import matplotlib.pyplot as plt
from matplotlib import colors

import common_functions

on_sample = sys.argv[1]
on_obsID = sys.argv[2]
job = sys.argv[3]

measure_cxb = True
#measure_cxb = False

include_cxb = True
#include_cxb = False

# background study
find_extended_src = True
select_mask_events = False

# all events
#find_extended_src = False
#select_mask_events = False

# select events in RoI
#find_extended_src = True
#select_mask_events = True

show_log_spectrum = False

ana_ccd_bins = [0]
#ana_ccd_bins = [1,2,3,4,5,6,7]

energy_cut_lower = 2000
energy_cut_upper = 12000

energy_array = [2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000]
cxb_energy_array = [2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,10500,11000,11500,12000]
cxb_radial_array = [0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28]

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
            delta_x = image_map_xry.xaxis[idx_x] - sky_ra_center
            for idx_y in range(0,len(image_map_xry.yaxis)):
                delta_y = image_map_xry.yaxis[idx_y] - sky_dec_center
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

def find_point_sources_with_gravity(mode,image_data,image_mask,threshold=4.0):

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

def draw_stacked_histogram(fig,hist_data,hist_bkg,bkg_colors,bkg_labels,xlabel,ylabel,plot_name,log_scale=False):

    n_bins = len(hist_data.xaxis)
    stack_bkg = []
    for h1 in range(0,len(hist_bkg)):
        stack_bkg += [MyArray1D(bin_start=hist_data.xaxis[0],bin_end=hist_data.xaxis[n_bins-1],pixel_scale=hist_data.xaxis[1]-hist_data.xaxis[0])]
        for h2 in range(h1,len(hist_bkg)):
            stack_bkg[h1].add(hist_bkg[h2])
    
    max_qpb = max(hist_bkg[len(hist_bkg)-1].yaxis)

    fig.clf()
    axbig = fig.add_subplot()
    for h in range(0,len(hist_bkg)):
        axbig.fill_between(stack_bkg[h].xaxis,0.,stack_bkg[h].yaxis,color=bkg_colors[h],label=bkg_labels[h])
    axbig.fill_between(stack_bkg[0].xaxis,stack_bkg[0].yaxis-stack_bkg[0].yerr,stack_bkg[0].yaxis+stack_bkg[0].yerr,color='k',alpha=0.2)
    axbig.errorbar(hist_data.xaxis,hist_data.yaxis,yerr=hist_data.yerr,color='k',label='Data')
    if log_scale:
        axbig.set_yscale('log')
        axbig.set_ylim(bottom=0.5*max_qpb)
    axbig.set_xlabel(xlabel)
    axbig.set_ylabel(ylabel)
    axbig.legend(loc='best')
    fig.savefig("%s"%(plot_name),bbox_inches='tight')
    axbig.remove()

def get_cxb_spectrum(cxb_output_dir,hist_cxb,obs_sky_l,obs_sky_b,detector,obsID=''):

    cxb_file_name = "%s/cxb_spectrum_%s.txt"%(cxb_output_dir,on_sample)
    print (f'read {cxb_file_name}')
    cxb_file = open(cxb_file_name)
    cxb_measurement = []
    cxb_measurement_weight = []
    use_this_data = True
    for line in cxb_file:
        if '#' in line:
            line_split = line.split(';')
            line_obsid = line_split[1]
            line_detector = line_split[2]
            line_coord = line_split[3]
            line_coord = line_coord.strip('( )')
            distance = pow(obs_sky_l-float(line_coord[0]),2)+pow(obs_sky_b-float(line_coord[1]),2)
            distance = pow(distance,0.5)
            if not obsID=='':
                if not obsID in line_obsid:
                    use_this_data = False
            if not detector in line_detector:
                use_this_data = False
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
        cxb_measurement_weight += [1.]
    cxb_file.close()

    if not obsID=='':
        print (f'Load data for obsID = {obsID}, detector = {detector}')
    print (f'{len(cxb_measurement)} CXB measurements are used.')

    hist_cxb_measurement = MyArray1D(bin_start=cxb_energy_array[0],bin_end=cxb_energy_array[len(cxb_energy_array)-1],pixel_scale=cxb_energy_array[1]-cxb_energy_array[0])
    for ch in range(0,len(cxb_energy_array)-1):
        n_samples = 0.
        avg_cxb = 0.
        total_weight = 0.
        for m in range(0,len(cxb_measurement)):
            n_samples += 1.
            total_weight += cxb_measurement_weight[m]
            avg_cxb += cxb_measurement[m][ch]*cxb_measurement_weight[m]
        avg_cxb = avg_cxb/total_weight
        hist_cxb_measurement.yaxis[ch] = avg_cxb
    for ch in range(0,len(cxb_energy_array)-1):
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

def get_cxb_radial(cxb_output_dir,hist_cxb,obs_sky_l,obs_sky_b,detector,obsID=''):

    cxb_file_name = "%s/cxb_radial_%s.txt"%(cxb_output_dir,on_sample)
    print (f'read {cxb_file_name}')
    cxb_file = open(cxb_file_name)
    cxb_measurement = []
    cxb_measurement_weight = []
    use_this_data = True
    for line in cxb_file:
        if '#' in line:
            line_split = line.split(';')
            line_obsid = line_split[1]
            line_detector = line_split[2]
            line_coord = line_split[3]
            line_coord = line_coord.strip('( )')
            distance = pow(obs_sky_l-float(line_coord[0]),2)+pow(obs_sky_b-float(line_coord[1]),2)
            distance = pow(distance,0.5)
            if not obsID=='':
                if not obsID in line_obsid:
                    use_this_data = False
            if not detector in line_detector:
                use_this_data = False
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
        cxb_measurement_weight += [1.]
    cxb_file.close()

    hist_cxb_measurement = MyArray1D(bin_start=cxb_radial_array[0],bin_end=cxb_radial_array[len(cxb_radial_array)-1],pixel_scale=cxb_radial_array[1]-cxb_radial_array[0])
    for ch in range(0,len(cxb_radial_array)-1):
        n_samples = 0.
        avg_cxb = 0.
        total_weight = 0.
        for m in range(0,len(cxb_measurement)):
            n_samples += 1.
            total_weight += cxb_measurement_weight[m]
            avg_cxb += cxb_measurement[m][ch]*cxb_measurement_weight[m]
        avg_cxb = avg_cxb/total_weight
        hist_cxb_measurement.yaxis[ch] = avg_cxb
    for ch in range(0,len(cxb_radial_array)-1):
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
        dist = hist_cxb.xaxis[x]
        hist_cxb.yaxis[x] = hist_cxb_measurement.get_bin_content(dist)
        hist_cxb.yerr[x] = hist_cxb_measurement.get_bin_error(dist)

def analyze_one_observation(obsID,detector):

    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print (f'analyze obsID {obsID} detector {detector}')

    input_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/ID'+obsID
    cxb_output_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample
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
    fig.savefig("%s/image_det_diff_lofreq_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
    axbig.remove()

    if find_extended_src:
        hist_gravity = find_point_sources_with_gravity(fft_lofreq_mode,image_det_diff_lofreq,image_det_mask)
        hist_gravity = find_point_sources_with_gravity(fft_lofreq_mode,image_det_diff_lofreq,image_det_mask)
        fig.clf()
        axbig = fig.add_subplot()
        axbig.plot(hist_gravity.xaxis,hist_gravity.yaxis)
        axbig.set_yscale('log')
        axbig.set_xlabel('Gravity')
        axbig.legend(loc='best')
        fig.savefig("%s/gravity_lofreq_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
        axbig.remove()

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
    fig.savefig("%s/image_det_mask_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
    axbig.remove()


    image_det_all = MyArray2D(pixel_scale=detx_scale)
    image_det_sci = MyArray2D(pixel_scale=detx_scale)
    image_det_qpb = MyArray2D(pixel_scale=detx_scale)
    image_det_spf = MyArray2D(pixel_scale=detx_scale)
    image_det_cxb = MyArray2D(pixel_scale=detx_scale)
    image_det_xry = MyArray2D(pixel_scale=detx_scale)
    detx_sci = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=map_rebin*detx_scale)
    detx_bkg = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=map_rebin*detx_scale)
    detx_qpb = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=map_rebin*detx_scale)
    detx_spf = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=map_rebin*detx_scale)
    detx_cxb = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=map_rebin*detx_scale)
    dety_sci = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=map_rebin*dety_scale)
    dety_bkg = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=map_rebin*dety_scale)
    dety_qpb = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=map_rebin*dety_scale)
    dety_spf = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=map_rebin*dety_scale)
    dety_cxb = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=map_rebin*dety_scale)
    spectrum_sci = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    spectrum_qpb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    spectrum_spf = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    spectrum_cxb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    spectrum_xry = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
    image_icrs_sci = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_xry = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_qpb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_spf = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    image_icrs_cxb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=image_icrs_size,pixel_scale=detx_scale*0.05/(60.*60.))
    radial_sci = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
    radial_qpb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
    radial_spf = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
    radial_cxb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
    radial_sci = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))

    area_curve = MyArray1D(bin_start=2000,bin_end=12000,pixel_scale=2000)
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

    for ccd in range(0,len(ana_ccd_bins)):
    
        arf_filename = '%s/sci_area_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
        arf_hdu_list = fits.open(arf_filename)
        arf_table = Table.read(arf_filename, hdu=1)
        for entry in range(0,len(arf_table)):
            entry_energy = arf_table[entry]['ENERGY']
            entry_area = arf_table[entry]['AREA']
            print (f'energy = {entry_energy}, area = {entry_area}')
            area_curve.fill(entry_energy,weight=entry_area)

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
        exposure = sci_hdu_list[1].header['obs_duration']
        print (f'exposure = {exposure}')
        sci_table = Table.read(sci_filename, hdu=1)
        for entry in range(0,len(sci_table)):
            evt_pi = sci_table[entry]['PI']
            evt_detx = sci_table[entry]['DETX']
            evt_dety = sci_table[entry]['DETY']
            evt_ra = sci_table[entry]['RA']
            evt_dec = sci_table[entry]['DEC']
            if evt_pi<energy_cut_lower: continue
            if evt_pi>energy_cut_upper: continue
            spectral_density = ch_scale/1000.*fov_size*exposure*area_curve.yaxis[0]
            spacial_density = (energy_cut_upper-energy_cut_lower)/1000.*pix_size*exposure*area_curve.yaxis[0]
            image_det_all.fill(evt_detx,evt_dety,weight=1./spacial_density)
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if not select_mask_events:
                if mask==1: continue
            else:
                if mask!=1: continue
            delta_ra = evt_ra - sky_ra_center
            delta_dec = evt_dec - sky_dec_center
            spectrum_sci.fill(evt_pi,weight=1./spectral_density)
            spectrum_comb_sci.fill(evt_pi,weight=1./spectral_density*1./len_run_list)
            image_det_sci.fill(evt_detx,evt_dety,weight=1./spacial_density)
            image_icrs_sci.fill(evt_ra,evt_dec,weight=1./spacial_density)
            image_icrs_comb_sci.fill(evt_ra,evt_dec,weight=1./spacial_density*1./len_run_list)
            detx_sci.fill(evt_detx,weight=1./spacial_density)
            dety_sci.fill(evt_dety,weight=1./spacial_density)

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
            spectral_density = ch_scale/1000.*fov_size*exposure*area_curve.yaxis[0]
            spacial_density = (energy_cut_upper-energy_cut_lower)/1000.*pix_size*exposure*area_curve.yaxis[0]
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if not select_mask_events:
                if mask==1: continue
            else:
                if mask!=1: continue
            delta_ra = evt_ra - sky_ra_center
            delta_dec = evt_dec - sky_dec_center
            spectrum_qpb.fill(evt_pi,weight=1./spectral_density*1./sample_scale)
            spectrum_comb_qpb.fill(evt_pi,weight=1./spectral_density*1./sample_scale*1./len_run_list)
            image_det_qpb.fill(evt_detx,evt_dety,weight=1./spacial_density*1./sample_scale)
            image_icrs_qpb.fill(evt_ra,evt_dec,weight=1./spacial_density*1./sample_scale)
            image_icrs_comb_qpb.fill(evt_ra,evt_dec,weight=1./spacial_density*1./sample_scale*1./len_run_list)
            detx_qpb.fill(evt_detx,weight=1./spacial_density*1./sample_scale)
            dety_qpb.fill(evt_dety,weight=1./spacial_density*1./sample_scale)

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
            spectral_density = ch_scale/1000.*fov_size*exposure*area_curve.yaxis[0]
            spacial_density = (energy_cut_upper-energy_cut_lower)/1000.*pix_size*exposure*area_curve.yaxis[0]
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if not select_mask_events:
                if mask==1: continue
            else:
                if mask!=1: continue
            delta_ra = evt_ra - sky_ra_center
            delta_dec = evt_dec - sky_dec_center
            spectrum_spf.fill(evt_pi,weight=1./spectral_density*1./sample_scale)
            spectrum_comb_spf.fill(evt_pi,weight=1./spectral_density*1./sample_scale*1./len_run_list)
            image_det_spf.fill(evt_detx,evt_dety,weight=1./spacial_density*1./sample_scale)
            image_icrs_spf.fill(evt_ra,evt_dec,weight=1./spacial_density*1./sample_scale)
            image_icrs_comb_spf.fill(evt_ra,evt_dec,weight=1./spacial_density*1./sample_scale*1./len_run_list)
            detx_spf.fill(evt_detx,weight=1./spacial_density*1./sample_scale)
            dety_spf.fill(evt_dety,weight=1./spacial_density*1./sample_scale)

    qpb_cnt = spectrum_qpb.integral()
    spf_cnt = spectrum_spf.integral()
    spf_frac = spf_cnt/qpb_cnt
    print (f'spf_cnt/qpb_cnt = {spf_cnt/qpb_cnt}')

    cxb_spectrum_measurement_tmp = []
    for ch in range(0,len(cxb_energy_array)-1):
        sci_cnt = spectrum_sci.integral(integral_range=[cxb_energy_array[ch],cxb_energy_array[ch+1]])
        spf_cnt = spectrum_spf.integral(integral_range=[cxb_energy_array[ch],cxb_energy_array[ch+1]])
        qpb_cnt = spectrum_qpb.integral(integral_range=[cxb_energy_array[ch],cxb_energy_array[ch+1]])
        bkg_cnt = spf_cnt+qpb_cnt
        nbins = (cxb_energy_array[ch+1]-cxb_energy_array[ch])/ch_scale
        cxb_spectrum_measurement_tmp += [(sci_cnt-bkg_cnt)/nbins]

    MakeRadialProjection(image_icrs_sci, radial_sci)
    MakeRadialProjection(image_icrs_qpb, radial_qpb)
    MakeRadialProjection(image_icrs_spf, radial_spf)

    cxb_radial_measurement_tmp = []
    for ch in range(0,len(cxb_radial_array)-1):
        sci_cnt = radial_sci.integral(integral_range=[cxb_radial_array[ch],cxb_radial_array[ch+1]])
        spf_cnt = radial_spf.integral(integral_range=[cxb_radial_array[ch],cxb_radial_array[ch+1]])
        qpb_cnt = radial_qpb.integral(integral_range=[cxb_radial_array[ch],cxb_radial_array[ch+1]])
        bkg_cnt = spf_cnt+qpb_cnt
        nbins = (cxb_radial_array[ch+1]-cxb_radial_array[ch])/(2.*detx_scale*0.05/(60.*60.))
        cxb_radial_measurement_tmp += [(sci_cnt-bkg_cnt)/nbins]

    if measure_cxb:
        edit_mode = 'w'
        if job!='0': edit_mode = 'a'

        cxb_file = open("%s/cxb_spectrum_%s.txt"%(cxb_output_dir,on_sample),edit_mode)
        cxb_file.write('# job %s; %s; %s; (%0.1f,%0.1f); %0.2f SPF \n'%(job,obsID,detector,sky_l_center,sky_b_center,spf_frac))
        for entry in range(0,len(cxb_spectrum_measurement_tmp)):
            cxb_file.write('%0.2e'%(cxb_spectrum_measurement_tmp[entry]))
            if entry!=len(cxb_spectrum_measurement_tmp)-1:
                cxb_file.write(' ')
            else:
                cxb_file.write('\n')
        cxb_file.close()

        cxb_file = open("%s/cxb_radial_%s.txt"%(cxb_output_dir,on_sample),edit_mode)
        cxb_file.write('# job %s; %s; %s; (%0.1f,%0.1f); %0.2f SPF \n'%(job,obsID,detector,sky_l_center,sky_b_center,spf_frac))
        for entry in range(0,len(cxb_radial_measurement_tmp)):
            cxb_file.write('%0.2e'%(cxb_radial_measurement_tmp[entry]))
            if entry!=len(cxb_radial_measurement_tmp)-1:
                cxb_file.write(' ')
            else:
                cxb_file.write('\n')
        cxb_file.close()

        get_cxb_spectrum(cxb_output_dir,spectrum_cxb,sky_l_center,sky_b_center,detector,obsID=obsID)
        get_cxb_radial(cxb_output_dir,radial_cxb,sky_l_center,sky_b_center,detector,obsID=obsID)

    elif include_cxb:
        get_cxb_spectrum(cxb_output_dir,spectrum_cxb,sky_l_center,sky_b_center,detector)
        spectrum_comb_cxb.add(spectrum_cxb,factor=1./len_run_list)
        get_cxb_radial(cxb_output_dir,radial_cxb,sky_l_center,sky_b_center,detector)
        radial_comb_cxb.add(radial_cxb,factor=1./len_run_list)

    cxb_sum = spectrum_cxb.integral()
    qpb_sum = spectrum_qpb.integral()
    cxb_2_qpb_ratio = cxb_sum/qpb_sum
    image_icrs_comb_cxb.add(image_icrs_comb_qpb,factor=cxb_2_qpb_ratio)
    detx_cxb.add(detx_qpb,factor=cxb_2_qpb_ratio)
    dety_cxb.add(dety_qpb,factor=cxb_2_qpb_ratio)

    fig.clf()
    axbig = fig.add_subplot()
    axbig.plot(area_curve.xaxis,area_curve.yaxis,color='k')
    axbig.set_xlabel('Energy [eV]')
    axbig.set_ylabel('Area [cm2]')
    axbig.legend(loc='best')
    fig.savefig("%s/area_curve_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    axbig = fig.add_subplot()
    axbig.plot(lightcurve_all.xaxis,lightcurve_all.yaxis,color='k')
    axbig.plot(lightcurve_sci.xaxis,lightcurve_sci.yaxis,color='r')
    axbig.set_yscale('log')
    axbig.set_xlabel('Time')
    axbig.set_ylabel('Count')
    axbig.legend(loc='best')
    fig.savefig("%s/lightcurve_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
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
    axbig.imshow(image_det_all.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax), norm=colors.LogNorm())
    axbig.contour(image_det_mask.zaxis[:,:], np.arange(0.0, 1.5, 1.0), linestyles='solid', colors='red', extent=(xmin,xmax,ymin,ymax))
    fig.savefig("%s/image_det_all_log_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
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
    fig.savefig("%s/image_det_all_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
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
    fig.savefig("%s/image_det_sci_job%s_%s_%s.png"%(output_dir,job,obsID,detector),bbox_inches='tight')
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
    save_name = "%s/spectrum_job%s_%s_%s.png"%(output_dir,job,obsID,detector)
    draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'Energy [eV]','Photons /cm2/s/sr/keV',save_name,show_log_spectrum)
    
    plot_data = detx_sci
    plot_bkg = []
    plot_color = []
    plot_label = []
    plot_bkg += [detx_cxb]
    plot_color += [color_list[2]]
    plot_label += ['CXB']
    plot_bkg += [detx_spf]
    plot_color += [color_list[1]]
    plot_label += ['SPF']
    plot_bkg += [detx_qpb]
    plot_color += [color_list[0]]
    plot_label += ['QPB']
    save_name = "%s/detx_job%s_%s_%s.png"%(output_dir,job,obsID,detector)
    draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'DETX','Photons /cm2/s/sr/keV',save_name)
    
    plot_data = dety_sci
    plot_bkg = []
    plot_color = []
    plot_label = []
    plot_bkg += [dety_cxb]
    plot_color += [color_list[2]]
    plot_label += ['CXB']
    plot_bkg += [dety_spf]
    plot_color += [color_list[1]]
    plot_label += ['SPF']
    plot_bkg += [dety_qpb]
    plot_color += [color_list[0]]
    plot_label += ['QPB']
    save_name = "%s/dety_job%s_%s_%s.png"%(output_dir,job,obsID,detector)
    draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'DETY','Photons /cm2/s/sr/keV',save_name)


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
        input_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/ID'+obsID
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

sky_ra_lower, sky_dec_lower, sky_ra_upper, sky_dec_upper, sky_ra_center, sky_dec_center = get_observation_pointings(on_run_list)
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
spectrum_comb_sci = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_comb_qpb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_comb_spf = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_comb_cxb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_comb_xry = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
radial_comb_sci = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_qpb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_spf = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_cxb = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
radial_comb_sci = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))

for run in on_run_list:
    obsID = run.split('_')[0]
    detector = run.split('_')[1]
    analyze_one_observation(obsID,detector)

output_dir = '/Users/rshang/xmm_analysis/output_plots/plot_%s/'%(on_sample)
common_functions.DrawSkyMap(fig,map_color,image_icrs_comb_sci,"%s/image_icrs_job%s_%s.png"%(output_dir,job,obsID))
common_functions.DrawSkyMap(fig,map_color,image_icrs_comb_sci,"%s/image_icrs_log_job%s_%s.png"%(output_dir,job,obsID),log_scale=True)

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
save_name = "%s/spectrum_comb_job%s_%s.png"%(output_dir,job,obsID)
draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'Energy [eV]','Photons /cm2/s/sr/keV',save_name,show_log_spectrum)
    
MakeRadialProjection(image_icrs_comb_sci, radial_comb_sci)
MakeRadialProjection(image_icrs_comb_qpb, radial_comb_qpb)
MakeRadialProjection(image_icrs_comb_spf, radial_comb_spf)

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
save_name = "%s/radial_comb_job%s_%s.png"%(output_dir,job,obsID)
draw_stacked_histogram(fig,plot_data,plot_bkg,plot_color,plot_label,'Angular distance [deg]','Photons /cm2/s/sr/keV',save_name,False)
    
