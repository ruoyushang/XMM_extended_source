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

#on_sample = 'extragalactic'
#on_obsID = 'ID0690900101' # 90% SPF
#on_obsID = 'ID0827241101' # 51% SPF
#on_obsID = 'ID0505460501' # 27% SPF
#on_obsID = 'ID0804860301' # 25% SPF
#on_obsID = 'ID0803160301' # 10% SPF
#on_obsID = 'ID0820560101' # 0% SPF
#on_obsID = 'ID0827200401' # 4% SPF
#on_obsID = 'ID0827200501' # 0% SPF
#on_obsID = 'ID0827251001' # 0% SPF
#on_obsID = 'ID0827251101' # 0% SPF

#on_sample = 'Cas_A'
#on_obsID = 'ID0412180101'
#on_obsID = 'ID0400210101' # Cas A Northern lobe
#on_obsID = 'ID0782961401' # angular distance to Cas A: 34.7 arcmin

on_sample = '3HWC_J1928_p178'
#on_obsID = 'ID0902120101'
#on_obsID = 'ID0503740101'
#on_obsID = 'ID0830192001'
#on_obsID = 'ID0406730101'
#on_obsID = 'ID0762980101'
#on_obsID = 'ID0742710301'
on_obsID = 'ID0822330301'
#on_obsID = 'ID0861880101' # 57% SPF
#on_obsID = 'ID0841190101' # bright source
#on_obsID = 'ID0851181701' # bright source
#on_obsID = 'ID0724270101' # SNR with escaped CRs?
#on_obsID = 'ID0724270201' # bright SNR
#on_obsID = 'ID0742710401' # broken chip
#on_obsID = 'ID0763240101' # strange ring patterns
#on_obsID = 'ID0822330201' # strange ring patterns

#detector = 'mos1'
detector = 'mos2'

ana_ccd_bins = [0]
#ana_ccd_bins = [1,2,3,4,5,6,7]

exclusion_inner = 0.
#exclusion_inner = 0.10
#exclusion_inner = 0.15
#exclusion_outer = 0.15
exclusion_outer = 1e10

point_source_cut = True
#point_source_cut = False

do_fit = False

energy_cut_lower = 2000
energy_cut_upper = 12000

energy_array = [2000,4000,6000,8000,10000,12000]

roi_ra = 350.85
roi_dec = 58.815

sample_scale = 10.

map_rebin = 3
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


def get_cxb_spectrum(hist_cxb):

    cxb_measurement = []
    cxb_measurement += [[1.2460764078351871, 0.2595672221132602, -0.06124044036042946, 0.010761918389746227, -0.21767261123909448]]
    cxb_measurement += [[1.5110231403649634, 0.2640095032122119, -0.3466999797149171, -0.4525475397740095, -0.5749894471585769]]
    cxb_measurement += [[1.5431756741827456, 0.6109999767024142, 0.1452067845868072, 0.22299613347240382, -0.0185044360227856]]
    cxb_measurement += [[1.873602737378438, 0.9950159236561099, 0.6393535129073321, 0.6287842997576867, 0.31489544574308825]]
    cxb_measurement += [[1.0822887673414292, 0.41326724730222963, 0.19304157013872897, 0.12473719995592134, 0.10555347885917783]]
    cxb_measurement += [[0.8917712578452764, 0.48677549677590903, 0.1328807248329376, -0.15747818141105654, -0.051373285693748015]]
    cxb_measurement += [[1.1162834709255116, 0.4128565130768429, 0.2193164617501438, 0.04738693633387427, -0.013712703188156183]]
    cxb_measurement += [[1.3413456972931797, 0.4922646398350615, 0.13633582524099525, 0.045889076594834675, -0.06825667292332441]]

    hist_cxb_measurement = MyArray1D(bin_start=energy_array[0],bin_end=energy_array[len(energy_array)-1],pixel_scale=energy_array[1]-energy_array[0])
    for ch in range(0,len(energy_array)-1):
        n_samples = 0.
        avg_cxb = 0.
        for m in range(0,len(cxb_measurement)):
            n_samples += 1.
            avg_cxb += cxb_measurement[m][ch]
        avg_cxb = avg_cxb/n_samples
        hist_cxb_measurement.yaxis[ch] = avg_cxb
    #print (f'hist_cxb_measurement.yaxis:')
    #print (f'{hist_cxb_measurement.yaxis}')

    for x in range(0,len(hist_cxb.xaxis)):
        energy = hist_cxb.xaxis[x]
        if energy<energy_cut_lower: continue
        if energy>energy_cut_upper: continue
        hist_cxb.yaxis[x] = hist_cxb_measurement.get_bin_content(energy)
        
def get_gxb_spectrum(hist_gxb):

    gxb_measurement = []
    gxb_measurement += [[1.3516506948299853, 0.5873110283798142, 0.20035191142546557, 0.09132174038363583, 0.14828325188954022]]
    gxb_measurement += [[2.096893323374876, 0.6825092676677877, 0.20085330998619902, -0.12528111058407812, -0.25872850490039534]]
    gxb_measurement += [[1.721944994819511, 0.6134754511413101, 0.10109597873251537, -0.01706314759658094, -0.15607266246362891]]
    gxb_measurement += [[1.764603682517303, 0.5002385196933191, -0.033018911413896256, -0.23136404597757654, -0.17491760224497677]]
    gxb_measurement += [[1.7929416062600456, 0.8932197038025826, 0.30368136538555723, 0.1469531272816006, -0.01657585506453758]]
    gxb_measurement += [[2.0574467391937805, 1.2317014735485052, 0.6749214331887166, 0.2727590099893934, 0.14625001792403572]]

    hist_gxb_measurement = MyArray1D(bin_start=energy_array[0],bin_end=energy_array[len(energy_array)-1],pixel_scale=energy_array[1]-energy_array[0])
    for ch in range(0,len(energy_array)-1):
        n_samples = 0.
        avg_gxb = 0.
        for m in range(0,len(gxb_measurement)):
            n_samples += 1.
            avg_gxb += gxb_measurement[m][ch]
        avg_gxb = avg_gxb/n_samples
        hist_gxb_measurement.yaxis[ch] = avg_gxb
    #print (f'hist_gxb_measurement.yaxis:')
    #print (f'{hist_gxb_measurement.yaxis}')

    for x in range(0,len(hist_gxb.xaxis)):
        energy = hist_gxb.xaxis[x]
        if energy<energy_cut_lower: continue
        if energy>energy_cut_upper: continue
        hist_gxb.yaxis[x] = hist_gxb_measurement.get_bin_content(energy)
        
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

def MakeGalacticProjection(image_map_xry, profile_x, profile_y):

    for idx_z in range(0,len(profile_x.xaxis)):
        profile_x.yaxis[idx_z] = 0.
    for idx_z in range(0,len(profile_y.xaxis)):
        profile_y.yaxis[idx_z] = 0.

    delta_pix_z = profile_x.xaxis[1]-profile_x.xaxis[0]
    for idx_z in range(0,len(profile_x.xaxis)):
        pix_z_lo = profile_x.xaxis[idx_z]
        pix_z_hi = profile_x.xaxis[idx_z]+delta_pix_z
        pix_integral_xry = 0.
        pix_count = 0.
        for idx_x in range(0,len(image_map_xry.xaxis)):
            pix_x = image_map_xry.xaxis[idx_x]
            if pix_x<pix_z_lo: continue
            if pix_x>=pix_z_hi: continue
            for idx_y in range(0,len(image_map_xry.yaxis)):
                pix_y = image_map_xry.yaxis[idx_y]
                pix_content_xry = image_map_xry.zaxis[idx_x,idx_y]
                pix_integral_xry += pix_content_xry
                pix_count += 1.
        profile_x.yaxis[idx_z] = pix_integral_xry

    delta_pix_z = profile_y.xaxis[1]-profile_y.xaxis[0]
    for idx_z in range(0,len(profile_y.xaxis)):
        pix_z_lo = profile_y.xaxis[idx_z]
        pix_z_hi = profile_y.xaxis[idx_z]+delta_pix_z
        pix_integral_xry = 0.
        pix_count = 0.
        for idx_y in range(0,len(image_map_xry.yaxis)):
            pix_y = image_map_xry.yaxis[idx_y]
            if pix_y<pix_z_lo: continue
            if pix_y>=pix_z_hi: continue
            for idx_x in range(0,len(image_map_xry.xaxis)):
                pix_x = image_map_xry.xaxis[idx_x]
                pix_content_xry = image_map_xry.zaxis[idx_x,idx_y]
                pix_integral_xry += pix_content_xry
                pix_count += 1.
        profile_y.yaxis[idx_z] = pix_integral_xry

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
            pix_x = image_map_xry.xaxis[idx_x]
            for idx_y in range(0,len(image_map_xry.yaxis)):
                pix_y = image_map_xry.yaxis[idx_y]
                distance = pow(pix_x*pix_x+pix_y*pix_y,0.5)*0.05/(60.*60.)
                if distance<pix_z_lo: continue
                if distance>=pix_z_hi: continue
                pix_content_xry = image_map_xry.zaxis[idx_x,idx_y]
                pix_integral_xry += pix_content_xry
                pix_count += 1.
        profile_r.yaxis[idx_z] = pix_integral_xry



sky_ra_center = 0.
sky_dec_center = 0.

input_dir = '/Users/rshang/xmm_analysis/output_plots/'+on_sample+'/'+on_obsID
output_dir = '/Users/rshang/xmm_analysis/output_plots/'

image_det_sci_fine = MyArray2D()
image_det_mask = MyArray2D(pixel_scale=map_rebin*detx_scale)
image_det_sci = MyArray2D(pixel_scale=map_rebin*detx_scale)
image_det_bkg = MyArray2D(pixel_scale=map_rebin*detx_scale)
image_det_diff = MyArray2D(pixel_scale=map_rebin*detx_scale)
image_det_diff_lofreq = MyArray2D(pixel_scale=map_rebin*detx_scale)
image_det_diff_hifreq = MyArray2D(pixel_scale=map_rebin*detx_scale)
for ccd in range(0,len(ana_ccd_bins)):

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    sci_hdu_list = fits.open(sci_filename)
    sci_table = Table.read(sci_filename, hdu=1)
    for entry in range(0,len(sci_table)):
        evt_pi = sci_table[entry]['PI']
        evt_detx = sci_table[entry]['DETX']
        evt_dety = sci_table[entry]['DETY']
        if evt_pi<energy_cut_lower: continue
        if evt_pi>energy_cut_upper: continue
        image_det_sci.fill(evt_detx,evt_dety)
        image_det_sci_fine.fill(evt_detx,evt_dety)

    sky_ra_center = sci_hdu_list[1].header['REF_RA']
    sky_dec_center = sci_hdu_list[1].header['REF_DEC']
    sky_l_center, sky_b_center = ConvertRaDecToGalactic(sky_ra_center, sky_dec_center)

    bkg_filename = '%s/qpb_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    bkg_hdu_list = fits.open(bkg_filename)
    bkg_table = Table.read(bkg_filename, hdu=1)
    for entry in range(0,len(bkg_table)):
        evt_pi = bkg_table[entry]['PI']
        evt_detx = bkg_table[entry]['DETX']
        evt_dety = bkg_table[entry]['DETY']
        if evt_pi<energy_cut_lower: continue
        if evt_pi>energy_cut_upper: continue
        image_det_bkg.fill(evt_detx,evt_dety,weight=1./sample_scale)
    bkg_filename = '%s/spf_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    bkg_hdu_list = fits.open(bkg_filename)
    bkg_table = Table.read(bkg_filename, hdu=1)
    for entry in range(0,len(bkg_table)):
        evt_pi = bkg_table[entry]['PI']
        evt_detx = bkg_table[entry]['DETX']
        evt_dety = bkg_table[entry]['DETY']
        if evt_pi<energy_cut_lower: continue
        if evt_pi>energy_cut_upper: continue
        image_det_bkg.fill(evt_detx,evt_dety,weight=1./sample_scale)

sky_ra_low = sky_ra_center-0.28
sky_ra_high = sky_ra_center+0.28
sky_dec_low = sky_dec_center-0.28
sky_dec_high = sky_dec_center+0.28
sky_l_low = sky_l_center-0.28
sky_l_high = sky_l_center+0.28
sky_b_low = sky_b_center-0.28
sky_b_high = sky_b_center+0.28

def smooth_map(image_data):

    print ('smoothing map...')
    kernel_radius = 4*detx_scale
    kernel_pix_size = int(kernel_radius/detx_scale)
    gaus_norm = 2.*np.pi*kernel_radius*kernel_radius
    image_kernel = MyArray2D(pixel_scale=detx_scale)
    central_bin_x = int(len(image_kernel.xaxis)/2)
    central_bin_y = int(len(image_kernel.yaxis)/2)
    for idx_x in range(0,len(image_kernel.xaxis)):
        for idx_y in range(0,len(image_kernel.yaxis)):
            pix_x = image_kernel.xaxis[idx_x]
            pix_y = image_kernel.yaxis[idx_y]
            distance = pow(pix_x*pix_x+pix_y*pix_y,0.5)
            pix_content = np.exp(-(distance*distance)/(2.*kernel_radius*kernel_radius))
            image_kernel.zaxis[idx_x,idx_y] = pix_content/gaus_norm

    image_smooth = MyArray2D(pixel_scale=detx_scale)
    for idx_x1 in range(0,len(image_data.xaxis)):
        print (f'smoothing pix x = {idx_x1}')
        for idx_y1 in range(0,len(image_data.yaxis)):
            old_content = image_data.zaxis[idx_x1,idx_y1]
            for idx_x2 in range(idx_x1-2*kernel_pix_size,idx_x1+2*kernel_pix_size):
                for idx_y2 in range(idx_y1-2*kernel_pix_size,idx_y1+2*kernel_pix_size):
                    if idx_x2<0: continue
                    if idx_y2<0: continue
                    if idx_x2>=len(image_data.xaxis): continue
                    if idx_y2>=len(image_data.yaxis): continue
                    scale = image_kernel.zaxis[central_bin_x+idx_x2-idx_x1,central_bin_y+idx_y2-idx_y1] 
                    image_smooth.zaxis[idx_x1,idx_y1] += old_content*scale

    print ('done smoothing map.')
    return image_smooth

def calculate_pix_gravity(image_data,src_x,src_y):

    total_gravity = 0.
    src_mass = image_data.get_bin_content(src_x,src_y)
    if src_mass==0.:
        return 0.
    min_distance = image_data.xaxis[1]-image_data.xaxis[0]
    for idx_x in range(0,len(image_data.xaxis)):
        for idx_y in range(0,len(image_data.yaxis)):
            pix_mass = image_data.zaxis[idx_x,idx_y]
            pix_x = image_data.xaxis[idx_x]
            pix_y = image_data.yaxis[idx_y]
            distance = pow(pow(src_x-pix_x,2)+pow(src_y-pix_y,2),0.5)
            distance = max(min_distance,distance)
            gravity = (src_mass*pix_mass)/distance
            total_gravity += gravity

    return total_gravity

def image_cleaning_with_svd(image_data,image_svd_cleaned):

    k = 5
    U_signal, S_signal, V_signal = np.linalg.svd(image_data.zaxis,full_matrices=False)
    image_svd_cleaned.zaxis = U_signal[:, :k] @ np.diag(S_signal[:k]) @ V_signal[:k, :]

def find_point_sources_with_gravity(image_data,image_mask):

    print ('calculating pixel gravity...')
    all_pix_gravity = []
    for idx_x in range(0,len(image_mask.xaxis)):
        for idx_y in range(0,len(image_mask.yaxis)):
            mask = image_mask.zaxis[idx_x,idx_y]
            if mask==1: continue
            pix_x = image_data.xaxis[idx_x]
            pix_y = image_data.yaxis[idx_y]
            gravity = calculate_pix_gravity(image_data,pix_x,pix_y)
            if gravity<=0.: continue
            all_pix_gravity += [np.log10(gravity)]

    mean_gravity = np.mean(all_pix_gravity)
    rms_gravity = 0.
    for entry in range(0,len(all_pix_gravity)):
        rms_gravity += pow(all_pix_gravity[entry]-mean_gravity,2)
    rms_gravity = pow(rms_gravity/len(all_pix_gravity),0.5)
    print (f'mean_gravity = {mean_gravity}')
    print (f'rms_gravity = {rms_gravity}')

    for idx_x in range(0,len(image_mask.xaxis)):
        for idx_y in range(0,len(image_mask.yaxis)):
            mask = image_mask.zaxis[idx_x,idx_y]
            if mask==1: continue
            pix_x = image_data.xaxis[idx_x]
            pix_y = image_data.yaxis[idx_y]
            gravity = calculate_pix_gravity(image_data,pix_x,pix_y)
            if gravity<=0.: continue
            if np.log10(gravity)>2.0*rms_gravity+mean_gravity:
                image_mask.zaxis[idx_x,idx_y] = 1

    max_gravity = max(all_pix_gravity)
    hist_gravity = MyArray1D(bin_start=mean_gravity-2.*rms_gravity,bin_end=mean_gravity+4.*rms_gravity,pixel_scale=6.*rms_gravity/100.)
    for entry in range(0,len(all_pix_gravity)):
        hist_gravity.fill(all_pix_gravity[entry])

    return hist_gravity


def find_point_sources(image_data,image_mask,threshold):

    point_source_threshold = threshold

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
            if significance>point_source_threshold:
                image_mask.zaxis[idx_x,idx_y] = 1

image_det_diff.zaxis = image_det_sci.zaxis - image_det_bkg.zaxis

fft_lofreq_mode = 4

fft_filter(image_det_mask,image_det_diff,image_det_diff_lofreq,0,fft_lofreq_mode)
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
fig.savefig("%s/%s_%s_image_det_diff_lofreq.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

hist_gravity = find_point_sources_with_gravity(image_det_diff_lofreq,image_det_mask)
fig.clf()
axbig = fig.add_subplot()
axbig.plot(hist_gravity.xaxis,hist_gravity.yaxis)
axbig.set_yscale('log')
axbig.set_xlabel('Gravity')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_gravity_lofreq.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fft_filter(image_det_mask,image_det_diff,image_det_diff_hifreq,fft_lofreq_mode,0)
fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_det_diff_hifreq.xaxis.min()
xmax = image_det_diff_hifreq.xaxis.max()
ymin = image_det_diff_hifreq.yaxis.min()
ymax = image_det_diff_hifreq.yaxis.max()
axbig.imshow(image_det_diff_hifreq.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_det_diff_hifreq.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

hist_gravity = find_point_sources_with_gravity(image_det_diff_hifreq,image_det_mask)
fig.clf()
axbig = fig.add_subplot()
axbig.plot(hist_gravity.xaxis,hist_gravity.yaxis)
axbig.set_yscale('log')
axbig.set_xlabel('Gravity')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_gravity_hifreq.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
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
fig.savefig("%s/%s_%s_image_det_mask.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()



image_det_xry = MyArray2D()
image_det_qpb = MyArray2D()
image_det_spf = MyArray2D()
image_det_cxb = MyArray2D()
image_det_gxb = MyArray2D()
image_icrs_xry = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=0.56,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
image_icrs_qpb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=0.56,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
image_icrs_spf = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=0.56,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
image_icrs_cxb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=0.56,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
image_icrs_gxb = MyArray2D(start_x=sky_ra_low,start_y=sky_dec_low,image_size=0.56,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
image_galactic_xry = MyArray2D(start_x=sky_l_low,start_y=sky_b_low,image_size=0.56,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
spectrum_sci = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_det_bkg = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_all_bkg = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_qpb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_spf = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_cxb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
spectrum_gxb = MyArray1D(bin_start=ch_low,bin_end=ch_high,pixel_scale=ch_scale)
detx_sci = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_bkg = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_qpb = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_spf = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_cxb = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
detx_gxb = MyArray1D(bin_start=detx_low,bin_end=detx_high,pixel_scale=detx_scale)
dety_sci = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=dety_scale)
dety_bkg = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=dety_scale)
dety_qpb = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=dety_scale)
dety_spf = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=dety_scale)
dety_cxb = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=dety_scale)
dety_gxb = MyArray1D(bin_start=dety_low,bin_end=dety_high,pixel_scale=dety_scale)
galx_xry = MyArray1D(bin_start=sky_l_low,bin_end=sky_l_high,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
galy_xry = MyArray1D(bin_start=sky_b_low,bin_end=sky_b_high,pixel_scale=map_rebin*detx_scale*0.05/(60.*60.))
skyr_xry = MyArray1D(bin_start=0,bin_end=0.28,pixel_scale=2.*detx_scale*0.05/(60.*60.))
area_curve = MyArray1D(bin_start=2000,bin_end=12000,pixel_scale=2000)
for ccd in range(0,len(ana_ccd_bins)):

    arf_filename = '%s/sci_area_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    arf_hdu_list = fits.open(arf_filename)
    arf_table = Table.read(arf_filename, hdu=1)
    for entry in range(0,len(arf_table)):
        entry_energy = arf_table[entry]['ENERGY']
        entry_area = arf_table[entry]['AREA']
        print (f'energy = {entry_energy}, area = {entry_area}')
        area_curve.fill(entry_energy,weight=entry_area)

    sci_filename = '%s/sci_events_%s_ccd%s.fits'%(input_dir,detector,ana_ccd_bins[ccd])
    sci_hdu_list = fits.open(sci_filename)
    exposure = sci_hdu_list[1].header['obs_duration']
    print (f'exposure = {exposure}')

    fov_size = 0.
    for idx_x in range(0,len(image_det_mask.xaxis)):
        for idx_y in range(0,len(image_det_mask.yaxis)):
            evt_cnt = image_det_sci.zaxis[idx_x,idx_y]
            mask = image_det_mask.zaxis[idx_x,idx_y]
            if evt_cnt==0.: continue
            if mask==1: continue
            fov_size += pow(map_rebin*detx_scale*0.05/(60.*60.),2)/3282.8 # steradian

    sci_table = Table.read(sci_filename, hdu=1)
    for entry in range(0,len(sci_table)):
        evt_pi = sci_table[entry]['PI']
        evt_detx = sci_table[entry]['DETX']
        evt_dety = sci_table[entry]['DETY']
        evt_ra = sci_table[entry]['RA']
        evt_dec = sci_table[entry]['DEC']
        distance_2_roi = DistanceToROI(evt_ra,evt_dec)
        if distance_2_roi<exclusion_inner or distance_2_roi>exclusion_outer: continue
        if evt_pi<energy_cut_lower: continue
        if evt_pi>energy_cut_upper: continue
        if point_source_cut:
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if mask==1: continue
        areatime = ch_scale/1000.*fov_size*exposure*area_curve.yaxis[0]
        spectrum_sci.fill(evt_pi,weight=1./areatime)
        image_det_xry.fill(evt_detx,evt_dety,weight=1./areatime)
        image_icrs_xry.fill(evt_ra,evt_dec,weight=1./areatime)
        detx_sci.fill(evt_detx,weight=1./areatime)
        dety_sci.fill(evt_dety,weight=1./areatime)

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
        if evt_pi<energy_cut_lower: continue
        if evt_pi>energy_cut_upper: continue
        if point_source_cut:
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if mask==1: continue
        areatime = ch_scale/1000.*fov_size*exposure*area_curve.yaxis[0]
        spectrum_qpb.fill(evt_pi,weight=1./areatime*1./sample_scale)
        spectrum_det_bkg.fill(evt_pi,weight=1./areatime*1./sample_scale)
        detx_qpb.fill(evt_detx,weight=1./areatime*1./sample_scale)
        dety_qpb.fill(evt_dety,weight=1./areatime*1./sample_scale)

        image_det_qpb.fill(evt_detx,evt_dety,weight=1./areatime*1./sample_scale)
        image_icrs_qpb.fill(evt_ra,evt_dec,weight=1./areatime*1./sample_scale)

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
        if evt_pi<energy_cut_lower: continue
        if evt_pi>energy_cut_upper: continue
        if point_source_cut:
            mask = image_det_mask.get_bin_content(evt_detx,evt_dety)
            if mask==1: continue
        areatime = ch_scale/1000.*fov_size*exposure*area_curve.yaxis[0]
        spectrum_spf.fill(evt_pi,weight=1./areatime*1./sample_scale)
        spectrum_det_bkg.fill(evt_pi,weight=1./areatime*1./sample_scale)
        detx_spf.fill(evt_detx,weight=1./areatime*1./sample_scale)
        dety_spf.fill(evt_dety,weight=1./areatime*1./sample_scale)

        image_det_spf.fill(evt_detx,evt_dety,weight=1./areatime*1./sample_scale)
        image_icrs_spf.fill(evt_ra,evt_dec,weight=1./areatime*1./sample_scale)

sci_bkg_ratio = []
for ch in range(0,len(energy_array)-1):
    sci_cnt = spectrum_sci.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    bkg_cnt = spectrum_det_bkg.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    if bkg_cnt>0.:
        sci_bkg_ratio += [sci_cnt/bkg_cnt]
    else:
        sci_bkg_ratio += [0.]
print ('sci_bkg_ratio = %s'%(sci_bkg_ratio))

bkg_cnt = spectrum_det_bkg.integral()
spf_cnt = spectrum_spf.integral()
print (f'spf_cnt/bkg_cnt = {spf_cnt/bkg_cnt}')

cxb_measurement_tmp = []
qpb_measurement_tmp = []
for ch in range(0,len(energy_array)-1):
    sci_cnt = spectrum_sci.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    bkg_cnt = spectrum_det_bkg.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    qpb_cnt = spectrum_qpb.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    nbins = (energy_array[ch+1]-energy_array[ch])/ch_scale
    cxb_measurement_tmp += [(sci_cnt-bkg_cnt)/nbins]
    qpb_measurement_tmp += [qpb_cnt/nbins]
print ('cxb_measurement += [%s]'%(cxb_measurement_tmp))

get_cxb_spectrum(spectrum_cxb)

gxb_measurement_tmp = []
for ch in range(0,len(energy_array)-1):
    sci_cnt = spectrum_sci.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    bkg_cnt = spectrum_det_bkg.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    cxb_cnt = spectrum_cxb.integral(integral_range=[energy_array[ch],energy_array[ch+1]])
    nbins = (energy_array[ch+1]-energy_array[ch])/ch_scale
    gxb_measurement_tmp += [(sci_cnt-bkg_cnt-cxb_cnt)/nbins]
print ('gxb_measurement += [%s]'%(gxb_measurement_tmp))

get_gxb_spectrum(spectrum_gxb)

spectrum_all_bkg.add(spectrum_det_bkg)
spectrum_all_bkg.add(spectrum_cxb)
spectrum_all_bkg.add(spectrum_gxb)

cxb_sum = spectrum_cxb.integral()
gxb_sum = spectrum_gxb.integral()
qpb_sum = spectrum_qpb.integral()
cxb_2_qpb_ratio = cxb_sum/qpb_sum
gxb_2_qpb_ratio = gxb_sum/qpb_sum
image_det_cxb.add(image_det_qpb,factor=cxb_2_qpb_ratio)
image_icrs_cxb.add(image_icrs_qpb,factor=cxb_2_qpb_ratio)
detx_cxb.add(detx_qpb,factor=cxb_2_qpb_ratio)
dety_cxb.add(dety_qpb,factor=cxb_2_qpb_ratio)
image_det_gxb.add(image_det_qpb,factor=gxb_2_qpb_ratio)
image_icrs_gxb.add(image_icrs_qpb,factor=gxb_2_qpb_ratio)
detx_gxb.add(detx_qpb,factor=gxb_2_qpb_ratio)
dety_gxb.add(dety_qpb,factor=gxb_2_qpb_ratio)

image_det_xry.add(image_det_qpb,factor=-1.)
image_det_xry.add(image_det_spf,factor=-1.)
image_det_xry.add(image_det_cxb,factor=-1.)
image_det_xry.add(image_det_gxb,factor=-1.)
image_icrs_xry.add(image_icrs_qpb,factor=-1.)
image_icrs_xry.add(image_icrs_spf,factor=-1.)
image_icrs_xry.add(image_icrs_cxb,factor=-1.)
image_icrs_xry.add(image_icrs_gxb,factor=-1.)
detx_bkg.add(detx_qpb)
detx_bkg.add(detx_spf)
detx_bkg.add(detx_cxb)
detx_bkg.add(detx_gxb)
dety_bkg.add(dety_qpb)
dety_bkg.add(dety_spf)
dety_bkg.add(dety_cxb)
dety_bkg.add(dety_gxb)

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_det_sci_fine.xaxis.min()
xmax = image_det_sci_fine.xaxis.max()
ymin = image_det_sci_fine.yaxis.min()
ymax = image_det_sci_fine.yaxis.max()
#im = axbig.imshow(image_det_sci_fine.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax), norm=colors.LogNorm())
im = axbig.imshow(image_det_sci_fine.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
#cbar = fig.colorbar(im)
#cbar.set_label('Count')
fig.savefig("%s/%s_%s_image_det_sci.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()


fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_det_xry.xaxis.min()
xmax = image_det_xry.xaxis.max()
ymin = image_det_xry.yaxis.min()
ymax = image_det_xry.yaxis.max()
axbig.imshow(image_det_xry.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_det_xry.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'RA'
label_y = 'DEC'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_icrs_xry.xaxis.min()
xmax = image_icrs_xry.xaxis.max()
ymin = image_icrs_xry.yaxis.min()
ymax = image_icrs_xry.yaxis.max()
axbig.imshow(image_icrs_xry.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_icrs_xry.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

ConvertRaDecMapToGalacticMap(image_icrs_xry,image_galactic_xry)
fig.clf()
axbig = fig.add_subplot()
label_x = 'Gal. l'
label_y = 'Gal. b'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_galactic_xry.xaxis.min()
xmax = image_galactic_xry.xaxis.max()
ymin = image_galactic_xry.yaxis.min()
ymax = image_galactic_xry.yaxis.max()
axbig.imshow(image_galactic_xry.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_galactic_xry.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
label_x = 'DETY'
label_y = 'DETX'
axbig.set_xlabel(label_x)
axbig.set_ylabel(label_y)
xmin = image_det_xry.xaxis.min()
xmax = image_det_xry.xaxis.max()
ymin = image_det_xry.yaxis.min()
ymax = image_det_xry.yaxis.max()
axbig.imshow(image_det_spf.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))
fig.savefig("%s/%s_%s_image_det_spf.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(spectrum_sci.xaxis,spectrum_sci.yaxis,yerr=spectrum_sci.yerr,color='k',label='Data')
axbig.plot(spectrum_qpb.xaxis,spectrum_qpb.yaxis,color='purple',label='QPB')
axbig.plot(spectrum_spf.xaxis,spectrum_spf.yaxis,color='blue',label='SPF')
axbig.plot(spectrum_cxb.xaxis,spectrum_cxb.yaxis,color='green',label='CXB')
axbig.plot(spectrum_gxb.xaxis,spectrum_gxb.yaxis,color='orange',label='GDE')
axbig.plot(spectrum_all_bkg.xaxis,spectrum_all_bkg.yaxis,color='red',label='Bkg')
#axbig.set_yscale('log')
#axbig.set_ylim(bottom=1)
axbig.set_xlabel('Energy [eV]')
axbig.set_ylabel('Photons /cm2/s/sr/keV')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_spectrum_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(detx_sci.xaxis,detx_sci.yaxis,yerr=detx_sci.yerr,color='k',label='Data')
axbig.plot(detx_qpb.xaxis,detx_qpb.yaxis,color='purple',label='QPB')
axbig.plot(detx_spf.xaxis,detx_spf.yaxis,color='blue',label='SPF')
axbig.plot(detx_cxb.xaxis,detx_cxb.yaxis,color='green',label='CXB')
axbig.plot(detx_gxb.xaxis,detx_gxb.yaxis,color='orange',label='GDE')
axbig.plot(detx_bkg.xaxis,detx_bkg.yaxis,color='red',label='Bkg')
#axbig.set_yscale('log')
#axbig.set_ylim(bottom=1)
axbig.set_xlabel('DETX')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_detx_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(dety_sci.xaxis,dety_sci.yaxis,yerr=dety_sci.yerr,color='k',label='Data')
axbig.plot(dety_qpb.xaxis,dety_qpb.yaxis,color='purple',label='QPB')
axbig.plot(dety_spf.xaxis,dety_spf.yaxis,color='blue',label='SPF')
axbig.plot(dety_cxb.xaxis,dety_cxb.yaxis,color='green',label='CXB')
axbig.plot(dety_gxb.xaxis,dety_gxb.yaxis,color='orange',label='GDE')
axbig.plot(dety_bkg.xaxis,dety_bkg.yaxis,color='red',label='Bkg')
#axbig.set_yscale('log')
#axbig.set_ylim(bottom=1)
axbig.set_xlabel('DETY')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_dety_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()


MakeGalacticProjection(image_galactic_xry, galx_xry, galy_xry)

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(galx_xry.xaxis,galx_xry.yaxis,yerr=galx_xry.yerr,color='k',label='Data')
axbig.set_xlabel('Gal. l')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_galx_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(galy_xry.xaxis,galy_xry.yaxis,yerr=galy_xry.yerr,color='k',label='Data')
axbig.set_xlabel('Gal. b')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_galy_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

MakeRadialProjection(image_det_xry, skyr_xry)

fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(skyr_xry.xaxis,skyr_xry.yaxis,yerr=skyr_xry.yerr,color='k',label='Data')
axbig.set_xlabel('Angular distance [deg]')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_skyr_model.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
axbig.plot(area_curve.xaxis,area_curve.yaxis,color='k')
axbig.set_xlabel('Energy [eV]')
axbig.set_ylabel('Area [cm2]')
axbig.legend(loc='best')
fig.savefig("%s/%s_%s_area_curve.png"%(output_dir,on_obsID,detector),bbox_inches='tight')
axbig.remove()


energy_threshold = 2000
energy_axis = []
res_count_axis = []
res_count_error = []
for entry in range(0,len(spectrum_sci.xaxis)):
    sci_entry_energy = spectrum_sci.xaxis[entry]
    sci_entry_count = spectrum_sci.yaxis[entry]
    bkg_entry_count = spectrum_all_bkg.yaxis[entry]
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

if do_fit:
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
