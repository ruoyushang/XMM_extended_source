
import numpy as np
import csv
from numpy.linalg import inv
import math

pattern_low = 0
pattern_high = 6
pattern_scale = 1
ch_low = 200
ch_high = 12000
ch_scale = 100
t_low = 0
t_high = 1
t_scale = 0.01
detx_low = -20000
detx_high = 20000
detx_scale = 200
dety_low = -20000
dety_high = 20000
dety_scale = 200
detr_low = 0
detr_high = 20000
detr_scale = 200
sky_ra_low = -0.28
sky_ra_high = 0.28
sky_dec_low = -0.28
sky_dec_high = 0.28
sky_scale = 500.*0.05/(60.*60.)


class MyArray2D:

    def __init__(self,start_x=-19474.5,start_y=-19474.5,image_size=39000,pixel_scale=200):
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
    def fill(self, value_x, value_y, weight=1.):
        key_idx_x = -1
        key_idx_y = -1
        for idx_x in range(0,len(self.xaxis)-1):
            if self.xaxis[idx_x]<=value_x and self.xaxis[idx_x+1]>value_x:
                key_idx_x = idx_x
        for idx_y in range(0,len(self.yaxis)-1):
            if self.yaxis[idx_y]<=value_y and self.yaxis[idx_y+1]>value_y:
                key_idx_y = idx_y
        if key_idx_x>=0 and key_idx_y>=0:
            self.zaxis[key_idx_x,key_idx_y] += 1.*weight
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
    def calc_ratio(self, src_array, bkg_array):
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                src = src_array.zaxis[idx_x,idx_y]
                bkg = bkg_array.zaxis[idx_x,idx_y]
                if not src>0.: continue
                if not bkg>0.: continue
                ratio = (src-bkg)/bkg
                self.zaxis[idx_x,idx_y] = ratio
    def flattening(self):
        pix_on_cnt = 0
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                if self.zaxis[idx_x,idx_y]!=0:
                    pix_on_cnt += 1
        avg_cnt = self.integral()/float(pix_on_cnt)
        for idx_x in range(0,len(self.xaxis)):
            for idx_y in range(0,len(self.yaxis)):
                if self.zaxis[idx_x,idx_y]!=0:
                    self.zaxis[idx_x,idx_y] = avg_cnt

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
    def fill(self, value, weight=1.):
        for entry in range(0,len(self.xaxis)-1):
            if self.xaxis[entry]<=value and self.xaxis[entry+1]>value:
                self.yaxis[entry] += 1.*weight
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
            if value_x>=self.xaxis[idx_x]:
                key_idx_x = idx_x
        return self.yaxis[key_idx_x]
    def normalize(self):
        reference = self.yaxis[0]
        for entry in range(0,len(self.yaxis)):
            self.yaxis[entry] = self.yaxis[entry]/reference
            self.yerr[entry] = self.yerr[entry]/reference
    def normalize2(self):
        for entry in range(0,len(self.yaxis)):
            if self.yaxis[entry]<0.:
                self.yaxis[entry] = 0.
        reference = self.integral()
        for entry in range(0,len(self.yaxis)):
            if reference>0.:
                self.yaxis[entry] = self.yaxis[entry]/reference
                self.yerr[entry] = self.yerr[entry]/reference
            else:
                self.yaxis[entry] = 0.
                self.yerr[entry] = 0.
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

def find_nearest_ref_det_idx(target_det,ref_det_list):

    min_dist = 1e10
    min_dist_idx_x = 0
    for idx_x in range(0,10):
        dist = abs(target_det[0]-ref_det_list[idx_x][0])
        if min_dist>dist:
            min_dist = dist
            min_dist_idx_x = idx_x

    min_dist = 1e10
    min_dist_idx_y = 0
    for idx_y in range(0,10):
        dist = abs(target_det[1]-ref_det_list[idx_y][1])
        if min_dist>dist:
            min_dist = dist
            min_dist_idx_y = idx_y

    return min_dist_idx_x, min_dist_idx_y

def find_nearest_ref_sky_idx(target_sky,ref_sky_list):

    min_dist = 1e10
    min_dist_idx_x = 0
    for idx_x in range(0,len(ref_sky_list)):
        dist = abs(target_sky[0]-ref_sky_list[idx_x][0])
        if min_dist>dist:
            min_dist = dist
            min_dist_idx_x = idx_x

    min_dist = 1e10
    min_dist_idx_y = 0
    for idx_y in range(0,len(ref_sky_list)):
        dist = abs(target_sky[1]-ref_sky_list[idx_y][1])
        if min_dist>dist:
            min_dist = dist
            min_dist_idx_y = idx_y

    return min_dist_idx_x, min_dist_idx_y


def LoadCoordinateMatrix(idx_ra,idx_dec,on_sample,on_obsID,detector):

    idx_ra_offset1 = idx_ra+1
    idx_dec_offset1 = idx_dec+0
    idx_ra_offset2 = idx_ra+0
    idx_dec_offset2 = idx_dec+1

    origin_filename = '../%s/%s/analysis/sky_det_ref/%s_esky2det_ra_%s_dec_%s.txt'%(on_sample,on_obsID,detector,idx_ra,idx_dec)
    offset1_filename = '../%s/%s/analysis/sky_det_ref/%s_esky2det_ra_%s_dec_%s.txt'%(on_sample,on_obsID,detector,idx_ra_offset1,idx_dec_offset1)
    offset2_filename = '../%s/%s/analysis/sky_det_ref/%s_esky2det_ra_%s_dec_%s.txt'%(on_sample,on_obsID,detector,idx_ra_offset2,idx_dec_offset2)
    #print ('read file: %s'%(origin_filename))

    origin_det = [0,0]
    origin_sky = [0,0]
    offset1_det = [0,0]
    offset1_sky = [0,0]
    offset2_det = [0,0]
    offset2_sky = [0,0]

    origin_file = open(origin_filename)
    last_line = False
    for line in origin_file:
        if last_line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            origin_det[0] = float(new_text_list[0])
            origin_det[1] = float(new_text_list[1].strip('\n'))
            break
        if 'Source RA' in line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            origin_sky[0] = float(new_text_list[4])
        if 'Source dec' in line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            origin_sky[1] = float(new_text_list[4])
        if '# detX       detY' in line:
            last_line = True

    #print ('origin_det = %s'%(origin_det))
    #print ('origin_sky = %s'%(origin_sky))

    offset1_file = open(offset1_filename)
    last_line = False
    for line in offset1_file:
        if last_line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            offset1_det[0] = float(new_text_list[0])
            offset1_det[1] = float(new_text_list[1].strip('\n'))
            break
        if 'Source RA' in line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            offset1_sky[0] = float(new_text_list[4])
        if 'Source dec' in line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            offset1_sky[1] = float(new_text_list[4])
        if '# detX       detY' in line:
            last_line = True

    #print ('offset1_det = %s'%(offset1_det))
    #print ('offset1_sky = %s'%(offset1_sky))

    offset2_file = open(offset2_filename)
    last_line = False
    for line in offset2_file:
        if last_line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            offset2_det[0] = float(new_text_list[0])
            offset2_det[1] = float(new_text_list[1].strip('\n'))
            break
        if 'Source RA' in line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            offset2_sky[0] = float(new_text_list[4])
        if 'Source dec' in line:
            text_list = line.split(' ')
            new_text_list = []
            for entry in text_list:
                if entry!='':
                    new_text_list += [entry]
            offset2_sky[1] = float(new_text_list[4])
        if '# detX       detY' in line:
            last_line = True

    #print ('offset2_det = %s'%(offset2_det))
    #print ('offset2_sky = %s'%(offset2_sky))

    mtx_delta_sky = np.matrix([ [offset1_sky[0]-origin_sky[0], offset1_sky[1]-origin_sky[1]] , [offset2_sky[0]-origin_sky[0], offset2_sky[1]-origin_sky[1]] ])
    mtx_delta_det = np.matrix([ [offset1_det[0]-origin_det[0], offset1_det[1]-origin_det[1]] , [offset2_det[0]-origin_det[0], offset2_det[1]-origin_det[1]] ])

    inv_mtx_delta_sky = inv(mtx_delta_sky)
    mtx_conv_sky_2_det = np.matmul(mtx_delta_det,inv_mtx_delta_sky)
    inv_mtx_delta_det = inv(mtx_delta_det)
    mtx_conv_det_2_sky = np.matmul(mtx_delta_sky,inv_mtx_delta_det)

    return origin_sky, origin_det, mtx_conv_det_2_sky, mtx_conv_sky_2_det

def HMS2deg(ra='', dec=''):
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split(':')]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    if ra:
        H, M, S = [float(i) for i in ra.split(':')]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC

def ReadHAWCTargetListFromFile():
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open('Cat_3HWC.txt')
    for line in inputFile:
        if line[0]=="#": continue
        if '- name:' in line:
            target_name = line.lstrip('   - name: ')
            target_name = target_name.strip('\n')
        if 'RA:' in line:
            target_ra = line.lstrip('     RA: ')
        if 'Dec:' in line:
            target_dec = line.lstrip('     Dec: ')
        if 'flux measurements:' in line:
            source_name += [target_name]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
            target_ra = ''
            target_dec = ''
    return source_name, source_ra, source_dec

def ReadSNRTargetListFromCSVFile():
    source_name = []
    source_ra = []
    source_dec = []
    source_size = []
    with open('SNRcat20221001-SNR.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            if len(row)==0: continue
            if '#' in row[0]: continue
            target_name = row[0]
            target_min_dist = row[13]
            if target_min_dist=='':
                target_min_dist = '1000'
            #if float(target_min_dist)>6.: continue
            target_size = row[15]
            if target_size=='':
                target_size = 0.
            target_ra = row[19]
            target_dec = row[20]
            source_name += [target_name]
            source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
            source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
            source_size += [0.5*float(target_size)/60.]
            #print('target_min_dist = %s'%(target_min_dist))
            #print('source_name = %s'%(source_name[len(source_name)-1]))
            #print('source_ra = %0.2f'%(source_ra[len(source_ra)-1]))
            #print('source_dec = %0.2f'%(source_dec[len(source_dec)-1]))
            #print(row)
    return source_name, source_ra, source_dec, source_size

def ReadATNFTargetListFromFile():
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_age = []
    source_edot = []
    inputFile = open('ATNF_pulsar_full_list.txt')
    for line in inputFile:
        if line[0]=="#": continue
        target_name = line.split(',')[0].strip(" ")
        if target_name=="\n": continue
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        target_dist = line.split(',')[3].strip(" ")
        target_age = line.split(',')[4].strip(" ")
        target_edot = line.split(',')[5].strip(" ")
        if target_dist=='*': continue
        if target_age=='*': continue
        if target_edot=='*': continue
        target_brightness = float(target_edot)/pow(float(target_dist),2)

        if float(target_edot)<1e35: continue
        #if float(target_dist)>6.: continue

        #ra_deg = float(HMS2deg(target_ra,target_dec)[0])
        #dec_deg = float(HMS2deg(target_ra,target_dec)[1])
        #gal_l, gal_b = ConvertRaDecToGalactic(ra_deg,dec_deg)
        #if abs(gal_b)<5.: continue

        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
        source_dist += [float(target_dist)]
        source_age += [float(target_age)]
        source_edot += [float(target_edot)]
    return source_name, source_ra, source_dec, source_dist, source_age

def DrawSkyMap(fig,map_color,image_data,save_name):

    xmin = image_data.xaxis.min()
    xmax = image_data.xaxis.max()
    ymin = image_data.yaxis.min()
    ymax = image_data.yaxis.max()

    target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age = ReadATNFTargetListFromFile()
    target_snr_name, target_snr_ra, target_snr_dec, target_snr_size = ReadSNRTargetListFromCSVFile()
    target_hwc_name, target_hwc_ra, target_hwc_dec = ReadHAWCTargetListFromFile()

    star_range = 0.8*(xmax-xmin)/2.
    source_ra = (xmax+xmin)/2.
    source_dec = (ymax+ymin)/2.
    psr_markers = []
    for star in range(0,len(target_psr_name)):
        if abs(source_ra-target_psr_ra[star])>star_range: continue
        if abs(source_dec-target_psr_dec[star])>star_range: continue
        psr_markers += [[target_psr_ra[star],target_psr_dec[star],target_psr_name[star]]]
    snr_markers = []
    for star in range(0,len(target_snr_name)):
        if abs(source_ra-target_snr_ra[star])>star_range: continue
        if abs(source_dec-target_snr_dec[star])>star_range: continue
        snr_markers += [[target_snr_ra[star],target_snr_dec[star],target_snr_name[star]]]
    hwc_markers = []
    for star in range(0,len(target_hwc_name)):
        if abs(source_ra-target_hwc_ra[star])>star_range: continue
        if abs(source_dec-target_hwc_dec[star])>star_range: continue
        hwc_markers += [[target_hwc_ra[star],target_hwc_dec[star],target_hwc_name[star]]]

    fig.clf()
    axbig = fig.add_subplot()
    label_x = 'RA'
    label_y = 'DEC'
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    axbig.imshow(image_data.zaxis[:,:],origin='lower',cmap=map_color,extent=(xmin,xmax,ymin,ymax))

    arrowprops = dict(arrowstyle="->")
    for star in range(0,len(psr_markers)):
        marker_size = 60
        axbig.scatter(psr_markers[star][0], psr_markers[star][1], s=1.5*marker_size, c='k', marker='+')
        axbig.annotate(psr_markers[star][2], xy=(psr_markers[star][0]+0.01, psr_markers[star][1]+0.01), xytext=(psr_markers[star][0]+0.05, psr_markers[star][1]+0.05), arrowprops=arrowprops)
    for star in range(0,len(snr_markers)):
        marker_size = 60
        axbig.scatter(snr_markers[star][0], snr_markers[star][1], s=1.5*marker_size, c='k', marker='+')
        axbig.annotate(snr_markers[star][2], xy=(snr_markers[star][0]+0.01, snr_markers[star][1]+0.01), xytext=(snr_markers[star][0]+0.05, snr_markers[star][1]+0.05), arrowprops=arrowprops)
    for star in range(0,len(hwc_markers)):
        marker_size = 60
        axbig.scatter(hwc_markers[star][0], hwc_markers[star][1], s=1.5*marker_size, c='k', marker='+')
        axbig.annotate(hwc_markers[star][2], xy=(hwc_markers[star][0]+0.01, hwc_markers[star][1]+0.01), xytext=(hwc_markers[star][0]+0.05, hwc_markers[star][1]+0.05), arrowprops=arrowprops)

    fig.savefig("%s"%(save_name),bbox_inches='tight')
    axbig.remove()

