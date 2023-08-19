
import numpy as np
import math

pattern_low = 0
pattern_high = 6
pattern_scale = 1
ch_low = 200
ch_high = 12000
ch_scale = 100
t_low = 0
t_high = 1
t_scale = 0.05
detx_low = -20000
detx_high = 20000
detx_scale = 500
detr_low = 0
detr_high = 20000
detr_scale = 500
sky_ra_low = 0
sky_ra_high = 0.5
sky_dec_low = 0
sky_dec_high = 0.5
sky_scale = 0.001


class MyArray2D:

    def __init__(self,coord='det',pixel_scale=200):
        start_x = -19474.5
        start_y = -19474.5
        image_size = 39000
        if coord=='radec':
            start_x = sky_ra_low
            start_y = sky_dec_low
            pixel_scale = sky_scale
            image_size = 0.5
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
        key_idx_x = 0
        key_idx_y = 0
        for idx_x in range(0,len(self.xaxis)-1):
            if self.xaxis[idx_x]<=value_x and self.xaxis[idx_x+1]>value_x:
                key_idx_x = idx_x
        for idx_y in range(0,len(self.yaxis)-1):
            if self.yaxis[idx_y]<=value_y and self.yaxis[idx_y+1]>value_y:
                key_idx_y = idx_y
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
            if self.xaxis[idx_x]<=value_x and self.xaxis[idx_x+1]>value_x:
                key_idx_x = idx_x
        return self.yaxis[key_idx_x]
    def normalize(self):
        reference = self.yaxis[0]
        for entry in range(0,len(self.yaxis)):
            self.yaxis[entry] = self.yaxis[entry]/reference
            self.yerr[entry] = self.yerr[entry]/reference
    def normalize2(self):
        reference = self.integral()
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

