from __future__ import division
import numpy as np
import matplotlib as plt
import scipy
from scipy import integrate, interpolate, signal
import scipy.stats as stats
import gdal
from gdalconst import *
gdal.UseExceptions()
import osr


#for future use, don't forget to download more than desired time period of interest b/c of int nan problem

fl = gdal.Open("/data/emily/WF/kate/fl.tif", GA_ReadOnly).ReadAsArray()

nanvalue = -9999
#sg settings
ws = 11
o = 6 
#data size
r = 604
c = 508
z515 = 253
z4 = 23
z16 = 2
z = z515 + z4 + z16

m4 = np.load('/home/emily/cv/2004/MOD13A2.005.npy').reshape((z4, r, c))
m16 = np.load('/home/emily/cv/2016/MOD13A2.005.npy').reshape((z16, r, c))
m515 = np.load('/data/emily/WF/kate/cv/MOD13A2.005.npy').reshape((z515, r, c))  

evi = np.concatenate((m4, m515, m16), axis =0)
evi2d = np.transpose(evi.reshape((z, r*c)))
del m4, m16, m515

#drop poor quality pixels
qualityCount = []  
evim = np.ma.masked_where(evi2d == evi2d.max(), evi2d)
evif = evim.filled(nanvalue) 

for rw in range(evif.shape[0]):
    pix = evif[rw,:]
    x = np.where(pix != nanvalue)   
    cnt = len(x[0])
    qualityCount.append(cnt)
    
q = np.asarray(qualityCount)
index_q = np.where(q>140) #more than 50% of pixels are good quality
qf_index = np.where(q<140)#less than 50% are good quality
evidc = evif[index_q,:].reshape((len(index_q[0]), evif.shape[1]))  
del evif

evi_mask = np.ma.masked_where(evidc == nanvalue, evidc)  #masked adaptive pixel subset
            

##interpolation
#evi_int = np.empty((evidc.shape[0], evidc.shape[1])) 
#for r in range(evidc.shape[0]): #loop through adaptive pixels
    #pixel = evi_mask[r,:]
    #x_ = (np.arange(len(pixel)))[~pixel.mask] #index of pixels with good data
    #y_ = pixel[~pixel.mask]  #pulls out good values
    #f = scipy.interpolate.interp1d(x_, y_, kind = 'linear', bounds_error = False)
    #x = np.arange(0,z,1)
    #ynew = f(x)
    #evi_int[r,:] = ynew
    
##drop edges
#evi_edge = evi_int[:,23:276] #drop 2004 and 16
#col_mean = stats.nanmean(evi_edge, axis = 0)
#ix = np.where(np.isnan(evi_edge))
#evi_edge[ix] = np.take(col_mean,ix[1])


##sg filter
#evi_sg = np.empty((evi_edge.shape[0], evi_edge.shape[1]))

#for rw in range(evidc.shape[0]):
    #pixel = evi_edge[rw,:]
    #sg = scipy.signal.savgol_filter(pixel, ws, o)
    #evi_sg[r,:] = sg

evi_sg = np.load('/data/emily/WF/kate/evisg_cv.npy')

#drop negatives for integration, no tvp when negative
evi_sg[evi_sg < 0] = 0

#integral
tvp = []
ix = range(0,253, 23)

for i in range(len(ix)):
    y = evi_sg[:,ix[i]:ix[i]+23]
    iy = []
    for rw in range(y.shape[0]):
        fxn = y[rw,:]
        out = integrate.simps(fxn, dx=.1, axis=0)
        iy.append(out)

    ifull = np.asarray(iy)
    tvp.append(ifull)

#from data to tiff
tvp = np.asarray(tvp)
tvp_im = []
for i in range(11):
    full = np.zeros((r*c))
    full[index_q] = tvp[i, :]  
    full[qf_index] = nanvalue
    full = full.reshape((r,c))
    tvp_im.append(full)
      
tvp_im = np.asarray(tvp_im)

#extract mean/std by tiff
huc = gdal.Open('/data/emily/WF/kate/cv_huc.tif').ReadAsArray()
#huc = np.tile(huc, (11,1,1))
ix = np.unique(huc) #65535 is nan

tvp_mn = np.zeros((11, r, c))
tvp_sd = np.zeros((11, r, c))

tvpm = np.ma.masked_where(tvp_im == nanvalue, tvp_im)

for y in range(11):
    for h in range(len(np.unique(huc))-1):
        data = tvpm[y,:,:][huc == ix[h]]
        if all(data == True):
            mn = nanvalue
            sd = nanvalue
        else:
            mn = np.nanmean(data)
            sd = np.nanstd(data)
        tvp_mn[y,:,:][huc == ix[h]] = mn
        tvp_sd[y,:,:][huc == ix[h]] = sd

tvp_mn[tvp_mn == 0] = nanvalue
tvp_sd[tvp_sd == 0] = nanvalue

#note, 492, 283 and 1095 hucs dropped b/c too much missing data

#fl only
tvp_mn_fl = np.zeros((11, r, c))
tvp_sd_fl = np.zeros((11, r, c))

tvpm = np.ma.masked_where(tvp_im == nanvalue, tvp_im)

for y in range(11):
    for h in range(len(np.unique(huc))-1):
        evi_sub = tvpm[y,:,:][huc == ix[h]]
        fl_sub = fl[4,:,:][huc == ix[h]]
        data = evi_sub[fl_sub == 1]
        if all(data == True):
            mn = nanvalue
            sd = nanvalue
        else:
            mn = np.nanmean(data)
            sd = np.nanstd(data)
        tvp_mn_fl[y,:,:][huc == ix[h]] = mn
        tvp_sd_fl[y,:,:][huc == ix[h]] = sd

tvp_mn_fl[tvp_mn_fl == 0] = nanvalue
tvp_sd_fl[tvp_sd_fl == 0] = nanvalue


#to tiffs

#drop water, which has 0 tvp

#tvp_im[tvp_im == 0] = nanvalue

infile = '/data/emily/WF/kate/cv.tif'
im = gdal.Open(infile).ReadAsArray()

#tvp_im
for i in range(11):
    arr = tvp_sd_fl[i,:,:]
    nodatav = -9999.0
    data = gdal.Open(infile)

    [cols, rows] = arr.shape
    trans = data.GetGeoTransform()
    proj = data.GetProjection()
    outfile = '/data/emily/WF/kate/cv_tvp/sd_tvp_fl' + str(i+5) + '.tif'

    outdriver = gdal.GetDriverByName('GTiff')
    outdata = outdriver.Create(str(outfile), rows, cols, 1, gdal.GDT_Float32)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outdata.SetProjection(srs.ExportToWkt())
    outdata.SetGeoTransform(trans)
    outdata.SetProjection(proj)
     
    #write array to file
    outdata.GetRasterBand(1).WriteArray(arr)
    outdata.GetRasterBand(1).SetNoDataValue(nodatav)
    outdata = None

    test = gdal.Open(outfile, GA_ReadOnly).ReadAsArray()
