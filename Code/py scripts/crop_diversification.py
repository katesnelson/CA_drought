import numpy as np
import matplotlib as plt
from __future__ import division
import gdal
from gdalconst import *
gdal.UseExceptions()
import numpy as np

d = '/data/emily/WF/kate'
#huc12
huc12 = gdal.Open(d + '/cv_huc.tif').ReadAsArray()
#landuse
#make sure reprojected with proper rows and columns
lulc = gdal.Open(d + '/cv_lulc.tif').ReadAsArray()
fl = gdal.Open("/data/emily/WF/kate/fl.tif", GA_ReadOnly).ReadAsArray()
ag_lulc = gdal.Open(d + '/cv_lulc/cv_ag_lulc.tif', GA_ReadOnly).ReadAsArray()




#crop diversification index
huc_idx = np.unique(huc12)
lulc_idx = np.unique(lulc)

cd_im = np.zeros(shape = lulc.shape)

for i in range(lulc.shape[0]):
    lc = lulc[i,:,:]    

    for h in range(len(np.unique(huc12))-1):
        data = lc[huc12 == huc_idx[h]]
        cdt = 0
        for l in range(len(np.unique(lulc))):
            total_area = data.shape[0]
            lulc_count = len(data[data == lulc_idx[l]])
            if lulc_count == 0:
                out == 0
            else:
                out = -(lulc_count/total_area)*np.log(lulc_count/total_area)
            cdt += out
        cd_im[i,:,:][huc12 == huc_idx[h]] = np.repeat(cdt, len(data))
        
cd_im[cd_im == 0] = -9999.0

#farmland
cd_im_fl = np.zeros(shape = lulc.shape)

for i in range(lulc.shape[0]):
    lc = lulc[i,:,:]    

    for h in range(len(np.unique(huc12))-1):
        lc_sub = lc[huc12 == huc_idx[h]]
        fl_sub = fl[4,:,:][huc12 == huc_idx[h]]
        data = lc_sub[fl_sub == 1]
        
        cdt = 0
        for l in range(len(np.unique(lulc))):
            total_area = data.shape[0]
            lulc_count = len(data[data == lulc_idx[l]])
            if lulc_count == 0:
                out == 0
            else:
                out = -(lulc_count/total_area)*np.log(lulc_count/total_area)
            cdt += out
        cd_im_fl[i,:,:][huc12 == huc_idx[h]] = np.repeat(cdt, len(lc_sub))
        
cd_im_fl[cd_im_fl == 0] = -9999.0

infile = '/data/emily/WF/kate/cv.tif'
im = gdal.Open(infile).ReadAsArray()

#tvp_im
for i in range(cd_im.shape[0]):
    arr = cd_im[i,:,:]
    nodatav = -9999.0
    data = gdal.Open(infile)

    [cols, rows] = arr.shape
    trans = data.GetGeoTransform()
    proj = data.GetProjection()
    outfile = '/data/emily/WF/kate/cd_cv/fl_cd' + str(i+7) + '.tif'

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
      
        
