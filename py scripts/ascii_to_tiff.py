# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 15:35:31 2016

@author: Emily Burchfield
"""
import numpy as np
from osgeo import osr, gdal
from gdal import *
import matplotlib.mlab as ml
import glob
import subprocess


#XYZ (long, lat, data) format ascii file without header information
#first load xyz into GIS to construct raster from point centroids

d = 'E:\\SPI_NLDAS'
file_list = sorted(glob.glob(d + '/*.asc'))
refIm_cv = 'C:\\Users\\Emily Burchfield\\Box Sync\\GIS\\USA\\CA\\cv\\cv.tif'
refIm = 'E:\\SPI_NLDAS\\refIm.tif'
extent = 'C:\\Users\\Emily Burchfield\\Box Sync\\GIS\\USA\\CA\\cv\\cv_extent.shp'

def ascii_to_tiff(infile, outfile, refIm):
    im = gdal.Open(refIm)
    ima = gdal.Open(refIm).ReadAsArray()
    row = ima.shape[0]; col = ima.shape[1]
     
    indata = np.genfromtxt(infile, delimiter=",", skip_header = True, dtype=None)
    lon = indata[:,0] #x
    lat = indata[:,1] #y
    data = indata[:,2]

    #create grid
    xmin, xmax, ymin, ymax = [min(lon), max(lon), min(lat), max(lat)]
    xi = np.linspace(xmin, xmax, col)
    yi = np.linspace(ymin, ymax, row)
    xi, yi = np.meshgrid(xi, yi)

    #linear interpolation
    zi = ml.griddata(lon, lat, data, xi, yi, interp = 'linear')
    final_array = np.asarray(np.rot90(np.transpose(zi)))

    #projection
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(outfile, col, row, 1, gdal.GDT_Float32)
    dst_ds.GetRasterBand(1).WriteArray(final_array)
    prj = im.GetProjection()
    dst_ds.SetProjection(prj)
    
    gt = im.GetGeoTransform()
    dst_ds.SetGeoTransform(gt)
    dst_ds = None

    final_tif = gdal.Open(outfile, GA_ReadOnly).ReadAsArray()
    return final_tif    

  
for i in range(len(file_list)):
    infile = file_list[i]
    outfile = infile[:-3] + "tif"
    ascii_to_tiff(infile, outfile, refIm)
    final_tiff = outfile[:-4] + "_ca.tif"
    #crop to extent
    subprocess.call(['gdalwarp', '-r', 'near', '-cutline', extent, '-crop_to_cutline', outfile, final_tiff, '-dstnodata', '9999'])

#create datacube
dc = []
ca_list = sorted(glob.glob(d + '\\*ca.tif'))

for i in range(len(ca_list)):
    x = gdal.Open(ca_list[i]).ReadAsArray()
    dc.append(x)
    
pdsi = np.asarray(dc)
pdsim = np.ma.masked_outside(pdsi, -100, 100)
#create annual values
year = np.arange(0,spi.shape[0], 12)

pdsi_y = []
for i in range(len(year)):
    y = pdsim[i:i+12, :,:]
    ysum = np.ma.sum(y, axis=0)
    pdsi_y.append(ysum)

final = np.asarray(pdsi_y)  
final[final == 0] = -9999.0  
