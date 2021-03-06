####
#### Surface differencing
#### @author Chris 01/05/14
####

#
#NB/
#Rows and columns and x origin are a digit out once read 
#in with gdal compared to the dem_par file.....
#

import os
import time as time 
import numpy as np 
from osgeo import gdal, gdalconst 
from osgeo.gdalconst import * 
from matplotlib import pyplot as plt 
from scipy import signal 
from scipy import ndimage 
from matplotlib import cm 

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

#~~~~~~~~~~~~~
#Functions
#~~~~~~~~~~~~~


def open_file_gdal(file_name):
	inds = gdal.Open(file_name, GA_ReadOnly)

	if inds is None:
		print "Really sorry Sir but I couldn't open this blasted file: " + file_name
		print '\nPerhaps you need an ENVI .hdr file? If so, just open the binary up in ENVI and one will be created for you!'
		os._exit(1)
	else:
		print "%s opened successfully" %file_name
	return inds
	

def image_size(inds):
	cols = inds.RasterXSize
	rows = inds.RasterYSize
	bands = inds.RasterCount
	print "columns: %i" %cols
	print "rows: %i" %rows
	print "bands: %i" %bands
	return cols, rows, bands


def georef_info(inds):
	geotransform = inds.GetGeoTransform()
	originX = geotransform[0]
	originY = geotransform[3]
	pixelWidth = geotransform[1]
	pixelHeight = geotransform[5]

	print "origin x: %i" %originX
	print "origin y: %i" %originY
	print "width: %2.2f" %pixelWidth
	print "height: %2.2f" %pixelHeight
	
	projection = inds.GetProjection()
	return geotransform, projection

def convert_to_2darray(inds, cols, rows, file_name):
	band = inds.GetRasterBand(1)
	datatype = band.DataType
	image_array = band.ReadAsArray(0, 0, cols, rows)
	image_array_name = file_name
	print type(image_array)
	print image_array.shape
	return image_array, image_array_name, band, datatype
	
################

def array_to_raster_format(array, cols, rows, bands, GDALdatatype):
	raster = driver.Create(array, cols, rows, bands, GDALdatatype)
	return raster


def project_raster(raster_projected, geotransform, projection):
	raster_projected.SetGeoTransform(geotransform)
	raster_projected.SetProjection(projection)		
	geotransform_raster = raster_projected.GetGeoTransform()
	
	originX = geotransform_raster[0]
	originY = geotransform_raster[3]
	pixelWidth = geotransform_raster[1]
	pixelHeight = geotransform_raster[5]

	print "origin x: %i" %originX
	print "origin y: %i" %originY
	print "width: %2.2f" %pixelWidth
	print "height: %2.2f" %pixelHeight
	
	return raster_projected
	
	
def data_to_raster(raster, array):
	raster_data = raster.GetRasterBand(1)
	print "raster_data: %s" %(raster_data)
	
	raster_data.WriteArray(array, 0, 0)	
	
	return raster_data
	
	
## Set file location

#file_name_1 = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/222a_lidar/bin/222a.helheim_post_0.5m.bin"
file_name_1 = r"/geog/data/sirius/epsilon/ggwillc/Gaussian_surface_filtering/Helheim/222/HELHEIM_222a_dem_gaussian_filter_sigma_120_order_0_crevasse_surface"
file_name_2 = r"/geog/data/sirius/epsilon/ggwillc/Gaussian_surface_filtering/Helheim/222/HELHEIM_222a_dem_gaussian_filter_sigma_50_order_0_crevasse_surface"

print '~~~~~~~~~~~~~~'
print 'Open file'
print '~~~~~~~~~~~~~~'

file_1 = open_file_gdal(file_name_1)
file_2 = open_file_gdal(file_name_2)

print '~~~~~~~~~~~~~~'
print 'Get image size'
print '~~~~~~~~~~~~~~'

cols_1, rows_1, bands_1 = image_size(file_1)
cols_2, rows_2, bands_2 = image_size(file_2)

print '~~~~~~~~~~~~~~'
print 'Get georeference information'
print '~~~~~~~~~~~~~~'

geotransform_1, projection_1 = georef_info(file_1)
geotransform_2, projection_2 = georef_info(file_2)

print '~~~~~~~~~~~~~~' 
print 'Convert image to 2D array'
print '~~~~~~~~~~~~~~'

image_array_1, image_array_name_1, bands_1, datatype_1 = convert_to_2darray(file_1, cols_1, rows_1, file_name_1)
image_array_2, image_array_name_2, bands_2, datatype_2 = convert_to_2darray(file_2, cols_2, rows_2, file_name_2)
print "data_type: %s" %(datatype_1)

print '~~~~~~~~~~~~~~' 
print 'Surface differencing'
print '~~~~~~~~~~~~~~'

difference = image_array_1 - image_array_2
#difference = image_array_1

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 
print 'Difference surface to binary'
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#difference_raster_format = array_to_raster_format('difference_surface.bin', cols_1, rows_1,  bands_1, 'GDT_Float32') # datatypemust be a GDALDataType..
difference_raster_format = array_to_raster_format('/geog/data/sirius/epsilon/ggwillc/Gaussian_surface_filtering/Helheim/difference_surfaces/difference_surface.bin',14001,17201,1,6) # datatypemust be a GDALDataType..

difference_raster_format_projected = project_raster(difference_raster_format, geotransform_1, projection_1)
difference_ouput = data_to_raster(difference_raster_format_projected, difference)

# SAVE RASTER FILE......

print "complete"
