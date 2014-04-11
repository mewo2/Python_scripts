####
# PYTHON
# How to:
#	(i)  open raster binary
#	(ii) convert binary to array
# 	(iii) display array
#	(iv) display a snip cross profile
# 	(v) save array back out as binary
####

import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
from matplotlib import cm # colour mapping

print '~~~~~~~~~~~~~~'
print 'Register driver'
print '~~~~~~~~~~~~~~'
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

print '~~~~~~~~~~~~~~'
print 'Set file location'
print '~~~~~~~~~~~~~~'
file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/201a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname

print '~~~~~~~~~~~~~~'
print 'open file'
print '~~~~~~~~~~~~~~'
inds = gdal.Open(file_name, GA_ReadOnly)

if inds is None:	
	print 'Could not open ' + file_name
	os._exit(1)
else:
	print "%s opened successfully" %file_name
	
print '~~~~~~~~~~~~~~'
print 'Get image size'
print '~~~~~~~~~~~~~~'
cols = inds.RasterXSize
rows = inds.RasterYSize
bands = inds.RasterCount

print "columns: %i" %cols
print "rows: %i" %rows
print "bands: %i" %bands

print '~~~~~~~~~~~~~~'
print 'Get georeference information'
print '~~~~~~~~~~~~~~'
geotransform = inds.GetGeoTransform()
originX = geotransform[0]
originY = geotransform[3]
pixelWidth = geotransform[1]
pixelHeight = geotransform[5]

print "origin x: %i" %originX
print "origin y: %i" %originY
print "width: %2.2f" %pixelWidth
print "height: %2.2f" %pixelHeight

print '~~~~~~~~~~~~~~' 
print 'Convert image to 2D array'
print '~~~~~~~~~~~~~~'
band = inds.GetRasterBand(1)
image_array = band.ReadAsArray(0, 0, cols, rows)
image_array_name = file_name
print type(image_array)
print image_array.shape

print '~~~~~~~~~~~~~~' 
print 'Resample array to make more manageable'
print '~~~~~~~~~~~~~~'
#print "%s resampled by a factor of 0.25" % (image_array_name)
#dem_025_refactor = ndimage.zoom(image_array, 0.25, order=0) ## essentially refactors the array by a factor the order value denotes the interpolation algorithm used
#print dem_025_refactor

print '~~~~~~~~~~~~~~' 
print 'Display resampled array'
print '~~~~~~~~~~~~~~'
#display = plt.imshow(dem_025_refactor)
#plt.title(image_array_name + ' (resampled)')
#plt.colorbar(display)
#image_output_original =  r'/home/staff/ggwillc/Desktop/filtering_output_images/'+ 'test_array_image.png'
#plt.savefig(image_output_original)
##plt.show()

print '~~~~~~~~~~~~~~'
print 'Create x and y grids'
print '~~~~~~~~~~~~~~'
x, y = np.mgrid[0:image_array.shape[0], 0:image_array.shape[1]]
print ' '
print 'x:'
print type(x)
print x.shape
print x.size
print x

print ' '
print 'y:'
print type(y)
print y.shape
print y.size
print y

print '~~~~~~~~~~~~~~'
print 'Create a slice profile of the three layers'
print '~~~~~~~~~~~~~~'

fig = plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(x[:,x.shape[0]/2],image_array[:,x.shape[0]/2],'r') # creates a plot with the x axis taken halfway along the x array - the y axis values are from the image-array at the x positions
ax1.set_title('Surface')
image_output_slice =  r'/home/staff/ggwillc/Desktop/filtering_output_images/'+ 'test_array_image_slice.png'
plt.savefig(image_output_slice)

print '~~~~~~~~~~~~~~'
print 'Export the array subtractions as a new ENVI binary file'
print '~~~~~~~~~~~~~~'

opath = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/"

if os.path.isdir(opath):
	print "output_path exists"	
else:
	print "output_path DOESN'T exist...\n"
	os.makedirs(r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/") 
	print "...but it does now"
output_path = opath + 'test_output_binary.bin'

# Creates a new raster data source
outDs = driver.Create(output_path, cols, rows, bands, gdal.GDT_Float32)

# Write metadata
outDs.SetGeoTransform(inds.GetGeoTransform())
outDs.SetProjection(inds.GetProjection())

#Write raster datasets

#for i in range(1):
outBand = outDs.GetRasterBand(1)
outBand.WriteArray(image_array)

print '~~~~~~~~~~~~~~'
print 'Clear variables'
print '~~~~~~~~~~~~~~'
 
band = None
image_array = None
