####
#### Surface trend removal
#### @author Chris 29/11/13
####

#
#NB/
#Rows and columns and x origin are a digit out once read 
#in with gdal compared to the dem_par file.....
#

import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
from matplotlib import cm # colour mapping

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

# Set file location
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/201a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/201a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname
file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/222a_lidar/bin/222a.helheim_post_0.5m.bin"

# open file
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

# Set pixel offset.....
print '~~~~~~~~~~~~~~' 
print 'Convert image to 2D array'
print '~~~~~~~~~~~~~~'
band = inds.GetRasterBand(1)
image_array = band.ReadAsArray(0, 0, cols, rows)
image_array_name = file_name
print type(image_array)
print image_array.shape

print "%s resampled by a factor of 0.25" % (image_array_name)
dem_025_refactor = ndimage.zoom(image_array, 0.25, order=0) ## essentially refactors the array by a factor the order value denotes the interpolation algorithm used
print dem_025_refactor

#~~~This works to display ~~~#
#display = plt.imshow(dem_025_refactor)
#plt.title(image_array_name + ' (resampled)')
#plt.colorbar(display)
#image_output_original =  r'/home/staff/ggwillc/Desktop/filtering_output_images/'+ '201a_original.png'
#plt.savefig(image_output_original)
#plt.show()

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
print 'Create filter'
print '~~~~~~~~~~~~~~'
kernel = 119 #59 #119
sizex = kernel 
sizey = kernel 

filter1 = np.ones([kernel,kernel])
#filter1 = np.ones([sizex,sizey])
filter1 = filter1/(kernel*kernel)
print type(filter1)
print filter1.shape

#print filter1
#plt.imshow(filter1)
#plt.show()

#plt.imshow(filter1)
#plt.show()

print '~~~~~~~~'
print 'Gaussian'
print '~~~~~~~~'

print "Calculating...."
method = "_gaussian"
order_value = 0
dem_gaussian_filter = ndimage.filters.gaussian_filter(image_array,sigma,order=order_value,mode='reflect')
print dem_gaussian_filter.shape

filtered_image = dem_maximum_filter
filtered_image_name = "HELHEIM_222a_dem%s_filter_sigma_order_%i" %(method, order_value)

print '~~~~~~~~~~~~~~'
print 'Surface differencing'
print '~~~~~~~~~~~~~~'

# original array - filtered array = crevasse array
crevasse_array = image_array - filtered_image
crevasse_array_name = filtered_image_name + "_crevasse_surface"
print crevasse_array.shape
print crevasse_array.size

print '~~~~~~~~~~~~~~'
print 'Export the array subtractions as a new ENVI binary file'
print '~~~~~~~~~~~~~~'

opath = r"/geog/data/sirius/epsilon/ggwillc/Gaussian_surface_filtering/Helheim/222/"

if os.path.isdir(opath):
	print "output_path exists"	
else:
	print "output_path DOESN'T exist...\n"
	os.makedirs(opath)
	print "...but it does now"
output_path = opath + crevasse_array_name

# Creates a new raster data source
outDs = driver.Create(output_path, cols, rows, bands, gdal.GDT_Float32)

# Write metadata
outDs.SetGeoTransform(inds.GetGeoTransform())
outDs.SetProjection(inds.GetProjection())

#Write raster datasets

#for i in range(1):
outBand = outDs.GetRasterBand(1)
outBand.WriteArray(crevasse_array)

# Close raster file
#outDs = None

'''
print '~~~~~~~~~~~~~~'
print 'Create a slice profile of the three layers'
print '~~~~~~~~~~~~~~'

# original array slice
fig = plt.figure()
ax1=fig.add_subplot(311)
ax1.plot(x[:,x.shape[0]/2],image_array[:,x.shape[0]/2],'r') # creates a plot with the x axis taken halfway along the x array - the y axis values are from the image-array at the x positions
ax1.set_title('Unfiltered surface')
#ax1.plot(x[:,x.shape[0]/2],filtered_image[:,x.shape[0]/2],'b')

# filtered array slice
ax2 = fig.add_subplot(312)
ax2.plot(x[:,x.shape[0]/2],filtered_image[:,x.shape[0]/2],'b') # creates a plot with the x axis taken halfway along the x array - the y axis values are from the image-array at the x positions
ax2.set_title('Filtered surface')

# difference slice 
ax3 = fig.add_subplot(313)
ax3.plot(x[:,x.shape[0]/2],crevasse_array[:,x.shape[0]/2],'g') # creates a plot with the x axis taken halfway along the x array - the y axis values are from the image-array at the x positions
ax3.set_title('Difference surface')

plt.show()
image_output = r'/home/staff/ggwillc/Desktop/filtering_output_images/' + filtered_image_name + 'slices.pdf'
plt.savefig(image_output)
plt.show()
'''

print '~~~~~~~~~~~~~~'
print 'Clear variables'
print '~~~~~~~~~~~~~~'
 
band = None
image_array = None

print '~~~~~~~~~~~~~~'
print '~~~~~~~~~~~~~~'
# end timing
endTime = time.time()
print 'The script took ' + str(endTime - startTime) + ' seconds'

plt.clf() # clears figures