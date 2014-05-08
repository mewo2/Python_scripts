#~~~~~~~~~~~~~~~~
#Raster Functions
#~~~~~~~~~~~~~~~~


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

def create_output_directory(opath):
	if os.path.isdir(opath):
		print "output_path exists"	
	else:
		print "output_path DOESN'T exist...\n"
		os.makedirs(opath)
		print "...but it does now"
	
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
	
	
def negative_anomaly_array(array):
	negative_anomaly_array = array.copy()
	'''
	ii,jj = array.shape
	for i in xrange(ii):
		for j in xrange(jj):
			if array[i,j] > 0:
				negative_anomaly_array[i,j] = 0
			else:
				negative_anomaly_array[i,j] = negative_anomaly_array[i,j]
	'''
	negative_anomaly_array[array > 0] = 0
	return negative_anomaly_array


def positive_anomaly_array(array):
	positive_anomaly_array = array.copy()
	positive_anomaly_array[array < 0] = 0
	return positive_anomaly_array
