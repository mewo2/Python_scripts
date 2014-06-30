import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

## Remember - subplot(131) means subplot frame of 1 row, 3 col and image position 1 
## In this instance 131, 132 and 133 == top left, top centre and top right
## In this instance 134, 135 and 136 == left, centre and right
## In this instance 137, 138 and 139 == bottom left, bottom centre and bottom right

filter_percentage =  50

img_1 = '/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/ROI_small/ROI_1.png'
img_2 = 'helheim_222a_sample_ELEVATION_ROI_1_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_2D_AXIS_FREQ.png' %(filter_percentage)
img_3 = 'helheim_222a_sample_ELEVATION_ROI_1_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_divided_by_brown_noise.png' %(filter_percentage)

img_4 = '/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/ROI_small/ROI_2.png'
img_5 = 'helheim_222a_sample_ELEVATION_ROI_2_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_2D_AXIS_FREQ.png' %(filter_percentage)
img_6 = 'helheim_222a_sample_ELEVATION_ROI_2_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_divided_by_brown_noise.png' %(filter_percentage)

img_7 = '/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/ROI_small/ROI_3.png'
img_8 = 'helheim_222a_sample_ELEVATION_ROI_3_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_2D_AXIS_FREQ.png' %(filter_percentage)
img_9 = 'helheim_222a_sample_ELEVATION_ROI_3_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_divided_by_brown_noise.png' %(filter_percentage)

img_10 = '/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/ROI_small/ROI_4.png'
img_11 = 'helheim_222a_sample_ELEVATION_ROI_4_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_2D_AXIS_FREQ.png' %(filter_percentage)
img_12 = 'helheim_222a_sample_ELEVATION_ROI_4_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_divided_by_brown_noise.png' %(filter_percentage)

img_13 = '/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/ROI_small/ROI_5.png'
img_14 = 'helheim_222a_sample_ELEVATION_ROI_5_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_2D_AXIS_FREQ.png' %(filter_percentage)
img_15 = 'helheim_222a_sample_ELEVATION_ROI_5_butterworth_1_low_pass_fft_output_butterworth_1_filter_%i_percent_divided_by_brown_noise.png' %(filter_percentage)

rows = 5
cols = 3
######

img=mpimg.imread(img_1)
plt.subplot(rows,cols,1)
plt.imshow(img)
plt.title("Elevation surface")
plt.axis('off')
plt.ylabel('ROI 1')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir)
#plt.show()

img=mpimg.imread(img_2)
plt.subplot(rows,cols,2)
plt.title("FFT (low pass)")
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir)
#plt.show()

img=mpimg.imread(img_3)
plt.subplot(rows,cols,3)
plt.title("FFT (brown)")
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

#######

img=mpimg.imread(img_4)
plt.subplot(rows,cols,4)
plt.imshow(img)
plt.axis('off')
plt.ylabel('ROI 2')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_5)
plt.subplot(rows,cols,5)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_6)
plt.subplot(rows,cols,6)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

#######

img=mpimg.imread(img_7)
plt.subplot(rows,cols,7)
plt.imshow(img)
plt.axis('off')
plt.ylabel('ROI 3')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_8)
plt.subplot(rows,cols,8)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_9)
plt.subplot(rows,cols,9)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.savefig(odir)
#plt.show()

#######

img=mpimg.imread(img_10)
plt.subplot(rows,cols,10)
plt.imshow(img)
plt.axis('off')
plt.ylabel('ROI 4')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_11)
plt.subplot(rows,cols,11)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_12)
plt.subplot(rows,cols,12)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.savefig(odir)
#plt.show()

#######

img=mpimg.imread(img_13)
plt.subplot(rows,cols,13)
plt.imshow(img)
plt.axis('off')
plt.ylabel('ROI 5')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_14)
plt.subplot(rows,cols,14)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test.pdf'
#plt.savefig(odir, dpi=400)
#plt.show()

img=mpimg.imread(img_15)
plt.subplot(rows,cols,15)
plt.imshow(img)
plt.axis('off')
odir = r'/home/staff/ggwillc/Desktop/IMAGE_test_filter_percentage_%i.pdf' %(filter_percentage)
plt.savefig(odir, dpi=400)
#plt.savefig(odir)
#plt.show()

#######
