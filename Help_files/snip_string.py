import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

a1 = 'helheim_222a_sample_MAX_ROI_2_butterworth_1_high_pass_fft_output_butterworth_1_filter_50_percent_2D_AXIS_FREQ.png'
a2 = a1.split('_pass')[0]
a3 = a2.split('MAX_')[1]
print a3

b1 = 'helheim_222a_sample_MAX_ROI_2_butterworth_1_high_pass_fft_output_butterworth_1_filter_50_percent_2D_AXIS_FREQ.png'
b2 = a1.split('_sample')[0]
print b2
