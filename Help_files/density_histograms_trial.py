"""
Go through density list directories
Display list values in histograms
Author: Chris 3/12/13
"""

import glob
import os as os
import matplotlib as matplotlib
import numpy as NP
import matplotlib.pyplot as plt

fldr201a = "/home/staff/ggwillc/Desktop/AL_ARSF_GRNLND_2013/LiDAR/201a/Count_density_lists/"
fldr201b = "/home/staff/ggwillc/Desktop/AL_ARSF_GRNLND_2013/LiDAR/201b/Count_density_lists/"
fldr203a = "/home/staff/ggwillc/Desktop/AL_ARSF_GRNLND_2013/LiDAR/203a/Count_density_lists/"
test = "/home/staff/ggwillc/Desktop/Test_count_density_hist/"

active_folder = ""
array = (fldr201a, fldr201b, fldr203a)
array_test = (test, fldr201a, fldr201b, fldr203a)

for i in range(0,1):

	active_folder = array[i];
	print "Active folder set to: %s" %(active_folder)
	
	path = "%s*" %(active_folder)
	print "Path variable set to: %s" %(path)
	
	### Need to exract these
	#factor = 1.375
	#cell = 2.75
	
        for fname in glob.glob(path):
            print "fname is: %s" %(fname)
	
            ############## If file size is less than a given value do this
	    ############## Otherwise flag the file up and its assoc. size
	    ############## Also, if histogram files doesn't exist
 
	    if(os.path.getsize(fname) < 6e8):

                x = NP.loadtxt(fname)
	    
                count_density_instance = fname.split('density_list_file_')[1]
                print "count_density_instance: %s" %(count_density_instance)
                factor_instance = count_density_instance.split('_')
                print "factor_instance: %s" %(factor_instance[0])
                
                temp_snip_1 = factor_instance[1]
                cell_size_instance = temp_snip_1.split('c_')[0]
                cell_size_alone = cell_size_instance.split('c')
                print "cell_size: %s" %(cell_size_alone[0])
                
	        max_list = max(x)
	        min_list = min(x)
	        mean_list = NP.mean(x)
	        median_list = NP.median(x)
	    
	        print "max value = %f" %(max_list)
	        print "min value = %f" %(min_list)
	        print "mean value = %f" %(mean_list)
	        print "median value = %f" %(median_list)
	    
	        # Histogram specificiation
                num_bins = 100
	        print "num_bins: %i" %(num_bins)
                n, bins, patches = plt.hist(x, num_bins, range=(0, 300), normed=False, facecolor='green', alpha=0.5)
	        plt.xlabel('Hits per m^2')
	        plt.ylabel('Frequency')
                #plt.title('Histogram of count density: \n factor %.2f with a cell size of %.2f m^2' %(factor, cell))
	        plt.title('Factor: %s Cell size: %s' %(factor_instance[0], cell_size_alone[0]))
	    #	
	    
	    ##
	    ### Plot and/or save chart 
            ##
	    
	    #plt.show()
                plt.subplots_adjust(left=0.15) # prevents clipping of y axis label
	        plt.savefig("%s_histogram.png" %(fname))
               
            else:
                print "Filename: %s \n File is larger than 600 MB - are you sure that this is to be processed? \n If you do want to plot it, you'll have to do it separately....." %(fname)




