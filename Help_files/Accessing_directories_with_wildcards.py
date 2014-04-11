"""
Go through directories using wildcards
Author: Chris 29/11/13
"""

import glob # package to let you use wild cards (* etc.)

path = "path/to/dir/*.csv"

for fname in glob.glob(path):
    print(fname)
	
	
	
	