'''
How to import functions saved along specific paths
'''

# The file FFT_functions.py contains various functions and is saved: /home/staff/ggwillc/Desktop/Python_scripts/functions
# To import these functions for use in another script, you must do the following at the top of the script in which you want to import the functions:

import sys
sys.path.insert(0, '/home/staff/ggwillc/Desktop/Python_scripts/functions')
import FFT_functions as FFT_functions

# The functions can then be used as FFT_functions.FunctionToUse()

