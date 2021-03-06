import os
import subprocess

#cor xyz file development

# loop through all cor files
# get file name
# select x, y and z cols using cut tool
# write to new file using old name + "_xyz.cor"

# change directory and list contents
path = "/home/staff/ggwillc/Desktop/Iceland_Bernd/Bernd_preprocess_FRI060913_SkeidararDepression"
os.chdir(path)
subprocess.call(["ls","-l"])

# open the file
txt = open('Depression1.cor','r')
print txt

# cut the xyz columns

print subprocess.call(['cut','-d', '\t', '-f', '4'], stdin=txt, stdout=x) 

'''
x_list[] = subprocess.call(["cut","-d'	'","-f4"], stdout=subprocess.PIPE) 
y_list[] = subprocess.call(["cut","-d"	"","-f6"], stdout=subprocess.PIPE)
z_list[] = subprocess.call(["cut","-d"	"","-f8"], stdout=subprocess.PIPE)
'''
