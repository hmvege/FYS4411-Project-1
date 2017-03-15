import numpy as np, matplotlib.pyplot as plt, os, sys

output_folder = "output"
omega10 = {} # For omega=1.0
omega01 = {} # For omega=0.1

for file in os.listdir(output_folder):
	if not file.split('.')[-1] == 'txt': continue # Ensuring we have a data file
	print file
	if 'omega' in file:
		electron_number = file.split('_')[-1].split('electrons')[0]
		d = np.loadtxt(file.replace(' ','\ '))
		print d

	# omega = file.split('omega')
	# print omega