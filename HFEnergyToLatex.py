import numpy as np, matplotlib.pyplot as plt, os, sys

output_folder = 'output'
omega10 = [] # For omega=1.0
omega01 = [] # For omega=0.1

for file in os.listdir(output_folder):
	if not file.split('.')[-1] == 'txt': continue # Ensuring we have a data file
	if file == 'HF_results1.txt': continue # ensuring only good files are checked
	if 'omega' in file:
		electron_number = file.split('_')[-1].split('electrons')[0]
		for line in open(output_folder + '/' + file,'r'):
			line_elements = line.split()
			electrons_stats = {}
			electrons_stats['electrons'] = electron_number
			electrons_stats['maxShell'] = line_elements[3]
			electrons_stats['HFIterations'] = line_elements[5]
			electrons_stats['HFEnergy'] = line_elements[7]
			omega01.append(electrons_stats)

	else:
		electron_number = file.split('_')[-1].split('electrons')[0]
		for line in open(output_folder + '/' + file,'r'):
			line_elements = line.split()
			electrons_stats = {}
			electrons_stats['electrons'] = electron_number
			electrons_stats['maxShell'] = line_elements[3]
			electrons_stats['HFIterations'] = line_elements[5]
			electrons_stats['HFEnergy'] = line_elements[7]
			omega10.append(electrons_stats)


def convert_to_latex(list_dictionary):
	'Shells  & $A = 2$	& $A = 6$ 	& $A = 12$ 	& $A = 20$ \\ \hline'
	string_to_build = '%10.10s & %10.10s & %10.10s & %10.10s & %10.10s  \\\ \n'
	latex_string = ''
	start_shell = 3
	stop_shell = 11
	for shell in xrange(start_shell, stop_shell+1):
		HFE2 = 0
		HFE6 = 0
		HFE12 = 0
		HFE20 = 0
		for stat_list in list_dictionary:
			if int(stat_list['electrons']) == 2 and int(stat_list['maxShell']) == shell: HFE2 = stat_list['HFEnergy']
			if int(stat_list['electrons']) == 6 and int(stat_list['maxShell']) == shell: HFE6 = stat_list['HFEnergy']
			if int(stat_list['electrons']) == 12 and int(stat_list['maxShell']) == shell: HFE12 = stat_list['HFEnergy']
			if int(stat_list['electrons']) == 20 and int(stat_list['maxShell']) == shell: HFE20 = stat_list['HFEnergy']

		latex_string += string_to_build % (shell, HFE2, HFE6, HFE12, HFE20)
			# latex_string += stat_list['HFEnergy'] + ' & '
		
	print latex_string

convert_to_latex(omega10)
convert_to_latex(omega01)