import numpy as np, matplotlib.pyplot as plt, os, sys, subprocess

# output_folder = 'output'
# omega10 = [] # For omega=1.0
# omega01 = [] # For omega=0.1

def get_data(output_folder, target_omega = None):
	data_list = []
	for file in os.listdir(output_folder):		
		if not file.split('.')[-1] == 'txt': continue # Ensuring we have a data file
		if file == 'HF_results1.txt': continue # ensuring only good files are checked
		if not 'omega' in file: continue # Ensuring we have omega value in file
		file_omega = float(file.split('omega')[-1].split('_')[0])
		if file_omega == target_omega:
			electron_number = file.split('_')[-1].split('electrons')[0]
			for line in open(output_folder + '/' + file,'r'):
				line_elements = line.split()
				electrons_stats = {}
				electrons_stats['electrons'] = electron_number
				electrons_stats['maxShell'] = line_elements[3]
				electrons_stats['HFIterations'] = line_elements[5]
				electrons_stats['HFEnergy'] = line_elements[7]
				data_list.append(electrons_stats)

		# else:
		# 	electron_number = file.split('_')[-1].split('electrons')[0]
		# 	for line in open(output_folder + '/' + file,'r'):
		# 		line_elements = line.split()
		# 		electrons_stats = {}
		# 		electrons_stats['electrons'] = electron_number
		# 		electrons_stats['maxShell'] = line_elements[3]
		# 		electrons_stats['HFIterations'] = line_elements[5]
		# 		electrons_stats['HFEnergy'] = line_elements[7]
		# 		data_list.append(electrons_stats)
	return data_list


def convert_to_latex(list_dictionary):
	#'Shells  & $A = 2$	& $A = 6$ 	& $A = 12$ 	& $A = 20$ \\ \hline'
	string_to_build = '%10.d & %10.6f & %10.6f & %10.6f & %10.6f  \\\ \n'
	latex_string = ''
	start_shell = 4
	stop_shell = 9
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

		latex_string += string_to_build % (int(shell), float(HFE2), float(HFE6), float(HFE12), float(HFE20))
			# latex_string += stat_list['HFEnergy'] + ' & '
		
	return latex_string

def runInTerminal(cmd):
	# Function for running commands in terminal
	print 'Running command: ', cmd
	subprocess.call(cmd, shell=True)

def main():
	# runInTerminal('./project1_src/project1') # For automating the data gathering

	output_folder = 'output2'
	omega10_data = get_data(output_folder, 1.0)
	omega028_data = get_data(output_folder, 0.5)
	omega05_data = get_data(output_folder, 0.28)
	print 'Omega = 0.28'
	print convert_to_latex(omega028_data)
	print 'Omega = 0.5'
	print convert_to_latex(omega05_data)
	print 'Omega = 1.0'
	print convert_to_latex(omega10_data)

if __name__ == '__main__':
	main()