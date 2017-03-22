import numpy as np, matplotlib.pyplot as plt, os, sys, subprocess

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
	return data_list

def check_iterations(iterations, max_iterations):
	if iterations == max_iterations:
		return '^*'
	elif iterations == 1:
		return '^\dagger'
	else:
		return ''

def convert_to_latex(list_dictionary, max_iter):
	string_to_build = '%10.d & $%s$ & $%s$ & $%s$ & $%s$  \\\ \n'
	latex_string = ''
	start_shell = 4
	stop_shell = 9
	bare_string = '%10.6f%s'
	for shell in xrange(start_shell, stop_shell+1):
		HFE2 = '%10.6f%s'
		HFE6 = '%10.6f%s'
		HFE12 = '%10.6f%s'
		HFE20 = '%10.6f%s'
		for stat_list in list_dictionary:
			if int(stat_list['electrons']) == 2 and int(stat_list['maxShell']) == shell and HFE2 == bare_string:
				HFE2 = HFE2 % (float(stat_list['HFEnergy']), check_iterations(int(stat_list['HFIterations']),max_iter))
			if int(stat_list['electrons']) == 6 and int(stat_list['maxShell']) == shell and HFE6 == bare_string:
				HFE6 = HFE6 % (float(stat_list['HFEnergy']), check_iterations(int(stat_list['HFIterations']),max_iter))
			if int(stat_list['electrons']) == 12 and int(stat_list['maxShell']) == shell and HFE12 == bare_string:
				HFE12 = HFE12 % (float(stat_list['HFEnergy']), check_iterations(int(stat_list['HFIterations']),max_iter));
			if int(stat_list['electrons']) == 20 and int(stat_list['maxShell']) == shell and HFE20 == bare_string:
				HFE20 = HFE20 % (float(stat_list['HFEnergy']), check_iterations(int(stat_list['HFIterations']),max_iter))		
		latex_string += string_to_build % (int(shell), HFE2, HFE6, HFE12, HFE20)
	return latex_string

def runInTerminal(cmd):
	# Function for running commands in terminal
	print 'Running command: ', cmd
	subprocess.call(cmd, shell=True)

def main():
	# runInTerminal('./project1_src/project1') # For automating the data gathering
	max_iter = 500
	output_folder = 'output'
	omega10_data = get_data(output_folder, 1.0)
	omega01_data = get_data(output_folder, 0.1)
	omega028_data = get_data(output_folder, 0.5)
	omega05_data = get_data(output_folder, 0.28)
	print 'Omega = 0.1'
	print convert_to_latex(omega01_data, max_iter)
	print 'Omega = 0.28'
	print convert_to_latex(omega028_data, max_iter)
	print 'Omega = 0.5'
	print convert_to_latex(omega05_data, max_iter)
	print 'Omega = 1.0'
	print convert_to_latex(omega10_data, max_iter)

if __name__ == '__main__':
	main()