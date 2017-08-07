from argparse import ArgumentParser
from decimal import *
import pandas as pd
import numpy as np
import os
import sys

def fort36Ave_condensed(fort36_file, out_location):
	df = pd.read_csv(fort36_file,
		sep = "\s+", header = None, usecols = [2, 6, 9],
		names = ['pH', 'E_type', 'E'])
	#print(df)
	AverE_df = df[ df['E_type'] == 'Ave.' ]

	# the agg.function is used on groupby objects and allows for users
	# to combine multiple statistics. For example, we can get values for 
	# both mean and std, grouped by pH, in one line.
	_res = AverE_df.groupby(['pH'], sort=False).agg( ['mean', 'std'] )
	_res.columns = ['Mean', 'Std']
	
	for type in ['Mean', 'Std']:
		_res[type] = ['%.3f' % round(val,3) for val in 
			_res[type]]

	file_location = out_location
	if not os.path.exists(file_location):
		os.makedirs(file_location)
	
	file_name = 'fort98'
	
	to_write = os.path.join(file_location, file_name)
	_res.to_csv(to_write, sep = ' ')


def fort36Conformers_std(fort36_file, out_location):
	#getcontext().prec = 2
	master_df = pd.DataFrame()

	df = pd.read_csv(fort36_file,
		sep = "\s+", header = None, usecols = [0, 1, 2],
		names = ['pH', 'Conformer', 'Occ'])
	#conformers = set([true_conformer if true_conformer != '='  true_conformer in df['Conformer'].values])	
	conformers = set(df['Conformer'].values)
	conformers.remove('=')
	
	df['Occ'] = [float(string[4:len(string)]) if string[0:3] == 'occ'
		else np.nan for string in df['Occ'].values]
	df = df.dropna()
 
	for conformer in sorted(conformers):
		df_conformer = df[ df['Conformer'] == str(conformer) ]
		to_write = df_conformer.groupby(['pH'], 
			sort = False).agg('std')
		to_write.columns = [str(conformer)]
		to_write[conformer] = ['%.3f' % round(val,3) for val in 
			to_write[conformer]]
	
		master_df = master_df.join(pd.DataFrame(to_write), 
			how = 'right')

	master_df.index.name = 'pH'	
	
	#transposed = master_df.transpose()
	
	file_location = out_location
	if not os.path.exists(file_location):
		os.makedirs(file_location)
	file_name = 'fort99'
	to_write_loc = os.path.join(file_location, file_name)

	master_df.transpose().to_csv(to_write_loc, sep = ' ', header = ['pH            0.00  1.00  \
2.00  3.00  4.00  5.00  6.00  7.00  8.00  9.00 10.00 11.00 12.00 13.00 14.00','', '', '',
		'', '', '', '', '', '', '','' , '', '', ''])
	
	std_cap = 0.100
	high_std = []
	analysis_file = open(to_write_loc, 'r').readlines()
	for line_index in range(len(analysis_file)):
		if(line_index == 0): pass
		else:
			for word in analysis_file[line_index].split():
				try:
					compare = float(word)
					if(compare >= std_cap):
						high_std.append(line_index)
				except ValueError:
					#print(word +  'is not a float')
					pass

	with open(to_write_loc, 'a') as append:
		append.write('''-------------------------------------------------------------------------------------------------
------------------------------------------SUMMARY------------------------------------------------
-------------------------------------HIGH STD CONFORMERS-----------------------------------------

''')
		if(len(high_std) == 0):
			summary = 'NO HIGH STD CONFORMERS FOUND'
		else:
			summary = ''.join([analysis_file[line] for
						line in high_std])
		append.write(summary)
		append.close()	


def ask_arguments():
	parser = ArgumentParser(description = '''Analyze values in fort36, 
		calculates mean of energies at each pH by default. Can
		process more values by if required by user''')
	required = parser.add_argument_group('''Required arguments
		please enter absolute paths where possible''')
	# required arguments below:
	required.add_argument('-f', '--fort36_location', required = True, 
		type = str, help = '''Location of fort36 file to be 
		analyzed''')
	required.add_argument('-out', '--output_location', required = True,
		type = str, help = '''Location of FOLDER containing output
		csv files''')
	# optional arguments below:
	optional = parser.add_argument_group('''Optional arguments''')
	optional.add_argument('-stdev', '--standard_deviation', 
		required = False, type = bool, help = '''Flag indicating
		whether stdev should be calculated. False by default''')
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	args = ask_arguments()
	fort36Ave_condensed(args.fort36_location, args.output_location)
	fort36Conformers_std(args.fort36_location, args.output_location)



