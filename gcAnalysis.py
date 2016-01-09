import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

number_of_samples = 4

class Data:

	# This set of tasks gets called when initializing the class
	def __init__(self):
		self.create_calibration()
		print self.control[1].name
		#self.aggregate_data()


	def create_calibration(self):
		self.control = [[] for x in range(len(glob.glob("./gc/controls/*.txt")))]
		current_list = 0

		for full_file_name in sorted(glob.glob("./gc/controls/*.txt")):

			if re.findall('internal_ref*',full_file_name):
				cols = ['internal_ref_time','internal_ref_intensity']
				self.control[current_list] = self.get_csv(full_file_name, cols)
				self.control[current_list].name = str(re.findall('ref\_([a-zA-Z0-9]*)', full_file_name)[0])
				self.control[current_list].type = 'control'
				current_list = current_list + 1
			else:
				self.control[current_list] = self.get_csv(full_file_name)
				self.control[current_list].name = str(re.findall('\_([a-zA-Z0-9]*)', full_file_name)[0])
				self.control[current_list].type = int(re.findall('([0-9]*)\_', full_file_name)[-1])
				current_list = current_list + 1




	# Construct the data array
	def aggregate_data(self):
		# Setting up a 2D empty list
		self.sample = [[] for x in range(number_of_samples)]

		# Goes through each of the files in ascending order and converts a list of times and intensities into a dataframe
		for full_file_name in sorted(glob.glob("./gc/*.txt")):
			# First get the name of the file and convert that to a number
			current_number = int(re.findall(r'\d+', full_file_name)[-1])

			# Create column headers based on the number it pulls from the file name
			cols = ['t' + str(int((current_number - 1) / number_of_samples)),
			        'i' + str(int((current_number - 1) / number_of_samples))]

			#Create a data frame for each sample and then concat onto it for future samplings of that sample
			if current_number <= number_of_samples:
				self.sample[current_number - 1] = self.get_csv(full_file_name, cols)
			else:
				self.sample[(current_number - 1) % number_of_samples] = \
					pd.concat([self.sample[(current_number - 1) % number_of_samples], self.get_csv(full_file_name, cols)], axis=1, join='inner')

		self.iterationNumber = int(current_number / number_of_samples)

	# Get data from *.csv and convert it to a dataframe
	def get_csv(self, full_file_name, cols= 'empty'):
		header_length = 22
		if cols == 'empty':
			return pd.read_csv(full_file_name, header= header_length, delimiter=r"\s+")
		else:
			return pd.read_csv(full_file_name, header= header_length, delimiter=r"\s+", names=cols)

	# Saves one sample to *.csv file
	def save_sample(self, value):
		self.sample[value].to_csv(str(value) + '.csv', sep='\t', encoding='utf-8')

	# Saves all samples to *.csv files
	def save_all(self):
		for current_sample in range(number_of_samples):
			self.save_sample(current_sample)

	# Returns a 2xN dataframe of one sample
	def get_sample(self, value):
		return self.sample[value]

	# Display entire 2D List
	def get_all(self):
		return self.sample

	# Get the number of sampling
	def get_iteration_number(self):
		return self.iterationNumber

	def plot_all(self):
		plt.close('all')
		f, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=True)

		for i in range(self.get_iteration_number()):
			ax0.plot(self.sample[0]['t' + str(i)], self.sample[0]['i' + str(i)])
		for i in range(data_set.get_iteration_number()):
			ax1.plot(self.sample[1]['t' + str(i)], self.sample[1]['i' + str(i)])
		for i in range(data_set.get_iteration_number()):
			ax2.plot(self.sample[2]['t' + str(i)], self.sample[2]['i' + str(i)])
		for i in range(data_set.get_iteration_number()):
			ax3.plot(self.sample[3]['t' + str(i)], self.sample[3]['i' + str(i)])

		plt.show()


data_set = Data()
#data_set.save_all()
#data_set.plot_all()

