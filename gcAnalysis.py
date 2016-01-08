import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re


class Data:
	# Get data from *.csv and convert it to a dataframe
	def get_csv(self, full_file_name, cols):
		header_length = 22
		return pd.read_csv(full_file_name, header= header_length, delimiter=r"\s+", names=cols)

	# This set of tasks gets called when initializing the class
	def __init__(self):
		# Setting up a 2D empty list
		self.sample = [[] for i in range(4)]

		# Goes through each of the files in ascending order and converts a list of times and intensities into a dataframe
		for full_file_name in sorted(glob.glob("./gc/*.txt")):
			# First get the name of the file and convert that to a number
			current_number = int(re.findall(r'\d+', full_file_name)[-1])

			# Create column headers based on the number it pulls from the file name
			cols = ['t' + str(int((current_number - 1) / 4)), 'i' + str(int((current_number - 1) / 4))]

			# Create the data frame if it is the first sample in the list
			# Otherwise concat the current data as a data frame with the already existing data frame for that sample
			if current_number <= 4:
				self.sample[current_number - 1] = self.get_csv(full_file_name, cols)
			elif current_number % 4 == 1:
				self.sample[0] = pd.concat([self.sample[0], self.get_csv(full_file_name, cols)], axis=1, join='inner')
			elif current_number % 4 == 2:
				self.sample[1] = pd.concat([self.sample[1], self.get_csv(full_file_name, cols)], axis=1, join='inner')
			elif current_number % 4 == 3:
				self.sample[2] = pd.concat([self.sample[2], self.get_csv(full_file_name, cols)], axis=1, join='inner')
			elif current_number % 4 == 0:
				self.sample[3] = pd.concat([self.sample[3], self.get_csv(full_file_name, cols)], axis=1, join='inner')

		self.iterationNumber = int(current_number / 4)

	# Display entire 2D List
	def get_all(self):
		return self.sample

	# Save One _sample, either 0, 1, 2, or 3 to csv file
	def save_sample(self, value):
		self.sample[value].to_csv(str(value) + '.csv', sep='\t', encoding='utf-8')

	# Returns a 2xN dataframe of the sample data
	def get_sample(self, value):
		return self.sample[value]

	# Get the number of sampling
	def get_iteration_number(self):
		return self.iterationNumber


data_set = Data()
data_set.save_sample(0)
data_set.save_sample(1)
data_set.save_sample(2)
data_set.save_sample(3)

df = data_set.get_all()
print len(df)
print df[0].shape
print df[1].shape
print df[2].shape
print df[3].shape

plt.close('all')
f, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=True)

for i in range(data_set.get_iteration_number()):
	ax0.plot(df[0]['t' + str(i)], df[0]['i' + str(i)])
for i in range(data_set.get_iteration_number()):
	ax1.plot(df[1]['t' + str(i)], df[1]['i' + str(i)])
for i in range(data_set.get_iteration_number()):
	ax2.plot(df[2]['t' + str(i)], df[2]['i' + str(i)])
for i in range(data_set.get_iteration_number()):
	ax3.plot(df[3]['t' + str(i)], df[3]['i' + str(i)])

plt.show()
