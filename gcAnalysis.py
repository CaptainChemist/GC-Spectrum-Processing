import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import glob
import re
import peakdetect

number_of_samples = 4
integration_range = 0.2
poly_calibration_fit_order = 2


# 2) Find the reactant peak elution time for the smallest concentration sample (min, max, center)
# 3) Calculate the reactant peak


class Data:
	# This set of tasks gets called when initializing the class
	def __init__(self):
		self.aggregate_calibration()
		self.create_calibration()

	# self.aggregate_data()

	def create_calibration(self):
		self.get_times()
		self.create_integration_table()

		for molecule in self.standard_integration_table:
			print molecule
			self.create_calibration_curve(self.standard_integration_table[molecule].keys(),
			                              self.standard_integration_table[molecule].values,
			                              poly_calibration_fit_order)

	def create_calibration_curve(self, time, intensity, order):
		print time
		print intensity
		x_new = np.arange(0, 5000, 1)
		coefs = poly.polyfit(time, intensity, order)
		ffit = poly.polyval(x_new, coefs)

		plt.scatter(intensity, time)
		plt.plot(ffit, x_new)
		plt.show()

	def create_integration_table(self):
		# Find the concentrations and molecules present and create an empty calibration dataframe
		concentrations = []
		molecules = []
		for run in self.control:
			if run.type not in concentrations and run.type != 'control':
				concentrations.append(run.type)
			if run.name not in molecules and run.type != 'control':
				molecules.append(run.name)
		self.standard_integration_table = pd.DataFrame(index=concentrations, columns=molecules).fillna(0)

		# Loop over the concentration calibration data and extract integrated intensities
		for molecule in self.control:
			if molecule.type == 'control':
				background = []
			else:
				self.standard_integration_table[molecule.name][molecule.type] \
					= self.calc_area(molecule.time, molecule.intensity,
					                 self.standard_times[molecule.name], integration_range)

	# Calculates the area of a section of a X,Y 1D plot based on peak center and width
	def calc_area(self, time, intensity, peak_center, integration_range):
		start = np.argmin(abs(np.array(time) - (peak_center - integration_range)))
		end = np.argmin(abs(np.array(time) - (peak_center + integration_range)))
		xs = np.array(time[start:end])
		ys = np.array(intensity[start:end])
		slope = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
		baseline = slope * (xs - xs[0]) + ys[0]
		integrated_intensity = np.trapz(ys - baseline, xs)

		return integrated_intensity

	# Gets the elution times for each standard and puts it in a table
	def get_times(self):
		self.standard_times = {}

		for run in self.control:
			if run.type == 'control':
				peaks = peakdetect.peakdetect(run.intensity, run.time)
				self.standard_times['control'] = float([peaks[0][1][0]][0])
			elif run.type == 4000:
				peaks = peakdetect.peakdetect(run.intensity, run.time)[0]

				for peak in peaks:
					if peak[0] > 3.5:
						self.standard_times[run.name] = float(peak[0])
						break

	# Construct the control list of data frames
	def aggregate_calibration(self):
		self.control = [[] for x in range(len(glob.glob("./gc/controls/*.txt")))]
		current_list = 0
		cols = ['time', 'intensity']

		for full_file_name in sorted(glob.glob("./gc/controls/*.txt")):
			self.control[current_list] = self.get_csv(full_file_name, cols)

			if re.findall('internal_ref*', full_file_name):
				self.control[current_list].name = str(re.findall('ref\_([a-zA-Z0-9]*)', full_file_name)[0])
				self.control[current_list].type = 'control'
			else:
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

			# Create a data frame for each sample and then concat onto it for future samplings of that sample
			if current_number <= number_of_samples:
				self.sample[current_number - 1] = self.get_csv(full_file_name, cols)
			else:
				self.sample[(current_number - 1) % number_of_samples] = \
					pd.concat([self.sample[(current_number - 1) % number_of_samples], self.get_csv(full_file_name, cols)], axis=1,
					          join='inner')

		self.iterationNumber = int(current_number / number_of_samples)

	# Get data from *.csv and convert it to a dataframe
	def get_csv(self, full_file_name, cols='empty'):
		header_length = 22
		if cols == 'empty':
			return pd.read_csv(full_file_name, header=header_length, delimiter=r"\s+")
		else:
			return pd.read_csv(full_file_name, header=header_length, delimiter=r"\s+", names=cols)

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
# data_set.save_all()
# data_set.plot_all()
