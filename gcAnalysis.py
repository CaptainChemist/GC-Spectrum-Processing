import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import glob
import re
import peakdetect

sample_names = ["4 mm benzil", "4 mm benzil, 40 nm CdS", "4 mm benzil, 400 nm CdS", "4 mm benzil, 4000 nm CdS"]
controls_location = "./gc/controls/*.txt"
samples_location = "./samples/2015-12-21/*.txt"
number_of_samples = 4
integration_range = 0.25
poly_calibration_fit_order = 1
points_in_poly_fit = 5000
concentration_to_find_time = 5000
min_elution_time = 3.5


class Data:
    # This set of tasks gets called when initializing the class
    def __init__(self):
        self.aggregate_calibration()
        self.create_calibration()
        self.aggregate_data()
        self.convert_data_to_concentration()

    def convert_data_to_concentration(self):
        molecules = self.standard_times.keys()
        runs = np.arange(0, self.iterationNumber-1, 1)
        intensity_table = [pd.DataFrame(index=runs, columns=molecules).fillna(0) for x in range(number_of_samples)]

        for molecule in molecules:
            i = -1
            for sample in self.sample:
                i = i + 1
                for run in range(self.iterationNumber):
                    intensity_table[i][molecule][run] = \
                                self.calc_area(sample['t'+str(run)], sample['i'+str(run)],
                          self.standard_times[molecule], integration_range)

        intensity_table = self.correct_control_peak(intensity_table)

        self.concentration_table = self.calc_concentration(intensity_table)

    def calc_concentration(self, intensity_table):
        i=-1
        for sample in intensity_table:
            intensity_table[i] = intensity_table[i].drop('control', 1)
            i = i + 1
            for molecule in sample:
                if molecule != 'control':
                    for run in sample.index:
                        intensity_index = np.argmin(abs(np.array(self.calibration_curves[molecule]) - intensity_table[i][molecule][run]))
                        intensity_table[i][molecule][run] = self.calibration_curves.iloc[intensity_index].name
        return intensity_table

    # Corrects the integrated intensities by the control peak intensity
    def correct_control_peak(self, intensity_table):
        control = intensity_table[0]['control'][0]
        for sample in range(len(intensity_table)):
            intensity_table[sample] = intensity_table[sample].multiply(control/intensity_table[sample]['control'], axis = 0)[:-1]
        return intensity_table

    def create_calibration(self):
        self.get_times()
        self.create_integration_table()
        self.calibration_curves = pd.DataFrame()
        for molecule in self.standard_integration_table:
            one_run = self.create_calibration_curve(
                        self.standard_integration_table[molecule].keys(),
                  self.standard_integration_table[molecule].values,
                  poly_calibration_fit_order)
            one_run.columns = [molecule]
            self.calibration_curves = pd.concat([self.calibration_curves, one_run], axis=1)

    # This creates a single calibration curve given a time and intensity array and a fitting order
    def create_calibration_curve(self, time, intensity, order):
        x_new = np.arange(0, points_in_poly_fit, .1)
        coefs = poly.polyfit(time, intensity, order)
        ffit = poly.polyval(x_new, coefs)
        return pd.DataFrame(data=ffit, index=x_new)

    # This creates a table of integrated intensities for each molecule and concentration
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
            if molecule.type != 'control':
                self.standard_integration_table[molecule.name][molecule.type] \
                    = self.calc_area(molecule.time, molecule.intensity,
                                     self.standard_times[molecule.name], integration_range)

        # Subtract zero concentration
        subtracted_calibration = pd.DataFrame(self.standard_integration_table.values-self.standard_integration_table.iloc[0, :].values,
                                              columns=self.standard_integration_table.columns, index=self.standard_integration_table.index)
        self.standard_integration_table = subtracted_calibration[1:]


    # Calculates the area of a section of a X,Y 1D plot based on peak center and width
    def calc_area(self, time, intensity, peak_center, int_range):
        start = np.argmin(abs(np.array(time) - (peak_center - int_range)))
        end = np.argmin(abs(np.array(time) - (peak_center + int_range)))
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
            elif run.type == concentration_to_find_time:
                peaks = peakdetect.peakdetect(run.intensity, run.time)[0]

                for peak in peaks:
                    if peak[0] > min_elution_time:
                        self.standard_times[run.name] = float(peak[0])
                        break

    # Construct the control list of data frames
    def aggregate_calibration(self):
        self.control = [[] for x in range(len(glob.glob(controls_location)))]
        current_list = 0
        cols = ['time', 'intensity']

        for full_file_name in sorted(glob.glob(controls_location)):
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
        self.sample = [pd.DataFrame() for x in range(number_of_samples)]

        # Goes through each of the files in ascending order and converts a list of times and intensities into a dataframe
        for full_file_name in sorted(glob.glob(samples_location)):
            # First get the name of the file and convert that to a number
            current_number = int(re.findall(r'\d+', full_file_name)[-1])

            # Create column headers based on the number it pulls from the file name
            cols = ['t' + str(int((current_number - 1) / number_of_samples)),
                    'i' + str(int((current_number - 1) / number_of_samples))]

            # Create a data frame for each sample and then concat onto it for future samplings of that sample
            self.sample[(current_number - 1) % number_of_samples] = \
                pd.concat([self.sample[(current_number - 1) % number_of_samples],
                           self.get_csv(full_file_name, cols)], axis=1, join='inner')

        self.iterationNumber = int(current_number / number_of_samples)

    # Get data from *.csv and convert it to a dataframe
    def get_csv(self, full_file_name, cols='empty'):
        header_length = 22
        if cols == 'empty':
            return pd.read_csv(full_file_name, header=header_length, delimiter=r"\s+")
        else:
            return pd.read_csv(full_file_name, header=header_length, delimiter=r"\s+", names=cols)

    ####################################################################################################
    ######################################PUBLIC API FUNCTION CALLS#####################################

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

    # Get a dataframe of the fitted poly calibration curves
    def get_calibration_curves(self):
        return self.calibration_curves

    # Get a dataframe of integrated intensities and calibration samples
    def get_integration_table(self):
        return self.standard_integration_table

    # Get a dictionary of standard elution times
    def get_standard_times(self):
        return self.standard_times

    # Get a list of dataframes of the concentrations over each run
    def get_concentration_table(self):
        return self.concentration_table

    # Plot the calibration line fits and integrated intensities
    def plot_integration_table(self, log=False):
        f, ax = plt.subplots(1)
        for molecule in self.standard_integration_table:
            ax.scatter(self.standard_integration_table.index, self.standard_integration_table[molecule])
            ax.plot(self.calibration_curves.index, self.calibration_curves[molecule])

        if log == True:
            ax.set_yscale('log')
            ax.set_xscale('log')

        plt.show()

    # Plot the intensity data for each sample over time
    def plot_calibrated_samples(self):
        f, axs = plt.subplots(number_of_samples, sharex=True, sharey=True)
        i = -1
        for sample in self.concentration_table:
            i = i + 1
            for molecule in sample:
                axs[i].plot(sample.index, sample[molecule]/1000)
                axs[i].set_ylabel('Concentration (mM)', fontsize = 10)
                axs[i].set_xlabel('Iteration', fontsize = 10)
                axs[i].set_title(sample_names[i], fontsize = 10)

        plt.show()

    # Plot the GC Elution Data
    def plot_all(self):
        f, axs = plt.subplots(number_of_samples, sharex=True, sharey=True)

        for sample in range(number_of_samples):
            for i in range(self.get_iteration_number()):
                axs[sample].plot(self.sample[sample]['t' + str(i)], self.sample[sample]['i' + str(i)])

        plt.show()


data_set = Data()

#print(data_set.get_standard_times())
print(data_set.get_integration_table())
#print(data_set.get_concentration_table())
#data_set.plot_integration_table(True)
data_set.plot_calibrated_samples()
#data_set.plot_all()
#print data_set.get_calibration_curves()
