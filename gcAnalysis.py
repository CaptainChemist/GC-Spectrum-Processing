import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import glob
import re
import peakdetect
from scipy.optimize import curve_fit
import scipy.special

sample_names = ["4 mm benzil", "4 mm benzil, 40 nm CdS", "4 mm benzil, 400 nm CdS", "4 mm benzil, 4000 nm CdS"]
controls_location = "./gc/controls/*.txt"
samples_location = "./samples/2015-12-21/*.txt"
header_length = 22
number_of_samples = 4
integration_range = .25
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

    def gauss(self, x, *p):
        A1, mu1, sigma1 = p
        return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))

    def two_gauss(self, x, *p):
        A1, mu1, sigma1, A2, mu2, sigma2 = p
        return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))+A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))

    def convert_data_to_concentration(self):
        _molecules = self.standard_times.keys()
        _runs = np.arange(0, self.iteration_number - 1, 1)
        _intensity_table = [pd.DataFrame(index=_runs, columns=_molecules).fillna(0) for x in range(number_of_samples)]

        for _molecule in _molecules:
            _i = -1
            for _sample in self.sample:
                _i = _i + 1
                for _run in range(self.iteration_number):
                    print _i
                    print _run

                    _intensity_table[_i][_molecule][_run] = \
                        self.calc_area(_sample['t' + str(_run)], _sample['i' + str(_run)],
                                       self.standard_times[_molecule], integration_range,
                                       _molecule)

        _intensity_table = self.correct_control_peak(_intensity_table)

        self.concentration_table = self.calc_concentration(_intensity_table)

    def calc_concentration(self, intensity_table):
        _i = -1
        for _sample in intensity_table:
            intensity_table[_i] = intensity_table[_i].drop('control', 1)
            _i = _i + 1
            for _molecule in _sample:
                if _molecule != 'control':
                    for _run in _sample.index:
                        intensity_index = np.argmin(
                                abs(np.array(self.calibration_curves[_molecule]) - intensity_table[_i][_molecule][_run]))
                        intensity_table[_i][_molecule][_run] = self.calibration_curves.iloc[intensity_index].name
        return intensity_table

    # Corrects the integrated intensities by the control peak intensity
    def correct_control_peak(self, intensity_table):
        _control = intensity_table[0]['control'][0]
        for _sample in range(len(intensity_table)):
            intensity_table[_sample] = intensity_table[_sample].multiply(
                    _control / intensity_table[_sample]['control'], axis=0)[:-1]
        return intensity_table

    def create_calibration(self):
        self.get_times()
        self.create_integration_table()
        self.calibration_curves = pd.DataFrame()
        for _molecule in self.standard_integration_table:
            _one_run = self.create_calibration_curve(
                    self.standard_integration_table[_molecule].keys(),
                    self.standard_integration_table[_molecule].values,
                    poly_calibration_fit_order)
            _one_run.columns = [_molecule]
            self.calibration_curves = pd.concat([self.calibration_curves, _one_run], axis=1)

    # This creates a single calibration curve given a time and intensity array and a fitting order
    def create_calibration_curve(self, time, intensity, order):
        _x_new = np.arange(0, points_in_poly_fit, .1)
        _coefs = poly.polyfit(time, intensity, order)
        _ffit = poly.polyval(_x_new, _coefs)
        return pd.DataFrame(data=_ffit, index=_x_new)

    # This creates a table of integrated intensities for each molecule and concentration
    def create_integration_table(self):
        # Find the concentrations and molecules present and create an empty calibration dataframe
        _concentrations = []
        _molecules = []
        for _run in self.control:
            if _run.type not in _concentrations and _run.type != 'control':
                _concentrations.append(_run.type)
            if _run.name not in _molecules and _run.type != 'control':
                _molecules.append(_run.name)
        self.standard_integration_table = pd.DataFrame(index=_concentrations, columns=_molecules).fillna(0)

        # Loop over the concentration calibration data and extract integrated intensities
        for _molecule in self.control:
            if _molecule.type != 'control':
                self.standard_integration_table[_molecule.name][_molecule.type] \
                    = self.calc_area(_molecule.time, _molecule.intensity,
                                     self.standard_times[_molecule.name], integration_range,
                                     _molecule.name)

        # Subtract zero concentration
        _subtracted_calibration = pd.DataFrame(
                self.standard_integration_table.values - self.standard_integration_table.iloc[0, :].values,
                columns=self.standard_integration_table.columns, index=self.standard_integration_table.index)
        self.standard_integration_table = _subtracted_calibration[1:]

    # Calculates the area of a section of a X,Y 1D plot based on peak center and width
    def calc_area(self, time, intensity, peak_center, int_range, two_gauss):
        _start = np.argmin(abs(np.array(time) - (peak_center - int_range)))
        _end = np.argmin(abs(np.array(time) - (peak_center + int_range)))
        _xs = np.array(time[_start:_end])
        _ys = np.array(intensity[_start:_end])
        print peak_center
        print two_gauss
        if two_gauss == 'b2enzil':
            try:
                popt, pcov = curve_fit(self.two_gauss, time, intensity, [1000,peak_center,0.1, 1000, peak_center+0.2, 0.1])

            except:
                popt, pcov = curve_fit(self.gauss, time, intensity, [1000,peak_center,0.1])

        else:
            popt, pcov = curve_fit(self.gauss, time, intensity, [1000,peak_center,0.1])
            #print popt
            #print scipy.integrate.quad(self.gauss, popt, popt[1]-int_range, popt[1]+int_range)

        print popt

        _slope = (_ys[-1] - _ys[0]) / (_xs[-1] - _xs[0])
        _baseline = _slope * (_xs - _xs[0]) + _ys[0]
        _integrated_intensity = np.trapz(_ys - _baseline, _xs)
        return _integrated_intensity

    # Gets the elution times for each standard and puts it in a table
    def get_times(self):
        self.standard_times = {}

        for _run in self.control:
            if _run.type == 'control':
                peaks = peakdetect.peakdetect(_run.intensity, _run.time)
                self.standard_times['control'] = float([peaks[0][1][0]][0])
            elif _run.type == concentration_to_find_time:
                peaks = peakdetect.peakdetect(_run.intensity, _run.time)[0]

                for peak in peaks:
                    if peak[0] > min_elution_time:
                        self.standard_times[_run.name] = float(peak[0])
                        break

    # Construct the control list of data frames
    def aggregate_calibration(self):
        self.control = [[] for x in range(len(glob.glob(controls_location)))]
        _current_list = 0
        _cols = ['time', 'intensity']

        for _full_file_name in sorted(glob.glob(controls_location)):
            self.control[_current_list] = self.get_csv(_full_file_name, _cols)
            if re.findall('internal_ref*', _full_file_name):
                self.control[_current_list].name = str(re.findall('ref\_([a-zA-Z0-9]*)', _full_file_name)[0])
                self.control[_current_list].type = 'control'
            else:
                self.control[_current_list].name = str(re.findall('\_([a-zA-Z0-9]*)', _full_file_name)[0])
                self.control[_current_list].type = int(re.findall('([0-9]*)\_', _full_file_name)[-1])

            _current_list = _current_list + 1

    # Construct the data array
    def aggregate_data(self):
        # Setting up a 2D empty list
        self.sample = [pd.DataFrame() for x in range(number_of_samples)]

        # Goes through the files in ascending order and converts a list of times and intensities into a dataframe
        for _full_file_name in sorted(glob.glob(samples_location)):
            # First get the name of the file and convert that to a number
            _current_number = int(re.findall(r'\d+', _full_file_name)[-1])

            # Create column headers based on the number it pulls from the file name
            _cols = ['t' + str(int((_current_number - 1) / number_of_samples)),
                    'i' + str(int((_current_number - 1) / number_of_samples))]

            # Create a data frame for each sample and then concat onto it for future samplings of that sample
            self.sample[(_current_number - 1) % number_of_samples] = \
                pd.concat([self.sample[(_current_number - 1) % number_of_samples],
                           self.get_csv(_full_file_name, _cols)], axis=1, join='inner')

        self.iteration_number = int(_current_number / number_of_samples)

    # Get data from *.csv and convert it to a dataframe
    def get_csv(self, full_file_name, cols='empty'):
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
        return self.iteration_number

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
        _f, _ax = plt.subplots(1)
        for _molecule in self.standard_integration_table:
            _ax.scatter(self.standard_integration_table.index, self.standard_integration_table[_molecule])
            _ax.plot(self.calibration_curves.index, self.calibration_curves[_molecule])

        if log == True:
            _ax.set_yscale('log')
            _ax.set_xscale('log')

        plt.show()

    # Plot the intensity data for each sample over time
    def plot_calibrated_samples(self):
        _f, _axs = plt.subplots(number_of_samples, sharex=True, sharey=True)
        _i = -1
        for _sample in self.concentration_table:
            _i = _i + 1
            for _molecule in _sample:
                _axs[_i].plot(_sample.index, _sample[_molecule] / 1000)
                _axs[_i].set_ylabel('Concentration (mM)', fontsize=10)
                _axs[_i].set_xlabel('Iteration', fontsize=10)
                _axs[_i].set_title(sample_names[_i], fontsize=10)

        plt.show()

    # Plot the GC Elution Data
    def plot_all(self):
        _f, _axs = plt.subplots(number_of_samples, sharex=True, sharey=True)

        for _sample in range(number_of_samples):
            for _i in range(self.get_iteration_number()):
                _axs[_sample].plot(self.sample[_sample]['t' + str(_i)], self.sample[_sample]['i' + str(_i)])

        plt.show()


data_set = Data()
data_set.save_all()

print(data_set.get_standard_times())
#print(data_set.get_integration_table())
#print(data_set.get_concentration_table())
#data_set.plot_integration_table(True)
#data_set.plot_calibrated_samples()
#data_set.plot_all()
#print data_set.get_calibration_curves()
