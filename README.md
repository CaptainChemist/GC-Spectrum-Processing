# GC-Spectrum-Processing

Motivation
=========
This script is for importing a series of *.csv files taken on a gas chromatogram instrument. This instrument can analyze
 the chemical composition of a sample by evaporating it onto a column (essentially a really long tube) which is then
 linearly heated. A particular molecule will desorb and be detected at a characteristic time based on its boiling point-
  higher boiling point molecules desorb at higher temperatures and thus later times and lower boiling point molecules
  desorb at low temperature, and thus come off earlier in time on the spectrum. We characterize which molecules are
  present based on the elution time and their concentration, which is proportional to the area of the peak.

The raw csv data consists of a header and series of times and intensities. Four samples are present which get repeatedly
 sampled in order (0, 1, 2, 3, 0, 1, etc) while the sample is illuminated with a laser pointer to detect how the amount
 of molecules will change over time.

Data Processing Workflow
=========================
# Creating the Calibration
1. Loop through the controls folder and store each spectrum as a dataframe within a list called **self.control**. Each
spectrum has the properties **name** and **type** which facilitates further processing later. We extract the molecule
name and concentration based on the naming scheme of the file, i.e. 500_benzil.txt would be 500 µM benzil control
for example.
2. Loop through **self.control** and find each unique target molecule. Utilize a peak detect function to find the first
non-solvent peak and record its elution time in a dictionary **self.standard_times**.
3. Create an integration table by looping through **self.control** and integrating with a window that corresponds to its
 elution time and the **integration_range**. These integrations are entered in to a table called
self.standard_integration_table and we subtract off any intensity from baseline with 0 µM blanks. The 0 µM blanks are
not stored in the final table.
4. Create a calibration curve for each target molecule by fitting the concentration table data to a second order
polynomial. Store it as **self.calibration_curves**.
5. If desired, we may plot the calibration curve by calling the function **plot_integration_table()**.

# Aggregating and Processing the Sample Data
1. Loop through all of the GC spectrum in the folder, skip the header and import the data, which is a time vs. intensity
 array as a pandas data frame (N x 2). Sort the data frames based on the number of the sample. Since we have 4 samples,
 spectra that are N % 4 == 1 will all be from the first sample for example and we can use this to append the data frames
  to the right sample data frame (N x Number of Iterations per Sample). **number_of_samples** may be customized to suit
  the experiment.
2. This will create a data set that we can address individually by calling data_set.get_sample(N), where N is the sample
 number.
3. A list of dataframes called **intensity_table** is created which has a list of data frames for each sample. Each
dataframe has the integrated intensity for each target molecule for each run of the sample. One the intensity for a
target molecule and run is found, we use **np.argmin** to find the corresponding index value on the calibration curve
for that target molecule and then find the corresponding concentration at that index. We store the calculated
concentration for that sample, run and target molecule as a list of dataframes **self.concentration_table**.
4. We may plot the concentration versus time for all samples by calling **plot_calibrated_samples()**.


Public API Function Calls
===========

First call the class **data_set = Data()**

**data_set.save_sample(sample_number)** - saves one sample with all iterations to a csv file

**data_set.save_all()** - saves all samples as individual csv files

**data_set.get_sample(sample_number)** - returns one sample as a dataframe

**data_set.get_all()** - returns a list of dataframes for each sample

**data_set.get_iteration_number()** - returns how many iterations of each sample are in the folder

**data_set.get_integration_table()** - returns a dataframe of samples vs. concentration

**data_set.get_standard_times()** - returns a dictionary of elution times vs. molecule

**data_set.get_concentration_table()** - returns a list of dataframes for each sample and the corresponding
products as a function of iteration number.

**data_set.plot_integration_table()** - plots integrated intensity vs. concentration for target molecules

**data_set.plot_calibrated_samples()** - creates a plot for each sample which shows product concentrations as a
function of iteration number.

**data_set.plot_all()** - plots all raw GC plots
