# GC-Spectrum-Processing

Motivation
=========
This script is for importing a series of *.csv files taken on a gas chromatogram instrument. This instrument can analyze the chemical composition of a sample by evaporating it onto a column (essentially a really long tube) which is then linearly heated. A particular molecule will desorb and be detected at a characteristic time based on its boiling point- higher boiling point molecules desorb at higher temperatures and thus later times and lower boiling point molecules desorb at low temperature, and thus come off earlier in time on the spectrum. We characterize which molecules are present based on the elution time and their concentration, which is proportional to the area of the peak. 

The raw csv data consists of a header and series of times and intensities. Four samples are present which get repeatedly sampled in order (0, 1, 2, 3, 0, 1, etc) while the sample is illuminated with a laser pointer to detect how the amount of molecules will change over time.

Data Processing Workflow
=========================
1. Loop through all of the GC spectrum in the folder, skip the header and import the data, which is a time vs. intensity array as a pandas data frame (N x 2).

2. Sort the data frames based on the number of the sample. Since we have 4 samples, spectra that are N % 4 == 1 will all be from the first sample for example and we can use this to append the data frames to the right sample data frame (N x Number of Iterations per Sample). 

3. This will create a data set that we can address individually by calling data_set.get_sample(N), where N is the sample number.

4. The script will also load in calibration samples for both reactant and products to determine:

  a. What time a molecule will desorb from the column
  
  b. How the integrated intensity directly relates to concentration

5. Once the script determines the calibration curve, it will calculate an absolute concentration for each molecule in each sample as it is sampled over time. This allows us to create a plot to track molecule concentration over time for each sample.


Functions
===========

First call the class data_set = Data()

data_set.save_sample(0) this will save all spectrum from a particular sample as one csv in the folder of the script

