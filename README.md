# GC-Spectrum-Processing

Overview
=========
This script is for importing a series of *.csv files taken on a gas chromatogram instrument. This instrument can analyze the chemical composition of a sample by evaporating it onto a column (essentially a really long tube) which is then linearly heated. Molecule's will desorb and be detected at a characteristic time based on its boiling point, higher boiling point molecules desorb at higher temperatures and thus later times and lower boiling point molecules desorb at low temperature, and thus come off earlier in time on the spectrum. We characterize which molecule's are present based on the elution time and their concentration, which is proportional to the area of the peak. 

The data consists of a header and series of times and intensities. Four samples are present which get repeatedly sampled in order (0, 1, 2, 3, 0, 1, etc) while the sample is illuminated with a laser pointer to detect how the concentration of photo-products evolve over time.
