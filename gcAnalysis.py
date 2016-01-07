import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

class Data:
	#This set of tasks gets called when initializing the class 
	def __init__(self):
		#Setting up a 2D empty list
		self.sample= [[] for i in range(4)]

		#Goes through each of the files in ascending order and converts a list of times and intensities into a dataframe
		for fname in sorted(glob.glob("./gc/*.txt")):
			#First get the name of the file and convert that to a number
			currentNumber = int(re.findall(r'\d+', fname)[-1])

			#For searching you need to add leading zeros so that it will go in order 01,02,...,10,11,etc
			if currentNumber < 10:
				padding = '0'
			else:
				padding = ''

			#Create column headers based on the number it pulls from the file name
			cols = ['t'+str(int((currentNumber-1)/4)),'i'+str(int((currentNumber-1)/4))]

			#Create the data frame if it is the first sample in the list
			#Otherwise concat the current data as a data frame with the already existing data frame for that sample
			if currentNumber <= 4:
				self.sample[currentNumber-1] = pd.read_csv(fname, header=22, delimiter=r"\s+", names=cols)
			elif currentNumber % 4 == 1:
				self.sample[0] = pd.concat([self.sample[0],pd.read_csv(fname, header=22, delimiter=r"\s+", names=cols)], axis = 1, join='inner')
			elif currentNumber % 4 == 2:
				self.sample[1] = pd.concat([self.sample[1],pd.read_csv(fname, header=22, delimiter=r"\s+", names=cols)], axis = 1, join='inner')
			elif currentNumber % 4 == 3:
				self.sample[2] = pd.concat([self.sample[2],pd.read_csv(fname, header=22, delimiter=r"\s+", names=cols)], axis = 1, join='inner')
			elif currentNumber % 4 == 0:
				self.sample[3] = pd.concat([self.sample[3],pd.read_csv(fname, header=22, delimiter=r"\s+", names=cols)], axis = 1, join='inner')

		self.iterationNumber = int(currentNumber / 4)
	#Display entire 2D List
	def getAll(self):
		return self.sample

	#Save One Sample, either 0, 1, 2, or 3 to csv file
	def saveSample(self,value):
		self.sample[value].to_csv(str(value)+'.csv', sep='\t', encoding='utf-8')

	#Returns a 2xN dataframe of the sample data
	def getSample(self,value):
		return self.sample[value]

	#Get the number of sampling 
	def getIterationNumber(self):
		return self.iterationNumber

dataSet = Data()
dataSet.saveSample(0)
dataSet.saveSample(1)
dataSet.saveSample(2)
dataSet.saveSample(3)



df = dataSet.getAll()
print len(df)
print df[0].shape
print df[1].shape
print df[2].shape
print df[3].shape

plt.close('all')
f, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=True)

for i in range(dataSet.getIterationNumber()):
	ax0.plot(df[0]['t'+str(i)],df[0]['i'+str(i)])
for i in range(dataSet.getIterationNumber()):
	ax1.plot(df[1]['t'+str(i)],df[1]['i'+str(i)])
for i in range(dataSet.getIterationNumber()):
	ax2.plot(df[2]['t'+str(i)],df[2]['i'+str(i)])
for i in range(dataSet.getIterationNumber()):
	ax3.plot(df[3]['t'+str(i)],df[3]['i'+str(i)])

plt.show()