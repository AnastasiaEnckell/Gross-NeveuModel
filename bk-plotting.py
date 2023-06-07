# This program plots a heatmap and a 3D plot of the probability density evolution with the Bravyi-Kitaev transformation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import seaborn as sns
import pandas as pd

n = 4 #number of qubits is 2^n

def matrix(x): # Bravyi-Kitaev matrix for the decoding
	if x==1:
		n1 = 1
		n2 = 0
		n3 = 1
		n4 = 1
		block = np.block([[n1, n2],[n3, n4]])
	elif x>1:
		n1 = matrix(x-1)
		n2 = np.block([[np.zeros((1,2**(x-1)))],[np.zeros((2**(x-1)-1,2**(x-1)))]] )
		n3 = np.block([[np.zeros((2**(x-1)-1,2**(x-1)))],[np.ones((1,2**(x-1)))]] )
		n4 = n1
		block = np.block([[n1, n2],[n3, n4]])
	return(block)

trans = matrix(n)
inverseabs = np.absolute(np.linalg.inv(trans))


with open("bravyi-kitaev.txt") as f:
	lines = f.readlines()

def write(lines): #write the data into arrays
	data = []
	x = []
	y = []
	pos = []
	posx = []
	index = []
	prob = []
	for i in range(0,len(lines)):
		s = lines[i].split('\n')
		data.append(s)
	for i in range(0, len(data)):
		x.append(data[i][0])
	for i in range(0, len(x)):
		y.append(x[i].split('\t'))
		pos.append(y[i][0].translate({ord(i): None for i in '[]'}))
		pos[i] = np.array(pos[i].split(','))
		if len(pos[i]) == 1:
			index.append(i)
		if len(pos[i])!= 1:
			prob.append(float(y[i][1]))
			posx.append(pos[i].astype(np.float))
	finaldata = []
	helping = np.zeros(2**n)
	check = posx[1]
	for i in range(0, len(posx)):
		intarray =  list(posx[i].astype(int))
		posx[i] = np.matmul(inverseabs, intarray)
		for j in range(0, len(posx[i])):
			if posx[i][j]%2 == 0 :
				posx[i][j] = 0.0
			if posx[i][j]%2 != 0 :
				posx[i][j] = 1.0
		count = 0
		for element in posx[i]:
			if element == 1:
				count = count + 1
		posx[i] = [element * prob[i] for element in posx[i]]
	for i in range(0, len(index)-1):
		ind1 = index[i]-i
		ind2 = index[i+1] - i-1
		helping = np.zeros(2**n)
		for j in range(ind1, ind2):
			helping = np.add(helping, posx[j])
		finaldata.append(helping)
	return finaldata

dataf = write(lines)

data_array = np.array(dataf)
x_data, y_data = np.meshgrid( np.arange(data_array.shape[1]), np.arange(data_array.shape[0]) )


# Heat map
fig = plt.figure()
df = pd.DataFrame(data_array)

ax = sns.heatmap(df, cmap="Blues_r")
ax.set_xlabel('x', fontsize=13)
ax.set_ylabel('t', fontsize=13)

plt.savefig("bravyi-kitaev_heatmap.png")


# 3D plot 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
z_data = data_array.flatten()
x_data = x_data.flatten() # number of sites
y_data = y_data.flatten() #timesteps
x_data1 = x_data[::2]
y_data1 = y_data[::2]
x_data2 = x_data[1::2]
y_data2 = y_data[1::2]
z_data1 = z_data[::2]
z_data2 = z_data[1::2]

ax.bar3d(x_data1-0.4, y_data1, np.zeros(len(z_data1)), 1, 1, z_data1, color='magenta',alpha=0.5)
ax.bar3d(x_data2-0.4, y_data2, np.zeros(len(z_data2)), 1, 1, z_data2, color='cyan',alpha=0.5)

ax.set_xlabel('x', fontsize=13)
ax.set_ylabel("t", fontsize=13)
ax.set_zlabel('ρ(x,t)', fontsize=13)

plt.savefig("bravyi-kitaev_3D.png")

