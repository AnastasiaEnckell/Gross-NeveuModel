import numpy as np
import csv
"""
This is a code for finding the indexes for the Pauli gates of the Bravyi-Kitaev transformation. 
Two output files can be used in the program optimizedbk.jl 
"""

n = 4 #number of qubits is 2^n

def matrix(x): #recursive function for the Bravyi-Kitaev matrix
	if x==1:
		n1 = 1
		n2 = 1
		n3 = 0
		n4 = 1
		block = np.block([[n1, n2],[n3, n4]])
	elif x>1:
		n1 = matrix(x-1)
		n2 = np.block([[np.ones((1,2**(x-1)))],[np.zeros((2**(x-1)-1,2**(x-1)))]] )
		n3 = np.block([[np.zeros((1,2**(x-1)))],[np.zeros((2**(x-1)-1,2**(x-1)))]] )
		n4 = n1
		block = np.block([[n1, n2],[n3, n4]])
	return(block)


def parityandflipset(matrix, j): # find the elements of the parity set P(j) or flip set F(j)
	matrix = matrix - np.identity(2**n)
	j = len(matrix)-1-j
	st = []
	for i in range(0, len(matrix)):
		if matrix[j][i] == 1 :
			st.append(len(matrix)- 1 -i)
	return st


def updateset(matrix, j): # find the elements of the update set U(j)
	matrix = matrix - np.identity(2**n)
	j = len(matrix)-1-j
	st = []
	for i in range(0, len(matrix)):
		if matrix[i][j] == 1 :
			st.append(len(matrix)- 1 -i)
	return st


def remainderset(parity, flip): # find the elements of the remainder set R(j)
	c = [element for element in parity if element not in flip] 
	return c


def rhoset(parity, remainderset, j): #find the elements of the rho set rho(j)
	if j %2 ==0 :
		rho = parity
	else :
		rho = remainderset
	return rho


def overlap(set1, set2): # function for (set1 union set2)\(set1 intersection set2)
	x = set1
	y = set2
	c = [element for element in set1 if element not in set2] 
	d = [element for element in set2 if element not in set1] 
	return c+d


def cut(set1, set2): # function for (set1 intersection set2)
	x = []
	for i in set1:
		if i in set2:
			x.append(i)
	return x


def but(set1, set2): # function for (set1\set2)
	x = set1
	y = set2
	c = [element for element in x if element not in y] 
	return c


def writetheindices(n): # function writes the indices of the exitation operators for 2^n qubits
	ind1plus = [[ 0,  6], [ 1 , 7]] # array for the H_1 with plus sign
	ind1minus = [[ 2  ,4],[ 3 , 5]] # array for the H_1 with minus sign
	ind2minus = [[ 0 , 4],[ 1 , 5]] # array for the H_1 with plus sign
	ind2plus = [[ 2 , 6],[ 3 , 7]] # array for the H_1 with minus sign
	for i in range(0, int((2**n)/2)-2):
		help1 = (ind1plus[i][0]+4) % (2**n)
		help2 = (ind1plus[i][1]+4) % (2**n)
		if help1<help2:
			ind1plus.append([help1, help2])
			ind2plus.append([(ind2plus[i][0]+4) % (2**n), (ind2plus[i][1]+4) % (2**n)])
			ind1minus.append([(ind1minus[i][0]+4) % (2**n), (ind1minus[i][1]+4) % (2**n)])
			ind2minus.append([(ind2minus[i][0]+4) % (2**n), (ind2minus[i][1]+4) % (2**n)])
		else:
			ind1plus.append([help2, help1])
			ind2plus.append([(ind2plus[i][1]+4) % (2**n), (ind2plus[i][0]+4) % (2**n)])
			ind1minus.append([(ind1minus[i][1]+4) % (2**n), (ind1minus[i][0]+4) % (2**n)])
			ind2minus.append([(ind2minus[i][1]+4) % (2**n), (ind2minus[i][0]+4) % (2**n)])
	return ind1plus, ind1minus, ind2plus, ind2minus


# Function positions(arr) creates 2 dimensional list of lists where every row has the indices of [[Pauli X],[Pauli Y],[Pauli Z]]. 
# The addition -100 means that the coefficient of the exitation operator is negative. 
def positions(arr):
	xyz12 =[]
	for i in range(0, len(arr)):
		ii = arr[i][0]
		jj = arr[i][1]
		Ui = updateset(bkmatrix, ii)
		Uj = updateset(bkmatrix, jj)
		Uij = overlap(Ui, Uj)
		Pi = parityandflipset(parityinverse, ii)
		Pj = parityandflipset(parityinverse, jj)
		Fi = parityandflipset(bkinverseabs, ii)
		Fj = parityandflipset(bkinverseabs, jj)
		Ri = remainderset(Pi, Fi)
		Rj = remainderset(Pj, Fj)
		P0ij = overlap(Pi, Pj)
		P1ij = overlap(Pi, Rj)
		P2ij = overlap(Ri, Pj)
		P3ij = overlap(Ri, Rj)
		aij = cut(Ui, Pj)
		if arr[i][0]%2 == 0:
			
			x = but(Uij, aij) + [ii]
			y = aij + [jj]
			z = but(P0ij, aij)
			xyz12.append([x,y,z])
			x = but(Uij, aij) + [jj]
			y = aij + [ii]
			z = but(P0ij, aij) + [-100] 
			xyz12.append([x,y,z])
		
		if arr[i][0]%2 != 0:
			if ii in Pj and jj in Ui:

				x = but(Uij, [jj]) + [ii]
				y =  []
				z = but(P2ij, [ii])
				xyz12.append([x,y,z])

				x = but(Uij, [jj]) + [ii]
				y =  []
				z = P1ij+[jj]+ [-100]
				xyz12.append([x,y,z])


			if ii not in Pj and jj in Ui:

				x = but(Uij, [jj]+aij) 
				y = [ii]+ aij
				z = but(P2ij, aij)+ [-100]
				xyz12.append([x,y,z])


				x = but(Uij, [jj]) + [ii]
				y =  []
				z = P1ij+[jj]+ [-100]
				xyz12.append([x,y,z])

			if ii in Pj and jj not in Ui:	
		
				x = Uij
				y = [ii] + [jj]
				z = but(P1ij, [ii])
				xyz12.append([x,y,z])
				
				x = Uij+ [jj]+[ii]
				y = []
				z = but(P2ij, [ii])
				xyz12.append([x,y,z])

			if ii not in Pj and jj not in Ui:

				x = but(Uij, aij) + [ii]
				y = aij + [jj]
				z = but(P1ij, aij)
				xyz12.append([x,y,z])

				x = but(Uij, aij) + [jj] 
				y = aij + [ii]
				z = but(P2ij, aij)+ [-100]
				xyz12.append([x,y,z])
	return xyz12


def find_list_index(lst, sublist):
    for i, item in enumerate(lst):
        if item == sublist:
            return i
    return -1

def delete_list(lst, sublist):
    index = find_list_index(lst, sublist)
    if index != -1:
        del lst[index]

def numberoperator(n): #function returns the indices of the Pauli-Z corresponding to the number operators
	plus = []
	minus = []
	numbopplus = []
	numbopminus = []
	for i in range(0, int(2**n / 4)): #the term with the coefficient (ma+r)
		plus.append(0+4*i)
		plus.append(1+4*i)
		minus.append(2+4*i)
		minus.append(3+4*i)
	for i in range(0, len(plus)):
		variable = parityandflipset(bkinverseabs, plus[i]) + [plus[i]]
		variable = [n+1 for n in variable]
		numbopplus.append(variable)
	for i in range(0, len(minus)):
		variable = parityandflipset(bkinverseabs, minus[i]) + [minus[i]]
		variable = [n+1 for n in variable]
		numbopminus.append(variable)
	return numbopplus, numbopminus


def exchangeoperator(n): #function returns the indices of the Pauli-Z corresponding to the exchange operators of the form a_j a^dagger_j * a_i a^dagger_i
	plusc = [[0,2],[0,3],[1,2], [1,3]] #indices with plus sign. For example, the fist array [0,2] coorresponds to a_0 a^dagger_0 * a_2 a^dagger_2
	minusc = [[0,1], [2,3]]
	z1 = []
	z2 = []
	for j in range(0,int(2**n/4)):
		for i in range(0, len(plusc)):
			fi = parityandflipset(bkinverseabs, plusc[i][0])+[plusc[i][0]]
			fj = parityandflipset(bkinverseabs, plusc[i][1])+[plusc[i][1]]
			fij = overlap(fi, fj)
			#print(fi)
			z1.append([fi,fj,fij])
			plusc[i][0] = plusc[i][0]+4
			plusc[i][1] = plusc[i][1]+4


	for j in range(0,int(2**n/4)):
		for i in range(0, len(minusc)):
			#print(plusc[i][0])
			fi = parityandflipset(bkinverseabs, minusc[i][0])+[minusc[i][0]]
			fj = parityandflipset(bkinverseabs, minusc[i][1])+[minusc[i][1]]
			fij = overlap(fi, fj)
			z2.append([fi,fj,fij])
			minusc[i][0] = minusc[i][0]+4
			minusc[i][1] = minusc[i][1]+4

	for i in range(0, len(z1)):
		for j in range(0, 3):
			z1[i][j] = [n+1 for n in z1[i][j]]

	for i in range(0, len(z2)):
		for j in range(0, 3):
			z2[i][j] = [n+1 for n in z2[i][j]]
	return z1, z2


#function combinetheoperators combines the operators that only contribute the Pauli-Z gates and returns a list of lists 
#where the negative numbers correspond to the coefficients: -1 for (ag^2)/2, -2 for -(ag^2)/2, -3 for (m*a+r)+(g^2)/2 and -4 for -(m*a+r)+(g^2)/2
def combinetheoperators(n):
	selfint1, selfint2 = numberoperator(n)
	z1, z2 = exchangeoperator(n)
	thefull = [[], [], [], []] #the empty list of lists for the indices of the coefficients [[(ag^2)/2], [(m*a+r)], [-(ag^2)/2], [-(m*a+r)]]
	for i in range(0, len(z1)):
		thefull[0].append(z1[i][0])
		thefull[0].append(z1[i][1])
		thefull[2].append(z1[i][2])
	for element in selfint1:
		thefull[1].append(element)
	for i in range(0, len(z2)):
		thefull[2].append(z2[i][0])
		thefull[2].append(z2[i][1])
		thefull[0].append(z2[i][2])
	for element in selfint2:
		thefull[3].append(element)
	g2aplus = [[],[]]
	g2aminus = [[],[]]
	maplus = [[],[]]
	maminus = [[],[]]


	for j in range(0, len(thefull[0])):
		if thefull[0][j] in thefull[0][j+1:]:
			g2aplus[1].append(thefull[0][j])
		if thefull[0][j] not in g2aplus[1]:
			g2aplus[0].append(thefull[0][j])

	for j in range(0, len(thefull[2])):
		if thefull[2][j] in thefull[2][j+1:]:
			g2aminus[1].append(thefull[2][j])
		if thefull[2][j] not in g2aminus[1]:
			g2aminus[0].append(thefull[2][j])

	for element in g2aplus[0]:  
		if element in g2aminus[0]:
			delete_list(g2aminus[0], element)
			delete_list(g2aplus[0], element)
	todelete = []
	for i in range(0, len(g2aplus[1])): 
		if g2aplus[1][i] in g2aminus[0]: 
			delete_list(g2aminus[0], g2aplus[1][i]) 
			if g2aplus[1][i] not in g2aplus[0]: 
				g2aplus[0].append(g2aplus[1][i]) 
				todelete.append(g2aplus[1][i])

	for element in todelete:
		if element in g2aplus[1]:
			delete_list(g2aplus[1], element)
	for element in thefull[1]:
		if element not in g2aplus[0]:
			maplus[0].append(element)
	if element in g2aplus[0]:
		maplus[1].append(element)
		delete_list(g2aplus[0], element)

	for element in thefull[3]:
		if element not in g2aplus[0]:
			maminus[0].append(element)
		if element in g2aplus[0]:
			#print("g2aplus[0]", element)
			maminus[1].append(element)
			delete_list(g2aplus[0], element)

	combine = []
	for i in range(0, len(g2aplus[0])):
		g2aplus[0][i] = g2aplus[0][i] + [-1]
		combine.append(g2aplus[0][i])
	for i in range(0, len(g2aminus[0])):
		g2aminus[0][i] = g2aminus[0][i] + [-2]
		combine.append(g2aminus[0][i])
	for i in range(0, len(maplus[1])):
		maplus[1][i] = maplus[1][i] + [-3]
		combine.append(maplus[1][i])
	for i in range(0, len(maminus[1])):
		maminus[1][i] = maminus[1][i] + [-4]
		combine.append(maminus[1][i])
	

	return combine

def optimize(l): #function returns the list in the lexicographic ordering
	helping = np.full((len(l), 2**n), 4) 
	if isinstance(l[0][0], list) :
		for i in range(0, len(l)):
			#print(l[i])
			for j in range(0, 3):
				for element in l[i][j]:
					if element > 0:
						helping[i][element-1]=j+1 
	else:
		for i in range(0, len(l)):
			for element in l[i]:
				if element > 0:
					helping[i][element-1]=3
	s = []
	for element in helping:
		s.append(int(''.join(str(x) for x in element)))

	sort_index = np.argsort(s)
	s.sort()
	newresult= []
	for i in range(0, len(l)):
		newresult.append(l[sort_index[i]])
	return newresult




parity = np.triu(np.ones((2**n,2**n)), k = 0) #parity matrix 
bkmatrix = matrix(n) # Bravyi-Kitaev matrix
bkinverse = np.linalg.inv(bkmatrix) # inverse BK matrix
bkinverseabs = np.absolute(np.linalg.inv(bkmatrix)) # absolute value of the inverse BK matrix
parityinverse = np.matmul(parity-np.identity(2**n), bkinverse) # matrix for the parity set


H1plus, H1minus, H2plus, H2minus = writetheindices(n) #the indices of the exitation operators



res1 = positions(H1plus)
res2 = positions(H1minus)
res3r = positions(H2plus)
res4r = positions(H2minus)

result = []

for i in range(0, len(res1)): #put all the indices of exitation operators in one list
	res1[i][2] = res1[i][2] + [-2]
	res2[i][2] = res2[i][2] + [-3]
	res3r[i][2] = res3r[i][2] + [-4]
	res4r[i][2] = res4r[i][2] + [-5]
	for j in range(0, 3):
		res1[i][j] = [n+1 for n in res1[i][j]]
		res2[i][j] = [n+1 for n in res2[i][j]]
		res3r[i][j] = [n+1 for n in res3r[i][j]]
		res4r[i][j] = [n+1 for n in res4r[i][j]]
	result.append(res1[i])
	result.append(res2[i])
	result.append(res3r[i])
	result.append(res4r[i])




combined = combinetheoperators(n) # combined indices of the operators correspondind to the Pauli-Z

optimized_z = optimize(combined)
optimized_xyz = optimize(result)


leng = 0
for n in range(0,len(optimized_xyz)):
	lengnew = max(len(l) for l in  optimized_xyz[n])
	if lengnew > leng:
		leng = lengnew


for i in range(0, len(optimized_xyz)):
	for j in range(0, 3):
		if len(optimized_xyz[i][j])<leng:
			for k in range(len(optimized_xyz[i][j]), leng):
				optimized_xyz[i][j].append(0)



lengnew = max(len(l) for l in  optimized_z)
if lengnew > leng:
	leng = lengnew


for i in range(0, len(optimized_z)):
		if len(optimized_z[i])<leng:
			for k in range(len(optimized_z[i]), leng):
				optimized_z[i].append(0)



with open("xyzdata.csv", "w", newline='') as file: 
	writer = csv.writer(file)
	for i in range(0, len(optimized_xyz)):
		for j in range(0, 3):
			writer.writerow(optimized_xyz[i][j])


with open("zdata.csv", "w", newline='') as file: 
	writer = csv.writer(file)
	for i in range(0, len(optimized_z)):
		writer.writerow(optimized_z[i])







