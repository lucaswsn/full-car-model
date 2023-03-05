'''
Author: Lucas Nascimento;
Data: 05/03/2023;
Rev: 1;
UFF-RJ;

numeric code for calculate eingevalues and eingevector for problem:
full car model
'''
from math import *
import numpy as np
from scipy import linalg

date_nf = open('out_nf.dat', 'w') #output file
date_ev = open('out_ev.dat', 'w') #output file

pi = 4*atan(1) #defination of pi

#-----------------Parameters of system---------------- 
m = 840 #kg
mf = 53 #kg
mr = 76 #kg
Ix = 820 #kg m²
Iy = 1100 #kg m²

a1 = 1.4 #m
a2 = 1.47 #m
b1 = 0.7 #m
b2 = 0.75 #m

kf = 10000 #N/m
kr = 13000 #N/m
ktf = 200000 #N/m
ktr = ktf #N/m

cf = 0 #N.s/m
cr = 0 #N.s/m

#-----------Assembley of matrix K------------
k11 = 2*kf+2*kr
k12 = b1*kf-b2*kf-b1*kr+b2*kr; k21 = k12
k31 = 2*a2*kr-2*a1*kf; k13=k31
k22 = b1*b1*kf+b2*b2*kf+b1*b1*kr+b2*b2*kr
k32 = a1*b2*kf-a1*b1*kf-a2*b1*kr+a2*b2*kr; k23=k32
k42 = -b1*kf ; k24=k42
k52 = b2*kf; k25=k52
k33 = 2*kf*a1*a1+2*kr*a2*a2
k44 = kf+ktf; k55 = k44

M = np.diag([m, Ix, Iy, mf, mf, mr, mr])
K = np.array([[k11, k12, k13, -kf, -kf, -kr,-kr],
              [k21, k22, k23, k24, k25, b1*kr,-b2*kr],
              [k31, k32, k33, a1*kf, a1*kf, -a2*kr, -a2*kr],
              [-kf, k42, a1*kf, k44, 0, 0, 0],
              [-kf, k52, a1*kf, 0, k55, 0, 0],
              [-kr, b1*kr, -a2*kr, 0, 0, kr+ktr, 0],
              [-kr, -b2*kr, -a2*kr, 0, 0, 0, kr+ktr]])


#-----------Evalueted Eigenvalue and Eingevector------------
def eigenvalueAndeigenvector():
    M_inv = linalg.inv(M)
    eigenvalue, eigenvector = linalg.eig(M_inv @ K)
    #eigenvalue = np.sort(eigenvalue)
    frequency = (np.sqrt(eigenvalue)) / (2*pi)
    writeDate(frequency, eigenvector) #Hz


#-----------Write in the files------------
def writeDate(frequency, eigenvector):
    Dim = len(frequency)

    #date_nf.write(f'Natural Frequency:' + '\n')
    for i in range(0, Dim): 
        date_nf.write(f'{frequency[i].real: .15f}' + '\n') #Hz display real part only
    date_nf.write('\n') 

    #date_ev.write(f'Eigenvectors:' + '\n')
    for i in range(0, Dim): 
        for j in range (0, Dim):
            date_ev.write(f'{eigenvector[i][j]: .15f}' + '\t') #Hz
        date_ev.write('\n')    


if __name__=='__main__':
    print('Calculating Natural Frequency...')
    eigenvalueAndeigenvector()
    print('End!')