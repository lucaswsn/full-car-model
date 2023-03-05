'''
Author: Lucas Nascimento;
Date: 04/03/2023;
Rev: 1;
UFF-RJ;

Numeric code runge-kutta 4th order for ODE system with fixal step
for problem: full-car model. 

'''
from math import *
import matplotlib.pyplot as plt
import numpy

#-----------------Parameters system---------------- 
Dim = 14 #Dimension of system

#Definition of variations
date = open('out_rk.dat', 'w') #output file
y = numpy.zeros(Dim) 
xt = numpy.zeros(Dim)
pi = 4*atan(1) #Defination of Pi
g = 9.81 #Defination of g

#RK parameters
t = 0.0 #inicial time of integration [s]
tf = 100 # final time of integration [s]
h = 0.01 #Range Kutta step [s]
tran = 0 #transient[s]
hh = h*0.5
h6 = h/6

#-----------------Initiation Values and parameters of ODE---------------- 
m = 840 #kg
mf = 53 #kg
mr = 76 #kg
Ix = 820 #kg m²
Iy = 11100 #kg m²

a1 = 1.4 #m
a2 = 1.47 #m
b1 = 0.7 #m
b2 = 0.75 #m

kf = 10000 #N/m
kr = 13000 #N/m
ktf = 1300 #N/m

cf = 940 #N.s/m
cr = 1450 #N.s/m

#Initiation Values:
y[0] = 0 # x [m]
y[1] = 0 #dx [m/s]
y[2] = 0 #phi [rad]
y[3] = 0 #dphi [rad/s]
y[4] = 0 #theta [rad]
y[5] = 0 #dtheta [rad/s]
y[6] = 0 #x1 [m]
y[7] = 0 #dx1 [m/s]
y[8] = 0 #x2 [m]
y[9] = 0 #dx2 [m/s]
y[10] = 0 #x3 [m]
y[11] = 0 #dx3 [m/s]
y[12] = 0 #x4 [m]
y[13] = 0 #dx4 [m/s]


#-----------function--------------------
def funcao(t, y):
    dy = numpy.zeros(Dim) #differentiation
    F = numpy.zeros(Dim) #force 

    for i in range(0, Dim):
     F[i] = -0.2*cos(t)    

    dy[0] = y[1]
    dy[1] = (1/m)*(-cf*(y[1]-y[7]+b1*y[3]-a1*y[5])-cf*(y[1]-y[9]-b2*y[3]-a1*y[5])-cr*(y[1]-y[11]-b1*y[3]+a2*y[5])-cr*(y[1]-y[13]+b2*y[3]+a2*y[5])-kf*(y[0]-y[6]+b1*y[2]-a1*y[4])-kf*(y[0]-y[8]-b2*y[2]-a1*y[4])-kr*(y[0]-y[10]-b1*y[2]+a2*y[4])-kr*(y[0]-y[12]+b2*y[2]+a2*y[4]))
    dy[2] = y[3]
    dy[3] = (1/Ix)*(-b1*cf*(y[1]-y[7]+b1*y[3]-a1*y[5])+b2*cf*(y[1]-y[9]-b2*y[3]-a1*y[5])+b1*cr*(y[1]-y[11]-b1*y[3]+a2*y[5])-b2*cr*(y[1]-y[13]+b2*y[3]+a2*y[5])-b1*kf*(y[0]-y[6]+b1*y[2]-a1*y[4])+b2*kf*(y[0]-y[8]-b2*y[2]-a1*y[4])+b1*kr*(y[0]-y[10]-b1*y[2]+a2*y[4])-b2*kr*(y[0]-y[12]+b2*y[2]+a2*y[4]))
    dy[4] = y[5]
    dy[5] = (1/Iy)*(a1*cf*(y[1]-y[7]+b1*y[3]-a1*y[5])+a1*cf*(y[1]-y[9]-b2*y[3]-a1*y[5])-a2*cr*(y[1]-y[11]-b1*y[3]+a2*y[5])-a2*cr*(y[1]-y[13]+b2*y[3]+a2*y[5])+a1*kf*(y[0]-y[6]+b1*y[2]-a1*y[4])+a1*kf*(y[0]-y[8]-b2*y[2]-a1*y[4])-a2*kr*(y[0]-y[10]-b1*y[2]+a2*y[4])-a2*kr*(y[0]-y[12]+b2*y[2]+a2*y[4]))
    dy[6] = y[7]
    dy[7] = (1/mf)*(cf*(y[1]-y[7]+b1*y[3]-a1*y[5])+kf*(y[0]-y[6]+b1*y[2]-a1*y[4])-ktf*(y[6]-F[1]))
    dy[8] = y[9]
    dy[9] = (1/mf)*(cf*(y[1]-y[9]-b2*y[3]-a1*y[5])+kf*(y[0]-y[8]-b2*y[2]-a1*y[4])-ktf*(y[8]-F[2]))
    dy[10] = y[11]
    dy[11] = (1/mr)*(cr*(y[1]-y[11]-b1*y[3]+a2*y[5])+kr*(y[0]-y[10]-b1*y[2]+a2*y[4])-ktf*(y[10]-F[3]))
    dy[12] = y[13]
    dy[13] = (1/mr)*(cr*(y[1]-y[13]+b2*y[3]+a2*y[5])+kr*(y[0]-y[12]+b2*y[2]+a2*y[4])-ktf*(y[12]-F[4]))

    return dy

#-----------integrator---------------
def integrator():
    global t, tf, y, tran, hh, h, h6

    while (t <= tf):
        writeDate() #call funtion of Write in the file
        th = t + hh

        k1 = funcao(t, y) #Calculate of K1
        for i in range(0, Dim):
            xt[i] = y[i] + hh*k1[i]
        
        k2 = funcao(th, xt) #Calculate of K2
        for i in range(0, Dim):
            xt[i] = y[i] + hh*k2[i]
           
        k3 = funcao(th, xt) #Calculate of K3
        for i in range(0, Dim):
            xt[i] = y[i] + h*k3[i]
           
        k4 = funcao(t + h, xt) #Calculate of K4
        for i in range(0, Dim):
            y[i] = y[i] + h6 * (k1[i] + 2.0*(k2[i]+k3[i]) + k4[i] ) #Calculate new value

        t = t + h

def writeDate():
    if (t >= tran): #Write befor of transiant
        date.write(f'{t:.15f}' + '\t')
        for i in range(0, Dim): 
            date.write(f'{y[i]:.15f}' + '\t')
        date.write('\n')

if __name__=='__main__':
    print('Runge Kutta 4th Start...')
    integrator()
    print('...Runge Kutta 4th End!')