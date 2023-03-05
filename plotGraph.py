'''
Author: Lucas Nascimento;
Date: 04/03/2023;
Rev: 1;
UFF-RJ;
code for plot graphts
'''
from math import *
import matplotlib.pyplot as plt
import numpy as np

values_column_rk = []
values_column_ev = []
values_column_nf = []
'''
values_column_rk [0] = t
values_column_rk  [1] = y[0]
values_column_rk [2] = y[1]...
values_column_rk  [15] = y[14]
'''
#---------plots-------------------
def plotGraphsRK():
    #Read dates of file
    for i in range(0, 15):
        df = np.loadtxt("out_rk.dat")[:, i]
        values_column_rk.append(df)

    #plots all graphs in .png
    for i in range(1, 15):
        plt.figure(dpi=180)
        plt.plot (values_column_rk[0], values_column_rk[i], color='black')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$y[' + str(i) +']$')
        plt.savefig('./plots/' + 'y' + str(i-1) +'.png')
        #plt.show()
    
def plotGraphsNF():
    x = np.array([1, 2, 3, 4, 5, 6, 7])
    #Read dates of file
    for i in range(0, 7):
        df = np.loadtxt("out_ev.dat")[:, i]
        df2 = np.loadtxt("out_nf.dat")[:]
        values_column_ev.append(df)
        values_column_nf.append(df2)
    #plots all graphs in .png
    for i in range(0, 7):
        plt.close()
        plt.figure(dpi=360)
        plt.rcParams['text.usetex'] = True
        plt.text(5.94, 0.06, r'$ \omega_{n}='+ f' {values_column_nf[0][i]: .3f}  [Hz]' +' $')
        #plt.xlabel(r'$$')
        plt.ylabel(r'$ u[' + str(i) +']$')
        plt.plot (x, values_column_ev[i], color='black')
        plt.savefig('./plots/' + 'u' + str(i) +'.png')
        #plt.show()
        plt.close()

if __name__=='__main__':
    print('Write Datas Start...')
    plotGraphsRK()
    plotGraphsNF()
    print('...Write Datas End!')