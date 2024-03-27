import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math

def showcontourf(mat,D,cmap=plt.cm.get_cmap('jet'),fsize=(12,12),vmin=0, vmax=100, timestep1 =0, name = str(11)):
    plt.clf()
    levels = np.arange(vmin,vmax,1)
    x=np.linspace(D[0],D[1],mat.shape[1])
    y=np.linspace(D[2],D[3],mat.shape[0])
    X,Y=np.meshgrid(x,y)
    z_max = np.max(mat)
    i_max,j_max = np.where(mat==z_max)[0][0], np.where(mat==z_max)[1][0]
    show_max = "U_max: {:.1f}".format(z_max)
    plt.plot(x[j_max],y[i_max],'ro') 
    plt.contourf(X, Y, mat, 100, cmap = cmap, origin = 'lower', levels = levels)
    plt.annotate(show_max,xy=(x[j_max],y[i_max]),xytext=(x[j_max],y[i_max]),fontsize=14)
    plt.colorbar()
    plt.xlabel('x', fontsize=20)
    plt.ylabel('y', fontsize=20)
    plt.axis('equal')
    #plt.draw()
    #plt.pause(0.1)
    plt.savefig('./{}{}temp.jpg'.format(str(name),str(timestep1)))
    #plt.show()
    plt.clf()