import numpy as numpy
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
def initiate_state(ny, nx):
    # 初始化网格参数-此处的参数均可随机设置
    density = 0
    conductivity = 0
    heat_capacity = 314 #iniate-capacity
    length = 0
    width = 0
    lattice_temp = 0
    grid = numpy.empty((ny, nx), dtype=Grid)
    for i in range(ny):
        for j in range(nx):
            if i < 0.75*ny:
                temperature = 300  # y轴中部以下的网格温度为300摄氏度
            else:
                # 计算线性梯度
                gradient = (5000 - 300) / (ny / 2 )
                temperature = 300#300 + (i - ny / 2) * gradient  # y轴中部以上的网格具有线性梯度'''
            #temperature=ele_temp[i,j]
            grid[i, j] = Grid(density, conductivity, heat_capacity, length, width, temperature, lattice_temp)
    U_in = numpy.full((ny,nx), 300)
    return grid,U_in
  
    
   
