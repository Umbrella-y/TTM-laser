import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
#用于定义初始的电子温度场信息，作为外部热源的输入方式之一
def initiate_state(ny, nx):
    # 初始化网格参数-此处的参数均可随机设置
    density = 0
    conductivity = 0
    heat_capacity = 0
    length = 0
    width = 0
    lattice_temp = 0
    grid = np.empty((ny, nx), dtype=Grid)
    for i in range(ny):
        for j in range(nx):
            if i < ny / 2:
                temperature = 300  # y轴中部以下的网格温度为300摄氏度
            else:
                # 计算线性梯度
                gradient = (3000 - 300) / (ny / 2 )
                temperature = 300 + (i - ny / 2) * gradient  # y轴中部以上的网格具有线性梯度
            grid[i, j] = Grid(density, conductivity, heat_capacity, length, width, temperature, lattice_temp)
    return grid

