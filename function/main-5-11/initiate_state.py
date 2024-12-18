import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
def initiate_state(nx,ny):
    # 初始化网格参数-此处的参数均可随机设置
    density = 0
    conductivity = 0
    heat_capacity = 0
    length = 0
    width = 0
    lattice_temp = 0
    grid = np.empty((nx, ny), dtype=Grid)
    for i in range(nx):
        for j in range(ny):
                if i >= 0.3*nx and i <= 0.5*nx and j >= 0.3*ny and j <= 0.5*ny:
                    temperature =  300  # 中心网格温度为100摄氏度
                else:
                    temperature = 300  # 其他网格温度为20摄氏度
                grid[i, j] = Grid(density, conductivity, heat_capacity, length, width, temperature, lattice_temp)
    return grid