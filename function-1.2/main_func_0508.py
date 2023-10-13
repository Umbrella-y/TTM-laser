import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from showcontourf import showcontourf
from Electron_conduct_calculation import Electron_conduct_calculation
from reading_data import reading_data
from calculate_velocity import calculate_velocity
from data_output import data_output
from initiate_state import initiate_state
# import lammps
# lammps()


#设置网格内的温度场信息，可以在后续使用外部文件进行导入
# #外部的信息格式应该为网格的x，y坐标、温度值信息，对应grid[i,j]的信息输入
# # 第一次初始化温度场
# 定义网格数
nx, ny = 20, 20
# 模拟体系x方向单位长度
x  = 404E-10 # A
dx = x/nx
# 模拟体系y方向单位长度
y  = 404E-10 #A
dy = y/ny
# 设置初始的电子温度场，后续统一使用grid进行保存
grid = initiate_state(nx, ny)   

density_file = r'density-profile.txt'
lattice_temp_file = r'lattice-temperature-profile.txt'

##--------------------------------------------------------------
gridsss, Latti_temp = Electron_conduct_calculation(density_file, lattice_temp_file, grid, nx, ny)
## 将温度场传递给计算函数
# 读入data文件，并开始根据所计算出的温度场
folder_path = r'F:\飞秒激光建模\python-realization'
file_path = os.path.join(folder_path, "relaxed-500ps.data")
#读取已有的data文件将其保存到df中
df = reading_data(file_path)
#根据前面读取的df的结果，将速度数据根据高斯分布计算出来
df1 = calculate_velocity(df,x=240E-10,nx=20,y=240E-10,ny=20,Latti_temp = Latti_temp)
data_output(df1,filename= 'test')



