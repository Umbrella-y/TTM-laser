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
from save_file_with_timestamp import save_file_with_timestamp
from initiate_state import initiate_state
# 指明工作路径
workpath = r'F:\飞秒激光建模\2023-6-16-elec-注释\2023-6-15-rebuild/step1'
# 定义网格数
nx, ny = 100, 50
# 模拟体系x方向单位长度
x  = 2*224.33166413475473E-10 # A
dx = x/nx
# 模拟体系y方向单位长度
y  = 2*112.16583206737737E-10 #A
dy = y/ny
# 设置初始的电子温度场，后续统一使用grid进行保存
grid = initiate_state(ny, nx)   
ss_time = time.time()
for step in range(0,50000,1):
    print('Now at step = {}'.format(step))
    os.chdir(workpath)
    # 更改路径后开始执行lammps
    os.system('lmp.exe -in in.read_data.lmp -sf omp')
    #将生成的dump文件进行重命名
    data_file = workpath + '/' + 'dump0.atom'
    new_file_name = 'dump{}.atom'.format(step+1)
    new_file_path = os.path.join(os.path.dirname(data_file), new_file_name)
    os.rename(data_file, new_file_path)
    ##----------------------------------
    density_file = r'density-profile.txt'
    #save_file_with_timestamp(density_file,step)#------------------------------------------------save_stamp
    lattice_temp_file = r'lattice-temperature-profile.txt'
    #save_file_with_timestamp(lattice_temp_file,step)#------------------------------------------------save_stamp
    folder_path = r'F:\飞秒激光建模\2023-6-16-elec-注释\2023-6-15-rebuild/function'
    ## 将温度场传递给计算函数
    grid, Latti_temp = Electron_conduct_calculation(density_file, lattice_temp_file, grid, nx, ny,step,x,y)
    # 读入data文件，并开始根据所计算出的温度场
    file_path = os.path.join(folder_path, workpath + '/'+ "relaxed-500ps.data")
    #读取已有的data文件将其保存到df中
    df,df_velocity = reading_data(file_path)
    #根据前面读取的df的结果，将速度数据根据高斯分布计算出来
    df1,df2 = calculate_velocity(df,df_velocity,x=x,nx=nx,y=y,ny=ny,Latti_temp = Latti_temp)
    data_output(df1, df2, x, y, filename= 'test')

ff_time = time.time()
print('Total time is {}!!!'.format(ff_time-ss_time))
