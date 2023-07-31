# -*- coding: UTF-8 -*-
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

workpath = r'/media/zyy/data1/2023-7-13-1mW/step1'
# 定义网格数
nx, ny = 40, 20
# 模拟体系x方向单位长度
x  = 2*80.958E-10 # A
dx = x/nx
# 模拟体系y方向单位长度
y  = 2*40.479E-10 #A
dy = y/ny
z  = 2* 4.0479E-10 #A

# 设置初始的电子温度场，后续统一使用grid进行保存
grid = initiate_state(ny, nx)   
ss_time = time.time()
name_step = 0
for step in range(0,5000,1):
    print('Now at step = {}'.format(step))
    os.chdir(workpath)
    if step%1 == 0 :
        name_step+=1
        ######## houxu jiashang if    yong linagge lmp_Mpi 
        os.system('lmp_mpi -in in.read_data.lmp')
        data_file = workpath + '/' + 'dump0.atom'
        new_file_name = 'dump{}.atom'.format(name_step+1)
        new_file_path = os.path.join(os.path.dirname(data_file), new_file_name)
        os.rename(data_file, new_file_path)
    else:
        os.system('lmp_mpi -in in.read_data_no_dump.lmp')
    ##----------------------------------
    density_file = r'density-profile.txt'
    #save_file_with_timestamp(density_file,step)#------------------------------------------------save_stamp
    lattice_temp_file = r'lattice-temperature-profile.txt'
    #save_file_with_timestamp(lattice_temp_file,step)#------------------------------------------------save_stamp
    folder_path = r'/media/zyy/data1/2023-7-13-1mW/function'
    ## 将温度场传递给计算函数
    grid, Latti_temp, Original_Latti_temp = Electron_conduct_calculation(density_file, lattice_temp_file, grid, nx, ny,step,x,y,z)
    # 读入data文件，并开始根据所计算出的温度场
    file_path = os.path.join(folder_path, workpath + '/'+ "relaxed-500ps.data")
    #读取已有的data文件将其保存到df中
    df,df_velocity = reading_data(file_path)
    #根据前面读取的df的结果，将速度数据根据高斯分布计算出来
    df1,df2 = calculate_velocity(df,df_velocity,x=x,nx=nx,y=y,ny=ny,
                                 Latti_temp = Latti_temp,Original_Latti_temp = Original_Latti_temp)
    data_output(df1,df2,x,y,z,filename= 'test')
    #break

ff_time = time.time()
print('Total time is {}!!!'.format(ff_time-ss_time))
