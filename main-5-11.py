# -*- coding: UTF-8 -*-
import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import multiprocessing
from showcontourf import showcontourf
from Electron_conduct_calculation import Electron_conduct_calculation
from reading_data import reading_data
from calculate_velocity import calculate_velocity
from data_output import data_output
from save_file_with_timestamp import save_file_with_timestamp
from initiate_state import initiate_state

workpath = r'E:/2023-9-17-hfxy/step1/'
# 定义网格数
nx, ny = 50, 100
# 模拟体系x方向单位长度
x  = 2*80.958E-10 # A
dx = x/nx
# 模拟体系y方向单位长度
y  = 2*121.437E-10 #A
dy = y/ny
z  = 2* 4.0479E-10 #A

# 设置初始的电子温度场，后续统一使用grid进行保存
grid = initiate_state(ny, nx)   

name_step = 0
for step in range(0,5000,1):
    ss_time = time.time()
    print('Now at step = {}'.format(step))
    os.chdir(workpath)
    print('LAMMPS now calculating')
    if step%10 == 0 and step != 0:
        name_step+=1
        ######## houxu jiashang if    yong linagge lmp_Mpi 
        os.system('mpiexec -np 4 lmp -in in.read_data.lmp')
        data_file = workpath + '/' + 'dump0.atom'
        new_file_name = 'dump{}.atom'.format(name_step+1)
        new_file_path = os.path.join(os.path.dirname(data_file), new_file_name)
        os.rename(data_file, new_file_path)
    else:
        os.system('mpiexec -np 4 lmp -in in.read_data_no_dump.lmp')
    ##----------------------------------
    density_file = r'density-profile.txt'
    #save_file_with_timestamp(density_file,step)#------------------------------------------------save_stamp
    lattice_temp_file = r'lattice-temperature-profile.txt'
    #save_file_with_timestamp(lattice_temp_file,step)#------------------------------------------------save_stamp
    folder_path = r'E:/2023-9-17-hfxy/function'
    ## 将温度场传递给计算函数
    print("Calculate heat condction")
    grid, Latti_temp, Original_Latti_temp = Electron_conduct_calculation(density_file, lattice_temp_file, grid, nx, ny,step,x,y,z)
    # 读入data文件，并开始根据所计算出的温度场
    file_path = os.path.join(folder_path, workpath + '/'+ "relaxed-500ps.data")
    #读取已有的data文件将其保存到df中
    print("Reading the Old data file now")
    df,df_velocity = reading_data(file_path)
    #根据前面读取的df的结果，将速度数据根据高斯分布计算出来
    print("Calculate New Velocities")
    df1,df2 = calculate_velocity(df,df_velocity,x=x,nx=nx,y=y,ny=ny,
                                 Latti_temp = Latti_temp,Original_Latti_temp = Original_Latti_temp)
    print("Now out put files")

    #file_path = 'E:/2023-9-17-hfxy/step1'
    #os.chmod(file_path, 0o755)

    data_output(df1,df2,x,y,z,filename= 'test')
    #break
    ff_time = time.time()
    print('This Cycle timecost is {}!!!'.format(ff_time-ss_time))
