# -*- coding: UTF-8 -*-
import cupy as cp
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
from inputs import inputs

workpath = r'/media/zyy/data/2023-11-3-cupy/step1/'
savepath = r'/media/zyy/data/2023-11-3-cupy/step1/'
# 定义网格数
nx, ny = nx_and_ny_settings(workpath)
# 定义模拟体系大小，默认为中心对称体系，不可随意修改
xhi, yhi, zhi = x_y_z_size_settings(workpath)
#必要时可以手动设置，平时默认读取
# 模拟体系x方向单位长度
x  = 2* xhi*1E-10 # A
dx = x/nx
# 模拟体系y方向单位长度
y  = 2* yhi*1E-10 #A
dy = y/ny
z  = 2* zhi*1E-10 #A
# 设置初始的电子温度场，后续统一使用grid进行保存
grid,U_out = initiate_state(ny, nx)  
name_step = 0
for step in range(0,20000,1):
    ss_time = time.time()
    print('Now at step = {}'.format(step))
    os.chdir(workpath)
    print('LAMMPS now calculating')
    lammps_time0 = time.time()
    if step%10 == 0 and step != 0:
        name_step+=1
        os.system('mpiexec -n 4 lmp_mpi -in in.read_data.lmp -sf gpu -pk gpu 0 neigh no')
        data_file = savepath + '/' + 'dump0.atom'
        new_file_name = 'dump{}.atom'.format(name_step+1)
        new_file_path = os.path.join(os.path.dirname(data_file), new_file_name)
        os.rename(data_file, new_file_path)
    else:
        os.system('mpiexec -n 4 lmp_mpi  -in in.read_data_no_dump.lmp -sf gpu -pk gpu 0 neigh no -screen none')
    ##----------------------------------
    lammps_time1 = time.time()
    print("LAMMPS COST = {}".format(lammps_time1-lammps_time0))
    density_file = r'density-profile.txt'
    #save_file_with_timestamp(density_file,step)#------------------------------------------------save_stamp
    lattice_temp_file = r'lattice-temperature-profile.txt'
    #save_file_with_timestamp(lattice_temp_file,step)#------------------------------------------------save_stamp
    folder_path = r'/media/zyy/data/2023-11-3-cupy/function'
    ## 将温度场传递给计算函数
    print("Calculate heat condction")
    heat_time0= time.time()
    grid, Latti_temp, Original_Latti_temp, U_out = Electron_conduct_calculation(U_out, savepath, density_file, lattice_temp_file, grid, nx, ny,step,x,y,z)
    heat_time1= time.time()
    print("HEAT CONDUCTION COST = {}".format(heat_time1-heat_time0))
    # 读入data文件，并开始根据所计算出的温度场
    file_path = workpath + '/' + 'relaxed-500ps.data'
    #读取已有的data文件将其保存到df中
    print("Reading the Old data file now")
    reading_time0 = time.time()
    df,df_velocity = reading_data(file_path)
    reading_time1 = time.time()
    print("READING FILE COST = {}".format(reading_time1-reading_time0))
    #根据前面读取的df的结果，将速度数据根据高斯分布计算出来
    print("Calculate New Velocities")
    velocity_time0 = time.time()
    df1,df2 = calculate_velocity(df,df_velocity,x=x,nx=nx,y=y,ny=ny,
                                 Latti_temp = Latti_temp,Original_Latti_temp = Original_Latti_temp)
    velocity_time1 = time.time()
    print("VELOCITY COST = {}".format(velocity_time1-velocity_time0))
    print("Now out put files")
    output_time0 = time.time()
    data_output(df1,df2,x,y,z,filename= 'test')
    output_time1 = time.time()
    print("OUTPUT COST = {}".format(output_time1-output_time0))
    #break
    ff_time = time.time()
    print('This Cycle timecost is {}!!!'.format(ff_time-ss_time))
