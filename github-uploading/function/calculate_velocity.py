import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
# 该函数是用于计算由传热方程计算出的晶格温度所对应的原子速度并赋还给原子本身的函数
def calculate_velocity(df,df_velocity,x,nx,y,ny,Latti_temp):
    # 计算在输入的data文件中的每一个部分的长度，分为x和y两个方向
    part_sizex = x / nx
    part_sizey = y / ny
    # 将df中存储的值从字符串格式转化为数值格式
    df['x'] = pd.to_numeric(df['x'])
    df['y'] = pd.to_numeric(df['y'])
    # 使用 cut() 函数将DataFrame中的数据划分为网格
    df['x_grid'] = pd.cut(df['x'], nx)
    df['y_grid'] = pd.cut(df['y'], ny)
    # 获取 cut() 函数使用的标签
    x_labels = pd.cut(df['x'], nx).unique()
    y_labels = pd.cut(df['y'], ny).unique()
    # 按 x 和 y 网格对DataFrame进行分组
    # 需要注意的是，我们在热传导的计算过程中，所采用的grid矩阵是ny*nx的，而非nx*ny的
    for i in range(ny):
        for j in range(nx):
            # 定义每次迭代时的下阈值和上阈值
            lower_y  = -0.5*y*1E10 + i*part_sizey*1E10
            hihger_y = lower_y +part_sizey*1E10
            lower_x  = -0.5*x*1E10 + j*part_sizex*1E10
            hihger_x = lower_x +part_sizex*1E10
            # 根据刚刚所划分的阈值，从Dataframe中定位数据，然后将现在通过传热方程计算出的晶格温度附上
            df.loc[(df['x'] > lower_x) & (df['x'] <= hihger_x
                ) & (df['y'] > lower_y) & (df['y'] <= hihger_y), 'lattice_temp_here'] = Latti_temp[i, j] 
            # 这里的Lattice_temp和前面的grid矩阵一样，都是ny*nx的
    # 根据由data文件提供初始速度值，根据公式计算出此时这个原子所对应的原始温度
    df.loc[:,'original_temp'] = (0.5*26.9815384*0.001/6.02E23*(
                                np.sqrt((df_velocity.loc[:,'vx'])**2 + (df_velocity.loc[:,'vy'])**2 + (df_velocity.loc[:,'vz'])**2)*100)**2)*2/3/1.38E-23
    # 通过将计算出的原始温度和传热方程计算出的温度值进行比较，得到一个缩放系数，用于将现在data文件中的温度变为传热方程中计算出的温度
    df_velocity.loc[:,'scale_factors'] = np.sqrt((df.loc[:,'lattice_temp_here']/df.loc[:,'original_temp']))
    # 将df里面所包含的晶格温度赋给df_velocity中
    df_velocity.loc[:,'lattice_temp_here'] = df.loc[:,'lattice_temp_here']
    # 将vx，vy，vz都按照我们之前所计算出的缩放因子进行缩放，不改变其方向，只改变其大小
    df_velocity.loc[:,'vx'] = df_velocity.loc[:,'vx']*df_velocity.loc[:,'scale_factors']
    df_velocity.loc[:,'vy'] = df_velocity.loc[:,'vy']*df_velocity.loc[:,'scale_factors']
    df_velocity.loc[:,'vz'] = df_velocity.loc[:,'vz']*df_velocity.loc[:,'scale_factors']
    # 将df_velocity中的原始温度，赋给df中
    df_velocity.loc[:,'original_temp'] = df.loc[:,'original_temp']
    # 此时根据新得到的速度，计算出新的温度，便于在计算中观察到新的计算温度是多少
    df_velocity.loc[:,'new_temp'] = (0.5*26.9815384*0.001/6.02E23*(
                                np.sqrt((df_velocity.loc[:,'vx'])**2 + (df_velocity.loc[:,'vy'])**2 + (df_velocity.loc[:,'vz'])**2)*100)**2)*2/3/1.38E-23
    print(df_velocity)
    return df,df_velocity