import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math

def calculate_velocity(df,df_velocity,x,nx,y,ny,Latti_temp,Original_Latti_temp):


    # Calculate the size of each part
    part_sizex = x / nx
    part_sizey = y / ny
    df['x'] = pd.to_numeric(df['x'])
    df['y'] = pd.to_numeric(df['y'])
    # Use the cut() function to divide the dataframe into grids
    df['x_grid'] = pd.cut(df['x'], nx)
    df['y_grid'] = pd.cut(df['y'], ny)
    # Get the labels used by the cut() function
    x_labels = pd.cut(df['x'], nx).unique()
    y_labels = pd.cut(df['y'], ny).unique()
    #print(x_labels)
    # Group the dataframe by the x and y grids
    #Latti_temp = np.ones((nx,ny))




    for i in range(ny):
        for j in range(nx):

            lower_y  = -y/2*1E10 + i*part_sizey*1E10
            hihger_y = lower_y +part_sizey*1E10
            lower_x  = -x/2*1E10 + j*part_sizex*1E10
            hihger_x = lower_x +part_sizex*1E10
            if i % 10 == 0 and j % 10 == 0:
                print("Now Calculating New Velocities at line({}, {})".format(i, j), end='\r')
            # 在循环之前计算索引条件
            condition = (df['x'] > lower_x) & (df['x'] <= hihger_x) & (df['y'] > lower_y) & (df['y'] <= hihger_y)
            # 使用预先计算的条件来更新DataFrame
            df.loc[condition, 'lattice_temp_here'] = Latti_temp[i, j]
            df.loc[condition, 'Original_Latti_temp'] = Original_Latti_temp[i, j]





    #计算方法2
    #根据计算的标量值，计算此时的动能大小并计算温度
    #df.loc[:,'original_temp'] = (0.5*26.9815384*0.001/6.02214076E23*(np.sqrt((df_velocity.loc[:,'vx'])**2 + (df_velocity.loc[:,'vy'])**2 + (df_velocity.loc[:,'vz'])**2)*100)**2)*2/3/1.380649E-23
    df.loc[:,'original_temp'] = df.loc[:,'Original_Latti_temp']
    df_velocity.loc[:,'scale_factors'] = np.sqrt((df.loc[:,'lattice_temp_here']/df.loc[:,'original_temp']))
    df_velocity.loc[:,'lattice_temp_here'] = df.loc[:,'lattice_temp_here']
    #print(df_velocity)
    df_velocity.loc[:,'vx'] = df_velocity.loc[:,'vx']*df_velocity.loc[:,'scale_factors']
    df_velocity.loc[:,'vy'] = df_velocity.loc[:,'vy']*df_velocity.loc[:,'scale_factors']
    df_velocity.loc[:,'vz'] = df_velocity.loc[:,'vz']*df_velocity.loc[:,'scale_factors']
    df_velocity.loc[:,'original_temp'] = df.loc[:,'original_temp']
    df_velocity.loc[:,'new_temp'] = (0.5*26.9815384*0.001/6.02214076E23*(np.sqrt((df_velocity.loc[:,'vx'])**2 + (df_velocity.loc[:,'vy'])**2 + (df_velocity.loc[:,'vz'])**2)*100)**2)*2/3/1.380649E-23
    #df = pd.concat([df,df_velocity],axis=0,ignore_index=True)
    #print(df_velocity)
    return df,df_velocity