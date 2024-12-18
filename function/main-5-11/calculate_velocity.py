import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math

def calculate_velocity(df,df_velocity,x,nx,y,ny,Latti_temp):


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
    for ix in range(nx):
        for iy in range(ny):
            lower_x  = -202 + ix*part_sizex*1E10
            hihger_x = lower_x +part_sizex*1E10
            lower_y  = -202 + iy*part_sizex*1E10
            hihger_y = lower_y +part_sizey*1E10
            #print(lower_x,hihger_x)
            df.loc[(df['x'] > lower_x) & (df['x'] <= hihger_x
                ) & (df['y'] > lower_y) & (df['y'] <= hihger_y), 'lattice_temp_here'] = Latti_temp[ix, iy]  
    #print(df)                      

    # Group the dataframe by the unique values in the 'lattice_temp_here' column
    #grouped = df.groupby('lattice_temp_here')
    # Iterate over the groups and generate the corresponding values
    #df_list = []
    #for name, group in grouped:
        # Generate the 'vx' and 'vy' values using a Gaussian distribution
        #print(group['lattice_temp_here'])
        #mean = math.sqrt(np.mean(group['lattice_temp_here'])*1.38E-23*3/2*2/(26.9815384E-3/6.02E23))*0.01#*1.23
        #std = 0
        #num_values = len(group)
        #vx_values = np.random.normal(0, 1, num_values)
        #vy_values = np.random.normal(0, 1, num_values)
        #vz_values = np.random.normal(0, 1, num_values)
        #v_mag_values = np.random.normal(mean, std, num_values)
        #___Scale the 'vx' and 'vy' components such that their sum equals the total velocity magnitude
        #scale_factors = v_mag_values / np.sqrt(vx_values**2 + vy_values**2 + vz_values**2)
        #vx_values *= scale_factors
        #vy_values *= scale_factors
        #_____Set the 'vz' values to zero
        #vz_values *= scale_factors
        #____Add the new columns to the group dataframe
        #group['vx'] = vx_values
        #group['vy'] = vy_values
        #group['vz'] = vz_values
        #df_list.append(group)
    #计算方法2
    #根据计算的标量值，计算此时的动能大小并计算温度
    df.loc[:,'original_temp'] = (0.5*26.9815384*0.001/6.02214076E23*(np.sqrt((df_velocity.loc[:,'vx'])**2 + (df_velocity.loc[:,'vy'])**2 + (df_velocity.loc[:,'vz'])**2)*100)**2)*2/3/1.380649E-23
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