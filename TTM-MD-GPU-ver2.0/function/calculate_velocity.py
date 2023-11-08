import cupy as cp
import numpy as numpy  
import pandas as pd

def calculate_velocity(df, df_velocity, x, nx, y, ny, Latti_temp, Original_Latti_temp):
    # 将单位从Angstrom转换为m
    x = x * 1E-10
    y = y * 1E-10

    # 将数据框df中的'x'和'y'列转换为数值类型
    df['x'] = pd.to_numeric(df['x'])
    df['y'] = pd.to_numeric(df['y'])

    # 为x和y创建均匀间隔的网格
    x_bins = numpy.linspace(-x/2, x/2, nx+1)
    y_bins = numpy.linspace(-y/2, y/2, ny+1)

    # 将数据框df中的'x'和'y'值分别映射到对应的x_grid和y_grid
    df['x_grid'] = numpy.digitize(df['x'], x_bins)
    df['y_grid'] = numpy.digitize(df['y'], y_bins)

    # 从Latti_temp和Original_Latti_temp中提取温度值
    df['lattice_temp_here'] = Latti_temp[df['y_grid']-1, df['x_grid']-1]
    df['Original_Latti_temp'] = Original_Latti_temp[df['y_grid']-1, df['x_grid']-1]

    # 创建'original_temp'列
    df['original_temp'] = df['Original_Latti_temp']
    df_velocity['lattice_temp_here'] = df['lattice_temp_here']

    # 计算温度缩放因子并应用于速度分量
    scale_factors = numpy.sqrt(df['lattice_temp_here'] / df['original_temp'])
    scale_factors[numpy.isnan(scale_factors)] = 1
    df_velocity['vx'] *= scale_factors
    df_velocity['vy'] *= scale_factors
    df_velocity['vz'] *= scale_factors

    # 计算新的温度值并将其添加到数据框
    df_velocity['original_temp'] = df['original_temp']
    df_velocity['new_temp'] = 0.5 * 26.9815384 * 0.001 / (6.02214076E23 * (numpy.sqrt(df_velocity['vx']**2 + df_velocity['vy']**2 + df_velocity['vz']**2) * 100)**2) * 2 / 3 / 1.380649E-23

    return df, df_velocity
