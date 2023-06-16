import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
# 该函数用于得到前面所有的计算结果，将计算结果写为lammps的data文件
def data_output(df1, df_velocity, x, y, filename):
    output_df1 = df1[['id', 'type', 'x', 'y', 'z', 'image1', 'image2', 'image3']]
    # 清空test.data中存留的数据并写入我们新的数据
    with open('test.data', 'w') as f:
        # 编写文件头
        f.write('LAMMPS data file via Python\n\n')
        # 编写原子类型和原子数量
        f.write('{} atoms\n1 atom types\n\n'.format(len(df1)))
        # 编写模拟盒子的x, y, z值
        f.write('{} {} xlo xhi\n{} {} ylo yhi\n-2.039378764861385 2.039378764861385 zlo zhi\n\n'.format(-0.5*x*1E10,0.5*x*1E10,-0.5*y*1E10,0.5*y*1E10))
        # 编写原子的相对原子质量
        f.write('Masses\n\n1 26.9815384\n\n')
        # 编写 'Atoms # atomic' 数据行
        f.write('Atoms # atomic\n\n')
        # 迭代前面计算出的df中的所有的数据，并进行写入 'id', 'type', 'x', 'y', 'z', 'image1', 'image2', 'image3' ，主要是原子的坐标信息
        for index, row in output_df1.iterrows():
            f.write(f"{row['id']} {row['type']} {row['x']} {row['y']} {row['z']} {row['image1']} {row['image2']} {row['image3']}\n")
        # 编写一个空行
        f.write('\n')
        # 开始编写'Velocities' 数据
        f.write('Velocities\n\n')
        # 从df_velocity中选出我们所需要的数据
        output_df2 = df_velocity[['id', 'vx', 'vy', 'vz']]
        # 将这些数据写入到data文件当中
        for index, row in output_df2.iterrows():
            f.write(f"{row['id']} {row['vx']} {row['vy']} {row['vz']}\n")
    #输出所有的数据作为存根，每一帧会刷新该存根
    df1.to_csv('total_Data_{}.csv'.format(str(filename)))