import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math

def data_output(df1,df_velocity,x,y,z,filename):
    write_x_low = -x*1E10/2
    write_x_hi = x*1E10/2
    write_y_low = -y*1E10/2
    write_y_hi = y*1E10/2
    write_z_low = -z*1E10/2
    write_z_hi = z*1E10/2
    output_df1 = df1[['id', 'type', 'x', 'y', 'z', 'image1', 'image2', 'image3']].astype(float)
    output_df2 = df_velocity[['id', 'vx', 'vy', 'vz']].astype(float)
    # Write the columns to a new text file with the specified format
    with open('test.data', 'w') as f:
        # Write the head
        f.write('LAMMPS data file via Python\n\n')
        f.write('{} atoms\n1 atom types\n\n'.format(len(df1)))
        f.write('{} {} xlo xhi\n{} {} ylo yhi\n{} {} zlo zhi\n\n'.format(write_x_low,write_x_hi,write_y_low,write_y_hi,write_z_low,write_z_hi))
        f.write('Masses\n\n1 26.9815384\n\n')
        # Write the 'Atoms # atomic' line
        f.write('Atoms # atomic\n\n')
        # Write the 'id', 'type', 'x', 'y', 'z', 'image1', 'image2', and 'image3' columns
        print('Now writing NEW DATA file')
        np.savetxt(f, output_df1, delimiter=' ',fmt= '%d %d %.5f %.5f %.5f %d %d %d')
        #for index, row in output_df1.iterrows():
           # if int(index)% 1000 == 0:
                # print("Now writing corrdiantes: line({})".format(int(index)),end = '\r')
            #f.write(f"{row['id']} {row['type']} {row['x']} {row['y']} {row['z']} {row['image1']} {row['image2']} {row['image3']}\n")
            
        # Write an empty line
        f.write('\n')
        # Write the 'Velocities' line
        f.write('Velocities\n\n')
        # Select the 'id', 'vx', 'vy', and 'vz' columns from the dataframe
        # Write the 'id', 'vx', 'vy', and 'vz' columns
        np.savetxt(f, output_df2,delimiter=' ', fmt= '%d %.5f %.5f %.5f')
        # for index, row in output_df2.iterrows():
        #     if int(index)% 1000 == 0:
        #         print("Now writing velocities: line({})".format(int(index)),end = '\r')
        #     f.write(f"{row['id']} {row['vx']} {row['vy']} {row['vz']}\n")

    df1.to_csv('total_Data_{}.csv'.format(str(filename)))