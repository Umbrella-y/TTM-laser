#import cupy as cp
import numpy as numpy
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math

def reading_data(file_path):
    with open(file_path, "r") as f:
        # Skip the first two lines
        next(f)
        next(f)
        # Read the third line and extract the total number of atoms
        total_atoms = int(f.readline().strip().split()[0])
        # Skip the lines until the "atom" keyword is encountered
        for line in f:
            if "Atoms" in line:
                break
        next(f)
        # Read the next total_atoms lines
        data = numpy.loadtxt(f, max_rows= total_atoms)
        df = pd.DataFrame(data, columns=['id', 'type', 'x', 'y', 'z','image1', 'image2', 'image3'])
        next(f)
        # Skip the lines until the "atom" keyword is encountered
        for line in f:
            if "Velocities" in line:
                break
        next(f)
        new_data = numpy.loadtxt(f, max_rows= total_atoms)
        df_velocity = pd.DataFrame(new_data, columns=['id', 'vx', 'vy', 'vz'])
        df_velocity['id'] = df_velocity['id'].astype(str)
        df_velocity[['vx', 'vy', 'vz']] = df_velocity[['vx', 'vy', 'vz']].astype(float)
    return df, df_velocity