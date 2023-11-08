import numpy as numpy 
import os
import pandas as pd

def nx_and_ny_settings(workpath):
    lmp_file = workpath + '/' + 'in.read_data.lmp'
    with open(lmp_file, 'r') as f:
        target_line = None
        for line in f:
            if "compute binchunk all chunk/atom bin/2d x lower" in line:
                target_line = line
                break
    # 如果找到目标行，提取x lower和y lower后面的值并计算1除以该值
    if target_line is not None:
        values = target_line.split()
        x_lower_index = values.index('x') + 2
        y_lower_index = values.index('y') + 2
        x_lower_value = float(values[x_lower_index])
        y_lower_value = float(values[y_lower_index])

        # 计算1除以x lower和y lower的值
        x_lower_inverse = 1.0 / x_lower_value
        y_lower_inverse = 1.0 / y_lower_value

        # 将结果转换为整数
        x_lower_inverse_int = int(x_lower_inverse)
        y_lower_inverse_int = int(y_lower_inverse)
        nx, ny = x_lower_inverse_int, y_lower_inverse_int# 模拟网格数目的设置，现在直接通过lmp文件获得，无需再进行手动设置
    else:
        print("目标行未找到")    
    return nx, ny

def x_y_z_size_settings(workpath):
    datafile = workpath + '/' + 'test.data'
    #读取模拟体系的长宽高参数
    with open(datafile, "r") as f:
        for _ in range(5):
            next(f)
        # 读取 xlo, ylo, zlo 的值
        xlo, xhi = map(float, next(f).split()[:2])
        ylo, yhi = map(float, next(f).split()[:2])
        zlo, zhi = map(float, next(f).split()[:2])
    #必要时可以手动设置，平时默认读取
    return xhi, yhi, zhi