import cupy as cp
import numpy as numpy
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from showcontourf import showcontourf

def Electron_conduct_calculation(U_in, savepath, density_file, lattice_temp_file, grid, nx, ny,step,x,y,z):
    def write_temperature_to_grid(U,file_path,step,x,y,z):
        ny,nx =U.shape
        write_x_low = -x*1E10/2
        write_x_hi = x*1E10/2
        write_y_low = -y*1E10/2
        write_y_hi = y*1E10/2
        write_z_low = -z*1E10/2
        write_z_hi = z*1E10/2
        with open(file_path, 'w') as file:
            file.write('ITEM: TIMESTEP\n')
            file.write(f'{step}\n')
            file.write('ITEM: BOX BOUNDS pp pp pp\n')
            file.write('{} {}\n{} {}\n{} {}\n'.format(write_x_low,write_x_hi,write_y_low,write_y_hi,write_z_low,write_z_hi))
            file.write("ITEM: DIMENSION\n")
            file.write("3\n")
            file.write("ITEM: GRID SIZE nx ny nz\n")
            file.write(f"{nx} {ny} 1\n")
            file.write("ITEM: GRID CELLS c_1:grid:data[1] c_1:grid:data[2] c_1:grid:data[3] c_1:grid:data[4] f_ttm:grid:Electron_temperature\n")
            # 写入温度数据
            grid_id =1
            for y in range(ny):
                for x in range(nx):
                    temperature = U[y, x]
                    file.write(f"{grid_id} {x+1} {y+1} 1 {temperature}\n")
                    grid_id += 1
        print(f"温度数据已成功写入文件：{file_path}")
    density = 0
    conductivity = 0
    heat_capacity = 314# J/m^3.K
    length = 0
    width = 0
    lattice_temp = 0
    ##----------------------------------------------------------------------------
    # 初始化----编制XY方向离散的系数矩阵A和B
    # x方向二阶导系数矩阵A
    A = (-2)*cp.eye(nx,k=0) + (1)*cp.eye(nx,k=-1) + (1)*cp.eye(nx,k=1)
    # y方向二阶导系数矩阵B 
    B = (-2)*cp.eye(ny,k=0) + (1)*cp.eye(ny,k=-1) + (1)*cp.eye(ny,k=1)
    #   [-2 1 0 ........... 0 0]
    #   [1 -2 1 ............ 0]
    # 初始化各矩阵参数
    ini = 0
    # 初始化晶格温度 矩阵
    T_lattice = ini*numpy.ones((ny,nx))
    # 初始化电子温度矩阵
    U = ini*numpy.ones((ny,nx))
    # 初始化密度矩阵
    Den_U = ini*numpy.ones((ny,nx))
    # 初始化热容矩阵
    Capa_U = ini*cp.ones((ny,nx))
    couple = ini*cp.ones((ny,nx))
    couple_lattice = ini*cp.ones((ny,nx))
    # 初始化热导率矩阵
    k = ini*numpy.ones((ny,nx))
    # 初始化热扩散系数矩阵
    a = 0*numpy.ones((ny,nx))
    # a : a = k / (den*Cp)              Coefficient of function.
    #     k   is the Thermal Conductivity, Unit: W/(m°C)
    #     den is the Density, Unit: kg/m3
    #     Cp  is the Specific heat capacity, Unit: J/(kg°C)
    # 初始化常系数矩阵
    fq = 0*cp.ones((ny,nx))
    # 初始化标志矩阵-用于处理
    flag = 1*numpy.ones((ny,nx))
    ##----------------------------------------------------------------------------
    duqu_time0 = time.time()
    ## 读取密度数据信息
    with open(density_file, 'r') as f:
        # Chunk Coord1 Coord2 Ncount density/number
        # Ignore the first three lines
        for i in range(3):
            next(f)
        # Read in the time stamp, grid count, and atom count
        line = next(f).split()
        # Update the material density
        coord_density = numpy.loadtxt(f, usecols = (1,2,3,4))
    ## 读取晶格温度信息
    with open(lattice_temp_file, 'r') as f:
        # Row Lattice-temperature
        # Ignore the first three lines
        for i in range(3):
            next(f)
        # Read in the time stamp, grid count, and atom count
        line = next(f).split()
        coord_l_temp = numpy.loadtxt(f, usecols=(1,2,3,4))
        coord_density = coord_density[:,3]
        coord_l_temp = coord_l_temp[:,3]
    # 将数据重塑为[400, 200]的矩阵
    mapped_density = coord_density.reshape(nx, ny).T
    mapped_l_temp = coord_l_temp.reshape(nx, ny).T
    U = U_in
    Den_U = mapped_density*2.99671213E30/1.806E29# 1/m^3
    Capa_U = 100*U# J/m^3 K
    T_lattice = mapped_l_temp
    zero_den_indices = Den_U == 0
    # 对于Den_U为零的位置，将相应的值置零
    k[zero_den_indices] = 0
    flag[zero_den_indices] = 0
    a[zero_den_indices] = 0
    fq[zero_den_indices] = 0
    couple[zero_den_indices] = 0
    couple_lattice[zero_den_indices] = 0
    U[zero_den_indices] = 300
    # 计算非零Den_U的位置的其他值
    non_zero_den_indices = Den_U > 0
    k[non_zero_den_indices] = 236 * (U[non_zero_den_indices] / T_lattice[non_zero_den_indices])
    couple[non_zero_den_indices] = 1 / Capa_U[non_zero_den_indices]
    couple_lattice[non_zero_den_indices] = 1 / 2.376E6
    a = Den_U*k/Capa_U
    #fq[non_zero_den_indices] = -10 / (Den_U[non_zero_den_indices] * Capa_U[non_zero_den_indices])
    # 如果fq的计算需要非零Den_U的值，取消注释上面的行并根据需要进行计算
    # 对于Den_U小于0.6的位置，计算a的值
    a[Den_U < 0.6] = 0.01 * Den_U[Den_U < 0.6] * k[Den_U < 0.6] / Capa_U[Den_U < 0.6]
    Original_Latti_temp = cp.asnumpy(T_lattice)
    T_lattice = cp.asarray(T_lattice)
    U = cp.asarray(U)
    # 初始化密度矩阵
    Den_U = cp.asarray(Den_U)
    k = cp.asarray(k)
    a = cp.asarray(a)
    fq = 0*cp.ones((ny,nx))
    flag = cp.asarray(flag)

    duqu_time1 = time.time()
    print("INITIATE COST = {}".format(duqu_time1-duqu_time0))
    # 模拟体系x方向单位长度
    dx = x/nx
    # 模拟体系y方向单位长度
    dy = y/ny
    # 模拟总时长
    T  = 1E-15 #fs
    # 模拟步数step数目
    Nt = 500
    # 模拟的单位时间步长
    dt = T/Nt 
    # 模拟中传热方程的x方向扩散项，y方向扩散项，ft常系数项（热源项）
    rx,ry,ft = a*dt/dx**2, a*dt/dy**2, fq*dt
    # 初始化模拟计时器
    start = time.time()
    # 模拟出图的绘图帧数项
    Frame = 1000
    # 电子-晶格耦合参数项
    Grt = 2.45E17
    ##--------------------------------------------------------------------------------
    ## 模拟主程序
    # 使用欧拉法进行时间离散-按步长进行迭代
    lattice_temp = pd.DataFrame(columns=['time','temp'])
    electron_temp = pd.DataFrame(columns=['time','temp'])
    index = 0
    print("Preparation Complete, Now Calculate starts")
    for p in range(Nt+1):
        t0 = time.time()
        # 真实时间
        tt = p*dt
        # 扩展温度矩阵和其他的相关参数矩阵
        U_big  =  cp.concatenate((U, U, U),axis=1)# manipulate in x 
        T_lattice_big = cp.concatenate((T_lattice, T_lattice, T_lattice),axis = 1)
        Den_big = cp.concatenate((Den_U, Den_U, Den_U), axis= 1)
        # x方向二阶导系数矩阵A
        A_big = (-2)*cp.eye(3*nx,k=0) + (1)*cp.eye(3*nx,k=-1) + (1)*cp.eye(3*nx,k=1)
        # y方向二阶导系数矩阵B 
        B_big = (-2)*cp.eye(ny,k=0) + (1)*cp.eye(ny,k=-1) + (1)*cp.eye(ny,k=1)
        #定义热源矩阵
        heat = cp.zeros((ny,nx))
        #S(y,t)= 2/√（π/ln2) * J_abs/ (L_op*t_p) * exp{-y/L_op} * exp{-4ln2* (t-t_p)^2/t_p^2} 
        #激光的能量在光斑范围内平均分布
        #J_abs 激光能量密度     L_op 光学穿透深度     t_p 激光脉冲宽度      y为激光传播方向的距离      
        J_abs=30    #文献中是30J/m^2     
        L_op=8E-9   #   m
        t_p=1E-13   #单位是s
        t=1E-13 
        if step > 15000:
            heat = heat
        else:
            for j in range (int(1/4*ny)):
                heat[j,:] = 0.93944*J_abs/(L_op*t_p) * cp.exp(-(j*y/ny)/L_op) * cp.exp(-2.77259*((t-t_p)/t_p)**2)
        heat_big = cp.concatenate((heat, heat,  heat),axis = 1)
        flag_big = cp.concatenate((flag, flag, flag),axis =1)
        couple_big = cp.concatenate((couple, couple, couple),axis = 1)
        couple_lattice_big = cp.concatenate((couple_lattice, couple_lattice, couple_lattice),axis =1)
        rx_big = cp.concatenate((rx, rx, rx),axis =  1)
        ry_big = cp.concatenate((ry, ry, ry),axis = 1)
        # 强制变量转换
        U_big = U_big.astype(float)
        rx_big = rx_big.astype(float)
        ry_big = ry_big.astype(float)
        A_big = A_big.astype(float)
        B_big = B_big.astype(float)
        heat_big =heat_big.astype(float)
        T_lattice_big = T_lattice_big.astype(float)
        U_big = U_big +rx_big*cp.dot(U_big,A_big)+ ry_big*cp.dot(B_big,U_big) + heat_big * dt*couple_big - (U_big-T_lattice_big)*Grt*dt*flag_big*couple_big
        T_lattice_big  = T_lattice_big  + (U_big -T_lattice_big )*Grt*dt*flag_big*couple_lattice_big 
        #_______________________________________边界识别算法开始____________________________________________
        condition_points = []
        rows, cols = flag_big.shape
        # 对flag_big矩阵进行上下左右偏移并累加
        flag_big = flag_big.get().astype(int)
        #flag_big = flag_big.astype(int)
        shifted_flags = numpy.zeros((rows, cols), dtype=int)
        # 定义上下左右四个方向的偏移
        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        # 遍历原始矩阵的每个元素
        for shift in directions:
            shifted_flags += numpy.roll(flag_big, shift, axis=(0, 1))
        # 修改温度值
        modified_U_big = numpy.copy(cp.asnumpy(U_big))
        # modified_U_big = U_big
        # 找出满足条件的点的坐标
        condition = (shifted_flags > 0) & (flag_big == 0)

        # 获取满足条件的点的索引
        condition_points = numpy.argwhere(condition)
        # 添加有效性检查，排除超出边界的点
        rows, cols = flag_big.shape
        valid_points = (0 < condition_points[:, 0]) & (condition_points[:, 0] < rows-1) & (0 < condition_points[:, 1]) & (condition_points[:, 1] < cols-1)
        condition_points = condition_points[valid_points]
        i, j = condition_points.T
        neighbors = [
            modified_U_big[i-1, j],
            modified_U_big[i+1, j],
            modified_U_big[i, j-1],
            modified_U_big[i, j+1],
        ] 
        max_neighbor_values = numpy.maximum.reduce(neighbors)
        modified_U_big[i,j] = max_neighbor_values
        U_big = cp.asarray(modified_U_big)
        #_____________________________________边界识别算法结束____________________________________________________________________
        ########## 绝热边界条件
        nambda = 1E20
        h = 1E-40
        ch = nambda/h
        U_big[:,0]  = U_big[:,1]#(300)#+ch*U_big[:,1]/dx)/(1+ch/dx)
        U_big[:,-1] = U_big[:,-2]#(300)#+ch*U_big[:,-2]/dx)/(1+ch/dx)
        U_big[0,:]  = U_big[1,:]#(300)#+ch*U_big[1,:]/dy)/(1+ch/dy)
        U_big[-1,:] = U_big[-2,:]#(300)#+ch*U_big[-2,:]/dy)/(1+ch/dy)
        size_x = U.shape[0]
        size_y = U.shape[1]
        size_x_big = U_big.shape[0]
        size_y_big = U_big.shape[1]
        # 把大矩阵重新切回以前的样子
        U = U_big[(size_x_big-size_x)//2:(size_x_big-size_x)//2+size_x,(size_y_big-size_y)//2:(size_y_big-size_y)//2+size_y]
        T_lattice = T_lattice_big[(size_x_big-size_x)//2:(size_x_big-size_x)//2+size_x,(size_y_big-size_y)//2:(size_y_big-size_y)//2+size_y]
        # 每一步打印模拟的结果，输出图像
        t1 = time.time()
        index+=1
        if p%10 ==0:
            print("This PerCycle has {} points need to be check. Timecost = {}".format(len(condition_points),t1-t0),end='\r')
            print('T = {:.2f} fs   max_U= {:.5f} max_lattice = {:.1f} Heat_contribution = {:.5f}'.format(tt*1E15,cp.max(U),cp.max(T_lattice),cp.max(heat)*dt*cp.max(couple)),'Step ={}'.format(step),end = '\r')
        if step%1 == 0 and p%100000000 == 0:
            write_temperature_to_grid(U,savepath +'electron_t_{}.grid.dump'.format(step),step,x,y,z)
    with open('lattice_temp.dat', 'a') as f:
        f.write("{},{}\n".format(tt,cp.max(T_lattice)))
    with open('electron_temp.dat', 'a') as f:
        f.write("{},{}\n".format(tt,cp.max(U)))
    ## 保存计算完成后的电子温度信息-用于下一次运行时读入
    U_out = cp.asnumpy(U)
    T_lattice = cp.asnumpy(T_lattice)
    return grid, T_lattice, Original_Latti_temp, U_out