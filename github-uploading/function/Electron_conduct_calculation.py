import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from showcontourf import showcontourf
#此函数是模拟中所使用的主函数，一切有关热传导的计算都是来自于此函数
def Electron_conduct_calculation(density_file, lattice_temp_file, grid, nx, ny,step,x,y):
    #现在的网格数nx和ny均从外部传入，无需内部再进行定义
    # 初始化网格参数-所有的初始参数均设置为0
    density = 0#密度
    conductivity = 0#电子热导率
    heat_capacity = 0#电子热容
    length = 0#计算单元长度——目前已废弃
    width = 0#计算单元长度——目前已废弃
    lattice_temp = 0#晶格温度
    # 初始化----编制XY方向离散的系数矩阵A和B
    # x方向二阶导系数矩阵A
    A = (-2)*np.eye(nx,k=0) + (1)*np.eye(nx,k=-1) + (1)*np.eye(nx,k=1)
    # y方向二阶导系数矩阵B 
    B = (-2)*np.eye(ny,k=0) + (1)*np.eye(ny,k=-1) + (1)*np.eye(ny,k=1)
    #   系数矩阵的形状如下所示
    #   [-2 1 0 ........... 0 0]
    #   [1 -2 1 ............. 0]
    #   ........................
    #   [0................ -2 1]
    #   [0.................1 -2]
    # 初始化各矩阵参数
    ini = 0
    # 初始化晶格温度矩阵
    T_lattice = ini*np.ones((ny,nx))
    # 初始化电子温度矩阵
    U = ini*np.ones((ny,nx))
    # 初始化密度矩阵
    Den_U = ini*np.ones((ny,nx))
    # 初始化热容矩阵
    Capa_U = ini*np.ones((ny,nx))
    # 初始化热导率矩阵
    k = ini*np.ones((ny,nx))
    # 初始化热扩散系数矩阵
    a = 0*np.ones((ny,nx))
    # a : a = k / (den*Cp)              Coefficient of function.
    #     k   is the Thermal Conductivity, Unit: W/(m°C)
    #     den is the Density, Unit: kg/m3
    #     Cp  is the Specific heat capacity, Unit: J/(kg°C)
    # 初始化常系数矩阵——目前暂时已废弃，等待后续使用
    fq = 0*np.ones((ny,nx))
    # 初始化标志矩阵——主要用于判定是否为空域并关闭热传导计算方程
    flag = 0*np.ones((ny,nx))
    ## 读取密度数据信息，并将信息保存到grid类中
    with open(density_file, 'r') as f:
        # Chunk Coord1 Coord2 Ncount density/number
        # 忽略前三行
        for i in range(3):
            next(f)
        # 读取目前的时间戳，读取目前所有的总grid数量，读取目前整个体系内的所有原子数量
        line = next(f).split()
        time_stamp = float(line[0])
        grid_count = int(line[1])
        atom_count = int(line[2])
        # 创建与grid数量一致的存储坐标和密度信息的矩阵
        coords = np.zeros((grid_count, 3))
        rho_r = np.zeros(grid_count)
        # 导入密度参数的初始化
        rho = 1
        # 遍历文件内的每一行，读取坐标和密度信息，并更新到对应的列表当中
        for i in range(grid_count):
            line = next(f).split()
            coords[i] = [float(line[1]), float(line[2]), float(line[3])]
            rho_r[i] = float(line[4])#密度信息保存
        #更新密度信息
        bi = 0#——bi只是一个用于调试的过程量
        for i in range(nx):
            for j in range(ny):
                    #这里的grid不是nx,ny排布的，是按照ny,nx排布的，但是由于密度信息是按照xy的顺序遍历的，这里进行了一定的调整。
                    grid[j,i].density = rho_r[bi]
                    bi+=1
    ##读取晶格温度文件，读取由lammps计算出的初始晶格温度信息的矩阵
    with open(lattice_temp_file, 'r') as f:
        # 忽略文件的头三行
        for i in range(3):
            next(f)
        # 读取目前的时间戳，读取目前所有的总grid数量
        line = next(f).split()
        time_stamp = float(line[0])
        grid_count = int(line[1])
        # 创建初始的晶格温度矩阵并初始化
        lattice_temp = np.zeros(grid_count)
        # 读取所有的晶格温度信息
        for i in range(grid_count):
            line = next(f).split()
            #更新晶格温度矩阵中的信息
            lattice_temp[i] = float(line[1])
        # bi只是一个过程量，用于调试过程中
        bi = 0
        # 将晶格温度信息赋予到grid类中保存的晶格温度信息当中
        for i in range(nx):
            for j in range(ny):
                    grid[j,i].lattice_temp = lattice_temp[bi]
                    #这里的grid不是nx,ny排布的，是按照ny,nx排布的，但是由于密度信息是按照xy的顺序遍历的，这里进行了一定的调整。
                    bi+=1
    # 将温度-密度-热容-晶格温度参数赋予网格grid矩阵，即初始的，或者上一步存留的信息进行移交
    for i in range(ny):
        for j in range(nx):       
            U[i,j] = grid[i,j].temperature
            Den_U[i,j] = grid[i,j].density
            Capa_U[i,j] = grid[i,j].heat_capacity
            T_lattice[i,j] = grid[i,j].lattice_temp
    # 原子密度和电子密度是相关的，这里通过原子密度来构建电子密度信息
    Den_U = Den_U*1.806E29*1E-30# 1/m^3
    ##--------------------------------------------------------------------------------
    ## 初始化模拟参数
    for i in range(ny):
        for j in range(nx):   
            k[i,j] = 100#*diffusivity*Capa_U[i,j]##这里在我之前的修改中改变了前面的100，的大小
            if Den_U[i,j] == 0:#如果判断到矩阵此处的密度乘以热容，变成了0，即此处没有原子，直接将此处置零，以免运算出错
                flag[i,j] = 0#置零方法为，提供一个flag变量，将此处的判定条件置为0
                a[i,j] = 0
                fq[i,j] = 0
            else:
                flag[i,j] = 1#判断是否是空域，空域则不需要计算
                #fq[i,j] = -10/(Den_U[i,j]*Capa_U[i,j])
                ##--------------------------------------------------------------------------
                a[i,j] = 0.777*Den_U[i,j]#k[i,j]/(Den_U[i,j]*Capa_U[i,j])##以后需要进行修改
                #----------------------------------------------------------------------------
    # 现在x和y已经能够从主程序main-5-11.py中获取
    # 模拟体系x方向单位长度
    #x  = 2*224.33166413475473E-10 # A
    dx = x/nx
    # 模拟体系y方向单位长度
    #y  = 2*112.16583206737737E-10 #A
    dy = y/ny
    # 模拟总时长
    T  = 1E-15 #fs
    # 模拟步数step数目
    Nt = 10000
    # 模拟的单位时间步长
    dt = T/Nt
    # 模拟中传热方程的x方向扩散项，y方向扩散项，ft常系数项（热源项）
    rx,ry,ft = a*dt/dx**2, a*dt/dy**2, fq*dt
    # 初始化模拟计时器
    start = time.time()
    # 模拟出图的绘图帧数项
    Frame = 100
    #绘图用空矩阵
    D = np.array([0,nx,0,ny])
    ##--------------------------------------------------------------------------------
    # 模拟使用的手动输入参数项
    # 电子-晶格耦合参数项
    Grt = 0.5E13##-------------------------------此处之后需要进行修改
    ##--------------------------------------------------------------------------------
    ## 模拟主程序
    # 使用欧拉法进行时间离散-按步长进行迭代
    # max electron -- max lattice，最大电子温度和最大晶格温度的信息保存
    lattice_temp = pd.DataFrame(columns=['time','temp'])
    electron_temp = pd.DataFrame(columns=['time','temp'])
    index = 0
    # 进行时间离散迭代
    for p in range(Nt+1):
        # 真实时间标志
        tt = p*dt
        # 求解矩阵内部节点
        # 划分矩阵节点
        x = np.linspace(D[0],D[1],U.shape[1])
        y = np.linspace(D[2],D[3],U.shape[0])
        X,Y = np.meshgrid(x,y)
        # 热源的峰值位置
        peak = 1000
        # 热源的时间和空间分量
        time_component = -np.power(p - peak, 2)/ (2 * np.power(0.1, 2))
        xy_component   = -np.power(X - 25, 2)/ (2 * np.power(0.1, 2))-np.power(Y - 25, 2)/ (2 * np.power(0.1, 2))
        #heat = 1.7E32*np.exp(xy_component)*np.exp(time_component)
        # 热源作用时间
        if p <= 1000:
            heat = 0#1.7E32*np.exp(xy_component)*np.exp(time_component)
        else:
            heat = 0
        # 温度矩阵1 = 温度矩阵0 + 温度矩阵x方向传热 + 温度矩阵y方向传热 +电子-晶格温度耦合项
        U = U + rx*np.dot(U,A) + ry*np.dot(B,U) + heat * dt - (U-T_lattice)*Grt*dt*flag#/(200/a)
        # 晶格温度1 = 晶格温度0 + 电子-晶格温度耦合项
        T_lattice = T_lattice + (U-T_lattice)*Grt*dt*flag#/(200/a)
        # 设置边界条件
        # 边界条件1-恒温
        U[:,0]  = U[:,-1] = 300
        U[0,:]  = U[-1,:] = 300
        # 边界条件2-与前一网格点温度保持一致
        #U[:,0]  = U[:,1]
        #U[:,-1] = U[:,-2]
        #U[0,:]  = U[1,:]
        #U[-1,:] = U[-2,:]
        # 边界条件3-向外界的恒定热流
        # ch = 0.1 
        # U[:,0] = (u_env + ch*U[:,1]/dx)/(1+ch/dx)
        # U[:,-1] = (u_env + ch*U[:,-2]/dx)/(1+ch/dx) 
        # U[0,:] = (u_env + ch*U[1,:]/dy)/(1+ch/dy)
        # U[-1,:] = (u_env + ch*U[-2,:]/dy)/(1+ch/dy) 
        end = time.time()
        # 每一步打印模拟的结果，输出图像
        
        index+=1
        print('T = {:.19f} s   max_U= {:.19f} min_U = {:.19f} max_lattice = {:.1f} min_lattice = {:.1f}'.format(tt,np.max(U),np.min(U),np.max(T_lattice),np.min(T_lattice)),end = '\r')
        if p%100000 == 0:
            #用于输出电子温度和晶格温度图像
            lattice_temp.loc[index]     = [tt,np.max(T_lattice)]
            electron_temp.loc[index]    = [tt,np.max(U)]
            showcontourf(U,D,vmax=2000,timestep1=step, name = str('electron'))
            showcontourf(T_lattice,D,vmax=2000,timestep1=step,name = str('lattice'))
    lattice_temp.to_csv('lattice_temp.txt')
    electron_temp.to_csv('electron_temp.txt')
    ## 保存计算完成后的电子温度信息-用于下一次运行时读入
    for i in range(ny):
        for j in range(nx):
            grid[i, j].temperature = U[i,j]  
            grid[i, j].lattice_temp = T_lattice[i,j]  

    return grid, T_lattice