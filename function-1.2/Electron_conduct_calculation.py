import numpy as np
from grid import Grid
import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from showcontourf import showcontourf

def Electron_conduct_calculation(density_file, lattice_temp_file, grid, nx, ny,step,x,y,z):
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
    #x  = 2*50E-10 # A
    #y  = 2*50E-10 #A
    # 定义网格数
    #nx, ny = 100, 100
    # 初始化网格参数-此处的参数均可随机设置
    density = 0
    conductivity = 0
    heat_capacity = 314# J/m^3.K
    length = 0
    width = 0
    lattice_temp = 0
    ##----------------------------------------------------------------------------
    # 初始化----编制XY方向离散的系数矩阵A和B
    # x方向二阶导系数矩阵A
    A = (-2)*np.eye(nx,k=0) + (1)*np.eye(nx,k=-1) + (1)*np.eye(nx,k=1)
    # y方向二阶导系数矩阵B 
    B = (-2)*np.eye(ny,k=0) + (1)*np.eye(ny,k=-1) + (1)*np.eye(ny,k=1)
    #   [-2 1 0 ........... 0 0]
    #   [1 -2 1 ............ 0]
    # 初始化各矩阵参数
    ini = 0
    # 初始化晶格温度 矩阵
    T_lattice = ini*np.ones((ny,nx))
    # 初始化电子温度矩阵
    U = ini*np.ones((ny,nx))
    # 初始化密度矩阵
    Den_U = ini*np.ones((ny,nx))
    # 初始化热容矩阵
    Capa_U = ini*np.ones((ny,nx))
    couple = ini*np.ones((ny,nx))
    couple_lattice = ini*np.ones((ny,nx))
    # 初始化热导率矩阵
    k = ini*np.ones((ny,nx))
    # 初始化热扩散系数矩阵
    a = 0*np.ones((ny,nx))
    # a : a = k / (den*Cp)              Coefficient of function.
    #     k   is the Thermal Conductivity, Unit: W/(m°C)
    #     den is the Density, Unit: kg/m3
    #     Cp  is the Specific heat capacity, Unit: J/(kg°C)
    # 初始化常系数矩阵
    fq = 0*np.ones((ny,nx))
    # 初始化标志矩阵-用于处理
    flag = 0*np.ones((ny,nx))
    ##----------------------------------------------------------------------------

    ## 读取密度数据信息
    with open(density_file, 'r') as f:
        # Chunk Coord1 Coord2 Ncount density/number
        # Ignore the first three lines
        for i in range(3):
            next(f)
        # Read in the time stamp, grid count, and atom count
        line = next(f).split()
        time_stamp = float(line[0])
        grid_count = int(line[1])
        atom_count = int(line[2])
        # Read in the grid information
        coords = np.zeros((grid_count, 3))
        rho_r = np.zeros(grid_count)
        rho = 1
        # Update the material density
        for i in range(grid_count):
            line = next(f).split()
            coords[i] = [float(line[1]), float(line[2]), float(line[3])]
            rho_r[i] = float(line[4])#密度信息保存
        #print(rho_r)
        bi = 0
        for i in range(nx):
            for j in range(ny):
                    grid[j,i].density = rho_r[bi]
                    bi+=1

    ##读取晶格温度信息
    with open(lattice_temp_file, 'r') as f:
        # Row Lattice-temperature
        # Ignore the first three lines
        for i in range(3):
            next(f)
        # Read in the time stamp, grid count, and atom count
        line = next(f).split()
        time_stamp = float(line[0])
        grid_count = int(line[1])
        # Read in the grid information
        lattice_temp = np.zeros(grid_count)
        # Update the material density
        for i in range(grid_count):
            line = next(f).split()
            lattice_temp[i] = float(line[4])#密度信息保存
        #print(rho_r)
        bi = 0
        for i in range(nx):
            for j in range(ny):
                    grid[j,i].lattice_temp = lattice_temp[bi]
                    bi+=1
                    #print(grid[i,j].lattice_temp) 

    # 将温度-密度-热容-晶格温度参数赋予网格grid矩阵
    for i in range(ny):
        for j in range(nx):       
            U[i,j] = grid[i,j].temperature
            Den_U[i,j] = grid[i,j].density
            Capa_U[i,j] = 100*U[i,j]# J/m^3 K
            T_lattice[i,j] = grid[i,j].lattice_temp

    Original_Latti_temp = T_lattice
    ##密度矩阵转化为电子密度矩阵
    Den_lattice = Den_U*(26.98/6.02214076E23)#lattice-density,which is not used now
    Den_U = (Den_U*2.99671213E30/1.806E29)# 1/m^3 # wuliangganghua
    ##--------------------------------------------------------------------------------
    ## 初始化模拟参数
    # 模拟开始时
    for i in range(ny):
        for j in range(nx):   
            ##这里在我之前的修改中改变了前面的100，的大小
            
            if Den_U[i,j] == 0 :#如果判断到矩阵此处的密度乘以热容，变成了0，即此处没有原子，直接将此处置零，以免运算出错
                #U[i,j] = max(U[i-1,j],U[i+1,j],U[i,j-1],U[i,j+1])
                k[i,j] = 0   
                U[i,j] = 0
                flag[i,j] = 0#置零方法为，提供一个flag变量，将此处的判定条件置为0
                a[i,j] = 0
                fq[i,j] = 0
                couple[i,j]= 0
                couple_lattice[i,j] = 0
                U[i,j] = 300
            # else:
            #     if grid[i,j].lattice_temp == 0:
            #         print(i,j)
            #         #T_lattice[i,j] = max(T_lattice[i-1,j],T_lattice[i+1,j],T_lattice[i,j-1],T_lattice[i,j+1])
            #         #U[i,j] = max(U[i-1,j],U[i+1,j],U[i,j-1],U[i,j+1])
                    
            #         k[i,j] = 0
            #         a[i,j] = 0
            #         flag[i,j] = 0
            #         fq[i,j] = 0
            else:
                k[i,j] = 236*(U[i,j]/T_lattice[i,j])
                couple[i,j] = 1/Capa_U[i,j]
                couple_lattice[i,j] = 1/2.376E6
                flag[i,j] = 1#判断是否是空域，空域则不需要计算
                #fq[i,j] = -10/(Den_U[i,j]*Capa_U[i,j])
                ##--------------------------------------------------------------------------
                if Den_U[i,j] >=0.6:
                    a[i,j] = Den_U[i,j]*k[i,j]/(Capa_U[i,j])##以后需要进行修改
                else:
                    a[i,j] = 0.01*Den_U[i,j]*k[i,j]/(Capa_U[i,j])
                #----------------------------------------------------------------------------

    # 模拟体系x方向单位长度
    dx = x/nx
    # 模拟体系y方向单位长度
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
    Grt = 2.45E17##-------------------------------此处之后需要进行修改
    ##--------------------------------------------------------------------------------
    ## 模拟主程序
    # 使用欧拉法进行时间离散-按步长进行迭代
    # max electron -- max lattice
    lattice_temp = pd.DataFrame(columns=['time','temp'])
    electron_temp = pd.DataFrame(columns=['time','temp'])
    index = 0
    for p in range(Nt+1):
        # 真实时间标志
        tt = p*dt
        # 扩展温度矩阵和其他的相关参数矩阵
        U_big  =  np.concatenate((U, U, U),axis=1)#manipulate in x 
        T_lattice_big = np.concatenate((T_lattice, T_lattice, T_lattice),axis = 1)
        Den_big = np.concatenate((Den_U, Den_U, Den_U), axis= 1)
        # x方向二阶导系数矩阵A
        A_big = (-2)*np.eye(3*nx,k=0) + (1)*np.eye(3*nx,k=-1) + (1)*np.eye(3*nx,k=1)
        # y方向二阶导系数矩阵B 
        B_big = (-2)*np.eye(ny,k=0) + (1)*np.eye(ny,k=-1) + (1)*np.eye(ny,k=1)
        #heat = np.zeros((ny,nx))
        #heat = 1.7E32*np.exp(xy_component)*np.exp(time_component)
        #if step >= 10000:
        #    heat = heat
        #else:
        #    heat[int(4/5*ny):int(5/5*ny),int(0/4*nx):int(4/4*nx)] = 1E22#*(step/50)

        #定义热源矩阵
        heat = np.zeros((ny,nx))
        #heat = 1.7E32*np.exp(xy_component)*np.exp(time_component)
        #S(y,t)= 2/√（π/ln2) * J_abs/ (L_op*t_p) * exp{-y/L_op} * exp{-4ln2* (t-t_p)^2/t_p^2} 
        #激光的能量在光斑范围内平均分布
        #J_abs 激光能量密度     L_op 光学穿透深度     t_p 激光脉冲宽度      y为激光传播方向的距离      
        J_abs=30    #文献中是30J/m^2     
        L_op=8E-9   #   m
        t_p=1E-13   #单位是s
        t=1E-13 
        if step >= 1000:
            heat = heat
        else:
            #heat[int(0/4*ny):int(1/4*ny),int(0/4*nx):int(4/4*nx)] = 1E22
            for j in range (int(1/4*ny)):
                #heat[j:]
                heat[j,:] = 0.93944*J_abs/(L_op*t_p) * np.exp(-(j*y/ny)/L_op) * np.exp(-2.77259*((t-t_p)/t_p)**2)
        heat_big = np.concatenate((heat, heat,  heat),axis = 1)
        flag_big = np.concatenate((flag, flag, flag),axis =1)
        couple_big = np.concatenate((couple, couple, couple),axis = 1)
        couple_lattice_big = np.concatenate((couple_lattice, couple_lattice, couple_lattice),axis =1)
        rx_big = np.concatenate((rx, rx, rx),axis =  1)
        ry_big = np.concatenate((ry, ry, ry),axis = 1)
        # qiangzhi bianliang zhuanhuan
        U_big = U_big.astype(float)
        rx_big = rx_big.astype(float)
        ry_big = ry_big.astype(float)
        A_big = A_big.astype(float)
        B_big = B_big.astype(float)
        heat_big =heat_big.astype(float)
        T_lattice_big = T_lattice_big.astype(float)

        U_big = U_big +rx_big*np.dot(U_big,A_big)+ ry_big*np.dot(B_big,U_big) + heat_big * dt*couple_big - (U_big-T_lattice_big)*Grt*dt*flag_big*couple_big
        T_lattice_big  = T_lattice_big  + (U_big -T_lattice_big )*Grt*dt*flag_big*couple_lattice_big 

        # # 温度矩阵1 = 温度矩阵0 + 温度矩阵x方向传热 + 温度矩阵y方向传热 +电子-晶格温度耦合项
        # U = U + rx*np.dot(U,A) + ry*np.dot(B,U) + heat * dt - (U-T_lattice)*Grt*dt*flag*couple#/(200/a)
        # # 晶格温度1 = 晶格温度0 + 电子-晶格温度耦合项
        # T_lattice = T_lattice + (U-T_lattice)*Grt*dt*flag*couple_lattice#/(200/a)
        # 设置边界条件
        #边界条件1-恒温

        #stime = time.time()
        

        rows, cols = flag_big.shape
        # 对flag_big矩阵进行上下左右偏移并累加
        flag_big = flag_big.astype(np.int32)
        shifted_flags = np.zeros((rows, cols), dtype=int)
        # 定义上下左右四个方向的偏移
        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        # 遍历原始矩阵的每个元素
        shifted_flags = sum(np.roll(flag_big, shift, axis=(0, 1)) for shift in directions)

        # 修改温度值
        modified_U_big = np.copy(U_big)

        # 找出满足条件的点的坐标
        condition = (shifted_flags > 0) & (flag_big ==0)# & (shifted_flags < 4)
        condition_points = np.argwhere(condition)

        # 遍历满足条件的点，并更新温度值
        for point in condition_points:
            i, j = point
            if 0 < i < rows - 1 and 0 < j < cols - 1:
                neighbors = [
                    modified_U_big[i-1, j],  # 上方点
                    modified_U_big[i+1, j],  # 下方点
                    modified_U_big[i, j-1],  # 左方点
                    modified_U_big[i, j+1]   # 右方点
                ]
        
                max_neighbor_value = max(neighbors)  # 找到最大邻居值
                modified_U_big[i, j] = max_neighbor_value  # 更新温度值
        U_big[1:-1, 1:-1] = modified_U_big[1:-1, 1:-1]




        #etime =time.time()
        #print(etime-stime)




        # for i in range(1,U_big.shape[0]-1):
        #     for j in range(1,U_big.shape[1]-1):
        #         if Den_big[i,j] == 0 :
        #             U_update = U_big[i,j]

        #             if Den_big[i-1,j] != 0 or Den_big[i+1,j] != 0 or Den_big[i,j-1] != 0 or Den_big[i,j+1] != 0:
        #                 U_update = max(U_update, U_big[i-1,j], U_big[i+1,j],  U_big[i,j-1], U_big[i,j+1])
        #             else:
        #                 U_update = 300
        #             U_big[i,j] =  U_update

        ########## juere bianjietiaojian
        nambda = 1E20
        h = 1E-40
        ch = nambda/h
        U_big[:,0]  = (300)#+ch*U_big[:,1]/dx)/(1+ch/dx)
        U_big[:,-1] = (300)#+ch*U_big[:,-2]/dx)/(1+ch/dx)
        U_big[0,:]  = (300)#+ch*U_big[1,:]/dy)/(1+ch/dy)
        U_big[-1,:] = U_big[-2,:]#(300)#+ch*U_big[-2,:]/dy)/(1+ch/dy)
        

        size_x = U.shape[0]
        size_y = U.shape[1]
        size_x_big = U_big.shape[0]
        size_y_big = U_big.shape[1]
        U = U_big[(size_x_big-size_x)//2:(size_x_big-size_x)//2+size_x,(size_y_big-size_y)//2:(size_y_big-size_y)//2+size_y]
        
        
        T_lattice = T_lattice_big[(size_x_big-size_x)//2:(size_x_big-size_x)//2+size_x,(size_y_big-size_y)//2:(size_y_big-size_y)//2+size_y]
        end = time.time()
        # 每一步打印模拟的结果，输出图像
        index+=1
        if p%100 ==0:
            print('T = {:.2f} fs   max_U= {:.5f} max_lattice = {:.1f} Heat_contribution = {:.5f}'.format(tt*1E15,np.max(U),np.max(T_lattice),np.max(heat)*dt*np.max(couple)),'Step ={}'.format(step))
        if step%1 == 0 and p%100000000 == 0:
            #lattice_temp.loc[index]     = [tt,np.max(T_lattice)]
            #electron_temp.loc[index]    = [tt,np.max(U)]
            write_temperature_to_grid(U,'electron_t_{}.grid.dump'.format(step),step,x,y,z)


            #showcontourf(U,D,vmin = np.min(U),vmax=np.max(U)+1,timestep1=step, name = str('electron'))
            #showcontourf(T_lattice,D,vmax=2000,timestep1=step,name = str('lattice'))
        #break
    with open('lattice_temp.dat', 'a') as f:
        f.write("{},{}\n".format(tt,np.max(T_lattice)))
    with open('electron_temp.dat', 'a') as f:
        f.write("{},{}\n".format(tt,np.max(U)))
    ## 保存计算完成后的电子温度信息-用于下一次运行时读入
    for i in range(ny):
        for j in range(nx):
            grid[i, j].temperature = U[i,j]  
            grid[i, j].lattice_temp = T_lattice[i,j]  

    return grid, T_lattice, Original_Latti_temp