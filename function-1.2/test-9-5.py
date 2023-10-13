import numpy as np

flag_big = np.array([[0, 0, 0, 1, 0],
                     [0, 1, 1, 1, 0],
                     [0, 1, 1, 1, 1],
                     [0, 1, 1, 1, 0],
                     [0, 0, 1, 1, 0]])

# 定义上下左右四个方向的偏移
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# 使用np.roll实现矩阵的偏移和相加
neighbor_sum = sum(np.roll(flag_big, shift, axis=(0, 1)) for shift in directions)

# 打印新矩阵
print(neighbor_sum)
