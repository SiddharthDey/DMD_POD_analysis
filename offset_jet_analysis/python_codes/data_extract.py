#from fluidfoam import readmesh
import os
import numpy as np
from sys import getsizeof

sol = '/home/aritra/pod/2d_offsetData/sets1'

dirs = os.listdir(sol)
dirs_float = []
for i in range(len(dirs)):
    # print(len(dirs[i]))
    dirs_float.append(float(dirs[i]))
dirs_sorted_index = np.argsort(dirs_float)

print(dirs[dirs_sorted_index[0]])

num_nodes = 217070
snapshot_matrix = np.zeros((num_nodes*2,len(dirs)/2))
count_temporal = 0
for i in range(len(dirs)):
    if i%2!=0:
	continue
    timename = dirs[dirs_sorted_index[i]]
    # print(time_stamp)
    # if str(time_stamp)[-1]=='0':
    #     time_stamp = int(time_stamp)
    # timename = str(time_stamp)
    # print(timename)
    file_name = sol + '/' + timename + "/" + 'U_surface_pod.vtk'
    vtk_lines = []
    with open(file_name) as f:
        vtk_lines = f.readlines()
    num_nodes = 217070
    count_spatial = 0
    for j in range(435084,len(vtk_lines)):
        vel_str = vtk_lines[j].split()
        u_vel = float(vel_str[0])
        v_vel = float(vel_str[1])
        # w_vel = float(vel_str[2])
  
        snapshot_matrix[count_spatial,count_temporal] = u_vel
        snapshot_matrix[count_spatial+num_nodes,count_temporal] = v_vel

        count_spatial += 1
        # snapshot_matrix[count1,time_stamp] = w_vel
        # count1 += 1
    if i%100==0:
        print("the snapshot number is = ",count_temporal)
        # print('size of the snapshot matrix is = ', getsizeof(snapshot_matrix)/(1024*1024*1024))
        # print('size of the vtk_lines is = ', getsizeof(vtk_lines)/(1024*1024*1024))
    count_temporal += 1

np.savetxt('/home/aritra/pod/2d_offsetData/snapshot_matrix_half.csv', snapshot_matrix, delimiter=',')