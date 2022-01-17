import mat73
import numpy as np


reconstructed_data_mat = mat73.loadmat('/home/aritra/pod/2d_offsetData/POD_reconstruction.mat')
reconstructed_data = reconstructed_data_mat['r_rank_reconstruction']
mesh_file = '/home/aritra/pod/2d_offsetData/sets1/1.4035981733264/U_surface_pod.vtk'

file_reconstruction = open('rank_5_reconstruction.vtk', 'w')


with open(mesh_file) as f:
    vtk_lines = f.readlines()
for k in range(435084):
    file_reconstruction.writelines(vtk_lines[k])

time_stamp = 2
reconstructed_data_u = np.real(reconstructed_data[:217070,time_stamp])
reconstructed_data_v = np.real(reconstructed_data[217070:,i])

for j in range(217070):
    L_real = [str(reconstructed_data_u[j]) + " " + str(reconstructed_data_v[j]) + " " + str(0.0)]

    file_real.writelines(L_real)
    file_real.write('\n')
    file_imag.writelines(L_imag)
    file_imag.write('\n')
file_real.close
file_imag.close
