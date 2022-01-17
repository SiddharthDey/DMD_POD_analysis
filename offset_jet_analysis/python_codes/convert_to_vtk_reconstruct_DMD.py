import mat73
import numpy as np


U1_mat = mat73.loadmat('/home/aritra/pod/2d_offsetData/reconstructed_data.mat')
U1 = U1_mat['X_hat_t']
mesh_file = '/home/aritra/pod/2d_offsetData/sets1/1.4035981733264/U_surface_pod.vtk'


file_real = open('DMD_mode_' + str(mode) + '_real' + '.vtk', 'w')
# file_imag = open('DMD_mode_' + str(mode) + '_imag' + '.vtk', 'w')

DMD_reconstructed_u_real = np.real(U1[:217070,0])
# DMD_mode_u_imag = np.imag(U1[:217070,mode])
DMD_reconstructed_v_real = np.real(U1[217070:,0])
# DMD_mode_v_imag = np.imag(U1[217070:,mode])

with open(mesh_file) as f:
    vtk_lines = f.readlines()
for k in range(435084):
    file_real.writelines(vtk_lines[k])
    # file_imag.writelines(vtk_lines[k])

for j in range(217070):
    L_real = [str(DMD_reconstructed_u_real[j]) + " " + str(DMD_reconstructed_v_real[j]) + " " + str(0.0)]
    # L_imag = [str(DMD_mode_u_imag[j]) + " " + str(DMD_mode_v_imag[j]) + " " + str(0.0)]
    file_real.writelines(L_real)
    file_real.write('\n')
    # file_imag.writelines(L_imag)
    # file_imag.write('\n')
file_real.close
file_imag.close
