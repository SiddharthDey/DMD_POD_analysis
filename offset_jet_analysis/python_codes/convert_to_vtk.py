import mat73
import numpy as np


U1_mat = mat73.loadmat('/home/aritra/pod/2d_offsetData/U1_SVD.mat')
U1 = U1_mat['U']
mesh_file = '/home/aritra/pod/2d_offsetData/sets1/1.4035981733264/U_surface_pod.vtk'

for i in range(5):
    file_real = open('POD_mode_' + str(i) + '_real' + '.txt', 'w')
    file_imag = open('POD_mode_' + str(i) + '_imag' + '.txt', 'w')

    POD_mode_u_real = np.real(U1[:217070,i])
    POD_mode_u_imag = np.imag(U1[:217070,i])
    POD_mode_v_real = np.real(U1[217070:,i])
    POD_mode_v_imag = np.imag(U1[217070:,i])
    with open(mesh_file) as f:
        vtk_lines = f.readlines()
    for k in range(len(vtk_lines)):
        file_real.writelines(vtk_lines[k])
        file_imag.writelines(vtk_lines[k])
        # file_real.write('\n')
        # file_imag.write('\n')

    for j in range(217070):
        L_real = [str(POD_mode_u_real[j]) + " " + str(POD_mode_v_real[j]) + " " + str(0.0)]
        L_imag = [str(POD_mode_u_imag[j]) + " " + str(POD_mode_v_imag[j]) + " " + str(0.0)]
        file_real.writelines(L_real)
        file_real.write('\n')
        file_imag.writelines(L_imag)
        file_imag.write('\n')
    file_real.close
    file_imag.close
    if i==0:
        break
    print("The file number is = ", i)