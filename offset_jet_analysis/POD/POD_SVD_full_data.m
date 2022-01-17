clc,clear;

Nt = 2848; 

D = csvread('snapshot_matrix.csv');
[U,Sigma,V] = svd(D,'econ');
mu = diag(Sigma);
save('U1_svd.mat','U');
save('V1_svd.mat','V');
save('eigenvalues1_svd.mat','mu');