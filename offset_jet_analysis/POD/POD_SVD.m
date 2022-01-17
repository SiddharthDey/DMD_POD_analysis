clc,clear;

Nt = 2848/2; 

D = csvread('snapshot_matrix_half.csv');
[U,Sigma,V] = svd(D,'econ');
mu = diag(Sigma);
reconstruction_error = zeros(50,1);
% for r = 1:50
%     r_rank_reconstruction = real(U(:,1:r)*Sigma(1:r,1:r)*V(:,1:r)');
%     error_matrix = abs(r_rank_reconstruction-D);
%     reconstruction_error(r) = sum(error_matrix.^2,'all');
% end

r = 5;
r_rank_reconstruction = real(U(:,1:r)*Sigma(1:r,1:r)*V(:,1:r)');
% save('U1_svd.mat','U');
% save('V1_svd.mat','V');
% save('eigenvalues1_svd.mat','mu');
% save('reconstruction_error.mat','reconstruction_error');
save('POD_reconstruction.mat','r_rank_reconstruction');