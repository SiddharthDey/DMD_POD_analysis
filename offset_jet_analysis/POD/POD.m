clc,clear;

Nt = 2848; 

M1 = csvread('snapshot_matrix_half.csv');
U = M1.';
clear M1;
%% 

C_s = (U*U')/(Nt-1);
[A_s, LAM_s] = eig(C_s,'vector');
[lambda_s,ilam_s] = sort(LAM_s,'descend');
A_s = A_s(:, ilam_s);
PHI_s = U'*A_s;

%%
norm_values = vecnorm(PHI_s);
phi_size = size(PHI_s);
PHI = zeros(phi_size(1),phi_size(2));
for i = 1:length(norm_values)
   PHI(:,i) = PHI_s(:,i)/norm_values(i); 
end
%PHI = normc(PHI_s); % Spatial modes
A = U*PHI;

% figure;
% plot(lambda_s(1:10)/sum(lambda_s),'o','LineWidth',2.5);
% title('Normalized eigen values');
%saveas(gcf,'Normalized_eigen_values.jpg');

save('eigenvalues.mat','lambda_s');
save('A_s.mat','A_s');
save('PHI_s.mat','PHI_s');
save('PHI.mat','PHI');
save('A.mat','A');

plot(mu(1:10)/sum(mu),'o','LineWidth',2.5);
title('Normalized eigen values');
saveas(gcf,'Normalized_eigen_values.jpg');




