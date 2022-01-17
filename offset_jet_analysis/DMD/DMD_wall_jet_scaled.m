clc,clear;

Nt = 2848/2; 
X = csvread('snapshot_matrix.csv');
dt = 2*5e-5;

r=100;
[U2,Sigma2,V2] = svd(X(:,1:end-1),'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
clear U2,
clear V2,
Atilde = U'*X(:,2:end)*V/Sigma;
A_hat = (Sigma^-0.5)*Atilde*(Sigma^0.5);
[W_hat,D_hat] = eig(A_hat);
W = (Sigma^0.5)*W_hat;
Phi = X(:,2:end)*V/Sigma*W;

z0 = Phi\X(:,1);
clear X;
clear U;
clear V;
lambda = diag(D_hat);
mu = log(lambda)/dt;
f = imag(log(lambda)/dt)/(2*pi);

power_values = zeros(r,1);
for i=1:r
    power_values(i) = norm(Phi(:,i))^2;
end
[B,I] = sort(power_values,'descend');

timestamp = 2;
X_hat_t = Phi*(D_hat^timestamp)*z0;

time_dynamics = zeros(r,Nt);
for i = 1:length(time_dynamics)
   for j = 1:length(lambda)
      time_dynamics(j,i) = lambda(j)^(i-1); 
   end
end

save('Sigma2_scaled.mat','Sigma2');
save('lambda_scaled.mat','lambda');
save('mu_scaled.mat','mu');
save('f_scaled.mat','f');
save('Phi_scaled.mat','lambda');
save('time_dynamics_scaled.mat','lambda');
save('reconstructed_data_scaled.mat','X_hat_t');
save('power_values_scaled.mat','power_values');

