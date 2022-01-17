clc,clear;

Nt = 2848/2; 
X = csvread('snapshot_matrix.csv');
dt = 2*5e-5;

r=100;
[U2,Sigma2,V2] = svd(X(:,1:end-1),'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
clear U2,
clear V2,
Atilde = U'*X(:,2:end)*V/Sigma;
[W,D] = eig(Atilde);

Phi = X(:,2:end)*V/Sigma*W;

z0 = Phi\X(:,1);
clear X;
clear U;
clear V;
lambda = diag(D);
mu = log(lambda)/dt;
f = imag(log(lambda)/dt)/(2*pi);

timestamp = 2;
X_hat_t = Phi*(D^timestamp)*z0;

time_dynamics = zeros(r,Nt);
for i = 1:length(time_dynamics)
   for j = 1:length(lambda)
      time_dynamics(j,i) = lambda(j)^(i-1); 
   end
end

save('Sigma2.mat','Sigma2');
save('lambda.mat','lambda');
save('mu.mat','mu');
save('f.mat','f');
save('Phi.mat','lambda');
save('time_dynamics.mat','time_dynamics');
save('reconstructed_data.mat','X_hat_t');

