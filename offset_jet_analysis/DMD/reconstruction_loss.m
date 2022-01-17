clc,clear;

Nt = 2848/2; 
X = csvread('snapshot_matrix_half.csv');
dt = 2*5e-5;

rank_array = zeros(50,1);
error_array = zeros(50,1);

for r=1:50
    [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
    clear U2,
    clear V2,
    Atilde = U'*X(:,2:end)*V/Sigma;
    [W,D] = eig(Atilde);

    Phi = X(:,2:end)*V/Sigma*W;

    z0 = Phi\X(:,1);
    clear U;
    clear V;

    timestamp = 2;
    X_hat_t = Phi*(D^timestamp)*z0;
    error_array(r) = sum((X_hat_t-X(3,:)).^2);
    clear X_hat_t;
end
save('error_array.mat','error_array');