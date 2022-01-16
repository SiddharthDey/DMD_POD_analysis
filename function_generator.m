clc,clear all;

x = linspace(0,3,300);
% x = linspace(0,3,150);
y = linspace(-0.5,0.5,100); 
% y = linspace(-3,3,300);
sim_time = linspace(0,20,2000);
a_n = [0.03,0.015,0.0075];
b_n = [0.05,0.035,0.02];
beta_n = [0.8,0.55,0.3];
gamma_n = [0.8,1.8,2.4];
q = zeros(length(y),length(x),length(sim_time));
q0_mat = zeros(length(y),length(x),length(sim_time));
q1_mat = zeros(length(y),length(x),length(sim_time));
q2_mat = zeros(length(y),length(x),length(sim_time));
q3_mat = zeros(length(y),length(x),length(sim_time));
i = 1;
for yy = y
    j = 1;
    for xx = x
        q0 = ones(1,2000)*exp(-(yy^2)/0.7);
        q1 = 0;
        q2 = 0;
        q3 = 0;
        for m = -10:10
            dn = (a_n(1)*xx) + b_n(1);
            q1 = q1+(((-1)^m)*exp(-(((xx-(beta_n(1)*m)-(gamma_n(1)*sim_time)).^2)/dn)-(yy^2/dn)));
            dn = (a_n(2)*xx) + b_n(2);
            q2 = q2+(((-1)^m)*exp(-(((xx-(beta_n(2)*m)-(gamma_n(2)*sim_time)).^2)/dn)-(yy^2/dn)));
            dn = (a_n(3)*xx) + b_n(3);
            q3 = q3+(((-1)^m)*exp(-(((xx-(beta_n(3)*m)-(gamma_n(3)*sim_time)).^2)/dn)-(yy^2/dn)));
        end
        q1 = q1*1;
        q2 = q2.*(exp(-sim_time/30)-0.1);
        q3 = q3.*(1-exp(-sim_time/20)+0.2);
        q0_mat(i,j,:) = q0;
        q1_mat(i,j,:) = q1;
        q2_mat(i,j,:) = q2;
        q3_mat(i,j,:) = q3;
%         q(i,j,:) = q0+q1+q2+q3;
        q(i,j,:) = q1+q2+q3;
        j = j+1;
    end
    disp(i);
    i = i+1;
end

U1 = zeros(length(sim_time),length(y)*length(x));
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        U1(:,count1) = q(i,j,:);
        count1 = count1+1;
    end
end

%% Data Visualization

clc;
timestamp_ind = 100;
x = linspace(0,3,300);
% x = linspace(0,3,150);
% y = linspace(-3,3,300);
y = linspace(-0.5,0.5,100);
[X,Y] = meshgrid(x,y);
Z0 = q0_mat(:,:,timestamp_ind);
Z1 = q1_mat(:,:,timestamp_ind);
Z2 = q2_mat(:,:,timestamp_ind);
Z3 = q3_mat(:,:,timestamp_ind);
Z = q(:,:,timestamp_ind);
% contour(q(:,:,timestamp_ind));
nexttile
contourf(X,Y,Z0)
title('q0')
nexttile
contourf(X,Y,Z1)
title('q1')
nexttile
contourf(X,Y,Z2)
title('q2')
nexttile
contourf(X,Y,Z3)
title('q3')
nexttile
contourf(X,Y,Z)
title('q=q0+q1+q2+q3')
% pcolor(q(:,:,timestamp_ind));
% surf(q(:,:,timestamp_ind), 'edgecolor', 'none'); view(2);

%% POD - Snapshot Method (with mean subtraction)

% clc;
% x = linspace(0,3,300);
% y = linspace(-0.5,0.5,100);
% sim_time = linspace(0,20,2000);
% Nt = length(sim_time);
% U = U1 - repmat(mean(U1,1), Nt, 1);
% C_s = (U*U')/(Nt-1);
% [A_s, LAM_s] = eig(C_s,'vector');
% [lambda_s,ilam_s] = sort(LAM_s,'descend');
% A_s = A_s(:, ilam_s);
% PHI_s = U'*A_s;
% k = 1;
% Utilde_k_s = A_s(:,k)*PHI_s(:,k)';

%% POD - SVD Method
clc,clear;
U1 = importdata('100y_300x_100m/U1.mat');
x = linspace(0,3,300);
y = linspace(-0.5,0.5,100);
sim_time = linspace(0,20,2000);
D = U1';
[U,Sigma,V] = svd(D,'econ');
mu = diag(Sigma);

%% POD - SVD Method (contd) - check different modes
mode = 1;
rank = 6;
A_mode = mu(mode)*U(:,mode)*V(:,mode)';

q_new = zeros(length(y),length(x),length(sim_time));
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        q_new(i,j,:) = A_mode(count1,:);
        count1 = count1+1;
    end
end

r = rank;
D_r = U(:,1:r)*Sigma(1:r,1:r)*V(:,1:r)';
q_D_r = zeros(length(y),length(x),length(sim_time));
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        q_D_r(i,j,:) = D_r(count1,:);
        count1 = count1+1;
    end
end

%visualize
timestamp_ind = 100;
[X,Y] = meshgrid(x,y);
Z = q_new(:,:,timestamp_ind);
Z_r = q_D_r(:,:,timestamp_ind);
figure;
plot(mu(1:10)/sum(mu),'o');
title('Normalized eigen values');
figure;
contourf(X,Y,Z);
title('mode reconstruction');
figure;
contourf(X,Y,Z_r);
title('rank r modes reconstruction');


%% POD - Snapshot Method (without mean subtraction)

clc,clear;
U1 = importdata('100y_300x_100m/U1.mat');
x = linspace(0,3,300);
y = linspace(-0.5,0.5,100);
sim_time = linspace(0,20,2000);
Nt = length(sim_time);
U = U1;
C_s = (U*U')/(Nt-1);
[A_s, LAM_s] = eig(C_s,'vector');
[lambda_s,ilam_s] = sort(LAM_s,'descend');
A_s = A_s(:, ilam_s);
PHI_s = U'*A_s;

PHI = normc(PHI_s); % Spatial modes
A = U*PHI;  % Time coefficients

%% POD - Snapshot Method - contd -  (check different modes)
k = 1;
Utilde_k_s = A_s(:,k)*PHI_s(:,k)';

A_t = zeros(2000,2);
A_t(:,1) = sim_time;
A_t(:,2) = A(:,k);

q_new = zeros(length(y),length(x),length(sim_time));
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        q_new(i,j,:) = Utilde_k_s(:,count1);
        count1 = count1+1;
    end
end


%frequency domain analysis
a_t = A_t(:,2);
L = length(a_t);
fft_a_t = fft(a_t);
PSD_a_t = abs(fft_a_t).^2;
Fs = 100;
f = Fs*(0:(L/2))/L;
P_a_t = PSD_a_t(1:L/2+1);

%visualize
timestamp_ind = 200;
[X,Y] = meshgrid(x,y);
Z = q_new(:,:,timestamp_ind);
figure;
plot(lambda_s(1:10)/sum(lambda_s(1:100)),'o','LineWidth',2.5);
% plot(lambda_s(1:6).^0.5/sum(lambda_s(1:100).^0.5),'o','LineWidth',2.5);
title('Normalized eigen values')
figure;
contourf(X,Y,Z);
title('q mode')
figure;
plot(f,P_a_t);
title('PSD')
figure;
plot(f,20*log(P_a_t));
title('20log(PSD)')

%% POD - Direct Method (with mean subtraction)

% clc;
% x = linspace(0,3,300);
% y = linspace(-0.5,0.5,100);
% sim_time = linspace(0,20,2000);
% Nt = length(sim_time);
% U = U1 - repmat(mean(U1,1), Nt, 1);
% C = (U'*U)/(Nt-1);
% [PHI LAM] = eig(C,'vector');
% [lambda,ilam] = sort(LAM,'descend');
% PHI = PHI(:, ilam);
% A = U*PHI;
% k = 1;
% Utilde_k = A(:,k)*PHI(:,k)';
% 
% plot(lambda(1:10)/sum(lambda))

%% POD - Direct Method (without mean subtraction)

clc;
x = linspace(0,3,300);
y = linspace(-0.5,0.5,100);
sim_time = linspace(0,20,2000);
Nt = length(sim_time);
U = U1;
C = (U'*U)/(Nt-1);
[PHI LAM] = eig(C,'vector');
[lambda,ilam] = sort(LAM,'descend');
PHI = PHI(:, ilam);
A = U*PHI;
k = 1;
Utilde_k = A(:,k)*PHI(:,k)';

q_new = zeros(length(y),length(x),length(sim_time));
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        q_new(i,j,:) = Utilde_k(:,count1);
        count1 = count1+1;
    end
end

timestamp_ind = 1800;
[X,Y] = meshgrid(x,y);
Z = q_new(:,:,timestamp_ind);
figure;
plot(lambda(1:10)/sum(lambda))
figure;
contourf(X,Y,Z);

%% DMD Method - Nathan Kutz Lecture

clc,clear;
U1 = importdata('100y_300x_100m/U1.mat');
x = linspace(0,3,300);
y = linspace(-0.5,0.5,100);
sim_time = linspace(0,20,2000);
X = U1.';
X1 = X(:,1:end-1);
X2 = X(:,2:end);
dt=sim_time(2)-sim_time(1);
r=50;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
% [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,r); Sigma=Sigma2(r,r); V=V2(:,r);

Atilde = U'*X2*V/Sigma;    
[W,D] = eig(Atilde);    
% Phi = U*W; %%%%%%% Projected DMD
Phi = X2*V/Sigma*W; %%%%%% Exact DMD
    
mu = diag(D);
omega = log(mu)/dt;

u0=U1(1,:).';
y0 = Phi\u0; 
u_modes = zeros(r,length(sim_time));
% u_modes = zeros(1,length(sim_time));
for iter = 1:length(sim_time)
     u_modes(:,iter) =(y0.*exp(omega*sim_time(iter)));
end
u_dmd = Phi*u_modes;

q_new = zeros(length(y),length(x),length(sim_time));
Utilde_k = u_dmd.';
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        q_new(i,j,:) = Utilde_k(:,count1);
        count1 = count1+1;
    end
end

timestamp_ind = 100;
[X,Y] = meshgrid(x,y);
Z = q_new(:,:,timestamp_ind);
figure;
plot(diag(Sigma2)/sum(diag(Sigma2)),'ko','Linewidth',[2])
title('normalized Singular Values');
figure;
contourf(X,Y,real(Z));
title('reconstructed data');
figure;
center = [0 0];
radius = 1;
plot(mu,'o');
viscircles(center,radius);
title('Ritz Values');
figure;
plot(omega,'o');
title('log(Ritz values)');

%%
norm_values = zeros(r,1);
for i=1:r
    norm_values(i) = norm(Phi(:,i));
end
[B,I] = sort(norm_values,'descend');

k = 31;
Phi_k = Phi(:,k);
u_mode_k = u_modes(k,:);
u_dmd_k = Phi_k*u_mode_k;
q_new = zeros(length(y),length(x),length(sim_time));
Utilde_k = u_dmd_k.';
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        q_new(i,j,:) = Utilde_k(:,count1);
        count1 = count1+1;
    end
end

% check_mat = zeros(30000,2000);
% for i = 1:50
%     Phi_k = Phi(:,i);
%     u_mode_k = u_modes(i,:);
%     u_dmd_k = Phi_k*u_mode_k;
%     check_mat = check_mat+u_dmd_k;
% end
% check_zeros = abs(check_mat-u_dmd);

timestamp_ind = 100;
[X,Y] = meshgrid(x,y);
Z = q_new(:,:,timestamp_ind);
% figure;
% plot(diag(Sigma2)/sum(diag(Sigma2)),'ko','Linewidth',[2])
% title('normalized Singular Values');
figure;
contourf(X,Y,real(Z));
title('reconstructed data');
% figure;
% center = [0 0];
% radius = 1;
% plot(mu,'o');
% viscircles(center,radius);
% title('Ritz Values');
% figure;
% plot(omega,'o');
% title('log(Ritz values)');

%% Nathan Kutz paper

clc,clear;
U1 = importdata('100y_300x_100m/U1.mat');
x = linspace(0,3,300);
y = linspace(-0.5,0.5,100);
sim_time = linspace(0,20,2000);
X = U1.'; 
X1 = X(:,1:end-1); 
X2 = X(:,2:end);
dt=sim_time(2)-sim_time(1); 
r=50;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);

Atilde = U'*X2*V/Sigma;
% [W,D] = eig(Atilde);
A_hat = (Sigma^-0.5)*Atilde*(Sigma^0.5);
[W_hat,D_hat] = eig(A_hat);
W = (Sigma^0.5)*W_hat;
Phi = X2*V/Sigma*W;

z0 = Phi\X(:,1);

lambda = diag(D_hat);
% lambda = diag(D);
mu = log(lambda)/dt;
f = imag(log(lambda)/dt)/(2*pi);

power_values = zeros(r,1);
for i=1:r
    power_values(i) = norm(Phi(:,i))^2;
end
[B,I] = sort(power_values,'descend');

%% visulalize

timestamp = 100;
% X_hat_t = Phi*(D^timestamp)*z0;
X_hat_t = Phi*(D_hat^timestamp)*z0;

error_matrix = abs(real(X_hat_t)-X(:,timestamp));
reconstruction_error = sum(error_matrix.^2,'all');

q_data = zeros(length(y),length(x));
Utilde_k_data = X(:,timestamp).';
count1_data = 1;
for i=1:length(y)
    for j =1:length(x)
        q_data(i,j) = Utilde_k_data(count1_data);
        count1_data = count1_data+1;
    end
end

q_new = zeros(length(y),length(x));
Utilde_k = X_hat_t.';
count1 = 1;
for i=1:length(y)
    for j =1:length(x)
        q_new(i,j) = Utilde_k(count1);
        count1 = count1+1;
    end
end


[X_grid,Y_grid] = meshgrid(x,y);
Z = q_new;
Z_data = q_data;
figure;
contourf(X_grid,Y_grid,real(Z));
title('reconstructed data');
figure;
contourf(X_grid,Y_grid,Z_data);
title('actual data');
figure;
center = [0 0];
radius = 1;
plot(lambda,'o','LineWidth',2.5);
viscircles(center,radius);
title('Ritz Values');
figure;
plot(mu,'ok','LineWidth',2.5);
title('lograthmic mapping vs frequency');
figure;
plot(power_values,'o','LineWidth',2.5);
title('Power of the DMD modes');

%% Nthan Kutz paper - iterate svd rank and get plot of reconstruction error ns rank

clc,clear;
U1 = importdata('100y_300x_100m/U1.mat');
x = linspace(0,3,300);
y = linspace(-0.5,0.5,100);
sim_time = linspace(0,20,2000);
X = U1.'; 
X1 = X(:,1:end-1); 
X2 = X(:,2:end);
dt=sim_time(2)-sim_time(1);
rank_array = zeros(50,1);
error_array = zeros(50,1);
for r=1:50
    [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
    Atilde = U'*X2*V/Sigma;
% [W,D] = eig(Atilde);
    A_hat = (Sigma^-0.5)*Atilde*(Sigma^0.5);
    [W_hat,D_hat] = eig(A_hat);
    W = (Sigma^0.5)*W_hat;
    Phi = X2*V/Sigma*W;

    z0 = Phi\X(:,1);

    lambda = diag(D_hat);
% lambda = diag(D);
    mu = log(lambda)/dt;
    f = imag(log(lambda)/dt)/(2*pi);

    power_values = zeros(r,1);
    for i=1:r
        power_values(i) = norm(Phi(:,i));
    end
    [B,I] = sort(power_values,'descend');

    timestamp = 100;
% X_hat_t = Phi*(D^timestamp)*z0;
    X_hat_t = Phi*(D_hat^timestamp)*z0;

    error_matrix = abs(real(X_hat_t)-X(:,timestamp));
    reconstruction_error = sum(error_matrix.^2,'all');
    rank_array(r) = r;
    error_array(r) = reconstruction_error;
    disp(r);
end
plot(rank_array,error_array,'LineWidth',2.5);