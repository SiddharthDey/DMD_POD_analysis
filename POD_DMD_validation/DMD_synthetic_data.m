%% Define time and space discretizations
clc,clear;
xi = linspace(-10,10,400) ;
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
[Xgrid ,T] = meshgrid(xi,t);
%% Create two spatiotemporal patterns
f1 = sech(Xgrid+3) .* (1*exp(1j*2.3*T));
f2 = (sech(Xgrid).*tanh (Xgrid)).*(2*exp(1j*2.8*T));

%% Combine signals and make data matrix
f = f1 + f2;
X = f'; % Data Matrix

%% FFT
clc,clear;
% % Fs = 1/dt;
% % x = f1(:,30);
% % N = length(x);

t = (0:0.01:10*pi);
x = sin(5*t)+sin(2*t);
dt = 0.01;
Fs = 1/dt;
N = length(t);

% % xdft = fft(x);
% % xdft = xdft(1:N/2+1);
% % psdx = (1/(Fs*N)) * abs(xdft).^2;
% % psdx(2:end-1) = 2*psdx(2:end-1);
% % freq = 0:Fs/length(x):Fs/2;
% % plot(freq,psdx)
% % grid on
% % title('Periodogram Using FFT')
% % xlabel('Frequency (Hz)')
% % ylabel('Power/Frequency (dB/Hz)')

% % y = fft(x);
% % f = (0:N-1)*(Fs/N);     % frequency range
% % pow = abs(y).^2/N;    % power of the DFT
% % plot(f,pow), title('2hz power'); 

xhat = fft(x);
xpower = abs(xhat(1:N/2+1))*2/N;
freqs = Fs*(0:(N/2))/N;
plot(freqs ,xpower ,'k','LineWidth',1.2)

%% Visualize f1, f2, and f
% figure;
% subplot(2,2,1);
% surfl(real(f1));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});

% surfl(real(f1));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% xlabel('x','FontSize',12,'FontWeight','bold');
% ylabel('t','FontSize',12,'FontWeight','bold');
% title('f_1(x,t)');
% 
% surfl(real(f2));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% xlabel('x','FontSize',12,'FontWeight','bold');
% ylabel('t','FontSize',12,'FontWeight','bold');
% title('f_2(x,t)');
% 
surfl(real(f));
shading interp; colormap(gray); view (-20,60);
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('t','FontSize',12,'FontWeight','bold');
title('f(x,t)');

% subplot(2,2,2);
% surfl(real(f2));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% subplot(2,2,3);
% surfl(real(f));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});

%% DMD Method - Nathan Kutz Lecture

clc,clear;
U1 = importdata('X.mat');
% x = linspace(0,3,300);
% y = linspace(-0.5,0.5,100);
% sim_time = linspace(0,20,2000);
X = U1;
X1 = X(:,1:end-1);
X2 = X(:,2:end);
xi = linspace(-10,10,400) ;
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
r=2;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
% [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,r); Sigma=Sigma2(r,r); V=V2(:,r);

Atilde = U'*X2*V/Sigma;    
[W,D] = eig(Atilde);    
%Phi = U*W; %%%%%%% Projected DMD
Phi = X2*V/Sigma*W; %%%%%% Exact DMD
    
mu = diag(D);
omega = log(mu)/dt;

u0=U1(:,1);
y0 = Phi\u0; 
time_dynamics = zeros(r,length(t));
% u_modes = zeros(1,length(sim_time));
for iter = 1:length(t)
     time_dynamics(:,iter) =(y0.*exp(omega*t(iter)));
end
u_dmd = Phi*time_dynamics;

% test that y0 or b are the mode amplitudes
% % % % b = [y0(1) 0;0 y0(2)];
% % % % td = zeros(2,length(t)-1);
% % % % for i = 1:length(t)-1
% % % %    td(1,i) = D(1,1)^(i-1);
% % % %    td(2,i) = D(2,2)^(i-1);
% % % % end
% % % % u_dmd_bb = Phi*b*td;
% % % % u_check = u_dmd_bb - u_dmd(:,1:199);
% % % % reconstruction_loss_bb = norm(u_dmd_bb - u_dmd(:,1:199),'fro');
% bb = W\(U'*u0);
% error_bb = norm(bb-y0);

% q_new = zeros(length(y),length(x),length(sim_time));
% Utilde_k = u_dmd.';
% count1 = 1;
% for i=1:length(y)
%     for j =1:length(x)
%         q_new(i,j,:) = Utilde_k(:,count1);
%         count1 = count1+1;
%     end
% end

% timestamp_ind = 100;
% [X,Y] = meshgrid(x,y);
% Z = q_new(:,:,timestamp_ind);
mode1 = real(Phi(:,1)*time_dynamics(1,:));
mode2 = real(Phi(:,2)*time_dynamics(2,:));
figure;
plot(diag(Sigma2(1:20,1:20))/sum(diag(Sigma2)),'o','Linewidth',3)
title('normalized Singular Values');




surfl(real(u_dmd.'));
shading interp; colormap(gray); view (-20,60);
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('t','FontSize',12,'FontWeight','bold');
title('f(x,t)_{DMD} rank = 2');

surfl(mode1.');
shading interp; colormap(gray); view (-20,60);
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('t','FontSize',12,'FontWeight','bold');
title('Reconstruction on mode 2');

surfl(mode2.');
shading interp; colormap(gray); view (-20,60);
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('t','FontSize',12,'FontWeight','bold');
title('Reconstruction on mode 1');



figure;
subplot(1,2,1);
surfl(real(X));
shading interp; colormap(gray); view (-20,60);
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});

subplot(1,2,2);
surfl(real(u_dmd));
shading interp; colormap(gray); view (-20,60);
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
title('Reconstructed Data');

figure;
center = [0 0];
radius = 1;
plot(mu,'o','LineWidth',2.5);
viscircles(center,radius);
title('Ritz Values');
% figure;
% plot(omega,'o');
% title('log(Ritz values)');

figure;
plot(t,real(time_dynamics(1,:)/y0(1)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.8*t)),'o','LineWidth',2);
title('Mode 1 Temporal');

figure;
plot(t,real(time_dynamics(2,:)/y0(2)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.3*t)),'o','LineWidth',2);
title('Mode 2 Temporal');

reconstruction_loss = norm(X-u_dmd,'fro');

figure;
plot(xi',real(Phi(:,1)),'LineWidth',2);
hold on;
plot(xi',(sech(xi).*tanh(xi))/norm(sech(xi).*tanh(xi)),'o','LineWidth',2);
title('Mode 1 Spatial');

figure;
plot(xi',real(Phi(:,2)),'LineWidth',2);
hold on;
plot(xi',sech(xi+3)/norm(sech(xi+3)),'o','LineWidth',2);
title('Mode 2 Spatial');

%% DMD - Nathan Kutz Paper

clc,clear;
U1 = importdata('X.mat');
X = U1;
X1 = X(:,1:end-1);
X2 = X(:,2:end);
xi = linspace(-10,10,400) ;
t = linspace(0,4*pi,200);
dt = t(2) - t(1);

r=2;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);

Atilde = U'*X2*V/Sigma;
[W,D] = eig(Atilde);
% A_hat = (Sigma^-0.5)*Atilde*(Sigma^0.5);
% [W_hat,D_hat] = eig(A_hat);
% W = (Sigma^0.5)*W_hat;

Phi = X2*V/Sigma*W;

z0 = Phi\X(:,1);

% lambda = diag(D_hat);
lambda = diag(D);
mu = log(lambda)/dt;
f = imag(log(lambda)/dt)/(2*pi);

power_values = zeros(r,1);
for i=1:r
    power_values(i) = norm(Phi(:,i))^2;
end
[B,I] = sort(power_values,'descend');

timestamp = 100;
X_hat_t = Phi*(D^timestamp)*z0;
% X_hat_t = Phi*(D_hat^timestamp)*z0;

time_dynamics = zeros(2,200);
% % % for i = 1:length(time_dynamics)
% % %    time_dynamics(1,i) = D(1,1)^(i-1);
% % %    time_dynamics(2,i) = D(2,2)^(i-1);
% % % end
for i = 1:length(time_dynamics)
   time_dynamics(:,i) = lambda.^(i-1);
end

x1_tilde = U'*X(:,1);
b = W\x1_tilde;
error_amp = abs(b-z0);
error_amp_norm = norm(error_amp);

mode_1_recon = Phi(:,1)*z0(1)*time_dynamics(1,:);  
mode_2_recon = Phi(:,2)*z0(2)*time_dynamics(2,:);
mode_recon = mode_1_recon + mode_2_recon;
error_mat = abs(real(mode_recon) - real(X));
error_mat_norm = norm(error_mat,'fro');
error_vec = abs(real(X_hat_t)-real(mode_recon(:,timestamp-1)));
error_vec_norm = norm(error_vec);
% error_vec_1 = 
%disp(error_vec_norm);

figure;
plot(diag(Sigma2(1:20,1:20))/sum(diag(Sigma2)),'ko','Linewidth',[2])
title('normalized Singular Values');

% figure;
% subplot(1,2,1);
% surfl(real(X));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% 
% subplot(1,2,2);
% surfl(real(u_dmd));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% title('Reconstructed Data');

figure;
center = [0 0];
radius = 1;
plot(lambda,'o','LineWidth',2.5);
viscircles(center,radius);
title('Ritz Values');
% figure;
% plot(omega,'o');
% title('log(Ritz values)');

figure;
plot(t,real(time_dynamics(1,:)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.8*t)),'o','LineWidth',2);
title('Mode 1 Temporal');

figure;
plot(t,real(time_dynamics(2,:)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.3*t)),'o','LineWidth',2);
title('Mode 2 Temporal');

% % % reconstruction_loss = norm(X-u_dmd,'fro');

figure;
plot(xi',real(Phi(:,1)),'LineWidth',2);
hold on;
plot(xi',(sech(xi).*tanh(xi))/norm(sech(xi).*tanh(xi)),'o','LineWidth',2);
title('Mode 1 Spatial');

figure;
plot(xi',real(Phi(:,2)),'LineWidth',2);
hold on;
plot(xi',sech(xi+3)/norm(sech(xi+3)),'o','LineWidth',2);
title('Mode 2 Spatial');

%% POD - SVD method

clc,clear;
xi = linspace(-10,10,400) ;
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
D = importdata('X.mat');
U = D';

[L,SIG,R] = svd(U);
PHI = R; % PHI are the spatial modes
A = U*PHI; % A are the time coefficients
lambda = diag(SIG).^2; % lambda are the eigenvalues

% % [U, Sigma, V] = svd(D);
% % 
% % disp(A(:,2)./V(:,2));


figure;
plot(xi',real(PHI(:,1)),'LineWidth',2);
hold on;
plot(xi',(sech(xi).*tanh(xi))/norm(sech(xi).*tanh(xi)),'o','LineWidth',2);
title('Mode 1 Spatial');

figure;
plot(xi',real(PHI(:,2)),'LineWidth',2);
hold on;
plot(xi',sech(xi+3)/norm(sech(xi+3)),'o','LineWidth',2);
title('Mode 2 Spatial');

figure;
plot(t,real(A(:,1)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.8*t)),'o','LineWidth',2);
title('Mode 1 Temporal');

figure;
plot(t,real(A(:,2)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.3*t)),'o','LineWidth',2);
title('Mode 2 Temporal');

%% POD - using correlation matrix

clc,clear;
M = importdata('X.mat');
xi = linspace(-10,10,400) ;
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
Nt = length(t);
U = M.';
C_s = (U*U')/(Nt-1);
[A_s, LAM_s] = eig(C_s,'vector');
[lambda_s,ilam_s] = sort(LAM_s,'descend');
A_s = A_s(:, ilam_s);
PHI_s = U'*A_s;

PHI = normc(PHI_s); % Spatial modes
A = U*PHI;

figure;
plot(xi',real(U(:,2)),'LineWidth',2);
hold on;
plot(xi',sech(xi+3)/norm(sech(xi+3)),'o','LineWidth',2);
title('Mode 2 Spatial');

figure;
plot(t,real(V(1,:)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.8*t)),'o','LineWidth',2);
title('Mode 1 Temporal');

figure;
plot(t,real(V(2,:)),'LineWidth',2);
hold on;
plot(t,real(exp(1j*2.3*t)),'o','LineWidth',2);
title('Mode 2 Temporal');

%% plotting POD and DMD modes and temporal dynamics together

clc,clear;
U1 = importdata('X.mat');
% x = linspace(0,3,300);
% y = linspace(-0.5,0.5,100);
% sim_time = linspace(0,20,2000);
X = U1;
X1 = X(:,1:end-1);
X2 = X(:,2:end);
xi = linspace(-10,10,400) ;
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
r=2;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
% [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,r); Sigma=Sigma2(r,r); V=V2(:,r);

Atilde = U'*X2*V/Sigma;    
[W,D] = eig(Atilde);    
% Phi = U*W; %%%%%%% Projected DMD
Phi = X2*V/Sigma*W; %%%%%% Exact DMD
    
mu = diag(D);
omega = log(mu)/dt;

u0=U1(:,1);
y0 = Phi\u0; 
time_dynamics = zeros(r,length(t));
% u_modes = zeros(1,length(sim_time));
for iter = 1:length(t)
     time_dynamics(:,iter) =(y0.*exp(omega*t(iter)));
end
u_dmd = Phi*time_dynamics;

xi_POD = linspace(-10,10,400) ;
t_POD = linspace(0,4*pi,200);
dt_POD = t(2) - t(1);
D_POD = importdata('X.mat');


[U_POD,Sigma_POD,V_POD] = svd(D_POD);


% q_new = zeros(length(y),length(x),length(sim_time));
% Utilde_k = u_dmd.';
% count1 = 1;
% for i=1:length(y)
%     for j =1:length(x)
%         q_new(i,j,:) = Utilde_k(:,count1);
%         count1 = count1+1;
%     end
% end

% timestamp_ind = 100;
% [X,Y] = meshgrid(x,y);
% Z = q_new(:,:,timestamp_ind);
mode1 = real(Phi(:,1)*time_dynamics(1,:));
mode2 = real(Phi(:,2)*time_dynamics(2,:));
figure;
plot(diag(Sigma2(1:20,1:20))/sum(diag(Sigma2)),'o','Linewidth',3)
title('normalized Singular Values');




% surfl(real(u_dmd.'));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% xlabel('x','FontSize',12,'FontWeight','bold');
% ylabel('t','FontSize',12,'FontWeight','bold');
% title('f(x,t)_{DMD} rank = 2');
% 
% surfl(mode1.');
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% xlabel('x','FontSize',12,'FontWeight','bold');
% ylabel('t','FontSize',12,'FontWeight','bold');
% title('Reconstruction on mode 2');
% 
% surfl(mode2.');
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% xlabel('x','FontSize',12,'FontWeight','bold');
% ylabel('t','FontSize',12,'FontWeight','bold');
% title('Reconstruction on mode 1');
% 
% 
% 
% figure;
% subplot(1,2,1);
% surfl(real(X));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% 
% subplot(1,2,2);
% surfl(real(u_dmd));
% shading interp; colormap(gray); view (-20,60);
% set(gca , 'YTick', numel(t)/4 * (0:4)),
% set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
% set(gca , 'XTick', linspace(1,numel(xi),3)),
% set(gca , 'Xticklabel',{'-10', '0', '10'});
% title('Reconstructed Data');

% figure;
% center = [0 0];
% radius = 1;
% plot(mu,'o','LineWidth',2.5);
% viscircles(center,radius);
% title('Ritz Values');
% figure;
% plot(omega,'o');
% title('log(Ritz values)');

% figure;
% plot(t,real(time_dynamics(1,:)/y0(1)),'LineWidth',2);
% hold on;
% plot(t,real(exp(1j*2.8*t)),'o','LineWidth',2);
% title('Mode 1 Temporal');
% 
% figure;
% plot(t,real(time_dynamics(2,:)/y0(2)),'LineWidth',2);
% hold on;
% plot(t,real(exp(1j*2.3*t)),'o','LineWidth',2);
% title('Mode 2 Temporal');

reconstruction_loss = norm(X-u_dmd,'fro');

figure;
plot(xi',-real(Phi(:,1)),'o',xi',(sech(xi).*tanh(xi))/norm(sech(xi).*tanh(xi)),xi',-U_POD(:,1),'LineWidth',2);
% title('Mode 1 Spatial');
title('Mode 1 Spatial','FontSize',12,'FontWeight','bold');
legend({'DMD','True','POD'},'Location','northeast','FontSize',12)
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
plot(xi',real(Phi(:,2)),'o',xi',sech(xi+3)/norm(sech(xi+3)),xi',U_POD(:,2),'LineWidth',2);
title('Mode 2 Spatial','FontSize',12,'FontWeight','bold');
legend({'DMD','True','POD'},'Location','northeast','FontSize',12)
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
plot(t,real(time_dynamics(1,:)/y0(1)),'o',t,real(exp(1j*2.8*t)),t,10.0826*real(V_POD(:,1)),'LineWidth',3);
title('Mode 1 Time Dynamics','FontSize',12,'FontWeight','bold');
legend({'DMD','True','POD'},'Location','northeast','FontSize',12)
xlabel('Time','FontSize',12,'FontWeight','bold')

figure;
plot(t,real(time_dynamics(2,:)/y0(2)),'o',t,real(exp(1j*2.3*t)),t,10.0972*real(V_POD(:,2)),'LineWidth',3);
title('Mode 2 Time Dynamics','FontSize',12,'FontWeight','bold');
legend({'DMD','True','POD'},'Location','northeast','FontSize',12)
xlabel('Time','FontSize',12,'FontWeight','bold')



% figure;
% plot(xi',real(Phi(:,2)),'LineWidth',2);
% hold on;
% plot(xi',sech(xi+3)/norm(sech(xi+3)),'o','LineWidth',2);
% title('Mode 2 Spatial');
