%% Data Generation
clc,clear;
TIME = 10; % t o t a l time o f the animation
W = 10; % Womersley number
Pa = 60; % dimens i onl e s s pr e s sur e c o e f f i c i e n t
modes = 5; % number o f modes f o r a l l approximat ions (max 7)
dt = 0.05; % time s t ep
dy = 0.05; % space s t ep

t = [0:dt:TIME]; 
n_t = length(t); 
y = [-1:dy:1]; 
n_y = length(y); % space d i s c r e t i z a t i o n
u_A_r = zeros(n_y,n_t); 

for j = 1:length(t)
Y = (1-cosh(W*sqrt(1i).*y)./(cosh(W*sqrt(1i))))*1i*Pa/W.^2 ;
u_A_r(:,j) = real(Y.*exp (1i*t(j)));
end

% Adding the mean f low component :
u_Mb = (1-y.^2 )*0.5 ; % mean f low
u_M = repmat (u_Mb, length(t),1) ; % r epe a t s o l u t i o n
u_A_R = u_M' + u_A_r ;

%% visualize data

sim_time = t(41);
figure;
plot(y,u_A_R(:,41),'LineWidth',2);
title('data at t=2sec');

%% POD - Snapshot Method (without mean subtraction)

clc,clear;
U1 = importdata('u_A_R.mat');
TIME = 10;
dt = 0.05;
dy = 0.05;
t = [0:dt:TIME]; 
Nt = length(t);
y = [-1:dy:1]; 
U = U1';
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

r = 1;
D_r_11 = A_s(:,1:r)*PHI_s(:,1:r)';
D_r_1 = D_r_11';
r = 2;
D_r_22 = A_s(:,1:r)*PHI_s(:,1:r)';
D_r_2 = D_r_22';
r = 3;
D_r_33 = A_s(:,1:r)*PHI_s(:,1:r)';
D_r_3 = D_r_33';

% A_t = zeros(2000,2);
% A_t(:,1) = sim_time;
% A_t(:,2) = A(:,k);

% q_new = zeros(length(y),length(x),length(sim_time));
% count1 = 1;
% for i=1:length(y)
%     for j =1:length(x)
%         q_new(i,j,:) = Utilde_k_s(:,count1);
%         count1 = count1+1;
%     end
% end


%frequency domain analysis
% % % a_t = A_t(:,2);
% % % L = length(a_t);
% % % fft_a_t = fft(a_t);
% % % PSD_a_t = abs(fft_a_t).^2;
% % % Fs = 100;
% % % f = Fs*(0:(L/2))/L;
% % % P_a_t = PSD_a_t(1:L/2+1);

%visualize
% timestamp_ind = 200;
% [X,Y] = meshgrid(x,y);
% Z = q_new(:,:,timestamp_ind);
% figure;
% plot(lambda_s(1:10)/sum(lambda_s(1:100)),'o','LineWidth',2.5);
% % plot(lambda_s(1:6).^0.5/sum(lambda_s(1:100).^0.5),'o','LineWidth',2.5);
% title('Normalized eigen values')
% figure;
% contourf(X,Y,Z);
% title('q mode')
% figure;
% plot(f,P_a_t);
% title('PSD')
% figure;
% plot(f,20*log(P_a_t));
% title('20log(PSD)')
figure;
plot(y,D_r_1(:,41),'LineWidth',2);
hold on;
plot(y,D_r_2(:,41),'LineWidth',2);
hold on;
plot(y,D_r_3(:,41),'o','LineWidth',2);
hold on;
plot(y,U1(:,41),'LineWidth',2);
title('U1-U2-U3-UR at t=2secs');

figure;
plot(lambda_s(1:10).^0.5/(lambda_s(1)^0.5),'o','LineWidth',2);
title('Normalized eigen values');


nexttile;
plot(y,PHI(:,1)/norm(PHI),'LineWidth',2);
title('\phi 1');
nexttile;
plot(t,A(:,1),'LineWidth',2);
title('\psi 1');
nexttile;
plot(y,PHI(:,2)/norm(PHI),'LineWidth',2);
title('\phi 2');
nexttile;
plot(t,A(:,2),'LineWidth',2);
title('\psi 2');
nexttile;
plot(y,PHI(:,3)/norm(PHI),'LineWidth',2);
title('\phi 3');
nexttile;
plot(t,A(:,3),'LineWidth',2);
title('\psi 3');
nexttile;
plot(y,PHI(:,4)/norm(PHI),'LineWidth',2);
title('\phi 4');
nexttile;
plot(t,A(:,4),'LineWidth',2);
title('\psi 4');
nexttile;
plot(y,PHI(:,5)/norm(PHI),'LineWidth',2);
title('\phi 5');
nexttile;
plot(t,A(:,5),'LineWidth',2);
title('\psi 5');

%% POD - SVD Method
clc,clear;
U1 = importdata('u_A_R.mat');
TIME = 10;
dt = 0.05;
dy = 0.05;
t = [0:dt:TIME]; 
y = [-1:dy:1]; 
D = U1;
[U,Sigma,V] = svd(D);
mu = diag(Sigma);

%% POD - SVD Method (contd) - check different modes
clc;
mode = 1;
A_mode = mu(mode)*U(:,mode)*V(:,mode)';

figure;
plot(mu(1:10)/mu(1),'o','LineWidth',2);
title('Normalized singular values');
ylabel('Normalized Singular Values','FontSize',12,'FontWeight','bold') 
xlabel('Number of modes','FontSize',12,'FontWeight','bold') 

r = 1;
D_r_1 = U(:,1:r)*Sigma(1:r,1:r)*V(:,1:r)';
r = 2;
D_r_2 = U(:,1:r)*Sigma(1:r,1:r)*V(:,1:r)';
r = 3;
D_r_3 = U(:,1:r)*Sigma(1:r,1:r)*V(:,1:r)';

figure;
set(gca,'FontSize',40)
plot(y,D_r_1(:,71),y,D_r_2(:,71),y,D_r_3(:,71),'o', y,U1(:,71),'LineWidth',2);
title('U1-U2-U3-UR at t=3.5secs','FontSize',12,'FontWeight','bold');
legend({'U_1','U_2','U_3','U_R'},'Location','northeast','FontSize',12)
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
set(gca,'FontSize',40)
plot(y,D_r_1(:,41),y,D_r_2(:,41),y,D_r_3(:,41),'o', y,U1(:,41),'LineWidth',2);
title('U1-U2-U3-UR at t = 2secs','FontSize',12,'FontWeight','bold');
legend({'U_1','U_2','U_3','U_R'},'Location','northeast','FontSize',12)
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold') 

 
%visualize

% figure;
% plot(y,A_mode(:,41),'LineWidth',2);
% title('mode reconstruction');
% figure;
% plot(y,D_r(:,41),'LineWidth',2);
% title('rank r modes reconstruction');
% figure;
% plot(y,D_r_1(:,41),'LineWidth',2);
% hold on;
% plot(y,D_r_2(:,41),'LineWidth',2);
% hold on;
% plot(y,D_r_3(:,41),'o','LineWidth',2);
% hold on;
% plot(y,U1(:,41),'LineWidth',2);
% title('U1-U2-U3-UR at t=2secs');
% legend({'y = sin(x)','y = cos(x)'},'Location','southwest')

% figure;
% plot(y,D_r_1(:,41),'LineWidth',2);
% hold on;
% plot(y,D_r_2(:,41),'LineWidth',2);
% hold on;
% plot(y,D_r_3(:,41),'o','LineWidth',2);
% hold on;
% plot(y,U1(:,41),'LineWidth',2);
% title('U1-U2-U3-UR at t=2secs');
% legend({'y = sin(x)','y = cos(x)'},'Location','southwest')

figure;
plot(y,U(:,1)/norm(U),'LineWidth',3);
title('\phi 1','FontSize',12,'FontWeight','bold');
ylabel('POD mode','FontSize',12,'FontWeight','bold') 
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold') 

figure;
plot(y,U(:,2)/norm(U),'LineWidth',3);
title('\phi 2','FontSize',12,'FontWeight','bold');
ylabel('POD mode','FontSize',12,'FontWeight','bold') 
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
plot(y,U(:,3)/norm(U),'LineWidth',3);
title('\phi 3','FontSize',12,'FontWeight','bold');
ylabel('POD mode','FontSize',12,'FontWeight','bold') 
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
plot(t,V(:,3),'LineWidth',3);
title('Time dynamics mode 3','FontSize',12,'FontWeight','bold');
ylabel('Time dynamics','FontSize',12,'FontWeight','bold') 
xlabel('time','FontSize',12,'FontWeight','bold')

figure;
plot(t,V(:,1),'LineWidth',3);
title('Time dynamics mode 1','FontSize',12,'FontWeight','bold');
ylabel('Time dynamics','FontSize',12,'FontWeight','bold') 
xlabel('time','FontSize',12,'FontWeight','bold')

figure;
plot(t,V(:,2),'LineWidth',3);
title('Time dynamics mode 2','FontSize',12,'FontWeight','bold');
ylabel('Time dynamics','FontSize',12,'FontWeight','bold') 
xlabel('time','FontSize',12,'FontWeight','bold')



nexttile;
plot(y,U(:,1)/norm(U),'LineWidth',3);
title('\phi 1','FontSize',12,'FontWeight','bold');
ylabel('POD mode','FontSize',12,'FontWeight','bold') 
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold') 
nexttile;
plot(t,V(:,1),'LineWidth',3);
title('\psi 1');
nexttile;
plot(y,U(:,2),'LineWidth',3);
title('\phi 2');
nexttile;
plot(t,V(:,2),'LineWidth',3);
title('\psi 2');
nexttile;
plot(y,U(:,3),'LineWidth',3);
title('\phi 3');
nexttile;
plot(t,V(:,3),'LineWidth',3);
title('\psi 3');
nexttile;
plot(y,U(:,4),'LineWidth',3);
title('\phi 4');
nexttile;
plot(t,V(:,4),'LineWidth',3);
title('\psi 4');
nexttile;
plot(y,U(:,5),'LineWidth',3);
title('\phi 5');
nexttile;
plot(t,V(:,5),'LineWidth',3);
title('\psi 5');


%% DMD Method - Nathan Kutz Lecture

clc,clear;
U1 = importdata('u_A_R.mat');
TIME = 10;
dt = 0.05;
dy = 0.05;
t = [0:dt:TIME]; 
n_t = length(t);
y = [-1:dy:1]; 
X = U1; 
X1 = X(:,1:end-1); 
X2 = X(:,2:end);
dt=t(2)-t(1);

r=3;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);

Atilde = U'*X2*V/Sigma;    
% Atilde = U'*X2*V*inv(Sigma);    
[W,D] = eig(Atilde);
Phi = U*W; %%%%%%% Projected DMD
% Phi = X2*V/Sigma*W; %%%%%% Exact DMD

mu = diag(D);
omega = log(mu)/dt;

u0=U1(:,1);
y0 = Phi\u0; 
u_modes = zeros(r,length(t));
% u_modes = zeros(1,length(sim_time));
for iter = 1:length(t)
     u_modes(:,iter) =(y0.*exp(omega*t(iter)));
end
u_dmd = Phi*u_modes;


Z = u_dmd(:,41);
% figure;
% plot(diag(Sigma2)/sum(diag(Sigma2)),'ko','Linewidth',[2])
% title('normalized Singular Values');
figure;
plot(y,Z)
title('reconstructed data');
figure;
center = [0 0];
radius = 1;
plot(mu,'o','LineWidth',2);
viscircles(center,radius);
title('Ritz Values');

figure;
theta = (0:1:100) *2*pi/100;

center = [0 0 0];
radius = 1;
plot(mu,'o','LineWidth',2);
amplitudes = abs(y0)/max(abs(y0));
scatter3(cos(theta),sin(theta),zeros(length(theta),1),'o','LineWidth',2.5);
hold on;
scatter3(real(mu),imag(mu),amplitudes,'o','LineWidth',2.5);
title('Ritz values with height as normalized amplitudes','FontSize',12,'FontWeight','bold');

figure;
plot(y,Phi(:,1),'LineWidth',3);
title('\phi 1','FontSize',12,'FontWeight','bold');
ylabel('POD mode','FontSize',12,'FontWeight','bold') 
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold') 

figure;
plot(y,Phi(:,2),'LineWidth',3);
title('\phi 2','FontSize',12,'FontWeight','bold');
ylabel('POD mode','FontSize',12,'FontWeight','bold') 
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
plot(y,Phi(:,3),'LineWidth',3);
title('\phi 3','FontSize',12,'FontWeight','bold');
ylabel('POD mode','FontSize',12,'FontWeight','bold') 
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
plot(t,real(u_modes(3,:)/y0(3)),'LineWidth',3);
title('Time dynamics mode 3','FontSize',12,'FontWeight','bold');
ylabel('Time dynamics','FontSize',12,'FontWeight','bold') 
xlabel('time','FontSize',12,'FontWeight','bold')
ytickformat('%.2f')

figure;
plot(t,real(u_modes(1,:)/y0(1)),'LineWidth',3);
title('Time dynamics mode 1','FontSize',12,'FontWeight','bold');
ylabel('Time dynamics','FontSize',12,'FontWeight','bold') 
xlabel('time','FontSize',12,'FontWeight','bold')

figure;
plot(t,real(u_modes(2,:)/y0(2)),'LineWidth',3);
title('Time dynamics mode 2','FontSize',12,'FontWeight','bold');
ylabel('Time dynamics','FontSize',12,'FontWeight','bold') 
xlabel('time','FontSize',12,'FontWeight','bold')

% figure;
% plot(omega,'o');
% title('log(Ritz values)');

Z_2sec = zeros(3,length(y));
Z_3_5sec = zeros(3,length(y));
for r = 1:3
    [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
    % [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,r); Sigma=Sigma2(r,r); V=V2(:,r);

    Atilde = U'*X2*V/Sigma;    
    [W,D] = eig(Atilde);    
    Phi = X2*V/Sigma*W;

    mu = diag(D);
    omega = log(mu)/dt;

    u0=U1(:,1);
    y0 = Phi\u0; 
    u_modes = zeros(r,length(t));
    % u_modes = zeros(1,length(sim_time));
    for iter = 1:length(t)
         u_modes(:,iter) =(y0.*exp(omega*t(iter)));
    end
    u_dmd = Phi*u_modes;
    Z_2sec(r,:) = u_dmd(:,41);
    Z_3_5sec(r,:) = u_dmd(:,71);
end

% figure;
% plot(y,Z_2sec(1,:),'LineWidth',2);
% hold on;
% plot(y,Z_2sec(2,:),'LineWidth',2);
% hold on;
% plot(y,Z_2sec(3,:),'o','LineWidth',2);
% hold on;
% plot(y,U1(:,41),'LineWidth',2);
% title('U1-U2-U3-UR at t=2secs');

figure;
set(gca,'FontSize',40)
plot(y,Z_2sec(1,:),y,Z_2sec(2,:),y,Z_2sec(3,:),'o', y,U1(:,41),'LineWidth',2);
title('U1-U2-U3-UR at t = 2 secs','FontSize',12,'FontWeight','bold');
legend({'U_1','U_2','U_3','U_R'},'Location','northeast','FontSize',12)
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

figure;
set(gca,'FontSize',40)
plot(y,Z_3_5sec(1,:),y,Z_3_5sec(2,:),y,Z_3_5sec(3,:),'o', y,U1(:,71),'LineWidth',2);
title('U1-U2-U3-UR at t = 3.5 secs','FontSize',12,'FontWeight','bold');
legend({'U_1','U_2','U_3','U_R'},'Location','northeast','FontSize',12)
xlabel('Non-dimensional y-coordinate','FontSize',12,'FontWeight','bold')

% figure;
% plot(y,Z_3_5sec(1,:),'LineWidth',2);
% hold on;
% plot(y,Z_3_5sec(2,:),'LineWidth',2);
% hold on;
% plot(y,Z_3_5sec(3,:),'o','LineWidth',2);
% hold on;
% plot(y,U1(:,71),'LineWidth',2);
% title('U1-U2-U3-UR at t=3.5secs');




