clc,clear;
load('r_700/lambda_700.mat');
load('r_700/mu_700.mat');
load('r_700/f_700.mat');
load('r_700/Sigma2_700.mat');
load('r_700/time_dynamics_700.mat');
load('r_700/power_values_scaled_700.mat');
load('r_700/DMD_rank_scaled_700.mat');
load('error_array_2_20.mat');
load('error_array_21_50.mat');
load('error_array_51_100.mat');

% load('lambda.mat');
% load('mu.mat');
% load('f.mat');
% load('Sigma2.mat');
% load('time_dynamics.mat');
% load('power_values_scaled_100.mat');
% load('DMD_ranks_scaled.mat');
% load('error_array_2_20.mat');
% load('error_array_21_50.mat');
%%
clc;
singular_values = diag(Sigma2);
[f_sorted,f_sorted_ind] = sort(f,'ascend');

power_sorted = power_values(f_sorted_ind);
figure;
center = [0 0];
radius = 1;
figure;
theta = (0:1:100) *2*pi/100;
plot(real(lambda),imag(lambda),'o',cos(theta),sin(theta),'--','LineWidth',2);
% viscircles(center,radius);
title('Ritz Values');
figure;
plot(mu,'o','LineWidth',1);
title('mu = log(lambda)/ Del t');

figure;
plot(singular_values(1:20)/sum(singular_values),'o','LineWidth',2);
title('Normalized singular values ');
ylabel('Normalized Singular Values','FontSize',12,'FontWeight','bold') 
xlabel('Number of modes','FontSize',12,'FontWeight','bold') 
figure;

plot(power_values,'ok','LineWidth',1);
title('Power of the DMD modes');
ylabel('Power','FontSize',12,'FontWeight','bold') 
xlabel('Mode number','FontSize',12,'FontWeight','bold')

figure;
plot(f_sorted(length(f_sorted)/2+1:end)*2*pi,power_sorted(length(f_sorted)/2+1:end),'o','LineWidth',1);
hold on;
plot(f_sorted(length(f_sorted)/2+1:end)*2*pi,power_sorted(length(f_sorted)/2+1:end),'LineWidth',2);
ylabel('Power','FontSize',12,'FontWeight','bold') 
xlabel('Frequency of mode (Hz)','FontSize',12,'FontWeight','bold')
title('DMD/Power spectrum');
%% time dynamics
figure;
mode = 679;
Nt = 2848/2; 
dt = 2*1e-5;
f_Hz = f*2*pi;
t_start = 1.3700476872162208;
t = [t_start:dt:(Nt-1)*dt + t_start]; 
disp('The growth/decay rate is');
disp(real(mu(mode)));
disp('The frequency of oscillation is');
disp(f_Hz(mode));
figure;
plot(t,real(time_dynamics(mode,:)),'LineWidth',2);
title('time dynamics mode 5');
disp('The magnitude of Ritz value is');
disp(abs(lambda(mode)));
%% reconstruction loss with svd rank
error_array_final = zeros(99,1);
rank = [2:100];
error_array_final(1:19) = error_array(2:20);
error_array_final(20:49) = error_array_21_50;
error_array_final(50:99) = error_array_51_100(51:100);
figure;
plot(rank,error_array_final,'LineWidth',2.5);
title('Reconstruction loss with svd rank');
ylabel('Reconstruction loss','FontSize',12,'FontWeight','bold') 
xlabel('Mode number','FontSize',12,'FontWeight','bold')


