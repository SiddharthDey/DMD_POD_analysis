clc,clear;
load('V1_SVD');
Nt = 2848/2; 
dt = 2*1e-5;
t_start = 1.3700476872162208;
t = [t_start:dt:(Nt-1)*dt + t_start]; 
plot(t,real(V(1,:)),'LineWidth',1);
title('time dynamics mode 1');
% disp('The magnitude of Ritz value is');
mode_time_dynamics = zeros(length(t),2);
mode_time_dynamics(:,1) = t;
mode_time_dynamics(:,2) = V(:,10);

%frequency domain analysis
% % % mode_time_dynamics = V(:,1);
% % % L = length(mode_time_dynamics);
% % % fft_a_t = fft(mode_time_dynamics);
% % % PSD_a_t = abs(fft_a_t).^2;
% % % Fs = 1/(5e-5);
% % % f = Fs*(0:(L/2))/L;
% % % P_a_t = PSD_a_t(1:L/2+1);
