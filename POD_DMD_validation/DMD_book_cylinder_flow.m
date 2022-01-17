%% load data
clear all, close all, clc
load DATA/FLUIDS/CYLINDER_ALL.mat
X1 = VORTALL(:,1:end-1);
X2 = VORTALL(:,2:end);

%% DMD
xi = linspace(-10,10,400) ;
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
r = 21;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
% [U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,r); Sigma=Sigma2(r,r); V=V2(:,r);

Atilde = U'*X2*V/Sigma;    
[W,D] = eig(Atilde);    
% Phi = U*W; %%%%%%% Projected DMD
Phi = X2*V/Sigma*W; %%%%%% Exact DMD
    
mu = diag(D);
omega = log(mu)/dt;

u0=X1(:,1);
y0 = Phi\u0; 
time_dynamics = zeros(r,length(t));
% u_modes = zeros(1,length(sim_time));
for iter = 1:length(t)
     time_dynamics(:,iter) =(y0.*exp(omega*t(iter)));
end
u_dmd = Phi*time_dynamics;
%%
plotCylinder(reshape(real(VORTALL(:,100)),nx,ny),nx,ny);

%% Plot DMD modes
for i=10:2:14
    plotCylinder(reshape(real(Phi(:,i)),nx,ny),nx,ny);
    plotCylinder(reshape(imag(Phi(:,i)),nx,ny),nx,ny);
end

%%  Plot DMD spectrum
figure;
center = [0 0];
radius = 1;
plot(mu,'o','LineWidth',3);
viscircles(center,radius);
title('Ritz Values','FontSize',12,'FontWeight','bold');
