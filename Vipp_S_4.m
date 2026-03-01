%% FULL-CAR 7-DOF + LQR + HARSH PSD + NORMAL BUMPS (PHYSICAL SEQUENCE)
clear; clc; close all;

%% ================== VEHICLE PARAMETERS ==================
ms   = 1500;  muf = 59;  mur = 59;
Ip   = 2160;  Ir  = 460;
kf   = 35000; kr  = 38000;
bf   = 1000;  br  = 1100;
ktf  = 190000; ktr = 190000;
Tf   = 0.505; Tr = 0.557;
a    = 1.4;   b  = 1.7;

%% ================== MASS, STIFFNESS, DAMPING ==================
M = diag([ms Ip Ir muf muf mur mur]);

G = [ 1 a  Tf;
      1 a -Tf;
      1 -b -Tr;
      1 -b  Tr];

E = [1 0 0 0 0 0 0;
     0 1 0 0 0 0 0;
     0 0 1 0 0 0 0];

Hs = G*E;

Eu = zeros(4,7);
Eu(1,4)=1; Eu(2,5)=1; Eu(3,6)=1; Eu(4,7)=1;

Rrel = Hs - Eu;

Kc = diag([kf kf kr kr]);
Cc = diag([bf bf br br]);

K = Rrel.'*Kc*Rrel;
C = Rrel.'*Cc*Rrel;

K(4,4)=K(4,4)+ktf; K(5,5)=K(5,5)+ktf;
K(6,6)=K(6,6)+ktr; K(7,7)=K(7,7)+ktr;

%% ================== INPUT MAPPING ==================
Fu = zeros(7,4);
Fu(1,:) = [1 1 1 1];
Fu(2,:) = [a a -b -b];
Fu(3,:) = [Tf -Tf Tr -Tr];
Fu(4,:) = [-1 0 0 0];
Fu(5,:) = [0 -1 0 0];
Fu(6,:) = [0 0 -1 0];
Fu(7,:) = [0 0 0 -1];

Fr = zeros(7,4);
Fr(4,1)=ktf; Fr(5,2)=ktf;
Fr(6,3)=ktr; Fr(7,4)=ktr;

%% ================== STATE SPACE ==================
Z7=zeros(7); I7=eye(7);
A=[Z7 I7; -M\K -M\C];
Bu=[zeros(7,4); M\Fu];
Br=[zeros(7,4); M\Fr];

%% ================== LQR ==================
Q=diag([ ...
    1e3 5e4 5e4 ...
    5e2*ones(1,4) ...
    5e2 5e3 5e3 ...
    2e2*ones(1,4)]);

R=0.1*eye(4);
K_lqr=lqr(A,Bu,Q,R);

Acl=A-Bu*K_lqr;
sys_cl=ss(Acl,Br,eye(14),zeros(14,4));

%% ================== TIME & SPEED ==================
v=80/3.6;                 % 80 km/h
wheelbase=a+b;
tau_rear=wheelbase/v;

fs=500; Tend=12;
t=(0:1/fs:Tend-1/fs).';

%% ================== HARSH PSD ROAD ==================
n0=0.1; w=2;
Gd0=8192e-6;
n_min=0.01; n_max=20;
rho=0.4;

z_FL_psd=road_psd(t,fs,v,Gd0,n0,w,n_min,n_max);
z_FR_psd=rho*z_FL_psd+sqrt(1-rho^2)*road_psd(t,fs,v,Gd0,n0,w,n_min,n_max);
z_RL_psd=delay(z_FL_psd,t,tau_rear);
z_RR_psd=delay(z_FR_psd,t,tau_rear);

% RMS boost
z_rms=0.03;
z_FL_psd=z_FL_psd*z_rms/std(z_FL_psd);
z_FR_psd=z_FR_psd*z_rms/std(z_FR_psd);
z_RL_psd=z_RL_psd*z_rms/std(z_RL_psd);
z_RR_psd=z_RR_psd*z_rms/std(z_RR_psd);

%% ================== NORMAL BUMPS / HOLES ==================
T_event=0.12;
h_bump=0.07; h_hole=-0.07;

z_FL_evt=bump(t,4.0,T_event,h_bump);
z_FR_evt=bump(t,4.5,T_event,h_hole);
z_RL_evt=bump(t,4.0+tau_rear,T_event,h_bump);
z_RR_evt=bump(t,4.5+tau_rear,T_event,h_hole);

%% ================== COMBINE & LIMIT ==================
z_FL=z_FL_psd+z_FL_evt;
z_FR=z_FR_psd+z_FR_evt;
z_RL=z_RL_psd+z_RL_evt;
z_RR=z_RR_psd+z_RR_evt;

z_FL=max(min(z_FL,0.2),-0.2);
z_FR=max(min(z_FR,0.2),-0.2);
z_RL=max(min(z_RL,0.2),-0.2);
z_RR=max(min(z_RR,0.2),-0.2);

zr=[z_FL.'; z_FR.'; z_RL.'; z_RR.'];

%% ================== SIMULATION ==================
[~,tSim,xCL]=lsim(sys_cl,zr.',t);
uCL=-(K_lqr*xCL.').';

%% ================== PLOTS ==================
figure;
plot(t,z_FL,t,z_FR,t,z_RL,t,z_RR,'LineWidth',1.3)
grid on; xlabel('Time [s]'); ylabel('Road [m]');
title('Road Inputs: Harsh PSD + Normal Bumps');
legend('FL','FR','RL','RR');

figure;
plot(tSim,xCL(:,2),tSim,xCL(:,3),'LineWidth',1.4)
grid on; xlabel('Time [s]'); ylabel('Angle [rad]');
title('Pitch and Roll'); legend('\theta','\phi');

figure;
plot(tSim,uCL,'LineWidth',1.3)
grid on; xlabel('Time [s]'); ylabel('Force [N]');
title('Actuator Forces'); legend('u1','u2','u3','u4');

%% ================== FUNCTIONS ==================
function z=bump(t,t0,T,h)
z=zeros(size(t));
idx=t>=t0 & t<=t0+T;
tau=(t(idx)-t0)/T;
z(idx)=h*(1-cos(pi*tau));
end

function z=road_psd(t,fs,v,Gd0,n0,w,nmin,nmax)
N=numel(t); df=fs/N;
f=(0:floor(N/2))'*df;
n=f/max(v,eps); n(1)=n0;
Gd=Gd0*(n/n0).^(-w);
Gd(n<nmin|n>nmax)=0;
Sz=Gd/max(v,eps);
A=sqrt(2*Sz*df).*exp(1j*2*pi*rand(size(f)));
A(1)=0;
Z=[A; conj(A(end-1:-1:2))];
z=real(ifft(Z)); z=z-mean(z);
end

function y=delay(x,t,tau)
y=interp1(t,x,t-tau,'linear',0);
end
