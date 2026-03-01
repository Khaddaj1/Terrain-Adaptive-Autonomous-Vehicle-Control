%% Vip_Wheel_NOzR_LoadEstimator_Recovery.m
% Full-car 7-DOF + baseline LQR
% FL: contact-recovery priority using NO-zr tire-force estimator
%   Ft_hat = muf*Zu1_ddot_hat - (kf*delta + bf*deldot) + uFL
%   Ft_hat clipped to >=0, drives smooth alpha blend (LQR <-> recovery)
%
% Controller does NOT use zr or zr_dot.
% Plant uses zr for simulation disturbance.

clear; clc; close all;

%% ================== VEHICLE PARAMETERS ==================
ms   = 1500;  muf = 59;  mur = 59;
Ip   = 2160;  Ir  = 460;
kf   = 35000; kr  = 38000;
bf   = 1000;  br  = 1100;
ktf  = 190000; ktr = 190000;
Tf   = 0.505; Tr  = 0.557;
a    = 1.4;   b   = 1.7;

%% ================== BUILD 2nd-ORDER MATRICES ==================
M = diag([ms Ip Ir muf muf mur mur]);

G = [ 1 a  Tf;     % FL
      1 a -Tf;     % FR
      1 -b -Tr;    % RL
      1 -b  Tr];   % RR

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

% Tire stiffness terms (linear plant)
K(4,4)=K(4,4)+ktf; K(5,5)=K(5,5)+ktf;
K(6,6)=K(6,6)+ktr; K(7,7)=K(7,7)+ktr;

%% ================== INPUT MAPS (ACTUATORS + ROAD) ==================
Fu = zeros(7,4);
Fu(1,:) = [ 1  1  1  1];        % heave
Fu(2,:) = [ a  a -b -b];        % pitch
Fu(3,:) = [ Tf -Tf  Tr -Tr];    % roll
Fu(4,:) = [-1  0  0  0];        % FL
Fu(5,:) = [ 0 -1  0  0];        % FR
Fu(6,:) = [ 0  0 -1  0];        % RL
Fu(7,:) = [ 0  0  0 -1];        % RR

Fr = zeros(7,4);
Fr(4,1)=ktf; Fr(5,2)=ktf;
Fr(6,3)=ktr; Fr(7,4)=ktr;

%% ================== STATE-SPACE ==================
Z7=zeros(7); I7=eye(7);
A  = [Z7 I7; -M\K -M\C];
Bu = [zeros(7,4); M\Fu];
Br = [zeros(7,4); M\Fr];

%% ================== BASELINE LQR ==================
Q = diag([ ...
    1e3, 2e4, 2e4, ...
    3e2*ones(1,4), ...
    3e2, 2e3, 2e3, ...
    1e2*ones(1,4) ]);

R = 1.0*eye(4);         % raise if you want less force
K_lqr = lqr(A, Bu, Q, R);

%% ================== ROAD (PLANT DISTURBANCE ONLY) ==================
v = 80/3.6;
wheelbase = a + b;
tau_rear  = wheelbase / v;

fs   = 500;
Tend = 12;
t    = (0:1/fs:Tend-1/fs).';

n0=0.1; wexp=2; Gd0=8192e-6;
n_min=0.01; n_max=20; rhoLR=0.4;

z_FL_psd = road_psd(t,fs,v,Gd0,n0,wexp,n_min,n_max);
z_FR_psd = rhoLR*z_FL_psd + sqrt(1-rhoLR^2)*road_psd(t,fs,v,Gd0,n0,wexp,n_min,n_max);
z_RL_psd = delay_sig(z_FL_psd,t,tau_rear);
z_RR_psd = delay_sig(z_FR_psd,t,tau_rear);

z_rms = 0.03; % harshness knob
z_FL_psd = z_FL_psd*z_rms/std(z_FL_psd);
z_FR_psd = z_FR_psd*z_rms/std(z_FR_psd);
z_RL_psd = z_RL_psd*z_rms/std(z_RL_psd);
z_RR_psd = z_RR_psd*z_rms/std(z_RR_psd);

T_event=0.12; h_bump=0.07; h_hole=-0.07;

z_FL_evt = bump_event(t,4.0, T_event,h_bump);
z_FR_evt = bump_event(t,4.5, T_event,h_hole);
z_RL_evt = bump_event(t,4.0+tau_rear, T_event,h_bump);
z_RR_evt = bump_event(t,4.5+tau_rear, T_event,h_hole);

z_FL = z_FL_psd + z_FL_evt;
z_FR = z_FR_psd + z_FR_evt;
z_RL = z_RL_psd + z_RL_evt;
z_RR = z_RR_psd + z_RR_evt;

zmax = 0.2;
z_FL = max(min(z_FL,zmax),-zmax);
z_FR = max(min(z_FR,zmax),-zmax);
z_RL = max(min(z_RL,zmax),-zmax);
z_RR = max(min(z_RR,zmax),-zmax);

zr = [z_FL.'; z_FR.'; z_RL.'; z_RR.'];  % 4 x N

%% ================== FL RECOVERY USING Ft_hat (NO zr) ==================
% alpha from Ft_hat:
% Ft_hat <= F_off -> alpha ~ 1 (recovery)
% Ft_hat >= F_on  -> alpha ~ 0 (LQR)
F_on  = 1200;     % [N] "good load"
F_off = 150;      % [N] "near unload"
tau_alpha = 0.10; % [s] blend smoothing

% Recovery law (no zr):
% u_rec = -Cw*Zu1dot - Kdel*delta - Cdel*deldot
Cw   = 7000;       % [N*s/m]
Kdel = 1.5e5;      % [N/m]
Cdel = 8e3;        % [N*s/m]

% Acceleration estimator (from Zu1_dot only)
tau_acc = 0.03;    % [s] filtering of a_hat (smaller = more responsive, noisier)

% Actuator shaping
u_max  = 15000;
tau_u  = 0.06;     % [s] command low-pass
du_max = 2.0e5;    % [N/s] slew limit

%% ================== SIMULATION ==================
x0 = zeros(14,1);

data.A = A; data.Bu = Bu; data.Br = Br;
data.K = K_lqr;

data.t = t; data.zr = zr;      % plant-only
data.a = a; data.Tf = Tf;

data.muf = muf; data.kf = kf; data.bf = bf;
data.F_on = F_on; data.F_off = F_off;
data.tau_alpha = tau_alpha;

data.Cw = Cw; data.Kdel = Kdel; data.Cdel = Cdel;
data.tau_acc = tau_acc;

data.u_max = u_max; data.tau_u = tau_u; data.du_max = du_max;

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[tSol,xSol] = ode45(@(tt,xx) dyn_closed(tt,xx,data), [t(1) t(end)], x0, opts);

%% ================== RECONSTRUCT HISTORIES ==================
Nsol = numel(tSol);
uSol     = zeros(Nsol,4);
alphaSol = zeros(Nsol,1);
FtHatSol = zeros(Nsol,1);
deltaSol = zeros(Nsol,1);

% IMPORTANT: reset controller persistents before replaying
clear control_FL_loadEstimator

for k=1:Nsol
    [uTmp, alphaSol(k), FtHatSol(k), deltaSol(k)] = control_FL_loadEstimator(tSol(k), xSol(k,:).', data);
    uSol(k,:) = uTmp.';
end

%% ================== PLOTS ==================
figure;
plot(t, z_FL, t, z_FR, t, z_RL, t, z_RR, 'LineWidth',1.1);
grid on; xlabel('Time [s]'); ylabel('Road [m]');
title('Road inputs (plant disturbance only)');
legend('FL','FR','RL','RR','Location','best');

figure;
plot(tSol, xSol(:,2), 'LineWidth',1.2); hold on;
plot(tSol, xSol(:,3), 'LineWidth',1.2);
grid on; xlabel('Time [s]'); ylabel('rad');
title('Pitch / Roll');
legend('\theta','\phi','Location','best');

figure;
plot(tSol, deltaSol, 'LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('\delta_{FL} = Z_{s1}-Z_{u1} [m]');
title('FL suspension deflection (used in F_{sus} and recovery)');

figure;
plot(tSol, FtHatSol, 'LineWidth',1.6); grid on;
xlabel('Time [s]'); ylabel('\hat F_{t,FL} [N] (clipped)');
title('Estimated FL tire force (NO z_r used)');
yline(F_off,'--','F_{off}','LineWidth',1.2);
yline(F_on,'--','F_{on}','LineWidth',1.2);

figure;
plot(tSol, alphaSol, 'LineWidth',1.6); grid on;
xlabel('Time [s]'); ylabel('\alpha');
title('Blend factor from \hat F_t: 0=LQR, 1=recovery');

figure;
ds=10;
tP = tSol(1:ds:end);
uFL = uSol(1:ds:end,1);
uFL_plot = smoothdata(uFL,'movmean',11);
plot(tP, uFL_plot, 'LineWidth',1.6); grid on;
xlabel('Time [s]'); ylabel('u_{FL} [N]');
title('FL actuator force (display downsampled + smoothed)');
yline(u_max,'--','LineWidth',1.1); yline(-u_max,'--','LineWidth',1.1);

disp('Done.');

%% =====================================================================
%% LOCAL FUNCTIONS
%% =====================================================================

function dx = dyn_closed(t, x, data)
    % Plant sees road input (disturbance)
    zr_t = interp1(data.t, data.zr.', t, 'linear', 'extrap').';  % 4x1

    % Controller does NOT use zr
    [u, ~, ~, ~] = control_FL_loadEstimator(t, x, data);

    dx = data.A*x + data.Bu*u + data.Br*zr_t;
end

function [u, alpha, Ft_hat, delta] = control_FL_loadEstimator(t, x, data)
    % States
    Zs   = x(1);  th  = x(2);  ph  = x(3);
    Zu1  = x(4);
    Zsd  = x(8);  thd = x(9);  phd = x(10);
    Zu1d = x(11);

    % FL sprung corner motion
    Zs1    = Zs  + data.a*th  + data.Tf*ph;
    Zs1dot = Zsd + data.a*thd + data.Tf*phd;

    % Suspension deflection and rate
    delta  = Zs1    - Zu1;
    deldot = Zs1dot - Zu1d;

    % Baseline LQR
    u_lqr = -(data.K * x);

    % --- Actuator command shaping persistents (u_applied) ---
    persistent u_applied prev_t_u
    if isempty(u_applied)
        u_applied = zeros(4,1);
        prev_t_u = t;
    end
    dt_u = max(0, t - prev_t_u);
    prev_t_u = t;

    % ==========================================================
    % 1) Estimate unsprung acceleration a_hat = d(Zu1dot)/dt filtered
    persistent prev_Zu1d a_hat prev_t_a
    if isempty(prev_Zu1d)
        prev_Zu1d = Zu1d;
        a_hat = 0;
        prev_t_a = t;
    end
    dt_a = max(0, t - prev_t_a);
    prev_t_a = t;

    if dt_a > 0
        a_raw = (Zu1d - prev_Zu1d) / dt_a;
        prev_Zu1d = Zu1d;

        % low-pass filter on acceleration estimate
        a_hat = a_hat + (dt_a/max(data.tau_acc,1e-6))*(a_raw - a_hat);
    end

    % 2) Suspension force on wheel (spring+damper)
    Fsus = data.kf*delta + data.bf*deldot;

    % 3) Tire force estimate (NO z_r)
    % Ft_hat = muf*a_u - Fsus + u_FL
    Ft_hat_raw = data.muf*a_hat - Fsus + u_applied(1);

    % Clip unilateral
    Ft_hat = max(Ft_hat_raw, 0);

    % ==========================================================
    % 4) Build alpha_target from Ft_hat (smooth blend)
    denom = max(data.F_on - data.F_off, 1e-6);
    alpha_target = (data.F_on - Ft_hat) / denom;  % low Ft => alpha->1
    alpha_target = min(max(alpha_target,0),1);

    persistent alpha_p prev_t_p
    if isempty(alpha_p)
        alpha_p = 0;
        prev_t_p = t;
    end
    dt = max(0, t - prev_t_p);
    prev_t_p = t;

    alpha_p = alpha_p + (dt/max(data.tau_alpha,1e-6))*(alpha_target - alpha_p);
    alpha_p = min(max(alpha_p,0),1);
    alpha = alpha_p;

    % 5) Recovery control (NO z_r)
    u_rec = -data.Cw*Zu1d - data.Kdel*delta - data.Cdel*deldot;

    % Blend ONLY FL
    u_cmd = u_lqr;
    u_cmd(1) = (1-alpha)*u_lqr(1) + alpha*u_rec;

    % Command saturation
    u_cmd = max(min(u_cmd, data.u_max), -data.u_max);

    % ==========================================================
    % 6) Apply actuator lag + slew rate limit (reduces chatter)
    if dt_u > 0
        u_filt = u_applied + (dt_u/max(data.tau_u,1e-6))*(u_cmd - u_applied);
    else
        u_filt = u_applied;
    end

    du = u_filt - u_applied;
    du_lim = data.du_max * max(dt_u,0);
    du = max(min(du, du_lim), -du_lim);

    u_applied = u_applied + du;

    % Final saturation
    u = max(min(u_applied, data.u_max), -data.u_max);
end

function z = bump_event(t, t0, T, h)
    z = zeros(size(t));
    idx = (t >= t0) & (t <= t0 + T);
    tau = (t(idx) - t0)/T;
    z(idx) = h*(1 - cos(pi*tau));
end

function z = road_psd(t,fs,v,Gd0,n0,w,nmin,nmax)
    N=numel(t); df=fs/N;
    f=(0:floor(N/2))'*df;
    n=f/max(v,eps); n(1)=n0;

    Gd = Gd0*(n/n0).^(-w);
    Gd(n<nmin | n>nmax) = 0;

    Sz = Gd/max(v,eps);
    Z1 = sqrt(2*Sz*df).*exp(1j*2*pi*rand(size(f)));
    Z1(1)=0;
    if mod(N,2)==0, Z1(end)=real(Z1(end)); end

    Zfull = [Z1; conj(Z1(end-1:-1:2))];
    z = real(ifft(Zfull));
    z = z - mean(z);
end

function y = delay_sig(x,t,tau)
    y = interp1(t,x,t-tau,'linear',0);
end
