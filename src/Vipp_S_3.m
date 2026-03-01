%% Full-car 7-DOF + LQR + Harsh road = (PSD roughness boosted) + (pothole/bump event)
% Road inputs: zr = [FL; FR; RL; RR] displacement [m]
% Control: u = -Kx (LQR)
%
% Key idea:
%   z_total(t) = z_PSD(t) + z_event(t)
%   PSD power boosted TWO ways:
%     (1) high ISO-level Gd0
%     (2) RMS scaling to target roughness (z_rms_target)

clear; clc; close all;

%% -------------------- Parameters (Table 2) --------------------
ms   = 1500;  muf  = 59;   mur  = 59;
Ip   = 2160;  Ir   = 460;
kf   = 35000; kr   = 38000;
Tf   = 0.505; Tr   = 0.557;
ktf  = 190000; ktr = 190000;
bf   = 1000;  br   = 1100;
a    = 1.4;   b    = 1.7;

%% -------------------- Build M, C, K --------------------
M = diag([ms, Ip, Ir, muf, muf, mur, mur]);

G = [ 1   a   Tf ;     % FL
      1   a  -Tf ;     % FR
      1  -b  -Tr ;     % RL
      1  -b   Tr ];    % RR

E = [1 0 0 0 0 0 0;     % Zs
     0 1 0 0 0 0 0;     % theta
     0 0 1 0 0 0 0];    % phi

Hs = G*E;

Eu       = zeros(4,7);
Eu(1,4)  = 1; Eu(2,5) = 1; Eu(3,6) = 1; Eu(4,7) = 1;

Rrel = Hs - Eu;

Kc = diag([kf kf kr kr]);
Cc = diag([bf bf br br]);

K = Rrel.' * Kc * Rrel;
C = Rrel.' * Cc * Rrel;

% Tire stiffness on unsprung
K(4,4) = K(4,4) + ktf;
K(5,5) = K(5,5) + ktf;
K(6,6) = K(6,6) + ktr;
K(7,7) = K(7,7) + ktr;

%% -------------------- Input mappings (CORRECT actuator mapping) --------------------
% Actuators act between sprung and unsprung: + on body DOFs, - on wheel DOFs
Fu = zeros(7,4);
Fu(1,:) = [ 1  1  1  1];          % heave
Fu(2,:) = [ a  a -b -b];          % pitch
Fu(3,:) = [ Tf -Tf  Tr -Tr];      % roll
Fu(4,:) = [-1  0  0  0];          % FL wheel
Fu(5,:) = [ 0 -1  0  0];          % FR wheel
Fu(6,:) = [ 0  0 -1  0];          % RL wheel
Fu(7,:) = [ 0  0  0 -1];          % RR wheel

% Road enters via tire stiffness: +k_t*zr in unsprung equations
Fr = zeros(7,4);
Fr(4,1) = ktf;  Fr(5,2) = ktf;
Fr(6,3) = ktr;  Fr(7,4) = ktr;

%% -------------------- First-order state-space --------------------
Z7 = zeros(7); I7 = eye(7);

A  = [ Z7     I7;
      -M\K  -M\C];       % 14x14

Bu = [zeros(7,4);
      M\Fu];             % 14x4

Br = [zeros(7,4);
      M\Fr];             % 14x4

%% -------------------- LQR --------------------
Q = diag([
    3e5,   ... % Zs
    5e4,   ... % theta
    5e6,   ... % phi
    5e2, 5e2, 5e2, 5e2, ... % Zu1..Zu4
    5e2,   ... % Zs_dot
    5e3,   ... % theta_dot
    5e3,   ... % phi_dot
    2e2, 2e2, 2e2, 2e2  ... % Zu_dot
]);

R =  0.0001*eye(4);
K_lqr = lqr(A, Bu, Q, R);

Acl = A - Bu*K_lqr;
sys_cl = ss(Acl, Br, eye(14), zeros(14,4)); % output all states

%% =====================================================================
%% Time + speed (80 km/h) + ROAD = PSD (boosted) + EVENT
%% =====================================================================
v = 80/3.6;             % [m/s]
wheelbase = a + b;
tau_rear  = wheelbase / v;

Tend = 12;              % [s]
fs   = 500;             % [Hz] (wideband)
dt   = 1/fs;
t    = (0:dt:Tend-dt).';
N    = numel(t);

%% -------------------- PSD rough road (ISO-like) --------------------
n0 = 0.1;               % [1/m]
w  = 2;

% (1) Boost PSD level (physically) -> VERY harsh
Gd0 = 8192e-6;          % [m^3] at n0 (extreme)

% Wide spatial bandwidth [1/m]
n_min = 0.01;
n_max = 20;             % large bandwidth (short wavelengths)

% Correlation
rho_LR = 0.35;

% Generate PSD roads
z_FL_psd = road_from_iso_psd_band(t, fs, v, Gd0, n0, w, n_min, n_max);
z_ind    = road_from_iso_psd_band(t, fs, v, Gd0, n0, w, n_min, n_max);
z_FR_psd = rho_LR*z_FL_psd + sqrt(max(0,1-rho_LR^2))*z_ind;

z_RL_psd = delay_signal(z_FL_psd, t, tau_rear);
z_RR_psd = delay_signal(z_FR_psd, t, tau_rear);

% (2) Boost PSD power by RMS scaling (guarantees "PSD matters")
z_rms_target = 0.07;    % [m] 3 cm RMS => VERY harsh broadband
z_FL_psd = z_FL_psd * (z_rms_target / std(z_FL_psd));
z_FR_psd = z_FR_psd * (z_rms_target / std(z_FR_psd));
z_RL_psd = z_RL_psd * (z_rms_target / std(z_RL_psd));
z_RR_psd = z_RR_psd * (z_rms_target / std(z_RR_psd));

%% -------------------- Deterministic harsh event (pothole/bump) --------------------
eventType = "pothole";  % "pothole" or "bump"

t_event = 6.0;          % [s] event start at front-left
T_event = 0.08;         % [s] shorter => harsher impact (more high freq)
depth   = 0.10;         % [m] 10 cm (harsh)

z_evt = zeros(size(t));
idx = (t >= t_event & t <= t_event + T_event);
tau = (t(idx) - t_event)/T_event;

switch eventType
    case "pothole"
        h = -abs(depth);
    case "bump"
        h =  abs(depth);
end

% Smooth, violent half-cosine event
z_evt(idx) = h*(1 - cos(pi*tau));
z_evt(idx) = -h*(1 - cos(pi*tau));

% Apply event on LEFT track only (you can change this)
z_FL_evt = z_evt;
z_FR_evt = zeros(size(t));
z_RL_evt = delay_signal(z_evt, t, tau_rear);
z_RR_evt = zeros(size(t));

%% -------------------- Combine road inputs --------------------
z_FL = z_FL_psd + z_FL_evt;
z_FR = z_FR_psd + z_FR_evt;
z_RL = z_RL_psd + z_RL_evt;
z_RR = z_RR_psd + z_RR_evt;

% Optional clamp (keeps things physical)
z_max = 0.20;  % 20 cm
z_FL = max(min(z_FL, z_max), -z_max);
z_FR = max(min(z_FR, z_max), -z_max);
z_RL = max(min(z_RL, z_max), -z_max);
z_RR = max(min(z_RR, z_max), -z_max);

zr = [z_FL.'; z_FR.'; z_RL.'; z_RR.'];    % 4xN

%% -------------------- Simulate closed-loop --------------------
Ucl = zr.';                                % Nx4
[~, tSim, xCL] = lsim(sys_cl, Ucl, t);

% Actuator forces
uCL = -(K_lqr * xCL.').';                  % Nx4

%% -------------------- Plots --------------------
% Road
figure;
plot(t, z_FL, 'LineWidth', 1.2); hold on;
plot(t, z_FR, 'LineWidth', 1.2);
plot(t, z_RL, 'LineWidth', 1.2);
plot(t, z_RR, 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('z_{ri}(t) [m]');
title('Road inputs = boosted PSD (harsh wideband) + event (pothole/bump)');
legend('FL','FR','RL','RR');

% Pitch/Roll and rates
figure;
plot(tSim, xCL(:,2), 'LineWidth', 1.5); hold on;
plot(tSim, xCL(:,3), 'LineWidth', 1.5);
plot(tSim, xCL(:,9), 'LineWidth', 1.5);
plot(tSim, xCL(:,10),'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('Amplitude');
title('Pitch/Roll and rates (LQR)');
legend('\theta_s','\phi_s','\dot\theta_s','\dot\phi_s');

% Actuator forces
figure;
plot(tSim, uCL(:,1), 'LineWidth', 1.5); hold on;
plot(tSim, uCL(:,2), 'LineWidth', 1.5);
plot(tSim, uCL(:,3), 'LineWidth', 1.5);
plot(tSim, uCL(:,4), 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('u_i [N]');
title('Actuator forces (LQR)');
legend('u_1','u_2','u_3','u_4');

% All states (separate figures)
stateNames = { ...
    'Z_s (m)', '\theta_s (rad)', '\phi_s (rad)', ...
    'Z_{u1} (m)', 'Z_{u2} (m)', 'Z_{u3} (m)', 'Z_{u4} (m)', ...
    '\dot Z_s (m/s)', '\dot\theta_s (rad/s)', '\dot\phi_s (rad/s)', ...
    '\dot Z_{u1} (m/s)', '\dot Z_{u2} (m/s)', '\dot Z_{u3} (m/s)', '\dot Z_{u4} (m/s)'};

for i = 1:14
    figure;
    plot(tSim, xCL(:,i), 'LineWidth', 1.2);
    grid on; xlabel('Time [s]'); ylabel(stateNames{i});
    title(['State response: ' stateNames{i}]);
end

%% -------------------- Quick diagnostics --------------------
fprintf('PSD RMS (FL): %.4f m\n', std(z_FL_psd));
fprintf('Total RMS (FL): %.4f m\n', std(z_FL));
fprintf('Total peak (FL): %.4f m\n', max(abs(z_FL)));

%% =====================================================================
%% Local functions (must be at end of script)
%% =====================================================================

function z = road_from_iso_psd_band(t, fs, v, Gd0, n0, w, n_min, n_max)
% Generate road displacement z(t) with ISO-like spatial PSD, band-limited.
% Gd(n) = Gd0*(n/n0)^(-w) [m^3], n in [1/m]
% Convert to temporal via f = v*n, and Sz(f) = (1/v) Gd(f/v) [m^2/Hz].

    N  = numel(t);
    df = fs/N;

    % One-sided temporal frequency grid
    f = (0:floor(N/2)).' * df;      % [Hz]
    n = f / max(v, eps);            % [1/m]

    % Spatial PSD (avoid n=0 singularity)
    n_safe = n;
    n_safe(1) = n0;
    Gd = Gd0 * (n_safe/n0).^(-w);

    % Band-limit
    band = (n >= n_min) & (n <= n_max);
    Gd(~band) = 0;

    % Temporal PSD
    Sz = (1/max(v, eps)) * Gd;      % [m^2/Hz]

    % Random phases
    phi = 2*pi*rand(size(f));

    % One-sided amplitude shaping
    A = sqrt(2*Sz*df);
    Z1 = A .* exp(1j*phi);

    % DC and Nyquist handling
    Z1(1) = 0;
    if mod(N,2)==0
        Z1(end) = real(Z1(end));
    end

    % Full spectrum (conjugate symmetry)
    Zfull = zeros(N,1);
    Zfull(1:numel(Z1)) = Z1;
    Zfull(numel(Z1)+1:end) = conj(Z1(end-1:-1:2));

    z = real(ifft(Zfull));
    z = z - mean(z);
end

function y = delay_signal(x, t, tau)
% Delay signal by tau seconds: y(t)=x(t-tau), zero-padded
    y = interp1(t, x, t - tau, 'linear', 0);
end
