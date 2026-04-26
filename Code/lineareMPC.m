%% lineareMPC_YALMIP_QP.m  (wie dein lineareMPC.m, aber mit Δu + Terminalkosten)
clear all; close all; clc;
yalmip('clear');

% ----------------------------
% 1) PARAMETER
% ----------------------------
Ts = 0.1;
N  = 30;
v_ref = 3.0;

% Modell: [X, Y, vX, vY], Input: [aX, aY]
A = [1, 0, Ts, 0;
     0, 1, 0, Ts;
     0, 0, 1, 0;
     0, 0, 0, 1];

B = [0.5*Ts^2, 0;
     0, 0.5*Ts^2;
     Ts, 0;
     0, Ts];

nx = 4; nu = 2;

% Gewichte
Q_pos = 10;
Q_vel = 1;
R     = 0.1;
Rdu   = 1.0;      % Gewicht für Eingangsänderung Δu (anpassen)

P_pos = Q_pos;    % Terminalgewicht (einfach)
P_vel = Q_vel;

% Bounds
u_min = [-5; -5];
u_max = [ 5;  5];

vx_min = -10; vx_max = 10;
vy_min = -10; vy_max = 10;

% Ratenbeschränkungen (optional, anpassen; hier Beispielwerte)
du_min = [-2; -2];
du_max = [ 2;  2];

% ----------------------------
% 2) YALMIP Variablen / Parameter
% ----------------------------
x = sdpvar(nx, N+1, 'full');       % Zustände
u = sdpvar(nu, N,   'full');       % Inputs

x0 = sdpvar(nx,1);                 % Parameter: aktueller Zustand
u_prev = sdpvar(nu,1);             % Parameter: u_{-1} (letzter angewandter Input)
x_ref_param = sdpvar(2, N+1,'full');% Parameter: Referenz nur für Positionen (X,Y)

% ----------------------------
% 3) MPC-Problem: Constraints + Objective
% ----------------------------
constraints = [];
objective   = 0;

constraints = [constraints, x(:,1) == x0];

for k = 1:N
    % Dynamik
    constraints = [constraints, x(:,k+1) == A*x(:,k) + B*u(:,k)];

    % Input-Box
    constraints = [constraints, u_min <= u(:,k) <= u_max];

    % Geschwindigkeits-Box
    constraints = [constraints, vx_min <= x(3,k) <= vx_max];
    constraints = [constraints, vy_min <= x(4,k) <= vy_max];

    % Tracking-Fehler
    pos_error = x(1:2,k) - x_ref_param(:,k);
    vel_error = x(3:4,k) - [v_ref; 0];

    % Δu definieren (k=1 relativ zu u_prev, sonst relativ zu u(:,k-1))
    if k == 1
        du = u(:,k) - u_prev;
    else
        du = u(:,k) - u(:,k-1);
    end

    % Kosten: Position + Velocity + Input + Δu
    objective = objective + ...
        Q_pos*(pos_error'*pos_error) + ...
        Q_vel*(vel_error'*vel_error) + ...
        R*(u(:,k)'*u(:,k)) + ...
        Rdu*(du'*du);

    % (optional) Ratenbeschränkung
    constraints = [constraints, du_min <= du <= du_max];
end

% Terminalkosten (auf x(:,N+1))
posN_error = x(1:2,N+1) - x_ref_param(:,N+1);
velN_error = x(3:4,N+1) - [v_ref; 0];
objective  = objective + ...
    P_pos*(posN_error'*posN_error) + P_vel*(velN_error'*velN_error);

% ----------------------------
% 4) Optimizer
% ----------------------------
options = sdpsettings('verbose', 0, 'solver', 'quadprog');

% Eingangsvektor: [x0; u_prev; x_ref(:)]
mpc_solver = optimizer(constraints, objective, options, ...
                       [x0; u_prev; x_ref_param(:)], ...
                       u(:,1));

% ----------------------------
% 5) Referenztrajektorie + Simulation
% ----------------------------
N_sim = 150;
time = (0:N_sim)*Ts;

ref_X = v_ref*time;
ref_Y = 2*sin(0.3*time);
ref_trajectory = [ref_X; ref_Y];

x_current = [0; 0; v_ref; 0];
u_last    = [0; 0];        % u_{-1} zu Beginn

X_traj = zeros(2, N_sim+1);
U_traj = zeros(nu, N_sim);
X_traj(:,1) = x_current(1:2);

for t = 1:N_sim
    % Horizont-Referenz extrahieren
    current_ref = zeros(2, N+1);
    for k = 0:N
        ref_idx = min(t+k, N_sim+1);
        current_ref(:,k+1) = ref_trajectory(:, ref_idx);
    end

    % MPC lösen
    input_vector = [x_current; u_last; current_ref(:)];
    [u_opt, errorcode] = mpc_solver(input_vector);

    if errorcode ~= 0 || isempty(u_opt)
        warning('MPC konnte nicht gelöst werden. Verwende Null-Eingang.');
        u_opt = zeros(nu,1);
    end

    % Anwenden + System
    x_current = A*x_current + B*u_opt;

    % Log
    X_traj(:,t+1) = x_current(1:2);
    U_traj(:,t)   = u_opt;

    % u_{-1} updaten für Δu im nächsten Schritt
    u_last = u_opt;
end

% ----------------------------
% 6) Plot
% ----------------------------
figure('Position',[100,100,800,500]);
plot(ref_X, ref_Y, 'r-', 'LineWidth', 2, 'DisplayName','Referenz'); hold on;
plot(X_traj(1,:), X_traj(2,:), 'b--', 'LineWidth', 2.5, 'DisplayName','MPC'); grid on; axis equal;
xlabel('X-Position [m]'); ylabel('Y-Position [m]');
title('YALMIP-QP-MPC (mit \Delta u und Terminalkosten)');
legend('Location','northeast');
