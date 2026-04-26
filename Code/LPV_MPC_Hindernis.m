%% LPV_MPC_RefShift.m
% LPV-MPC (QP via YALMIP/quadprog)
% Kreis-Referenz + 4 kleine Quadrat-Hindernisse
% Idee: Referenz wird im Hindernisbereich lokal radial verschoben (innen/außen),
% danach wieder zurück (smooth mit alpha/beta)

clear; close all; clc;
yalmip('clear');

%% 1) Parameter
Ts    = 0.1;
N     = 40;
N_sim = 200;

lf = 1.105; lr = 1.738; L = lf + lr;

a_min = -5;   a_max = 2.5;
w_min = -0.6; w_max = 0.6;
v_min = 0;    v_max = 10;
d_min = -0.5; d_max = 0.5;

v_ref = 3.0;

% Tracking (stark, damit er nach dem Offset wieder zurück kommt)
Q_pos = 80;
Q_vel = 1.0;
Q_del = 0.05;
R_u   = 0.05;
R_du  = 0.2;

P_pos = 200;

nx = 5; nu = 2;

%% 2) Basis-Referenz: Kreis
time = (0:N_sim)*Ts;
omega_path = 2*pi/(N_sim*Ts);
theta_all = omega_path*time;

Rad = 10;
ref_base = [Rad*cos(theta_all); Rad*sin(theta_all)];

%% 3) Hindernisse als kleine Quadrate (nur Visualisierung) + Referenz-Shift-Parameter
sq_half   = 0.60; % Quadrat Halbseite [m]
obs_theta = deg2rad([45, 135, 225, 315]);

% +1  => Referenz nach außen schieben (radial nach außen)
% -1  => Referenz nach innen schieben (radial nach innen)
shift_side = [+1, +1, +1, +1];

% Shift-Stärke (wie weit Referenz verschoben wird)
shift_mag = 1.00;       

% Hindernisbereich (Winkel)
alpha = deg2rad(14);    % Kernbereich (klein = nur kurz)
beta  = deg2rad(8);     % Übergang weich rein/raus

%% 4) QP Setup (YALMIP) – nur Tracking auf verschobene Referenz
x = sdpvar(nx, N+1, 'full');
u = sdpvar(nu, N,   'full');

x0     = sdpvar(nx,1);
u_prev = sdpvar(nu,1);

Apar = sdpvar(nx,nx,'full');
Bpar = sdpvar(nx,nu,'full');
cpar = sdpvar(nx,1,'full');

% Referenz als Parameter (verschoben!)
rpar = sdpvar(2, N+1, 'full');

constraints = [x(:,1) == x0];
objective   = 0;

for k = 1:N
    constraints = [constraints, x(:,k+1) == Apar*x(:,k) + Bpar*u(:,k) + cpar];

    constraints = [constraints, a_min <= u(1,k) <= a_max];
    constraints = [constraints, w_min <= u(2,k) <= w_max];
    constraints = [constraints, v_min <= x(4,k) <= v_max];
    constraints = [constraints, d_min <= x(5,k) <= d_max];

    if k == 1
        du = u(:,k) - u_prev;
    else
        du = u(:,k) - u(:,k-1);
    end

    e_pos = x(1:2,k) - rpar(:,k);
    e_vel = x(4,k) - v_ref;
    e_del = x(5,k);

    objective = objective + ...
        Q_pos*(e_pos'*e_pos) + ...
        Q_vel*(e_vel'*e_vel) + ...
        Q_del*(e_del'*e_del) + ...
        R_u*(u(:,k)'*u(:,k)) + ...
        R_du*(du'*du);
end

% Terminal: starkes zurückziehen auf Referenz
e_posN = x(1:2,N+1) - rpar(:,N+1);
objective = objective + P_pos*(e_posN'*e_posN);

opts = sdpsettings('verbose',0,'solver','quadprog');

in_vars = [x0; u_prev; Apar(:); Bpar(:); cpar(:); rpar(:)];
mpc_qp  = optimizer(constraints, objective, opts, in_vars, u(:,1));

%% 5) Simulation
x_current = [Rad; 0; pi/2; v_ref; 0];
u_last    = [0; 0];

X_log = zeros(nx, N_sim+1);
X_log(:,1) = x_current;

% Zum Plotten: shifted reference über Zeit (optional)
ref_shifted_log = zeros(2, N_sim+1);
ref_shifted_log(:,1) = ref_base(:,1);

for t = 1:N_sim
    % LPV linearization (frozen)
    [A, B, c] = lin_bicycle_euler(x_current, u_last, Ts, lr, L);

    % Horizon-Referenz (verschoben) vorberechnen
    r_h = zeros(2, N+1);
    for kk = 0:N
        idx = min(t+kk, N_sim+1);

        r0 = ref_base(:,idx);
        th = theta_all(idx);

        r_h(:,kk+1) = shifted_reference_point(r0, th, obs_theta, shift_side, shift_mag, alpha, beta);
    end

    % log für plot (aktueller Schritt)
    ref_shifted_log(:,t+1) = r_h(:,2);

    % QP lösen
    in = [x_current; u_last; A(:); B(:); c(:); r_h(:)];
    [u0, err] = mpc_qp(in);
    if err ~= 0 || isempty(u0)
        u0 = zeros(nu,1);
    end

    % Pflanze (nichtlinear)
    x_current = f_bicycle_euler(x_current, u0, Ts, lr, L);
    u_last    = u0;

    X_log(:,t+1) = x_current;
end

%% 6) 1 Figure (Basis-Referenz + gefahrene Trajektorie + Quadrate)
figure('Position',[120,120,900,520]);

plot(ref_base(1,:), ref_base(2,:), 'r-', 'LineWidth', 2); hold on;
plot(X_log(1,:),    X_log(2,:),    'b--',  'LineWidth', 2.5);

% optional: verschobene Referenz anzeigen 
plot(ref_shifted_log(1,:), ref_shifted_log(2,:), 'g-.', 'LineWidth', 1.5);

% Quadrate zeichnen
for i = 1:length(obs_theta)
    th = obs_theta(i);
    p0 = Rad*[cos(th); sin(th)];
    tvec = [-sin(th); cos(th)];
    nvec = [ cos(th); sin(th)];

    corners = [ ...
        p0 + sq_half*tvec + sq_half*nvec, ...
        p0 - sq_half*tvec + sq_half*nvec, ...
        p0 - sq_half*tvec - sq_half*nvec, ...
        p0 + sq_half*tvec - sq_half*nvec ];
    corners = [corners, corners(:,1)];
    plot(corners(1,:), corners(2,:), 'k-', 'LineWidth', 1.6);
end

grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]');
legend('Referenz (Kreis)', 'LPV-MPC Trajektorie', 'verschobene Referenz', 'Quadrat-Hindernisse', 'Location','best');
title(sprintf('LPV-MPC, N=%d, Ts=%.2f', N, Ts));

%% ------------------ Funktionen (nur 1x am Ende) ------------------

function r_shift = shifted_reference_point(r0, theta, obs_theta, shift_side, shift_mag, alpha, beta)
% r0: Basis-Referenzpunkt auf Kreis
% theta: aktueller Winkel (atan2(y,x) oder vorgegeben)
% shift_side(i): +1 nach außen, -1 nach innen
% shift_mag: maximale Verschiebung [m]
% alpha/beta: Sektor und Übergang

% Radialvektor (Normalenrichtung)
n = [cos(theta); sin(theta)];

% Weight w: 1 im Kern, 0 außerhalb, linear in Transition
w_total = 0;
s_total = 0;

for i = 1:length(obs_theta)
    d  = wrapToPi_num(theta - obs_theta(i));
    ad = abs(d);

    if ad <= alpha
        w = 1;
    elseif ad <= alpha + beta
        w = 1 - (ad - alpha)/beta;
    else
        w = 0;
    end

    % mehrere Hindernisse addieren (wenn mal nahe)
    w_total = w_total + w;
    s_total = s_total + w * shift_side(i);
end

% begrenzen, damit es nicht explodiert
w_total = min(1, w_total);
s_total = max(-1, min(1, s_total)); % [-1,1]

offset = shift_mag * w_total * s_total;  % >0 nach außen, <0 nach innen

r_shift = r0 + offset*n;
end

function x_next = f_bicycle_euler(x, u, Ts, lr, L)
X=x(1); Y=x(2); phi=x(3); v=x(4); del=x(5);
a=u(1); om=u(2);

k = lr/L;
    beta  = atan(k * tan(del));

x_next = [ ...
    X   + Ts*(v*cos(phi + beta));
    Y   + Ts*(v*sin(phi + beta));
    phi + Ts*((v/L)*cos(beta)*tan(del));
    v   + Ts*a;
    del + Ts*om ];
end

function [A,B,c] = lin_bicycle_euler(x,u,Ts,lr,L)
    % Linearisierung des diskreten Euler-Modells:
    % x+ = f(x,u)  =>  x+ ≈ A x + B u + c

    X   = x(1); Y   = x(2); phi = x(3); v = x(4); del = x(5);
    a   = u(1); om  = u(2);

    k = lr/L;
    tdel  = tan(del);
    beta  = atan(k * tdel);

    % d beta / d del
    % beta = atan(k*tan(del))
    db_dDel = (k * (1/cos(del))^2) / (1 + (k*tdel)^2);

    s = sin(phi + beta);
    cphi = cos(phi + beta);
    sb = sin(beta);
    cb = cos(beta);

    A = eye(5);
    B = zeros(5,2);

    % Xn = X + Ts*(v*cos(phi+beta))
    A(1,3) = Ts*(v * (-s));              % dXn/dphi
    A(1,4) = Ts*(cphi);                  % dXn/dv
    A(1,5) = Ts*(v * (-s) * db_dDel);    % dXn/ddel

    % Yn = Y + Ts*(v*sin(phi+beta))
    A(2,3) = Ts*(v * (cphi));            % dYn/dphi
    A(2,4) = Ts*(s);                     % dYn/dv
    A(2,5) = Ts*(v * (cphi) * db_dDel);  % dYn/ddel

    % phin = phi + Ts*(v/L)*cos(beta)*tan(del)
    A(3,4) = Ts*(1/L)*cb*tdel;               % dphin/dv
    A(3,5) = Ts*(v/L)*(-sb*db_dDel*tdel + cb*(1/(cos(del))^2));       % dphin/ddel

    % vn = v + Ts*a
    B(4,1) = Ts;

    % deln = del + Ts*omega
    B(5,2) = Ts;

    % Affiner Anteil c = f(x,u) - A x - B u
    xplus = [ ...
        X   + Ts*(v*cos(phi + beta));
        Y   + Ts*(v*sin(phi + beta));
        phi + Ts*((v/L)*cos(beta)*tan(del));
        v   + Ts*a;
        del + Ts*om ];

    c = xplus - A*x - B*u;
end

function a = wrapToPi_num(a)
    a = mod(a + pi, 2*pi) - pi;
end
