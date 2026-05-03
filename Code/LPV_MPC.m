clear; close all; clc;
yalmip('clear');

% ----------------------------
% 1) Parameter
% ----------------------------
Ts    = 0.1;
N     = 10;
% N_sim = 200;

lf = 1.105;
lr = 1.738;
L  = lf + lr;

% Grenzen
a_min = -5;   a_max = 2.5;
w_min = -0.6; w_max = 0.6;
v_min = 0;    v_max = 10;
d_min = -0.5; d_max = 0.5;

use_rate_constraints = true;
du_min = [-1.0; -0.4];
du_max = [ 1.0;  0.4];

v_ref = 3.0;

% Kosten 
Q_pos = 20;
Q_phi = 0.1;
Q_vel = 1.0;
Q_del = 0.5;
R_u   = 0.05;
R_du  = 0.2;

P_pos = Q_pos;
P_phi = Q_phi;
P_vel = Q_vel;
P_del = Q_del;

nx = 5;  % [X; Y; phi; v; delta]
nu = 2;  % [a; omega]

% ----------------------------
% 2) Referenz: voller Kreis
% ----------------------------
% time = (0:N_sim+N)*Ts;
% Tsim = (N_sim+N)*Ts;
% omega_path = 2*pi / Tsim;
% theta = omega_path * time;
% 
R = 10;
% ref_X = R * cos(theta);
% ref_Y = R * sin(theta);
x_current = [R; 0; pi/2; v_ref; 0];
[ref_X,ref_Y,phi_ref,delta_ref,n_seg]=generateCircularTrajectory(R,v_ref,Ts,L,pi/2,N);
N_sim=n_seg-N;
ref_traj = [ref_X; ref_Y];

% ----------------------------
% 3) YALMIP QP-Setup (einmalig)
%    Wir verwenden ein "frozen" LPV-Modell über den Horizont:
%    A,B,c werden pro Zeitschritt neu berechnet, aber dann als konstant im Horizont genutzt.
% ----------------------------
x = sdpvar(nx, N+1, 'full');
u = sdpvar(nu, N,   'full');

% Parameter fürs QPa
% x0     = sdpvar(nx,1);
 % u_prev = sdpvar(nu,1);
% x_ref  = sdpvar(2, N+1, 'full');   % nur X,Y Referenz

% Apar = sdpvar(nx,nx,'full');
% Bpar = sdpvar(nx,nu,'full');

u_prev = [0;0];


%x_current = [R; 0; pi/2; v_ref; 0];


X_log = zeros(nx, N_sim+1);
U_log = zeros(nu, N_sim);
X_log(:,1) = x_current;
xt = x_current;
Bpar = zeros(5,2);
Bpar(4,1) = 1;
Bpar(5,2) = 1;

Bpard = Ts * Bpar;

for t = 1:N_sim
    phit = xt(3);
    delt = xt(5);
    betat = atan((lr/L)*tan(delt));
    pt = [cos(phit+betat); sin(phit+betat); ((cos(betat)*tan(delt))/L)];
    constraints = [x(:,1) == xt];
    objective   = 0;
    Apar = zeros(5,5);
    Apar(1,4) = pt(1);
    Apar(2,4) = pt(2);
    Apar(3,4) = pt(3);
    
    Apard{1} = eye(5) + Ts * Apar;

    if t == 1 
        for j = 2:N
            Apard{j} = Apard{1};
        end
    else
        for j = 2:N
        
            xtn = s{j};
            phit = xtn(3);
            delt = xtn(5);
            betat = atan((lr/L)*tan(delt));
            pt = [cos(phit+betat); sin(phit+betat); ((cos(betat)*tan(delt))/L)];
            Apar = zeros(5,5);
            Apar(1,4) = pt(1);
            Apar(2,4) = pt(2);
            Apar(3,4) = pt(3);
        
            Apard{j} = eye(5) + Ts * Apar;
        end
    end

    for k = 1:N
        % Lineare affine Dynamik
        constraints = [constraints, x(:,k+1) == Apard{k}*x(:,k) + Bpard*u(:,k)];
    
        % Bounds (Inputs)
        constraints = [constraints, a_min <= u(1,k) <= a_max];
        constraints = [constraints, w_min <= u(2,k) <= w_max];
    
        % Bounds (States)
        constraints = [constraints, v_min <= x(4,k+1) <= v_max];
        constraints = [constraints, d_min <= x(5,k+1) <= d_max];
    
        % Delta-u
        if k == 1
            du = u(:,k) - u_prev;
        else
            du = u(:,k) - u(:,k-1);
        end
        if use_rate_constraints
            constraints = [constraints, du_min <= du <= du_max];
        end
    
        % Kosten
        e_pos = x(1:2,k) - ref_traj(:,t+k-1); 
        e_phi = x(3,k) - phi_ref(:,t+k-1);         
        e_vel = x(4,k) - v_ref;
        e_del = x(5,k) - delta_ref;
    
        objective = objective + ...
            Q_pos*(e_pos'*e_pos) + ...
            Q_phi*(e_phi'*e_phi) + ...
            Q_vel*(e_vel'*e_vel) + ...
            Q_del*(e_del'*e_del) + ...
            R_u*(u(:,k)'*u(:,k)) + ...
            R_du*(du'*du);
    end

    % Terminalkosten
    e_posN = x(1:2,N+1) - ref_traj(:,t+N); 
    e_phiN = x(3,N+1)-phi_ref(:,t+N);
    e_velN = x(4,N+1) - v_ref;
    e_delN = x(5,N+1)-delta_ref;
    
    objective = objective + ...
        P_pos*(e_posN'*e_posN) + ...
        P_phi*(e_phiN'*e_phiN) + ...
        P_vel*(e_velN'*e_velN) + ...
        P_del*(e_delN'*e_delN);
    
    % QP Solver
    opts = sdpsettings('verbose',0,'solver','quadprog');
    
    % Optimizer: Ausgabe ist u(:,1)
    diagnostics = optimize(constraints, objective, opts);
    
    if diagnostics.problem ~= 0
        % Fallback: sichere Eingabe (z.B. 0)
        warning("QP infeasible/failed at t=%d: %s", t, yalmiperror(diagnostics.problem));
        ut = u_prev;
    else
        ut = value(u(:,1));
        for j = 1:N
            s{j}= value(x(:,j+1));
            if j == 1
                s_NL{j} = f_bicycle_euler(xt, value(u(:,j)), Ts, lr, L);
            else
                s_NL{j} = f_bicycle_euler(s_NL{j-1}, value(u(:,j)), Ts, lr, L);
            end
        end
    end     
    if t == 100
        disp('Erg')
    end
    u_prev = ut;

% ----------------------------
% 4) Simulation
% ----------------------------

    % Horizon-Referenz
    % ref_h = zeros(2, N+1);
    % for kk = 0:N
    %     idx = min(t+kk, N_sim+1);
    %     ref_h(:,kk+1) = ref_traj(:,idx);
    % end

    % --- LPV Schritt: Linearisieren um (x_current, u_last)
    % [A, B] = lin_bicycle_euler(xt, ut, Ts, lr, L);

    % % QP lösen
    % in = [xt; ut; ref_h(:); A(:); B(:)];
    % [u0, err] = mpc_qp(in);
    % if err ~= 0 || isempty(u0)
    %     u0 = zeros(nu,1);
    % end

    % --- "echte" nichtlineare Pflanze weiter simulieren
    x_next = f_bicycle_euler(xt, ut, Ts, lr, L);

    X_log(:,t+1) = x_next;
    U_log(:,t)       = ut;
    xt    = x_next;

    % Fortschritt anzeigen
    if mod(t, 1) == 0
        fprintf('Simulationsschritt %d/%d\n', t, N_sim);
    end
end

% Plot
figure('Position',[120,120,900,520]);
plot(ref_traj(1,:), ref_traj(2,:), 'r-', 'LineWidth', 2); hold on;
plot(X_log(1,:), X_log(2,:), 'b--', 'LineWidth', 2.5);
grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]');
legend('Referenz (voller Kreis)','LPV-MPC (QP, linearisiert)','Location','best');
title(sprintf('LPV-MPC Einspurmodell: Kreisfolge, N=%d, Ts=%.2f', N, Ts));

% ----------------------------
% 5) Hilfsfunktionen
% ----------------------------
function x_next = f_bicycle_euler(x, u, Ts, lr, L)
    X   = x(1); Y   = x(2); phi = x(3); v = x(4); del = x(5);
    a   = u(1); om  = u(2);

    k = lr/L;
    beta  = atan(k * tan(del));

    x_next = [ ...
        X   + Ts*(v*cos(phi + beta));
        Y   + Ts*(v*sin(phi + beta));
        phi + Ts*((v/L)*cos(beta)*tan(del));
        v   + Ts*a;
        del + Ts*om ];
end

function [A,B] = lin_bicycle_euler(x,u,Ts,lr,L)
    % Linearisierung des diskreten Euler-Modells:
    % x+ = f(x,u)  =>  x+ ≈ A x + B u

    X   = x(1); Y   = x(2); phi = x(3); v = x(4); del = x(5);
    a   = u(1); om  = u(2);

    k = lr/L;
    tdel  = tan(del);
    beta  = atan(k * tdel);

   

    s = sin(phi + beta);
    cphi = cos(phi + beta);
    sb = sin(beta);
    cb = cos(beta);

    A = eye(5);
    B = zeros(5,2);

    % Xn = X + Ts*(v*cos(phi+beta))
   
    A(1,4) = Ts*(cphi);                  % dXn/dv
   

    % Yn = Y + Ts*(v*sin(phi+beta))
   
    A(2,4) = Ts*(s);                     % dYn/dv
   

    % phin = phi + Ts*(v/L)*cos(beta)*tan(del)
    A(3,4) = Ts*(1/L)*cb*tdel;               % dphin/dv
    

    % vn = v + Ts*a
    B(4,1) = Ts;

    % deln = del + Ts*omega
    B(5,2) = Ts;

    
  
end

function [X_ref,Y_ref,phi_ref,delta_ref,n_seg]=generateCircularTrajectory(R,v_ref,Ts,L,phi0,N)
% traj: Struct with X_ref, Y_ref, phi_ref, v_ref, delta_ref, t
% Number of segments
circum=2*pi*R;
n_seg=ceil(circum/(v_ref*Ts))+N;
% angles along circle
theta=linspace(0,2*pi,n_seg+1); 
% Global positions
X_ref=R*cos(theta);
Y_ref=R*sin(theta);
% Heading (tangent to circle, ccw)
phi_ref=theta+phi0; 
% Steering delta=atan(L/R)
delta_ref=atan(L/R);
end