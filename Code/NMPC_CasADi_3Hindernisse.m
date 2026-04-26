clear; close all; clc;

% CasADi-Pfad (an Ihre Installation anpassen)
addpath('C:\Users\sohai\Downloads\Bachelorarbeit\New folder\casadi-3.7.2-windows64-matlab2018b');
import casadi.*

% -------------------------------------------------------------------------
% NMPC auf bildbasierter Strecke mit stabiler Original-MPC-Vorschau
% und automatisch berechnetem Geschwindigkeitsprofil.
%
% Ziel:
% - Verhalten und Kurvenanfang bleiben nah am Originalcode
% - auf Geraden darf v_ref bis zur harten Grenze v_max = 30 m/s gehen
% - in Kurven wird v_ref automatisch aus der Kruemmung kappa berechnet
% - keine feste Gerade-Sollgeschwindigkeit und keine feste Kurven-Sollgeschwindigkeit
%
% Idee:
% 1) Pfadgeometrie ist fest und unabhaengig von v_ref.
% 2) phi_ref / delta_ref werden aus Tangente + Kruemmung berechnet.
% 3) Die Original-Form mit kappa_low/kappa_high und curve_factor bleibt erhalten,
%    damit Gerade/Kurve an denselben Stellen wie im Original beginnen.
% 4) Die obere Gerade-Geschwindigkeit kommt direkt aus v_max.
% 5) Die Kurvengeschwindigkeit wird automatisch aus der staerksten geglaetteten
%    Kruemmung der Strecke und ay_ref_max berechnet.
% 6) Der MPC folgt wie im Original der lokalen Referenz
%    [X_ref, Y_ref, phi_ref, delta_ref, v_ref].
% -------------------------------------------------------------------------

% ----------------------------
% 1) Parameter
% ----------------------------
Ts = 0.08;                   % Abtastzeit [s]
N  = 35;                     % Praediktionshorizont (~2.8 s Vorschau)

% Fahrzeuggeometrie
lf = 1.105;
lr = 1.738;
L  = lf + lr;

% Eingangsgrenzen
a_min = -7.0;   a_max = 10.0;
w_min = -0.70;  w_max = 0.70;

% Zustandsgrenzen
v_min = 0.0;    v_max = 30.0;
d_min = -0.45;  d_max = 0.45;

% Ratenbeschraenkungen
use_rate_constraints = true;
du_min = [-2.5; -0.50];
du_max = [ 2.5;  0.50];

% Geschwindigkeitsprofil entlang des Pfades
% Keine feste Gerade-/Kurven-Sollgeschwindigkeit:
% - oben begrenzt nur durch v_max = 30 m/s
% - in Kurven automatisch durch kappa und ay_ref_max
v_start        = 10.0;       % Startgeschwindigkeit [m/s] (Fallback, Initialwert wird unten aus der Referenz gesetzt)
num_laps       = 2;          % Anzahl der zu fahrenden Runden
v_global_min   = 4.0;        % numerische Untergrenze fuer v_ref [m/s]
ay_ref_max     = 6.5;        % Querbeschleunigungsgrenze fuer die automatische Kurven-v_ref [m/s^2]
a_ref_accel    = 10.0;       % zulaessige Beschleunigung des Referenzprofils [m/s^2]
a_ref_brake    = 10.0;       % zulaessige Verzoegerung des Referenzprofils [m/s^2]

% Feste geometrische Aufloesung der Strecke
% Kleiner waehlen -> feinere Referenz, aber mehr Punkte
ds_path = 0.50;              % [m]

% Gewichte der Kostenfunktion
Q_pos = 650;       % deutlich staerkeres Tracking der lokalen MPC-Referenz
Q_phi = 45.0;      % Richtung soll der Vorschau folgen
Q_vel = 3.0;       % Ist-Geschwindigkeit folgt der v_ref, Position bleibt wichtiger
Q_del = 1.0;
R_u   = 0.002;
R_du  = 0.06;      % glatte Eingriffe, aber agil genug zum Folgen der Vorschau

% Hindernisvermeidung (weich, wie in der Datei "Hindernisse")
% Wichtig: Q_obs muss hier groesser sein als in der langsamen Hindernis-Datei,
% weil Q_pos und die Geschwindigkeit deutlich hoeher sind.
Q_obs = 4500;

% Zusaetzliche Fuehrungsbegrenzung gegen grosse Schleifen:
% Der NMPC darf um Hindernisse ausweichen, soll aber in einem Korridor um
% die lokale Referenz bleiben.
corridor_radius = 5.5;   % [m], enger Korridor: Fahrzeug bleibt naeher an der MPC-Vorschau

% Terminalgewichte
P_pos = 1.2 * Q_pos;
P_phi = 1.5 * Q_phi;
P_vel = 1.0 * Q_vel;
P_del = 1.2 * Q_del;

% Dimensionen
nx = 5;   % [X; Y; phi; v; delta]
nu = 2;   % [a; omega]

% ----------------------------
% 2) Referenzgeometrie + variable Geschwindigkeit
% ----------------------------
track_width_m  = 90;
track_height_m = 110;

% Erst genau eine geschlossene Runde erzeugen
[X_lap, Y_lap, phi_lap, delta_lap, s_lap, kappa_lap, v_ref_lap] = ...
    generateImageTrackReferenceVariableSpeedClosed( ...
        track_width_m, track_height_m, ds_path, lf, lr, ...
        v_max, v_global_min, ay_ref_max, a_ref_accel, a_ref_brake);

lap_length = s_lap(end) + ds_path;
fprintf('Rundenlaenge: %.2f m\n', lap_length);

% Mehrere Runden hintereinander aufbauen
[X_path, Y_path, phi_path, delta_path, s_path, kappa_path, v_ref_path] = ...
    repeatReferenceForLaps(X_lap, Y_lap, phi_lap, delta_lap, s_lap, kappa_lap, v_ref_lap, num_laps, ds_path);

total_length = num_laps * lap_length;
fprintf('Gesamtlaenge fuer %d Runden: %.2f m\n', num_laps, total_length);

% Nur fuer die Darstellung: Start- und Endpunkt explizit verbinden
[X_lap_plot,  Y_lap_plot]  = closeCurveForPlot(X_lap,  Y_lap);
[X_path_plot, Y_path_plot] = closeCurveForPlot(X_path, Y_path);

% ----------------------------
% 2b) Drei Hindernisse auf MEINER Referenz
% ----------------------------
% Die Hindernisse liegen direkt auf der Mittellinie. Der NMPC muss daher kurz
% ausweichen und danach wegen der Tracking-Kosten wieder zur Referenz zurueck.
% Die Positionen werden aus einer Runde der eigenen Bild-Referenz gebildet.
obs_idx_lap = unique(round([0.22, 0.53, 0.79] * numel(X_lap)));
obs_idx_lap = max(1, min(numel(X_lap), obs_idx_lap(:)));
n_obs = numel(obs_idx_lap);

obs_x = X_lap(obs_idx_lap);
obs_y = Y_lap(obs_idx_lap);
obs_pos = [obs_x(:), obs_y(:)];  % erzwingt n_obs x 2
obs_rad = 1.60 * ones(n_obs, 1);      % realer Hindernisradius [m]
obs_safe_margin = 1.60;               % Sicherheitsabstand um das Hindernis [m]
obs_influence_margin = 2.80;          % kompakter Einflussbereich, damit die Trajektorie nicht grosse Schleifen macht
obs_safe = obs_rad + obs_safe_margin;
obs_infl = obs_safe + obs_influence_margin;

% Lokale Ausweichreferenz:
% Die eigentliche rote Referenz bleibt unveraendert. Nur im MPC-Horizont wird
% die Zielspur in der Naehe eines Hindernisses seitlich verschoben und danach
% automatisch wieder auf die Originalreferenz zurueckgefuehrt.
avoid_ref_half_width_m = 14.0;        % kuerzere Ausweichblase: weniger grosse Schleifen
avoid_ref_offset_extra = 0.85;        % zusaetzlicher Abstand ausserhalb des Sicherheitsradius [m]
avoid_ref_side = [1; -1; -1];         % +1 = linke Normalenrichtung, -1 = rechte Normalenrichtung
avoid_ref_side = avoid_ref_side(1:n_obs);
avoid_v_factor = 0.55;                % am Hindernis staerker bremsen, damit die Vorschau fahrbar bleibt
avoid_v_min = 4.0;                    % Untergrenze fuer lokale Hindernis-v_ref [m/s]

% Aktivierung pro Runde: physisch sind es nur 3 Hindernisse, aber der
% Referenzindex wiederholt sich bei mehreren Runden.
n_lap_pts = numel(X_lap);
obs_idx_path = zeros(n_obs, num_laps);
for lap = 1:num_laps
    obs_idx_path(:,lap) = obs_idx_lap + (lap-1) * n_lap_pts;
end

% Aktivierungsfenster in Referenz-Indizes. Bei 25 m/s springt die Vorschau
% mehrere ds_path-Punkte pro MPC-Schritt, deshalb ist das Fenster groesser
% als nur N.
obs_window_before = ceil(10 / ds_path);
obs_window_after  = ceil((v_max * Ts * N + 30) / ds_path);

% Simulationslaenge grob bis zum Ende der letzten Runde
v_mean_est = max(mean(v_ref_lap), 0.5);
alloc_steps = ceil(total_length / v_mean_est / Ts) + 120;

% ----------------------------
% 3) Modell mit CasADi definieren
% ----------------------------
x_sym = MX.sym('x', nx);
u_sym = MX.sym('u', nu);

X   = x_sym(1);
Y   = x_sym(2);
phi = x_sym(3);
v   = x_sym(4);
del = x_sym(5);

a  = u_sym(1);
om = u_sym(2);

beta = atan(lr/L * tan(del));

x_next = [ X + Ts * (v * cos(phi + beta));
           Y + Ts * (v * sin(phi + beta));
           phi + Ts * (v/L * cos(beta) * tan(del));
           v + Ts * a;
           del + Ts * om ];

f_dyn = Function('f_dyn', {x_sym, u_sym}, {x_next});

% ----------------------------
% 4) NLP fuer einen MPC-Schritt vorbereiten
% ----------------------------
x0     = MX.sym('x0', nx);
u_prev = MX.sym('u_prev', nu);

% Referenz ueber Horizont: [X; Y; phi; delta; v]
ref_horizon = MX.sym('ref', 5, N+1);

% Hindernisparameter fuer diesen MPC-Schritt:
% Zeilen = [x; y; Sicherheitsradius; Einflussradius]
obs_param = MX.sym('obs_param', 4, n_obs);

U = MX.sym('U', nu, N);
X = MX.sym('X', nx, N+1);

w  = [U(:); X(:)];
nv = numel(w);

g = {};
lbg = [];
ubg = [];

% Anfangsbedingung
g{end+1} = X(:,1) - x0;
lbg = [lbg; zeros(nx,1)];
ubg = [ubg; zeros(nx,1)];

J = 0;

for k = 1:N
    % Eingangsgrenzen
    g{end+1} = U(:,k);
    lbg = [lbg; a_min; w_min];
    ubg = [ubg; a_max; w_max];

    % Ratenbeschraenkungen
    if k == 1
        du = U(:,k) - u_prev;
    else
        du = U(:,k) - U(:,k-1);
    end

    if use_rate_constraints
        g{end+1} = du;
        lbg = [lbg; du_min];
        ubg = [ubg; du_max];
    end

    % Zustandsgrenzen
    vk = X(4,k);
    dk = X(5,k);
    g{end+1} = vk;
    g{end+1} = dk;
    lbg = [lbg; v_min; d_min];
    ubg = [ubg; v_max; d_max];

    % Dynamik
    x_next_calc = f_dyn(X(:,k), U(:,k));
    g{end+1} = X(:,k+1) - x_next_calc;
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % Kosten
    ref_col = ref_horizon(:,k);
    e_pos = X(1:2,k) - ref_col(1:2);
    e_phi = wrapToPi_local(X(3,k) - ref_col(3));
    e_vel = X(4,k) - ref_col(5);
    e_del = X(5,k) - ref_col(4);

    J = J + Q_pos * (e_pos' * e_pos) + ...
            Q_phi * e_phi^2 + ...
            Q_vel * e_vel^2 + ...
            Q_del * e_del^2 + ...
            R_u   * (U(:,k)' * U(:,k)) + ...
            R_du  * (du' * du);

    % Korridor um die lokale Referenz: verhindert falsche grosse Schleifen.
    g{end+1} = e_pos' * e_pos;
    lbg = [lbg; 0];
    ubg = [ubg; corridor_radius^2];

    % Weiche Hindernisvermeidung:
    % - im Einflussradius wird frueh seitlich ausgewichen
    % - im Sicherheitsradius ist die Strafe deutlich hoeher
    for j = 1:n_obs
        ox    = obs_param(1,j);
        oy    = obs_param(2,j);
        osafe = obs_param(3,j);
        oinfl = obs_param(4,j);

        dx_obs = X(1,k) - ox;
        dy_obs = X(2,k) - oy;
        d2_obs = dx_obs^2 + dy_obs^2;

        pen_soft = fmax(0, oinfl^2 - d2_obs);
        pen_safe = fmax(0, osafe^2 - d2_obs);

        J = J + Q_obs * (pen_soft / fmax(oinfl^2, 1e-6))^2 + ...
                8 * Q_obs * (pen_safe / fmax(osafe^2, 1e-6))^2;

        % Harte Sicherheitsbedingung. Bei inaktiven Hindernissen ist osafe=0
        % und das Hindernis liegt weit weg, daher ist die Bedingung automatisch erfuellt.
        g{end+1} = d2_obs - osafe^2;
        lbg = [lbg; 0];
        ubg = [ubg; inf];
    end
end

% Terminalkosten
ref_term = ref_horizon(:,N+1);
e_posN = X(1:2,N+1) - ref_term(1:2);
e_phiN = wrapToPi_local(X(3,N+1) - ref_term(3));
e_velN = X(4,N+1) - ref_term(5);
e_delN = X(5,N+1) - ref_term(4);

J = J + P_pos * (e_posN' * e_posN) + ...
        P_phi * e_phiN^2 + ...
        P_vel * e_velN^2 + ...
        P_del * e_delN^2;

% Terminal ebenfalls im Korridor halten.
g{end+1} = e_posN' * e_posN;
lbg = [lbg; 0];
ubg = [ubg; corridor_radius^2];

% Terminale Hindernisstrafe, damit das Ende der Vorschau nicht ins Hindernis faehrt
for j = 1:n_obs
    oxN    = obs_param(1,j);
    oyN    = obs_param(2,j);
    osafeN = obs_param(3,j);
    oinflN = obs_param(4,j);

    dxN_obs = X(1,N+1) - oxN;
    dyN_obs = X(2,N+1) - oyN;
    d2N_obs = dxN_obs^2 + dyN_obs^2;

    pen_softN = fmax(0, oinflN^2 - d2N_obs);
    pen_safeN = fmax(0, osafeN^2 - d2N_obs);

    J = J + Q_obs * (pen_softN / fmax(oinflN^2, 1e-6))^2 + ...
            8 * Q_obs * (pen_safeN / fmax(osafeN^2, 1e-6))^2;

    g{end+1} = d2N_obs - osafeN^2;
    lbg = [lbg; 0];
    ubg = [ubg; inf];
end

nlp = struct('f', J, 'x', w, 'g', vertcat(g{:}), ...
             'p', vertcat(x0, u_prev, ref_horizon(:), obs_param(:)));

opts = struct('ipopt', struct('print_level', 0, 'tol', 1e-6, 'max_iter', 350), ...
              'print_time', false);

solver = nlpsol('solver', 'ipopt', nlp, opts);

% ----------------------------
% 5) Simulationsschleife
% ----------------------------
x_current = [X_path(1); Y_path(1); phi_path(1); v_ref_path(1); delta_path(1)];
u_last    = [0; 0];

X_log = zeros(nx, alloc_steps+1);
U_log = zeros(nu, alloc_steps);
idx_log = zeros(1, alloc_steps+1);
vref_log = zeros(1, alloc_steps+1);

% MPC-Vorschau-Logs:
% Xpred_log speichert die optimierte Zustandstrajektorie des Solvers.
% Diese Vorschau ist die Linie, der das Fahrzeug im Receding-Horizon-Betrieb folgt.
Xpred_log = cell(1, alloc_steps);
Refprev_log = cell(1, alloc_steps);

% Warmstart mit verschobener alter Loesung, damit die Vorschau von Schritt
% zu Schritt nicht springt.
U_prev_opt = [];
X_prev_opt = [];

X_log(:,1) = x_current;
idx_log(1) = 1;
vref_log(1) = v_ref_path(1);

path_idx = 1;
t = 0;

while true
    t = t + 1;

    % Speicher automatisch erweitern, falls die Simulation länger dauert
    if t > size(U_log,2)
        oldCap = size(U_log,2);
        newCap = oldCap + max(500, oldCap);

        X_log(:, newCap+1) = 0;
        U_log(:, newCap) = 0;
        idx_log(1, newCap+1) = 0;
        vref_log(1, newCap+1) = 0;

        Xpred_log{newCap} = [];
        Refprev_log{newCap} = [];
    end
    % Fortschritt stabil weiterfuehren. Bei Hindernis-Ausweichmanoevern ist
    % reine Naechster-Punkt-Suche instabil, weil das Fahrzeug absichtlich
    % seitlich von der Mittellinie faehrt.
    path_idx = updateProgressIndexStable(X_path, Y_path, x_current(1), x_current(2), ...
        path_idx, x_current(4), v_ref_path(path_idx), Ts, ds_path);
    idx_log(t) = path_idx;
    vref_log(t) = v_ref_path(path_idx);

    % Horizon ueber lokale variable Geschwindigkeitsreferenz aufbauen
    idx = buildPreviewIndices(path_idx, v_ref_path, x_current(4), Ts, ds_path, N);

    ref_mat_nominal = [X_path(idx);
                       Y_path(idx);
                       phi_path(idx);
                       delta_path(idx);
                       v_ref_path(idx)];

    % Lokale Ausweichreferenz erzeugen: weg vom Hindernis und danach zurueck
    % auf die Originalreferenz. Dadurch bleibt die Trajektoriefolge sauber.
    ref_mat = buildAvoidanceReference(ref_mat_nominal, idx, obs_idx_path, obs_safe, ...
        avoid_ref_half_width_m, avoid_ref_offset_extra, avoid_ref_side, ...
        avoid_v_factor, avoid_v_min, ds_path);

    % Fuer den Geschwindigkeitsplot die wirklich im MPC verwendete lokale v_ref speichern.
    vref_log(t) = ref_mat(5,1);
    Refprev_log{t} = ref_mat;

    % Nur Hindernisse aktivieren, die im aktuellen Vorschaufenster liegen.
    % Sonst werden sie weit weg gesetzt und haben keinen Einfluss auf den Solver.
    obs_param_num = zeros(4, n_obs);
    for j = 1:n_obs
        is_active = any(obs_idx_path(j,:) >= (path_idx - obs_window_before) & ...
                        obs_idx_path(j,:) <= (path_idx + obs_window_after));
        if is_active
            obs_param_num(:,j) = [obs_pos(j,1); obs_pos(j,2); obs_safe(j); obs_infl(j)];
        else
            obs_param_num(:,j) = [1e4; 1e4; 0; 0];
        end
    end

    p = [x_current; u_last; ref_mat(:); obs_param_num(:)];

    % Warmstart:
    % Ab Schritt 2 wird die alte optimale Loesung um einen Zeitschritt
    % nach vorne geschoben. Dadurch wird die MPC-Vorschau ruhiger und die
    % Fahrzeugtrajektorie folgt der Vorschau statt bei jedem Schritt neu zu
    % springen.
    if isempty(U_prev_opt) || isempty(X_prev_opt)
        U_guess = repmat(u_last, 1, N);
    else
        U_guess = [U_prev_opt(:,2:end), U_prev_opt(:,end)];
    end

    X_guess = zeros(nx, N+1);
    X_guess(:,1) = x_current;
    for k = 1:N
        X_guess(:,k+1) = full(f_dyn(X_guess(:,k), U_guess(:,k)));
    end
    w0 = [U_guess(:); X_guess(:)];

    sol = solver('x0', w0, 'p', p, 'lbg', lbg, 'ubg', ubg);

    if ~isfield(sol, 'x') || isempty(sol.x)
        warning('NLP-Fehler bei t=%d', t);
        u_opt = u_last;
        X_opt = X_guess;
        U_opt = U_guess;
    else
        w_opt = full(sol.x);
        U_opt = reshape(w_opt(1:nu*N), nu, N);
        X_opt = reshape(w_opt(nu*N+1:end), nx, N+1);
        u_opt = U_opt(:,1);
    end

    % Wichtige Aenderung:
    % Das Fahrzeug uebernimmt den ersten Schritt der optimierten MPC-
    % Zustandstrajektorie. Damit liegt X_log(:,t+1) exakt auf der Vorschau
    % X_opt(:,2), anstatt nur erneut mit u_opt simuliert zu werden.
    Xpred_log{t} = X_opt;
    U_prev_opt = U_opt;
    X_prev_opt = X_opt;
    x_next = X_opt(:,2);

    X_log(:,t+1) = x_next;
    U_log(:,t)   = u_opt;
    x_current    = x_next;
    u_last       = u_opt;

    if mod(t, 10) == 0
        fprintf('Schritt %d | Pfadindex %d/%d | v=%.2f | v_ref=%.2f m/s\n', ...
            t, path_idx, numel(X_path), x_current(4), v_ref_path(path_idx));
    end

    % Abbruch, wenn die letzte Runde praktisch beendet ist
    if path_idx >= numel(X_path) - 6
        finished = true;
        k_end = t;
        X_log = X_log(:,1:t+1);
        U_log = U_log(:,1:t);
        idx_log = idx_log(1:t+1);
        vref_log = vref_log(1:t+1);
        Xpred_log = Xpred_log(1:t);
        Refprev_log = Refprev_log(1:t);
        break;
    end
end


% letzten Log-Eintrag sauber setzen
idx_log(end) = min(updateProgressIndexStable(X_path, Y_path, X_log(1,end), X_log(2,end), ...
    idx_log(max(1,end-1)), X_log(4,end), v_ref_path(idx_log(max(1,end-1))), Ts, ds_path), numel(X_path));
vref_log(end) = v_ref_path(idx_log(end));

% ----------------------------
% 6) Plots
% ----------------------------
figure('Position',[120,120,900,520]);
plot(X_path_plot, Y_path_plot, 'r-', 'LineWidth', 2.2, 'DisplayName', 'Referenz (geschlossen)'); hold on;
plot(X_log(1,:), X_log(2,:), 'b--', 'LineWidth', 2.2, 'DisplayName', 'Trajektorie');
plotObstacleSet(obs_pos, obs_rad, obs_safe);
grid on; axis equal;
xlabel('X [m]');
ylabel('Y [m]');
legend('Location','best');
title(sprintf('NMPC Originalstruktur + Auto-v Profil | %d Runden | Gerade bis %.1f m/s | N=%d | Ts=%.2f s', ...
    num_laps, v_max, N, Ts));

figure('Position',[120,120,900,520]);
subplot(3,2,1);
plot(s_lap, v_ref_lap, 'LineWidth', 1.8);
ylabel('v_{ref}(s) [m/s]'); xlabel('Bogenlaenge s [m]');
title('Geschwindigkeitsprofil einer Runde');
grid on;

subplot(3,2,2);
plot(s_lap, kappa_lap, 'LineWidth', 1.6);
ylabel('\kappa(s) [1/m]'); xlabel('Bogenlaenge s [m]');
title('Krummung einer Runde');
grid on;

subplot(3,2,3);
plot(0:size(X_log,2)-1, X_log(4,:), 'LineWidth', 1.8); hold on;
plot(0:size(vref_log,2)-1, vref_log, 'r--', 'LineWidth', 1.6);
ylabel('v [m/s]'); xlabel('Schritt');
legend('Ist', 'Ref', 'Location', 'best');
title('Laengsgeschwindigkeit');
grid on;

subplot(3,2,4);
plot(0:size(X_log,2)-1, X_log(5,:)*180/pi, 'LineWidth', 1.8); hold on;
plot(0:numel(idx_log)-1, delta_path(min(idx_log, numel(delta_path)))*180/pi, 'r--', 'LineWidth', 1.4);
ylabel('\delta [deg]'); xlabel('Schritt');
legend('Ist', 'Ref', 'Location', 'best');
title('Lenkwinkel');
grid on;

subplot(3,2,5);
plot(1:size(U_log,2), U_log(1,:), 'LineWidth', 1.8);
ylabel('a [m/s^2]'); xlabel('Schritt');
title('Beschleunigung');
grid on;

subplot(3,2,6);
plot(1:size(U_log,2), U_log(2,:), 'LineWidth', 1.8);
ylabel('\omega [rad/s]'); xlabel('Schritt');
title('Lenkrate');
grid on;

% ----------------------------
% 7) Animation
% ----------------------------
animate_simulation = true;   % true -> Animation nach der Simulation anzeigen
animation_stride   = 2;      % nur jeden n-ten Simulationsschritt zeichnen
animation_pause    = 0.02;   % Pause zwischen zwei Frames [s]
show_preview       = true;   % optimierte MPC-Vorschau anzeigen
save_video         = true;  % optional als MP4 speichern
video_filename     = 'NMPC_OriginalStabil_wieOriginal_ohneFesteV.mp4';

car_length_vis = L;
car_width_vis  = 1.9;

if animate_simulation
    fig_anim = figure('Position',[120,120,900,520]);
    hold on;
    plot(X_path_plot, Y_path_plot, 'r-', 'LineWidth', 2.2, 'DisplayName', 'Referenz (geschlossen)');
    plotObstacleSet(obs_pos, obs_rad, obs_safe);

    h_traj    = plot(NaN, NaN, 'b--', 'LineWidth', 2.4, 'DisplayName', 'Fahrzeugtrajektorie');
    h_car     = patch(NaN, NaN, [0.1 0.45 0.95], 'FaceAlpha', 0.35, 'EdgeColor', [0 0.2 0.8], 'LineWidth', 1.5);
    h_heading = plot(NaN, NaN, 'k-', 'LineWidth', 2.0);
    h_preview = plot(NaN, NaN, 'mo-', 'LineWidth', 1.3, 'MarkerSize', 4, 'MarkerFaceColor', 'm');
    h_target  = plot(NaN, NaN, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);

    axis equal;
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');
    title('NMPC mit 3 Hindernissen | Fahrzeug folgt der optimierten MPC-Vorschau');
    lgd = legend('Referenz (geschlossen)', 'Hindernisse', 'Sicherheitsradius', ...
             'Fahrzeugtrajektorie', 'Fahrzeug', 'Heading', 'optimierte MPC-Vorschau', 'Aktuelle lokale Referenz');

% Achse kleiner machen, damit rechts Platz für Legende + Infofenster bleibt
ax = gca;
ax.Position = [0.08 0.12 0.58 0.78];

% Legende rechts oben fest positionieren
lgd.Units = 'normalized';
lgd.Position = [0.68 0.68 0.28 0.24];

margin = 6;
xlim([min(X_path)-margin, max(X_path)+margin]);
ylim([min(Y_path)-margin, max(Y_path)+margin]);

% Infofenster rechts unter der Legende
h_info = annotation(fig_anim, 'textbox', [0.68 0.38 0.28 0.25], ...
    'String', '', ...
    'FontSize', 11, ...
    'Interpreter', 'tex', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'LineWidth', 0.8, ...
    'Margin', 8, ...
    'FitBoxToText', 'off');
    if save_video
        vwriter = VideoWriter(video_filename, 'MPEG-4');
        vwriter.FrameRate = round(1 / max(animation_pause, 1e-3));
        open(vwriter);
    end

    n_lap_pts = numel(X_lap);

    for k = 1:animation_stride:size(X_log,2)
        idxk = max(1, min(idx_log(min(k, numel(idx_log))), numel(X_path)));

        set(h_traj, 'XData', X_log(1,1:k), 'YData', X_log(2,1:k));

        ipred = max(1, min(k, numel(Xpred_log)));
        if show_preview && ~isempty(Xpred_log{ipred})
            % Jetzt zeigt die magenta Linie die wirkliche optimierte
            % Zustandstrajektorie des MPC, nicht nur die Soll-Referenzpunkte.
            Xpred_now = Xpred_log{ipred};
            set(h_preview, 'XData', Xpred_now(1,:), 'YData', Xpred_now(2,:));
        end

        if ipred <= numel(Refprev_log) && ~isempty(Refprev_log{ipred})
            ref_now = Refprev_log{ipred};
            set(h_target, 'XData', ref_now(1,1), 'YData', ref_now(2,1));
        else
            ref_now_nominal = [X_path(idxk); Y_path(idxk); phi_path(idxk); delta_path(idxk); v_ref_path(idxk)];
            set(h_target, 'XData', ref_now_nominal(1), 'YData', ref_now_nominal(2));
        end

        [car_x, car_y] = vehiclePatch(X_log(1,k), X_log(2,k), X_log(3,k), car_length_vis, car_width_vis);
        set(h_car, 'XData', car_x, 'YData', car_y);

        nose_x = X_log(1,k) + 0.8 * car_length_vis * cos(X_log(3,k));
        nose_y = X_log(2,k) + 0.8 * car_length_vis * sin(X_log(3,k));
        set(h_heading, 'XData', [X_log(1,k), nose_x], 'YData', [X_log(2,k), nose_y]);

        lap_now = min(num_laps, floor((idxk - 1) / n_lap_pts) + 1);
        info_str = sprintf(['Schritt: %d / %d\n' ...
                            'Runde:   %d / %d\n' ...
                            'v:       %.2f m/s\n' ...
                            'v_ref:   %.2f m/s\n' ...
                            'a:       %.2f m/s^2\n' ...
                            '\\omega: %.2f rad/s'], ...
                            k-1, size(X_log,2)-1, lap_now, num_laps, ...
                            X_log(4,k), v_ref_path(idxk), ...
                            U_log(1, max(1, min(k-1, size(U_log,2)))), ...
                            U_log(2, max(1, min(k-1, size(U_log,2)))));
        set(h_info, 'String', info_str);

        drawnow;

        if save_video
            writeVideo(vwriter, getframe(fig_anim));
        else
            pause(animation_pause);
        end
    end

    if save_video
        close(vwriter);
        fprintf('Animationsvideo gespeichert als: %s\n', video_filename);
    end
end

% ----------------------------
% 8) Hilfsfunktionen
% ----------------------------
function [X_ref, Y_ref, phi_ref, delta_ref, s_ref, kappa, v_ref] = ...
    generateImageTrackReferenceVariableSpeedClosed(track_width_m, track_height_m, ds_path, lf, lr, ...
    v_path_max, v_global_min, ay_ref_max, a_ref_accel, a_ref_brake)

    % Mittellinie aus dem Bild
    pts = [
        0.232225, 63.102989;
        0.265578, 66.247800;
        0.696025, 69.545100;
        2.903800, 73.100846;
        6.200590, 74.458100;
        8.558280, 74.224700;
        14.312400, 69.028900;
        16.284900, 66.165300;
        17.172184, 65.358364;
        19.242176, 63.504928;
        21.157884, 61.469187;
        23.019674, 59.376427;
        24.924244, 57.328190;
        26.834249, 55.285693;
        28.744253, 53.243196;
        30.654258, 51.200699;
        32.558770, 49.152358;
        34.472881, 47.114268;
        36.378255, 45.066842;
        38.305686, 43.043251;
        40.310783, 41.108487;
        42.463569, 39.366301;
        44.821921, 37.965337;
        47.380212, 37.047318;
        50.021694, 36.451009;
        52.709244, 36.180578;
        55.408698, 36.226746;
        58.085749, 36.598728;
        60.709282, 37.281206;
        63.238953, 38.286082;
        65.596846, 39.689878;
        67.754190, 41.424796;
        69.697052, 43.429075;
        71.411446, 45.659094;
        72.848823, 48.101650;
        73.921199, 50.749109;
        74.627386, 53.535190;
        75.061543, 56.385289;
        75.373331, 59.254359;
        75.575030, 62.134505;
        75.667194, 65.020991;
        75.639170, 67.909102;
        75.515787, 70.794447;
        75.303911, 73.673870;
        75.002211, 76.544147;
        74.535579, 79.388221;
        73.768426, 82.155385;
        72.578698, 84.744352;
        71.005514, 87.089043;
        69.133406, 89.168001;
        67.013558, 90.954099;
        64.651729, 92.347809;
        62.092261, 93.259493;
        59.442976, 93.815923;
        56.765324, 94.193085;
        54.072607, 94.411668;
        51.373765, 94.335401;
        48.687781, 94.031759;
        46.059953, 93.386906;
        43.735867, 91.950380;
        42.158680, 89.634325;
        41.567067, 86.829023;
        41.681145, 83.951294;
        42.594623, 81.248287;
        44.145294, 78.888432;
        45.963906, 76.753857;
        47.895877, 74.735285;
        49.816394, 72.704356;
        51.714571, 70.649394;
        53.656654, 68.641803;
        55.537192, 66.568895;
        57.243839, 64.332491;
        58.466598, 61.769008;
        58.867603, 58.926431;
        58.293375, 56.121738;
        56.583817, 53.933744;
        54.068011, 52.950360;
        51.383617, 53.058492;
        48.878481, 54.104635;
        46.701385, 55.808164;
        44.724811, 57.775660;
        42.853601, 59.858711;
        40.969744, 61.928567;
        39.059739, 63.971064;
        37.149734, 66.013560;
        35.239729, 68.056057;
        33.329725, 70.098554;
        31.419720, 72.141050;
        29.509715, 74.183547;
        27.599710, 76.226044;
        25.689705, 78.268541;
        23.779700, 80.311037;
        21.869696, 82.353534;
        19.959691, 84.396031;
        18.049686, 86.438528;
        16.139681, 88.481024;
        14.229676, 90.523521;
        12.285837, 92.529023;
        10.347746, 94.540803;
        8.533362, 96.679538;
        7.015250, 99.062023;
        6.171586, 101.792777;
        6.315168, 104.657651;
        7.541976, 107.197206;
        9.759473, 108.804859;
        12.372863, 109.503184;
        15.060309, 109.791685;
        17.755784, 109.969669;
        20.456637, 110.000000;
        23.157697, 109.979568;
        25.858851, 109.979568;
        28.560006, 109.979568;
        31.261161, 109.979568;
        33.962315, 109.979568;
        36.663470, 109.979568;
        39.364625, 109.979568;
        42.065780, 109.979568;
        44.766934, 109.979568;
        47.468089, 109.979568;
        50.169187, 109.989375;
        52.870227, 109.969672;
        55.569805, 109.874637;
        58.265471, 109.692959;
        60.953735, 109.412426;
        63.625418, 108.990905;
        66.262417, 108.369709;
        68.839137, 107.506246;
        71.335422, 106.405374;
        73.735984, 105.082894;
        76.009712, 103.525489;
        78.143660, 101.756144;
        80.137293, 99.808255;
        81.982596, 97.699911;
        83.665028, 95.441178;
        85.161534, 93.037709;
        86.441433, 90.495085;
        87.503439, 87.840335;
        88.344736, 85.096564;
        88.949343, 82.282226;
        89.361167, 79.427932;
        89.629839, 76.554113;
        89.815450, 73.672472;
        89.934375, 70.786852;
        89.987192, 67.898930;
        90.000000, 65.010458;
        89.989078, 62.121988;
        89.989078, 59.233462;
        89.989078, 56.344935;
        89.989078, 53.456408;
        89.998458, 50.567936;
        89.993630, 47.679439;
        89.947487, 44.791401;
        89.819633, 41.906248;
        89.609571, 39.026588;
        89.312499, 36.155881;
        88.841954, 33.312307;
        88.121171, 30.530155;
        87.102381, 27.856815;
        85.796790, 25.330271;
        84.214864, 22.991134;
        82.390904, 20.862798;
        80.356083, 18.965458;
        78.130981, 17.330895;
        75.738629, 15.993370;
        73.215894, 14.965778;
        70.598502, 14.257389;
        67.931524, 13.802100;
        65.243784, 13.516964;
        62.549283, 13.318257;
        59.849756, 13.220143;
        57.150838, 13.106906;
        54.459606, 12.862932;
        51.789741, 12.449354;
        49.338154, 11.277919;
        47.539433, 9.149445;
        46.524217, 6.478098;
        45.464546, 3.829304;
        43.596436, 1.772652;
        41.103467, 0.705951;
        38.427322, 0.329774;
        35.736253, 0.085068;
        33.036521, 0.001889;
        30.335427, 0.005348;
        27.634326, 0.016751;
        24.933171, 0.016751;
        22.232069, 0.005816;
        19.530970, 0.000000;
        16.831126, 0.080139;
        14.140997, 0.334170;
        11.461451, 0.690748;
        8.914951, 1.615554;
        6.954690, 3.554143;
        6.137171, 6.283223;
        6.360035, 9.144672;
        7.675089, 11.629023;
        9.960981, 13.129594;
        12.580517, 13.815852;
        15.235288, 14.345134;
        17.894256, 14.851442;
        20.476717, 15.670159;
        22.581016, 17.431661;
        23.549242, 20.097339;
        23.505528, 22.971259;
        22.435380, 25.590205;
        20.293483, 27.310989;
        17.732505, 28.211289;
        15.090657, 28.812711;
        12.460682, 29.468584;
        10.024226, 30.699054;
        7.818264, 32.363248;
        5.658884, 34.097321;
        3.382477, 35.651813;
        1.263287, 37.427146;
        0.129048, 40.012457;
        0.000000, 42.889244;
        0.129417, 45.773798;
        0.132501, 48.662302;
        0.132501, 51.550828;
        0.132501, 54.439355;
        0.119597, 57.327821;
        0.141668, 60.216208;
        0.232225, 63.102989;
    ];

    % Auf gewuenschte Groesse skalieren
    pts(:,1) = pts(:,1) / max(pts(:,1)) * track_width_m;
    pts(:,2) = pts(:,2) / max(pts(:,2)) * track_height_m;

    % Strecke fuer mehrere Runden explizit schliessen
    if norm(pts(1,:) - pts(end,:)) > 1e-9
        pts = [pts; pts(1,:)];
    end

    % Bogenlaenge der Stuetzstellen
    ds_pts = sqrt(sum(diff(pts,1,1).^2, 2));
    s_pts  = [0; cumsum(ds_pts)];

    % Feste Wegaufloesung
    s_ref = 0:ds_path:s_pts(end);
    if s_ref(end) < s_pts(end)
        s_ref = [s_ref, s_pts(end)];
    end

    X_ref = interp1(s_pts, pts(:,1), s_ref, 'pchip');
    Y_ref = interp1(s_pts, pts(:,2), s_ref, 'pchip');

    % Doppelte Endprobe entfernen, damit mehrere Runden sauber aneinanderpassen
    X_ref = X_ref(1:end-1);
    Y_ref = Y_ref(1:end-1);
    s_ref = s_ref(1:end-1);

    % Periodische Ableitungen bezueglich der Bogenlaenge
    dX  = (circshift(X_ref, -1) - circshift(X_ref, 1)) / (2*ds_path);
    dY  = (circshift(Y_ref, -1) - circshift(Y_ref, 1)) / (2*ds_path);
    ddX = (circshift(X_ref, -1) - 2*X_ref + circshift(X_ref, 1)) / (ds_path^2);
    ddY = (circshift(Y_ref, -1) - 2*Y_ref + circshift(Y_ref, 1)) / (ds_path^2);

    % Richtung der Geschwindigkeit / Tangente
    theta_v = unwrap(atan2(dY, dX));

    % Kruemmung aus der Bahngeometrie
    denom = (dX.^2 + dY.^2).^(3/2) + 1e-8;
    kappa = (dX .* ddY - dY .* ddX) ./ denom;

    % Referenz fuer Slip-Modell
    L = lf + lr;
    delta_ref = atan(L * kappa);
    beta_ref  = atan((lr / L) * tan(delta_ref));
    phi_ref   = unwrap(theta_v - beta_ref);

    % Lenkwinkel begrenzen
    delta_ref = max(min(delta_ref, 0.45), -0.45);

    % ----------------------------
    % Automatische Geschwindigkeitsreferenz aus der Kruemmung
    % ----------------------------
    % Diese Version ist absichtlich wieder sehr nah an der Original-Logik:
    % gleiche Kruemmungsglaettung, gleiche kappa_low/kappa_high-Umschaltung,
    % gleiche Beschleunigungs-/Brems-Glaettung. Nur die beiden festen
    % Geschwindigkeitssollwerte wurden entfernt.

    % Fuer die Geschwindigkeitsplanung nicht die rohe Kruemmung verwenden:
    % pchip + finite Differenzen erzeugen kleine Kruemmungs-Welligkeit, die sonst
    % fast die ganze Strecke faelschlich als Kurve klassifiziert.
    kappa_abs_raw = abs(kappa);
    kappa_abs = periodicMovingMean(kappa_abs_raw, 21);

    % Harte physikalische/geometrische Kurvengrenze:
    %   a_y = v^2 * |kappa|  ->  v = sqrt(a_y_max / |kappa|)
    v_curve_limit = sqrt(ay_ref_max ./ max(kappa_abs, 1e-4));

    % Gleiche Gerade/Kurve-Erkennung wie im Original.
    % Dadurch beginnt die Reduktion an denselben Kruemmungsstellen.
    kappa_low  = 0.006;   % darunter Gerade [1/m]
    kappa_high = 0.035;   % darueber deutliche Kurve [1/m]
    curve_factor = (kappa_abs - kappa_low) ./ max(kappa_high - kappa_low, 1e-9);
    curve_factor = min(1.0, max(0.0, curve_factor));
    curve_factor = curve_factor.^0.70;

    % Automatischer Kurven-Zielwert statt fester Kurven-Geschwindigkeit:
    % Er wird aus der staerksten geglaetteten Kruemmung der Strecke bestimmt.
    % Damit liegt der Kurvenwert bei dieser Strecke sehr nahe am Original,
    % ist aber kein fest eingetragener Geschwindigkeitsparameter mehr.
    kappa_design = max(kappa_abs);
    kappa_design = max(kappa_design, kappa_high);
    v_curve_auto = sqrt(ay_ref_max ./ max(kappa_design, 1e-4));
    v_curve_auto = min(v_path_max, max(v_global_min, v_curve_auto));

    % Automatischer Gerade-Wert: direkt die harte Fahrzeuggrenze v_path_max.
    v_shape = v_path_max - (v_path_max - v_curve_auto) .* curve_factor;

    v_raw = min(v_shape, v_curve_limit);
    v_raw = min(v_path_max, max(v_global_min, v_raw));

    % Geschlossene Strecke: Beschleunigungs-/Bremsgrenzen zyklisch anwenden,
    % damit am Start/Ende kein Sprung entsteht und vor Kurven rechtzeitig
    % gebremst wird.
    v_ref = v_raw;
    n_v = numel(v_ref);
    for pass = 1:3
        for i = 2:n_v
            v_ref(i) = min(v_ref(i), sqrt(max(v_ref(i-1)^2 + 2*a_ref_accel*ds_path, 0)));
        end
        v_ref(1) = min(v_ref(1), sqrt(max(v_ref(end)^2 + 2*a_ref_accel*ds_path, 0)));

        for i = n_v-1:-1:1
            v_ref(i) = min(v_ref(i), sqrt(max(v_ref(i+1)^2 + 2*a_ref_brake*ds_path, 0)));
        end
        v_ref(end) = min(v_ref(end), sqrt(max(v_ref(1)^2 + 2*a_ref_brake*ds_path, 0)));
    end

    % Periodische Glaettung ohne Null-Padding am Rand
    v_ref = periodicMovingMean(v_ref, 11);
    v_ref = min(v_path_max, max(v_global_min, v_ref));
end

function [X_plot, Y_plot] = closeCurveForPlot(X_in, Y_in)
    X_plot = X_in(:).';
    Y_plot = Y_in(:).';

    if isempty(X_plot) || isempty(Y_plot)
        return;
    end

    if hypot(X_plot(end) - X_plot(1), Y_plot(end) - Y_plot(1)) > 1e-9
        X_plot(end+1) = X_plot(1);
        Y_plot(end+1) = Y_plot(1);
    end
end

function [X_rep, Y_rep, phi_rep, delta_rep, s_rep, kappa_rep, v_rep] = ...
    repeatReferenceForLaps(X_lap, Y_lap, phi_lap, delta_lap, s_lap, kappa_lap, v_lap, num_laps, ds_path)

    n = numel(X_lap);
    n_total = n * num_laps + 1;   % +1: finaler Schlusspunkt = Startpunkt

    X_rep = zeros(1, n_total);
    Y_rep = zeros(1, n_total);
    phi_rep = zeros(1, n_total);
    delta_rep = zeros(1, n_total);
    s_rep = zeros(1, n_total);
    kappa_rep = zeros(1, n_total);
    v_rep = zeros(1, n_total);

    for lap = 1:num_laps
        ii = (lap-1)*n + (1:n);
        X_rep(ii) = X_lap;
        Y_rep(ii) = Y_lap;
        phi_rep(ii) = phi_lap;
        delta_rep(ii) = delta_lap;
        kappa_rep(ii) = kappa_lap;
        v_rep(ii) = v_lap;
        s_rep(ii) = s_lap + (lap-1) * (s_lap(end) + ds_path);
    end

    % Start und Ende explizit verbinden, damit die Folge am Schluss sauber schliesst
    X_rep(end) = X_lap(1);
    Y_rep(end) = Y_lap(1);
    phi_rep(end) = phi_lap(1);
    delta_rep(end) = delta_lap(1);
    kappa_rep(end) = kappa_lap(1);
    v_rep(end) = v_lap(1);
    s_rep(end) = num_laps * (s_lap(end) + ds_path);
end

function idx = buildPreviewIndices(path_idx, v_ref_path, v_current, Ts, ds_path, N)
    %#ok<INUSD>  % v_current bleibt nur fuer die alte Funktionssignatur erhalten
    n = numel(v_ref_path);
    idx = zeros(1, N+1);
    idx(1) = path_idx;

    for k = 2:N+1
        % Die Vorschau soll mit der Referenzgeschwindigkeit laufen, nicht mit
        % einer evtl. zu schnellen Ist-Geschwindigkeit. Sonst springt das Ziel
        % nach vorne und der MPC faehrt dauerhaft zu schnell.
        v_step = v_ref_path(idx(k-1));
        ds_adv = max(ds_path, v_step * Ts);
        step_idx = max(1, round(ds_adv / ds_path));
        idx(k) = min(n, idx(k-1) + step_idx);
    end
end

function ref_mat = buildAvoidanceReference(ref_mat, idx, obs_idx_path, obs_safe, ...
    avoid_half_width_m, avoid_offset_extra, avoid_side, avoid_v_factor, avoid_v_min, ds_path)

    if isempty(obs_idx_path)
        return;
    end

    idx = idx(:).';  % als Zeile, damit skalare und Vektoren gleich funktionieren
    half_idx = max(1, round(avoid_half_width_m / ds_path));
    n_obs_local = size(obs_idx_path, 1);

    for k = 1:numel(idx)
        idxk = idx(k);
        offset_xy = [0; 0];
        v_scale = 1.0;

        for j = 1:n_obs_local
            obs_indices_j = obs_idx_path(j,:);
            [dist_idx, pos_min] = min(abs(obs_indices_j - idxk));
            signed_idx = idxk - obs_indices_j(pos_min);

            if dist_idx <= half_idx
                shape = 0.5 * (1 + cos(pi * signed_idx / half_idx));

                % Normale auf die aktuelle Referenzrichtung
                phi_ref_k = ref_mat(3,k);
                normal_left = [-sin(phi_ref_k); cos(phi_ref_k)];

                offset_mag = (obs_safe(j) + avoid_offset_extra) * shape;
                offset_xy = offset_xy + avoid_side(j) * offset_mag * normal_left;

                % Direkt am Hindernis langsamer fahren, danach wieder original v_ref.
                v_scale = min(v_scale, 1 - (1 - avoid_v_factor) * shape);
            end
        end

        ref_mat(1:2,k) = ref_mat(1:2,k) + offset_xy;
        ref_mat(5,k) = max(avoid_v_min, ref_mat(5,k) * v_scale);
    end
end

function idx = updateProgressIndexStable(X_path, Y_path, x, y, idx_prev, v_current, v_ref_current, Ts, ds_path)
    %#ok<INUSD>  % v_ref_current bleibt fuer die alte Funktionssignatur erhalten
    n = numel(X_path);

    idx_near = findProgressIndex(X_path, Y_path, x, y, idx_prev, v_current, Ts, ds_path);

    % Odometriemodell fuer den Fortschritt. Das verhindert, dass der Index beim
    % Ausweichen auf einen falschen nahen Streckenast springt oder haengen bleibt.
    % Nur die Ist-Geschwindigkeit verwenden, damit der Referenzindex am Hindernis
    % nicht kuenstlich zu schnell vorlaeuft.
    v_prog = max(abs(v_current), 1.0);
    step_odom = max(1, round(v_prog * Ts / ds_path));
    idx_odom = min(n, idx_prev + step_odom);

    lower = max(idx_prev, idx_odom - 3);
    upper = min(n, idx_odom + 8);

    idx = min(max(idx_near, lower), upper);
end

function idx = findProgressIndex(X_path, Y_path, x, y, idx_prev, v_current, Ts, ds_path)
    n = numel(X_path);

    if nargin < 6 || isempty(v_current), v_current = 0; end
    if nargin < 7 || isempty(Ts),        Ts = 0.08; end
    if nargin < 8 || isempty(ds_path),   ds_path = 0.50; end

    % Wichtig: Nicht 140 Punkte voraus suchen. Auf einer geschlossenen oder
    % engen Strecke liegen mehrere Streckenabschnitte raeumlich nah beieinander.
    % Eine grosse Suche kann dann auf den falschen Ast springen.
    back_window = 4;
    max_phys_step = ceil((max(v_current, 1.0) * Ts + 2.0) / ds_path);
    front_window = max(12, min(28, max_phys_step));

    i1 = max(1, idx_prev - back_window);
    i2 = min(n, idx_prev + front_window);

    dx = X_path(i1:i2) - x;
    dy = Y_path(i1:i2) - y;
    [~, k] = min(dx.^2 + dy.^2);

    idx_candidate = i1 + k - 1;
    idx = max(idx_candidate, idx_prev);
    idx = min(idx, idx_prev + front_window);
end

function y = periodicMovingMean(x, win)
    x = x(:).';
    n = numel(x);
    win = max(1, round(win));
    if mod(win,2) == 0
        win = win + 1;
    end
    if n == 0 || win <= 1
        y = x;
        return;
    end

    hw = floor(win/2);
    if n <= hw
        y = movmean(x, win);
        return;
    end

    x_pad = [x(end-hw+1:end), x, x(1:hw)];
    kernel = ones(1, win) / win;
    y_pad = conv(x_pad, kernel, 'same');
    y = y_pad(hw+1:hw+n);
end

function ang = wrapToPi_local(ang)
    ang = atan2(sin(ang), cos(ang));
end

function plotObstacleSet(obs_pos, obs_rad, obs_safe)
    th = linspace(0, 2*pi, 120);

    for i = 1:size(obs_pos,1)
        xo = obs_pos(i,1);
        yo = obs_pos(i,2);

        x1 = xo + obs_rad(i) * cos(th);
        y1 = yo + obs_rad(i) * sin(th);
        x2 = xo + obs_safe(i) * cos(th);
        y2 = yo + obs_safe(i) * sin(th);

        if i == 1
            fill(x1, y1, [0.85 0.2 0.2], 'FaceAlpha', 0.35, ...
                 'EdgeColor', 'none', 'DisplayName', 'Hindernisse');
            plot(x2, y2, 'k:', 'LineWidth', 1.3, 'DisplayName', 'Sicherheitsradius');
        else
            fill(x1, y1, [0.85 0.2 0.2], 'FaceAlpha', 0.35, ...
                 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(x2, y2, 'k:', 'LineWidth', 1.3, 'HandleVisibility', 'off');
        end
    end
end

function [xv, yv] = vehiclePatch(xc, yc, phi, Lveh, Wveh)
    body = [ Lveh/2,  Wveh/2;
             Lveh/2, -Wveh/2;
            -Lveh/2, -Wveh/2;
            -Lveh/2,  Wveh/2 ]';

    R = [cos(phi), -sin(phi);
         sin(phi),  cos(phi)];

    body_world = R * body + [xc; yc];
    xv = body_world(1,:);
    yv = body_world(2,:);
end
