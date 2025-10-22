
clear; clc; close all

%% --- Parameters (C-17 @ MTOW, sea-level standard-day) ---
param = struct( ...
    'T0',   12.7e5, ...     % N, 4 × 40,400 lbf (F117-PW-100)
    'g',    9.81,   ...     % m/s^2
    'S',    353.0,  ...     % m^2, wing area = 3800 ft^2
    'rho',  1.225,  ...     % kg/m^3, sea-level
    'mu',   0.020,  ...     % rolling friction (keep model default/tuned)
    'm',    3.2925e5, ...   % kg, MTOW = 585,000 lb
    'vmax', 120     ...     % m/s, keep model default/tuned
);

% Takeoff high-lift (from PDF for takeoff CL; rest kept as in your model)
to = struct( ...
    'CLmax',  3.156, ...    % Takeoff CL at SL (PDF)
    'CLroll', 0.75,  ...    % keep model default/tuned
    'CD0',    0.050, ...    % keep model default/tuned
    'k',      0.045  ...    % keep model default/tuned
);

% Derived reference speeds
W   = param.m * param.g;
Vs  = sqrt( 2*W / (param.rho * param.S * to.CLmax) ); % stall in TO config
Vlof= 1.10 * Vs;                                      % liftoff target
Vr  = 1.05 * Vlof;                                    % rotation start
V2  = 1.20 * Vs;                                      % safety climb speed

%% --- Time/grid & ICs ---
t0 = 0; tf = 160;     % s
h  = 0.5;             % s step (for midpoint RK2)
y0 = [5; 0];          % [v0; x0] = [m/s; m]

N  = floor((tf - t0)/h);
t  = t0 + (0:N)*h;

%% --- Integrate ground roll with event at liftoff ---
opts = odeset('Events', @(t,y) ev_liftoff(t,y,param,to,Vlof), 'RelTol',1e-7, 'AbsTol',1e-8);
[t_ref, y_ref] = ode45(@(t,y) takeoff_ode(t,y,param,to,Vr,Vlof), [t0 tf], y0, opts);

% Also run your midpoint RK2 for comparison (no event handling)
[t_mid, y_mid] = midpoint_rk2(@(t,y) takeoff_ode(t,y,param,to,Vr,Vlof), t0, y0, h, N);

% Ground-roll distance and liftoff speed from event run
x_ground = y_ref(end,2);
v_lof    = y_ref(end,1);

%% --- Air segment to 35 ft at V2 (excess-thrust climb angle) ---
q2   = 0.5*param.rho*V2^2;
CD2  = to.CD0 + to.k*to.CLmax^2;
D2   = q2 * param.S * CD2;
T2   = param.T0 * max(0, 1 - (V2/param.vmax)^2);
gamma= asin( max(0, (T2 - D2) / (param.m*param.g) ) );     % climb angle (rad)

h_screen = 35/3.28084;                                    % 35 ft in meters
gamma_safe = max(deg2rad(0.5), gamma);                    % prevent tiny angles
s_air = h_screen / tan(gamma_safe);

s_TO = x_ground + s_air;

%% --- Diagnostics / Residuals (using instantaneous CL schedule) ---
CL_ref = arrayfun(@(v) cl_schedule(v,Vr,Vlof,to), y_ref(:,1));
residual_ode45 = 0.5 * param.rho * param.S .* CL_ref .* (y_ref(:,1).^2) - W;

CL_mid = arrayfun(@(v) cl_schedule(v,Vr,Vlof,to), y_mid(1,:));
residual_mid = 0.5 * param.rho * param.S .* CL_mid .* (y_mid(1,:).^2) - W;

%% --- Report distances ---
fprintf('C-17 Takeoff distances (sea-level ISA, MTOW):\n');
fprintf('  Ground roll to liftoff:  %.0f m  (%.0f ft)\n', x_ground, x_ground*3.28084);
fprintf('  Air segment to 35 ft:    %.0f m  (%.0f ft)\n', s_air,    s_air*3.28084);
fprintf('  TOTAL to 35 ft:          %.0f m  (%.0f ft)\n', s_TO,      s_TO*3.28084);

%% --- Plots (kept close to your originals) ---
figure('Color','w');

subplot(3,2,1); hold on; grid on
plot(t_ref, y_ref(:,1), 'k-', 'LineWidth', 1.4, 'DisplayName','ode45 (with event)');
plot(t_mid, y_mid(1,:), 'ro', 'DisplayName','Midpoint RK2');
title('Velocity vs Time'); xlabel('Time [s]'); ylabel('v(t) [m/s]'); legend

subplot(3,2,2); hold on; grid on
plot(t_ref, y_ref(:,2), 'k-', 'LineWidth', 1.4, 'DisplayName','ode45 (with event)');
plot(t_mid, y_mid(2,:), 'ro', 'DisplayName','Midpoint RK2');
xline(t_ref(end), '--', 'Liftoff (event)');
title('Position vs Time'); xlabel('Time [s]'); ylabel('x(t) [m]'); legend

subplot(3,2,3); hold on; grid on
plot(y_ref(:,2), y_ref(:,1), 'k-', 'LineWidth', 1.4, 'DisplayName','ode45');
plot(y_mid(2,:), y_mid(1,:), 'ro', 'DisplayName','Midpoint RK2');
title('Velocity vs Position'); xlabel('x(t) [m]'); ylabel('v(t) [m/s]'); legend

subplot(3,2,4); hold on; grid on
plot(t_ref, residual_ode45, 'k-', 'DisplayName','ode45');
plot(t_mid, residual_mid, 'ro', 'DisplayName','Midpoint RK2');
title('Residual vs Time'); xlabel('Time [s]');
ylabel('0.5\rho S C_L v^2 - mg [N]'); legend

subplot(3,2,5); hold on; grid on
plot(y_ref(:,2), residual_ode45, 'k-', 'DisplayName','ode45');
plot(y_mid(2,:), residual_mid, 'ro', 'DisplayName','Midpoint RK2');
xlabel('x(t) [m]'); ylabel('0.5\rho S C_L v^2 - mg [N]');
title('Residual vs Position'); legend

%% ===================== Local functions =====================

function dy = takeoff_ode(~, y, p, to, Vr, Vlof)
    v = y(1);  x = y(2); %#ok<NASGU>

    % Thrust lapse vs speed (simple parabola)
    T = p.T0 * max(0, 1 - (v/p.vmax)^2);

    % CL schedule: roll -> rotate -> CLmax
    CL = cl_schedule(v, Vr, Vlof, to);

    % Drag with induced term
    CD = to.CD0 + to.k * CL^2;

    % Aerodynamic forces
    q = 0.5 * p.rho * v^2;
    L = q * p.S * CL;
    D = q * p.S * CD;

    % Rolling resistance
    Nw   = max(0, p.m * p.g - L);
    Ffr  = p.mu * Nw;

    dvdt = (T - D - Ffr) / p.m;
    dxdt = v;
    dy   = [dvdt; dxdt];
end

function CL = cl_schedule(v, Vr, Vlof, to)
    if v < Vr
        CL = to.CLroll;
    elseif v <= Vlof
        s  = (v - Vr) / max(1e-6, (Vlof - Vr));   % 0→1 during rotation
        CL = to.CLroll + s * (to.CLmax - to.CLroll);
    else
        CL = to.CLmax;
    end
end

function [value, isterminal, direction] = ev_liftoff(~, y, p, to, Vlof)
    v = y(1);
    q = 0.5 * p.rho * v^2;
    L = q * p.S * to.CLmax;
    % Trigger when BOTH: speed >= Vlof AND L >= W
    value      = min([v - Vlof, L - p.m*p.g]);
    isterminal = 1;
    direction  = +1;
end

function [t, Y] = midpoint_rk2(f, t0, y0, h, N)

    Y = zeros(numel(y0), N+1);
    t = t0 + (0:N)*h;
    Y(:,1) = y0;
    for k = 1:N
        tk = t(k);    yk = Y(:,k);
        k1 = f(tk, yk);
        k2 = f(tk + 0.5*h, yk + 0.5*h*k1);
        Y(:,k+1) = yk + h * k2;
    end
end