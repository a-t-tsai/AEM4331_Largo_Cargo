%% --- Parameters (C-17 Landing, MLW, sea-level ISA) ---
param = struct( ...
    'T0',   7.18e5, ...   % [KNOWN] N, total static thrust (4 × 40,400 lbf F117-PW-100)
    'g',    9.81,   ...   % [KNOWN] m/s², standard gravity
    'S',    353.0,  ...   % [KNOWN] m², wing area (3,800 ft²)
    'rho',  1.225,  ...   % [KNOWN] kg/m³, ISA sea-level density
    'mu',   0.020,  ...   % [ASSUMED] rolling friction (taxi/ground roll baseline)
    'm',    2.277e5, ...  % [KNOWN] kg, MLW ≈ 502,100 lb
    'vmax', 120     ...   % [ASSUMED] m/s, thrust-lapse scaling factor
);

% --- C-17 Landing (50 ft to full stop) — MLW, dry runway (realistic) ---
% Tip: ensure your landing mass is set to MLW:
% param.m = 2.277e5;  % [KNOWN] kg, ~502,100 lb (full-flap MLW)

ldg = struct( ...
    'CLmax',     5.014,      ... % [KNOWN] from your PDF (landing CL at SL)
    'CD0',       0.070,      ... % [ASSUMED] landing parasitic drag (flaps+gear)
    'k',         0.050,      ... % [ASSUMED] induced drag factor
    'mu_brake',  0.22,       ... % [ASSUMED] effective braking (conservative dry)
    'CLdump',    0.40,       ... % [ASSUMED] more lift retained during rollout → longer stop
    'beta_rev',  0.05,       ... % [ASSUMED] idle/ reverse 
    'gamma_app', deg2rad(3)  ... % [KNOWN/CONVENTION] 3° approach from 50 ft
);

Vs_ldg = sqrt( 2*param.m*param.g / (param.rho*param.S*ldg.CLmax) );  % [PHYSICS]
Vapp   = 1.30 * Vs_ldg;                                             % [CONVENTION]
Vtd    = 1.25 * Vs_ldg;                                             % [ASSUMED/TUNED]

h_scr  = 50 / 3.28084;    % 50 ft screen height in meters

% --- Air distance from 50 ft to touchdown along the glidepath ---
s_air_ldg = h_scr / tan(ldg.gamma_app);

% --- Ground roll integration from touchdown to full stop ---
y0_ldg = [Vtd; 0];  % [v0; x0] at touchdown
optsL  = odeset('Events', @(t,y) ev_stop(t,y), 'RelTol',1e-7, 'AbsTol',1e-8);
[t_ldg, y_ldg] = ode45(@(t,y) landing_ode(t,y,param,ldg), [0 300], y0_ldg, optsL);

x_stop = y_ldg(end,2);
v_min  = y_ldg(end,1);

% --- Totals ---
LFL_total = s_air_ldg + x_stop;

fprintf('\n=== Landing (from 50 ft, SL ISA) ===\n');
fprintf('  Air segment (50 ft to TD):  %.0f m  (%.0f ft)\n', s_air_ldg, s_air_ldg*3.28084);
fprintf('  Ground roll to stop:        %.0f m  (%.0f ft)\n', x_stop,    x_stop*3.28084);
fprintf('  LFL (50 ft to full stop):   %.0f m  (%.0f ft)\n', LFL_total,  LFL_total*3.28084);

% --- Quick plots for landing (similar style to takeoff) ---
figure('Color','w');
subplot(1,2,1); hold on; grid on
plot(t_ldg, y_ldg(:,1), 'k-', 'LineWidth', 1.4);
title('Landing: Velocity vs Time'); xlabel('Time [s]'); ylabel('v(t) [m/s]');
subplot(1,2,2); hold on; grid on
plot(t_ldg, y_ldg(:,2), 'k-', 'LineWidth', 1.4);
title('Landing: Position vs Time'); xlabel('Time [s]'); ylabel('x(t) [m]');

%% ===================== Local functions (landing) =====================
function dv = T_rev_available(v, p, ldg)
    % Reverse thrust model (as a fraction of static T0, lapses with v)
    Tmag = ldg.beta_rev * p.T0 * max(0, 1 - (v/p.vmax)^2);
    dv   = -Tmag; % negative thrust (acts backward)
end

function dy = landing_ode(~, y, p, ldg)
    v = max(0, y(1));  %#ok<NASGU> keep non-negative numerically

    % Aerodynamics during rollout (spoilers/lift-dump active)
    CL = ldg.CLdump;
    CD = ldg.CD0 + ldg.k * CL^2;

    q = 0.5 * p.rho * v^2;
    L = q * p.S * CL;
    D = q * p.S * CD;

    % Normal load and braking friction (rolling resistance included in mu_brake)
    Nw   = max(0, p.m * p.g - L);
    Ffr  = ldg.mu_brake * Nw;

    % Reverse thrust (negative)
    Trev = T_rev_available(v, p, ldg);

    dvdt = (-D - Ffr + Trev) / p.m;
    dxdt = v;
    dy   = [dvdt; dxdt];
end

function [value,isterminal,direction] = ev_stop(~, y)
    % Stop when speed reaches zero
    value      = y(1);     % v = 0
    isterminal = 1;
    direction  = -1;
end

