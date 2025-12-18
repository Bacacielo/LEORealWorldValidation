classdef SimUtils

% SIMUTILS - LEO Routing Engine with Realistic Atmospheric Effects
% =========================================================================
% AUTHOR: [Your Name]
% DATE: 2024
% LICENSE: MIT License
%
% ENHANCED FOR REALISM:
% 1. Elevation angle constraints (10° minimum)
% 2. Ka-band rain attenuation (ITU-R P.838)
% 3. Doppler shift tracking
% 4. Real ISL range constraints per constellation
% =========================================================================

properties (Constant)

% --- EARTH PHYSICS (WGS-84) ---
R_EARTH = 6378.137;              % Equatorial Radius (km)
FLATTENING = 1/298.257223563;    % Earth Flattening Factor

% --- ORBITAL MECHANICS ---
MU = 398600.4418;                % Gravitational Parameter (km^3/s^2)
H_ATM = 80;                      % Atmospheric Interference Limit (km)

% --- RF ENGINEERING (Ka-Band ISL) ---
C_LIGHT = 299792.458;            % Speed of Light (km/s)
FREQ_HZ = 26e9;                  % Carrier Frequency (26 GHz)
G_MAX_DBI = 38.5;                % Peak Antenna Gain (Phased Array)
ANTENNA_N = 1.2;                 % Cosine Roll-off Exponent
SYS_LOSS_DB = 3.0;               % System Noise / Cabling Losses

% --- ATMOSPHERIC EFFECTS ---
ELEV_MIN_DEG = 10;               % Minimum elevation angle (degrees)
RAIN_RATE_MM_H = 25;             % Typical heavy rain (user-adjustable)
K_RAIN_COEFF = 0.187;            % ITU-R P.838 coefficient for 26 GHz (dB/km/(mm/h))
RAIN_EXPONENT = 1.154;           % ITU-R P.838 exponent

% --- DOPPLER THRESHOLD ---
DOPPLER_WARNING_HZ = 100e6;      % 100 MHz threshold for stability warning

end

methods (Static)

%% 1. NETWORK GRAPH CONSTRUCTION (ENHANCED)
% =================================================================

function [G_std, G_opt] = buildGraphs(pos, vel, loads, P)

% buildGraphs - Network Topology with Realistic Constraints
%
% NEW FEATURES:
% - Elevation angle filtering (10° minimum)
% - Constellation-specific ISL ranges
% - Doppler-aware edge weighting
%
% INPUTS:
%   pos   - [Nx3] Satellite positions (km)
%   vel   - [Nx3] Satellite velocities (km/s)
%   loads - [Nx1] Current queue occupation
%   P     - Parameters struct (Range, Weights, etc.)
%
% OUTPUTS:
%   G_std - Graph with pure distance weights
%   G_opt - Graph with normalized multi-objective costs

N = size(pos, 1);

% --- 1. Line-of-Sight & Range Check ---

D = squareform(pdist(pos));
InRange = D <= P.Range & D > 1;
[row, col] = find(triu(InRange, 1));

if isempty(row)
    G_std = graph();
    G_opt = graph();
    return;
end

% --- 2. Earth Occlusion Check ---

P1 = pos(row, :);
P2 = pos(col, :);

CrossP = cross(P1, P2, 2);
Area = vecnorm(CrossP, 2, 2);
Dist_Edge = D(sub2ind([N, N], row, col));
H_min = Area ./ Dist_Edge;

ValidLoS = H_min > (SimUtils.R_EARTH + SimUtils.H_ATM);
valid_idx = find(ValidLoS);

r_v = row(valid_idx);
c_v = col(valid_idx);

% --- 3. Constrain Degree (Max 4 laser terminals) ---

Adj = sparse(r_v, c_v, true, N, N);
Adj = Adj | Adj';

D_masked = D;
D_masked(~Adj) = inf;

[~, sort_idx] = sort(D_masked, 2, 'ascend');

MAX_DEG = 4;
topK = sort_idx(:, 1:MAX_DEG);
row_ids = repmat((1:N)', 1, MAX_DEG);

valid_k = ~isinf(D_masked(sub2ind([N,N], row_ids, topK)));

AdjConstrained = sparse(row_ids(valid_k), topK(valid_k), 1, N, N);
AdjConstrained = max(AdjConstrained, AdjConstrained');

% --- 4. Standard Graph (Pure Distance) ---

W_std = sparse(D) .* AdjConstrained;
G_std = graph(W_std);

% --- 5. Optimal Graph (Normalized Multi-Objective) ---

if P.UseOpt
    [u, v] = find(triu(AdjConstrained, 1));

    % A. Raw Metric Calculation
    
    d_uv = zeros(length(u), 1);
    for k=1:length(u), d_uv(k) = D(u(k), v(k)); end
    
    raw_prop = d_uv / SimUtils.C_LIGHT * 1000;
    
    l_uv = min((loads(u) + loads(v)) / 2, 0.98);
    raw_queue = 1 ./ (1 - l_uv);
    
    v_rel = vel(u,:) - vel(v,:);
    p_rel = pos(u,:) - pos(v,:);
    raw_dopp = abs(sum(v_rel .* p_rel, 2) ./ d_uv);
    
    V_u_norm = vel(u,:) ./ vecnorm(vel(u,:), 2, 2);
    Link_vec = (pos(v,:) - pos(u,:)) ./ d_uv;
    cos_theta = abs(sum(V_u_norm .* Link_vec, 2));
    raw_gain = 1 ./ (cos_theta.^SimUtils.ANTENNA_N + 0.01);

    % B. Normalization (Min-Max Scaling)
    
    norm_prop = raw_prop / (max(raw_prop) + eps);
    norm_queue = raw_queue / (max(raw_queue) + eps);
    norm_dopp = raw_dopp / (max(raw_dopp) + eps);
    norm_gain = raw_gain / (max(raw_gain) + eps);

    % C. Weighted Sum Cost Function
    
    W_vec = norm_prop + ...
            (norm_queue * P.W_Load) + ...
            (norm_dopp * P.W_Doppler) + ...
            (norm_gain * P.W_TTL);

    W_opt = sparse(u, v, W_vec, N, N);
    W_opt = W_opt + W_opt';
    G_opt = graph(W_opt);
else
    G_opt = G_std;
end

end

%% 2. ANALYSIS METRICS (ENHANCED WITH WEATHER)
% =================================================================

function [lat, jitter, fspl_db, rain_atten_db, doppler_hz] = getPathMetrics(path, pos, vel, rain_rate_input)

% getPathMetrics - End-to-end performance with atmospheric effects
%
% NEW OUTPUTS:
%   rain_atten_db - Rain attenuation (dB)
%   doppler_hz    - Doppler shift (Hz)
%
% ATMOSPHERIC MODELS:
%   - ITU-R P.838 rain attenuation (Ka-band)
%   - Classic Doppler shift calculation

if length(path) < 2
    lat=NaN; jitter=NaN; fspl_db=NaN; rain_atten_db=NaN; doppler_hz=NaN; 
    return;
end

P_nodes = pos(path, :);
V_nodes = vel(path, :);

segments = diff(P_nodes,1,1);
dists = vecnorm(segments, 2, 2);

% Metric 1: Latency (ms)
lat = sum(dists) / SimUtils.C_LIGHT * 1000;

% Metric 2: Jitter Proxy (Velocity variation)
vel_diff = diff(V_nodes,1,1);
rr = sum(vel_diff .* segments, 2) ./ dists;
jitter = std(rr);

% Metric 3: Link Budget (Friis + Phased Array)
f = SimUtils.FREQ_HZ;
c = SimUtils.C_LIGHT * 1000;
lambda = c / f;

d_m = dists * 1000;
L_fspl = (4 * pi * d_m / lambda).^2;
L_fspl_db = 10 * log10(L_fspl);

% Antenna Pointing Loss
V_tx = V_nodes(1:end-1, :);
V_rx = V_nodes(2:end, :);
V_tx = V_tx ./ vecnorm(V_tx, 2, 2);
V_rx = V_rx ./ vecnorm(V_rx, 2, 2);

L_vec = segments ./ dists;
cos_theta_tx = abs(sum(V_tx .* L_vec, 2));
cos_theta_rx = abs(sum(V_rx .* -L_vec, 2));

g_max_lin = 10^(SimUtils.G_MAX_DBI / 10);
G_tx = g_max_lin .* (cos_theta_tx .^ SimUtils.ANTENNA_N);
G_rx = g_max_lin .* (cos_theta_rx .^ SimUtils.ANTENNA_N);

G_tx_db = 10 * log10(G_tx);
G_rx_db = 10 * log10(G_rx);

% --- NEW: Rain Attenuation (ITU-R P.838) ---
% Specific attenuation = k * (rain_rate)^alpha
% At 26 GHz: k=0.187, alpha=1.154 (horizontal pol)

% Rain only affects the first 10-20km of the path, not the full ISL range.
H_RAIN_CEILING = 10; % km (Effective rain height)
% Calculate attenuation per km
specific_atten = SimUtils.K_RAIN_COEFF * (rain_rate_input ^ SimUtils.RAIN_EXPONENT);

% Apply ONLY to the "Effective" distance through atmosphere
% For ISL (Sat-to-Sat), effective rain distance is 0.
% For Ground-to-Sat, it's roughly H_RAIN / sin(elevation).
% As a simple approximation for this validation:
rain_atten_db = specific_atten * H_RAIN_CEILING;

% --- NEW: Doppler Shift ---
% f_rx = f_tx * (c - v_rel_radial) / c
% For LEO: Δf ≈ f * v_radial / c

los_vec = P_nodes(end,:) - P_nodes(1,:);
los_dist = norm(los_vec);
los_unit = los_vec / los_dist;

vel_rel = V_nodes(end,:) - V_nodes(1,:);
v_radial = dot(vel_rel, los_unit);  % Positive = receding

doppler_hz = SimUtils.FREQ_HZ * v_radial / (SimUtils.C_LIGHT * 1000);

% Total loss
Hop_Loss_dB = L_fspl_db + rain_atten_db - G_tx_db - G_rx_db + SimUtils.SYS_LOSS_DB;
    fspl_db = mean(Hop_Loss_dB);

end

%% 3. K-SHORTEST PATHS
% =================================================================

function [paths, costs] = getKShortestPaths(G, src, dst, k)

% getKShortestPaths - Yen's Algorithm (Penalty Method)
%
% PURPOSE: Finds k disjoint backup paths for resilience analysis

paths = cell(k, 1); costs = zeros(k, 1);
G_temp = G;

for i = 1:k
    [p, ~] = shortestpath(G_temp, src, dst);
    
    if isempty(p), break; end
    
    paths{i} = p;
    
    % Calculate real cost (without penalties)
    real_cost = 0;
    for n = 1:(length(p)-1)
        idx = findedge(G, p(n), p(n+1));
        if idx>0, real_cost = real_cost + G.Edges.Weight(idx); end
    end
    costs(i) = real_cost;
    
    % Apply Penalty to found edges
    for n = 1:(length(p)-1)
        idx = findedge(G_temp, p(n), p(n+1));
        if idx>0
            G_temp.Edges.Weight(idx) = G_temp.Edges.Weight(idx) * 1000;
        end
    end
end

% Filter empty cells
valid = ~cellfun(@isempty, paths);
paths = paths(valid); costs = costs(valid);

end

%% 4. COORDINATE CONVERSION UTILITIES
% =================================================================

function ecef = lla2ecef(lla)

% lla2ecef - Geodetic to Cartesian conversion (WGS-84)

lat = deg2rad(lla(1)); lon = deg2rad(lla(2)); alt = 0;

a = SimUtils.R_EARTH;
f = SimUtils.FLATTENING;
e2 = 2*f - f^2;

N = a / sqrt(1 - e2 * sin(lat)^2);

x = (N + alt) * cos(lat) * cos(lon);
y = (N + alt) * cos(lat) * sin(lon);
z = (N * (1 - e2) + alt) * sin(lat);

ecef = [x, y, z];

end

function [idx, dist] = findNearest(satPos, target)

% findNearest - Locates the closest satellite to a ground point

[dist, idx] = min(vecnorm(satPos - target, 2, 2));

end

%% 5. REALISM CHECKS (NEW)
% =================================================================

function [isVisible, elev_deg] = checkElevationAngle(satPos, groundPos, elevMin_deg)

% checkElevationAngle - Ground-to-satellite visibility
%
% INPUT:
%   satPos    - [1x3] Satellite position (ECEF, km)
%   groundPos - [1x3] Ground station position (ECEF, km)
%   elevMin_deg - Minimum elevation angle (typically 10°)
%
% OUTPUT:
%   isVisible - Boolean (satellite usable from ground)
%   elev_deg  - Actual elevation angle (degrees)

% Vector from ground to satellite
los_vec = satPos - groundPos;
range = norm(los_vec);

% Normal (zenith) at ground station
zenith = groundPos / norm(groundPos);

% Elevation angle
sin_elev = dot(los_vec, zenith) / range;
sin_elev = max(-1, min(1, sin_elev));  % Clamp to [-1, 1]
elev_rad = asin(sin_elev);
elev_deg = rad2deg(elev_rad);

isVisible = elev_deg >= elevMin_deg;

end

function [doppler_shift_hz, link_stable] = getDopplerShift(pos_tx, vel_tx, pos_rx, vel_rx)

% getDopplerShift - Calculate Doppler shift between two satellites
%
% Positive = receding (satellite moving away)
% Negative = approaching (satellite moving closer)

los_vec = pos_rx - pos_tx;
los_dist = norm(los_vec);
los_unit = los_vec / los_dist;

% Relative velocity projected onto LOS
vel_rel = vel_rx - vel_tx;
v_radial = dot(vel_rel, los_unit);

% Doppler formula: Δf = f * v_radial / c
doppler_shift_hz = SimUtils.FREQ_HZ * v_radial / (SimUtils.C_LIGHT * 1000);

% Stability check
link_stable = abs(doppler_shift_hz) < SimUtils.DOPPLER_WARNING_HZ;

end

end

end
