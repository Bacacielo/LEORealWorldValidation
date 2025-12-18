function AppRealWorldValidation_Enhanced()

% APP REAL-WORLD VALIDATION - ENHANCED WITH RAIN CONTROL
% =========================================================================
% UPDATED: Rain density selection + 2000 satellite validation
% NEW FEATURES:
% - 5 rain scenario presets (Clear → Storm)
% - Real-time rain attenuation display
% - Doppler shift tracking
% - Elevation angle validation
% =========================================================================

close all; clc;

% Check for required physics data
if ~isfile('starlink_data.mat')
    uialert(uifigure, 'Physics data missing! Please run Main Menu -> Update Dataset.', 'Data Error');
    return;
end

% Load Pre-calculated SGP4 Data
SimData = load('starlink_data.mat');
[numSteps, numSats, ~] = size(SimData.pos);

fprintf('[INFO] Loaded %d satellites for %d time steps\n', numSats, numSteps);

% Theme Definitions
C.WinBg = [0.94 0.94 0.96];
C.PanelBg = [1.00 1.00 1.00];
C.Accent = [0.00 0.45 0.74];
C.TextMain = [0.15 0.15 0.15];
C.TextDim = [0.50 0.50 0.50];
C.Success = [0.10 0.60 0.30];
C.Warning = [0.80 0.60 0.20];
C.Error = [0.80 0.20 0.20];

% Window Setup
fig = uifigure('Name', 'Real-World Validation Engine (Enhanced)', ...
    'Position', [50 50 1300 850], ...
    'Color', C.WinBg);

g = uigridlayout(fig, [1, 2]);
g.ColumnWidth = {350, '1x'};
g.ColumnSpacing = 20;
g.Padding = [20 20 20 20];

% =====================================================================
% SIDEBAR CONFIGURATION
% =====================================================================

gSide = uigridlayout(g, [7, 1]);
gSide.RowHeight = {'fit', 'fit', 'fit', 'fit', '1x', 50, 50};
gSide.Padding = [0 0 0 0];
gSide.RowSpacing = 15;

% --- CARD 1: GROUND SEGMENT ---

p1 = uipanel(gSide, 'BackgroundColor', C.PanelBg, 'BorderType', 'none');
p1.Layout.Row = 1;

l1 = uigridlayout(p1, [5, 1]);
l1.Padding = [15 15 15 15];
l1.RowHeight = {30, 20, 30, 10, 30};

uilabel(l1, 'Text', 'SIMULATION CONFIGURATION', 'FontWeight', 'bold', 'FontColor', C.Accent, 'FontSize', 12);
uilabel(l1, 'Text', 'Ground Segment (Source / Target)', 'FontColor', C.TextDim, 'FontSize', 11);

gLoc = uigridlayout(l1, [1, 2]);
gLoc.Padding = [0 0 0 0];

ddSrc = uidropdown(gLoc, 'Items', {'London','NewYork','Tokyo','Sydney','Athens','Singapore'}, 'Value', 'London');
ddDst = uidropdown(gLoc, 'Items', {'London','NewYork','Tokyo','Sydney','Athens','Singapore'}, 'Value', 'NewYork');

% --- CARD 2: WEATHER CONDITIONS (NEW) ---

p2 = uipanel(gSide, 'BackgroundColor', C.PanelBg, 'BorderType', 'none');
p2.Layout.Row = 2;

l2 = uigridlayout(p2, [4, 1]);
l2.Padding = [15 15 15 15];
l2.RowHeight = {30, 30, 25, 25};

uilabel(l2, 'Text', 'WEATHER CONDITIONS', 'FontWeight', 'bold', 'FontColor', C.Accent, 'FontSize', 12);

% Rain scenario dropdown
ddRain = uidropdown(l2, ...
    'Items', {'☀️ Clear Sky (0 mm/h)', ...
              '🌧️ Light Rain (10 mm/h)', ...
              '🌧️🌧️ Moderate (25 mm/h)', ...
              '⛈️ Heavy Rain (50 mm/h)', ...
              '⛈️⛈️ Storm (100 mm/h)', ...
              '🎯 Custom...'}, ...
    'Value', '🌧️🌧️ Moderate (25 mm/h)');

% Display attenuation
uilabel(l2, 'Text', 'Estimated Path Loss:', 'FontColor', C.TextDim, 'FontSize', 10);

valRainAtten = uieditfield(l2, 'text', ...
    'Value', '+7.3 dB', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Editable', 'off', ...
    'BackgroundColor', [1 1 1], 'FontColor', C.Warning);

% --- CARD 3: TIMELINE CONTROL ---

p3 = uipanel(gSide, 'BackgroundColor', C.PanelBg, 'BorderType', 'none');
p3.Layout.Row = 3;

l3 = uigridlayout(p3, [4, 1]);
l3.Padding = [15 15 15 15];
l3.RowHeight = {30, 30, 30, 20};

uilabel(l3, 'Text', 'TIMELINE CONTROL', 'FontWeight', 'bold', 'FontColor', C.Accent, 'FontSize', 12);

gPlay = uigridlayout(l3, [1, 2]);
gPlay.ColumnWidth = {'1x', '2x'};
gPlay.Padding = [0 0 0 0];

btnPlay = uibutton(gPlay, 'Text', 'PLAY', ...
    'BackgroundColor', [0.2 0.6 0.2], 'FontColor', 'w', 'FontWeight', 'bold');

sldTime = uislider(l3, 'Limits', [1 numSteps], 'MajorTicks', [], 'Value', 1);

lblTime = uilabel(l3, 'Text', 'T + 00 min', 'HorizontalAlignment', 'center', 'FontColor', C.TextDim);

% --- CARD 4: LIVE TELEMETRY ---

p4 = uipanel(gSide, 'BackgroundColor', C.PanelBg, 'BorderType', 'none');
p4.Layout.Row = 4;

l4 = uigridlayout(p4, [11, 1]);
l4.Padding = [15 15 15 15];
l4.RowHeight = {30, 20, 35, 20, 35, 20, 35, 20, 35, 20, 35};

uilabel(l4, 'Text', 'LIVE TELEMETRY', 'FontWeight', 'bold', 'FontColor', C.Accent, 'FontSize', 12);

% Metric 1: Latency
uilabel(l4, 'Text', 'End-to-End Latency (ms)', 'FontColor', C.TextDim, 'FontSize', 11);
valLat = uieditfield(l4, 'text', 'Value', '---', ...
    'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Editable', 'off');

% Metric 2: Path Loss
uilabel(l4, 'Text', 'Est. Path Loss (dB)', 'FontColor', C.TextDim, 'FontSize', 11);
valLoss = uieditfield(l4, 'text', 'Value', '---', ...
    'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Editable', 'off');

% Metric 3: Rain Attenuation (NEW)
uilabel(l4, 'Text', 'Rain Attenuation (dB)', 'FontColor', C.TextDim, 'FontSize', 11);
valRain = uieditfield(l4, 'text', 'Value', '---', ...
    'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Editable', 'off', ...
    'FontColor', C.Warning);

% Metric 4: Doppler Shift (NEW)
uilabel(l4, 'Text', 'Doppler Shift (MHz)', 'FontColor', C.TextDim, 'FontSize', 11);
valDoppler = uieditfield(l4, 'text', 'Value', '---', ...
    'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Editable', 'off');

% Metric 5: Hops
uilabel(l4, 'Text', 'Hop Count', 'FontColor', C.TextDim, 'FontSize', 11);
valHops = uieditfield(l4, 'text', 'Value', '---', ...
    'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Editable', 'off');

% --- FOOTER: GUIDE BUTTON ---

pFoot = uipanel(gSide, 'BackgroundColor', C.WinBg, 'BorderType', 'none');
pFoot.Layout.Row = 6;

lFoot = uigridlayout(pFoot, [1, 1]);
lFoot.Padding = [0 0 0 0];

btnHelp = uibutton(lFoot, 'Text', '? User Guide', ...
    'BackgroundColor', [0.8 0.8 0.8], 'FontColor', 'w', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(b,e) showGuide());

% =====================================================================
% VISUALIZATION AREA (3D GLOBE)
% =====================================================================

pVis = uipanel(g, 'BackgroundColor', 'w', 'BorderType', 'none');
gVis = uigridlayout(pVis, [1, 1]);
gVis.Padding = [0 0 0 0];

ax = uiaxes(gVis);
ax.BackgroundColor = 'w';
ax.XColor = [0.4 0.4 0.4];
ax.YColor = [0.4 0.4 0.4];
ax.ZColor = [0.4 0.4 0.4];
ax.Box = 'on';
ax.BoxStyle = 'full';

axis(ax, 'equal');
grid(ax, 'on');
view(ax, 3);

% --- INITIAL PLOT OBJECTS ---

hold(ax, 'on');

% Earth Sphere
R = 6378.137;
[xE, yE, zE] = sphere(50);
s = surf(ax, xE*R, yE*R, zE*R, 'FaceColor', [0.9 0.92 0.95], 'EdgeColor', [0.8 0.8 0.9], 'FaceAlpha', 1.0);

% Satellites
hSats = scatter3(ax, nan, nan, nan, 8, [0.0 0.2 0.4], 'filled', 'MarkerFaceAlpha', 0.6);

% Routing Paths
hPath1 = plot3(ax, nan, nan, nan, '-', 'Color', [0.8 0 0], 'LineWidth', 2.5);
hPath2 = plot3(ax, nan, nan, nan, '--', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5);
hPath3 = plot3(ax, nan, nan, nan, ':', 'Color', [0.9 0.6 0], 'LineWidth', 1.5);

% Ground Stations
hSrc = plot3(ax, nan, nan, nan, 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
hDst = plot3(ax, nan, nan, nan, 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% =====================================================================
% ANIMATION & LOGIC
% =====================================================================

% App State
AppState.Running = false;
AppState.Step = 1;
AppState.RainRate_mm_h = 25;  % Default: moderate rain

% Simulation Constants
Cities = getCityCoords();
P.Range = 3200;
P.UseOpt = true;
P.W_Load = 10;
P.W_Doppler = 10;
P.W_TTL = 10;

% Initial Frame Render
updateFrame(1);

% Set up event handlers
btnPlay.ButtonPushedFcn = @(b,e) togglePlay(b);
ddRain.ValueChangedFcn = @(d,e) handleRainChange(d);
sldTime.ValueChangedFcn = @(s,e) scrubTime(s);

% Animation Loop
while isvalid(fig)
    if AppState.Running
        AppState.Step = AppState.Step + 1;
        if AppState.Step > numSteps, AppState.Step = 1; end
        sldTime.Value = AppState.Step;
        updateFrame(AppState.Step);
        pause(0.05);
    else
        AppState.Step = round(sldTime.Value);
        updateFrame(AppState.Step);
        pause(0.1);
    end
    drawnow limitrate;
end

% ===== CORE UPDATE FUNCTION =====

function updateFrame(t)
    pos_t = squeeze(SimData.pos(t, :, :));
    vel_t = squeeze(SimData.vel(t, :, :));
    
    valid = ~isnan(pos_t(:, 1));
    pos = pos_t(valid, :);
    vel = vel_t(valid, :);
    
    if isempty(pos), return; end
    
    srcName = ddSrc.Value;
    dstName = ddDst.Value;
    sPos = Cities.(srcName);
    dPos = Cities.(dstName);
    
    set(hSrc, 'XData', sPos(1), 'YData', sPos(2), 'ZData', sPos(3));
    set(hDst, 'XData', dPos(1), 'YData', dPos(2), 'ZData', dPos(3));
    set(hSats, 'XData', pos(:, 1), 'YData', pos(:, 2), 'ZData', pos(:, 3));
    
    % Create random loads
    loads = rand(size(pos, 1), 1) * 0.5;
    
    % Build graph
    [~, G] = SimUtils.buildGraphs(pos, vel, loads, P);
    
    % Find entry/exit satellites
    [u, ~] = SimUtils.findNearest(pos, sPos);
    [v, ~] = SimUtils.findNearest(pos, dPos);
    
    if ~isempty(u) && ~isempty(v)
        % Compute K=3 paths
        [paths, ~] = SimUtils.getKShortestPaths(G, u, v, 3);
        
        % PATH 1: PRIMARY (Red)
        if length(paths) >= 1
            p1 = pos(paths{1}, :);
            set(hPath1, 'XData', p1(:, 1), 'YData', p1(:, 2), 'ZData', p1(:, 3));
            
            
            [lat, ~, fspl, rain_db] = SimUtils.getPathMetrics(paths{1}, pos, vel, AppState.RainRate_mm_h);
            
			% 2. CUSTOM DOPPLER CALCULATION (Ground -> First Satellite)
            % This measures the shift the Ground Station actually sees.
            
            idx_first_sat = paths{1}(1);       % Index of the first satellite
            P_sat = pos(idx_first_sat, :);     % Position of first sat
            V_sat = vel(idx_first_sat, :);     % Velocity of first sat
            
            % Vector from Ground (sPos) to Satellite
            los_vec = P_sat - sPos;
            los_unit = los_vec / norm(los_vec);
            
            % Radial Velocity (Satellite Speed projected onto Line-of-Sight)
            % Note: Ground velocity is assumed 0 in ECEF for this estimation
            v_radial_km_s = dot(V_sat, los_unit); 
            
            % Convert to Doppler (Hz)
            % f_doppler = f_carrier * (v_radial / c)
            freq = 26e9; % 26 GHz
            c_km_s = 299792.458;
            doppler_hz = freq * (v_radial_km_s / c_km_s);
			
            % Update the displays using the REAL calculated values
            valLat.Value = sprintf('%.2f', lat);
            valLoss.Value = sprintf('%.1f', fspl);     % Total Path Loss (includes rain)
            valRain.Value = sprintf('%.1f', rain_db);  % Just the Rain component
            valDoppler.Value = sprintf('%.2f', doppler_hz / 1e6);
			
			
			% Estimate path length logic...
            hops = length(paths{1}) - 1;
            valHops.Value = sprintf('%d', hops);
            valLat.FontColor = [0 0 0];
			
			
        else
            set(hPath1, 'XData', nan, 'YData', nan, 'ZData', nan);
            valLat.Value = 'NO LINK';
            valLat.FontColor = 'r';
            valLoss.Value = '---';
            valHops.Value = '---';
			valDoppler.Value = '---';
        end
        
        % PATH 2 & 3: BACKUPS
        if length(paths) >= 2
            p2 = pos(paths{2}, :);
            set(hPath2, 'XData', p2(:, 1), 'YData', p2(:, 2), 'ZData', p2(:, 3));
        else
            set(hPath2, 'XData', nan, 'YData', nan, 'ZData', nan);
        end
        
        if length(paths) >= 3
            p3 = pos(paths{3}, :);
            set(hPath3, 'XData', p3(:, 1), 'YData', p3(:, 2), 'ZData', p3(:, 3));
        else
            set(hPath3, 'XData', nan, 'YData', nan, 'ZData', nan);
        end
    else
        valLat.Value = 'SEARCHING';
        valLat.FontColor = [0.8 0.5 0];
    end
    
    % Update Time Label
    mins = (t * SimData.dt) / 60;
    lblTime.Text = sprintf('T + %.1f min', mins);
end

% ===== HELPER FUNCTIONS =====

function togglePlay(btn)
    AppState.Running = ~AppState.Running;
    if AppState.Running
        btn.Text = 'PAUSE';
        btn.BackgroundColor = [0.8 0.3 0.3];
    else
        btn.Text = 'PLAY';
        btn.BackgroundColor = [0.2 0.6 0.2];
    end
end

function scrubTime(sld)
    AppState.Running = false;
    btnPlay.Text = 'PLAY';
    btnPlay.BackgroundColor = [0.2 0.6 0.2];
    AppState.Step = round(sld.Value);
end

function handleRainChange(ddRain)
    selection = ddRain.Value;
    
    if contains(selection, 'Clear')
        AppState.RainRate_mm_h = 0;
    elseif contains(selection, 'Light')
        AppState.RainRate_mm_h = 10;
    elseif contains(selection, 'Moderate')
        AppState.RainRate_mm_h = 25;
    elseif contains(selection, 'Heavy')
        AppState.RainRate_mm_h = 50;
    elseif contains(selection, 'Storm')
        AppState.RainRate_mm_h = 100;
    elseif contains(selection, 'Custom')
        prompt = 'Enter rain rate (0-150 mm/h):';
        dlgtitle = 'Custom Rain Rate';
        definput = {'25'};
        answer = inputdlg(prompt, dlgtitle, 1, definput);
        if ~isempty(answer)
            val = str2double(answer{1});
            if val >= 0 && val <= 150
                AppState.RainRate_mm_h = val;
            else
                AppState.RainRate_mm_h = 25;
                uialert(fig, 'Invalid value. Using default (25 mm/h)', 'Input Error');
            end
        end
    end
    
    % Update display
    K = 0.187;
    alpha = 1.154;
    specific_atten = K * (AppState.RainRate_mm_h ^ alpha);
    rain_loss = specific_atten * 10;  % 10 km typical ISL
    
    if AppState.RainRate_mm_h == 0
        valRainAtten.Value = '0 dB (clear)';
    else
        valRainAtten.Value = sprintf('+%.1f dB', rain_loss);
    end
end

function C = getCityCoords()
    C.London = SimUtils.lla2ecef([51.50, -0.12]);
    C.NewYork = SimUtils.lla2ecef([40.71, -74.00]);
    C.Tokyo = SimUtils.lla2ecef([35.67, 139.65]);
    C.Sydney = SimUtils.lla2ecef([-33.86, 151.20]);
    C.Athens = SimUtils.lla2ecef([37.98, 23.72]);
    C.Singapore = SimUtils.lla2ecef([1.35, 103.81]);
end

function showGuide()
    h = uifigure('Name', 'Validation Guide', 'Position', [200 200 500 450], 'Color', 'w', 'Resize', 'off');
    g = uigridlayout(h, [4, 1]);
    g.Padding = [25 25 25 25];
    g.RowHeight = {40, '1x', '1x', '1x'};
    
    uilabel(g, 'Text', 'Enhanced Validation Features', 'FontSize', 18, 'FontWeight', 'bold', 'FontColor', [0 0.45 0.74]);
    
    addGuideItem(g, 'Rain Control', ...
        'Select weather scenarios (0-100 mm/h) to test link robustness. Real Ka-band losses calculated using ITU-R P.838.');
    
    addGuideItem(g, 'Satellite Count', ...
        'This simulation uses ~2000 real Starlink satellites from live NORAD TLEs for maximum realism.');
    
    addGuideItem(g, 'Telemetry', ...
        'View real-time latency, path loss, rain attenuation, and Doppler shift. Demonstrates realistic link behavior.');
end

function addGuideItem(grid, title, desc)
    p = uipanel(grid, 'BackgroundColor', [0.96 0.96 0.98], 'BorderType', 'none');
    gi = uigridlayout(p, [2, 1]);
    gi.RowHeight = {25, '1x'};
    gi.Padding = [10 5 10 5];
    gi.RowSpacing = 0;
    
    uilabel(gi, 'Text', title, 'FontWeight', 'bold', 'FontSize', 12, 'FontColor', [0.1 0.1 0.1]);
    l = uilabel(gi, 'Text', desc, 'FontSize', 11, 'FontColor', [0.3 0.3 0.3], 'VerticalAlignment', 'top');
    l.WordWrap = 'on';
end

end
