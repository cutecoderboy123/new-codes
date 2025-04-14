%% Countermovement Jump (CMJ) Trial 2 - Biomechanics Lab

% Clear everything
close all; 
clear; 
clc; 

disp('Starting to process CMJ Trial 2 data...');

%% Set up constants
mass = 75.5;          % kg
grav = 9.81;          % gravity m/s^2
L_LL = 0.840;         % left leg length m
R_LL = 0.840;         % right leg length m
Leg_L = (L_LL + R_LL) / 2; % avg leg length
R_mark = 0.012;       % some marker offset thing

%% Load all the data files for Trial 2
% Load trajectory data
opts2_traj = detectImportOptions('Group 3 S CMJ 02 Trial 2 trajectories.csv', 'NumHeaderLines', 5);
traj2 = readmatrix('Group 3 S CMJ 02 Trial 2 trajectories.csv', opts2_traj);

% Load force data 
opts2_force = detectImportOptions('Group 3 S CMJ 02 Trial 2 force.csv', 'NumHeaderLines', 5);
force2 = readmatrix('Group 3 S CMJ 02 Trial 2 force.csv', opts2_force);

%% Get all the basic stuff from the data
% Get bodyweight - average when standing still
BW2 = mean(force2(force2(:,2)==0, 5)); 
f_max2 = size(traj2,1); % number of frames
Freq2 = 100; % Hz sampling rate I think
t2 = (0:f_max2-1)/Freq2; % time array

% Filter out the force data we need
ReSamp2 = find(force2(:,2) == 0); 
combined_normal2 = -force2(ReSamp2,5); % negative because of direction I think?
t_force2 = (0:length(combined_normal2)-1)/1000; % time for force data

% Get the markers - column numbers from the manual
LASI2 = traj2(:,72:74); 
RASI2 = traj2(:,75:77); 
LPSI2 = traj2(:,78:80); 
RPSI2 = traj2(:,81:83);

%% Gold Standard Method - all the complicated stuff from lab manual
% Calculate pelvic width
pelvic_w2 = LASI2 - RASI2;
ASIS_Dist2 = mean(sqrt(sum(pelvic_w2.^2,2))) / 1000; % convert to m

% These are constants from the manual - not sure what they all mean
X_dist2 = 0.1288 * Leg_L - 0.04856;
C2 = Leg_L * 0.115 - 0.0153;
B2 = deg2rad(18); 
T2 = deg2rad(28.4);

% Calculate hip joint center locations in local coords
LHJC_L2 = [(C2*sin(T2) - 0.5*ASIS_Dist2); 
           (-X_dist2 - R_mark)*cos(B2) + (C2*cos(T2)*sin(B2)); 
           (-X_dist2 - R_mark)*sin(B2) - (C2*cos(T2)*cos(B2))] * 1000;
           
RHJC_L2 = [-(C2*sin(T2) - 0.5*ASIS_Dist2); 
           (-X_dist2 - R_mark)*cos(B2) + (C2*cos(T2)*sin(B2)); 
           (-X_dist2 - R_mark)*sin(B2) - (C2*cos(T2)*cos(B2))] * 1000;

% Get pelvis origin
pelvis_org2 = (RASI2 + LASI2) / 2;
LHJC2 = zeros(3, f_max2); 
RHJC2 = zeros(3, f_max2);

% This loop is super confusing - not sure I understand it all
for f = 1:f_max2
    % Calculate reference frame or something?
    i2 = (RASI2(f,:).' - pelvis_org2(f,:).') / norm(RASI2(f,:) - pelvis_org2(f,:));
    L22 = pelvis_org2(f,:).' - RPSI2(f,:).'; 
    L22 = L22 / norm(L22);
    k2 = cross(i2, L22); 
    j2 = cross(k2, i2); 
    GRL2 = [i2, j2, k2];
    
    % Calculate hip joint centers in global coords
    LHJC2(:,f) = GRL2 * LHJC_L2 + pelvis_org2(f,:).';
    RHJC2(:,f) = GRL2 * RHJC_L2 + pelvis_org2(f,:).';
end

% Calculate COM - just take average of hip centers for now
COM2 = (LHJC2 + RHJC2) / 2;
COM_Z2 = COM2(3,:) / 1000; % convert to meters (z component only)

% Get initial height from first second - seems to work ok
first_second_idx2 = find(t2 <= 1);
initial_height_gold_standard2 = mean(COM_Z2(first_second_idx2));

%% Kinematic Method - easier to understand than the gold standard one
% F = ma → a = F/m
net_force2 = combined_normal2 - BW2; % Subtract bodyweight
acceleration2 = net_force2 / mass; 

% v = ∫a dt
velocity2 = cumtrapz(t_force2, acceleration2);

% Find takeoff and landing - when force drops below threshold
foot_off_idx2 = find(combined_normal2 < 50, 1, 'first');
foot_land_idx2 = find(combined_normal2 < 50, 1, 'last');

% Get max velocity at takeoff
max_velocity2 = max(velocity2(1:foot_off_idx2));

% d = ∫v dt
displacement_kin2 = cumtrapz(t_force2, velocity2) + initial_height_gold_standard2;

%% Projectile Motion Calculation
% Time from takeoff to landing
flight_time2 = t_force2(foot_land_idx2) - t_force2(foot_off_idx2);
mid_flight_time2 = t_force2(foot_off_idx2) + flight_time2/2; % middle of flight

% Get just the flight phase times
proj_t2 = t_force2(t_force2 >= t_force2(foot_off_idx2) & t_force2 <= t_force2(foot_land_idx2));
proj_t_rel2 = proj_t2 - t_force2(foot_off_idx2); % times relative to takeoff

% Ballistic equation: h = v0*t - 0.5*g*t^2 + h0
proj_disp2 = max_velocity2 .* proj_t_rel2 - 0.5 * grav .* proj_t_rel2.^2 + initial_height_gold_standard2;

%% Energy Method - simplest way
% Conservation of energy: mgh = 0.5mv^2
% h = v^2/(2g)
energy_height2 = (max_velocity2^2)/(2*grav) + initial_height_gold_standard2;

%% Make plots for trial 2
fig1 = figure('Name', 'Height Methods - CMJ Trial 2');

% 1. Gold Standard Method plot
subplot(2, 2, 1);
plot(t2, COM_Z2, 'k', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Gold Standard Method');
grid on;

% 2. Kinematic Method plot
subplot(2, 2, 2);
plot(t_force2, displacement_kin2, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Kinematic Method');
grid on;

% 3. Projectile Motion plot
subplot(2, 2, 3);
plot(proj_t2, proj_disp2, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Projectile Motion');
grid on;

% 4. Energy Method plot - just a line
subplot(2, 2, 4);
plot(t_force2, ones(size(t_force2)) * energy_height2, 'm--', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Energy Method');
grid on;

%% Force-time graph showing jump phases
fig2 = figure('Name', 'Force vs Time - Trial 2');
plot(t_force2, combined_normal2, 'k', 'LineWidth', 1.5);
hold on;

% Add vertical lines to mark takeoff and landing
xline(t_force2(foot_off_idx2), 'r--', 'Foot Off', 'LineWidth', 1.5);
xline(t_force2(foot_land_idx2), 'b--', 'Foot Land', 'LineWidth', 1.5);

xlabel('Time (s)'); 
ylabel('Force (N)');
title('Force-Time Graph - Trial 2');
legend('Vertical GRF', 'Foot Off', 'Foot Land');
grid on;

%% Plot all the movement parameters in one figure
% Create the motion parameters figure
fig3 = figure('Name', 'Motion Parameters - Trial 2');

% Plot acceleration
subplot(2, 2, 1);
plot(t_force2, acceleration2, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Accel (m/s^2)');
title('Acceleration vs Time');
grid on;

% Plot velocity
subplot(2, 2, 2);
plot(t_force2, velocity2, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Vel (m/s)');
title('Velocity vs Time');
grid on;

% Plot displacement
subplot(2, 2, 3);
plot(t_force2, displacement_kin2, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Height vs Time');
grid on;

% Plot all three together
subplot(2, 2, 4);
plot(t_force2, acceleration2, 'r', 'LineWidth', 1.5); 
hold on;
plot(t_force2, velocity2, 'g', 'LineWidth', 1.5);
plot(t_force2, displacement_kin2, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Values');
title('All Parameters');
legend('Accel', 'Vel', 'Height');
grid on;