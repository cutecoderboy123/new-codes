%% Squat Jump (SJ) Analysis - Trial 2
% Lab Assignment code for analyzing squat jump data

% Start fresh
close all; 
clear; 
clc; 

disp('Starting SJ Trial 2 analysis...');

%% Constants for calculations
mass = 75.5;          % mass (kg)
grav = 9.81;          % gravity (m/s^2)
L_LL = 0.840;         % left leg length (m)
R_LL = 0.840;         % right leg length (m)
Leg_L = (L_LL + R_LL) / 2; % average leg length
R_mark = 0.012;       % some offset thing

%% Load the data files
% Get the trajectory data first
opts_traj = detectImportOptions('Group 3 S SJ 02 Trial 2  trajectories.csv', 'NumHeaderLines', 5);
traj = readmatrix('Group 3 S SJ 02 Trial 2  trajectories.csv', opts_traj);

% Now get the force plate data
opts_force = detectImportOptions('Group 3 S SJ 02 Trial 2 force.csv', 'NumHeaderLines', 5);
force = readmatrix('Group 3 S SJ 02 Trial 2 force.csv', opts_force);

%% Get the basic data from files
% Calculate bodyweight from initial standing
BW = mean(force(force(:,2)==0, 5)); 
f_max = size(traj,1); % number of frames
Freq = 100; % sampling freq (Hz)
t = (0:f_max-1)/Freq; % time array 

% Filter and get force data
ReSamp = find(force(:,2) == 0); % not sure what this filter does exactly
combined_normal = -force(ReSamp,5); % vertical force - negative bc of direction?
t_force = (0:length(combined_normal)-1)/1000; % time for force data

% Get marker positions - these are from the manual, not sure why these columns
LASI = traj(:,72:74); % left front
RASI = traj(:,75:77); % right front
LPSI = traj(:,78:80); % left back
RPSI = traj(:,81:83); % right back

%% Gold Standard Method - complicated stuff
% Get pelvic width
pelvic_w = LASI - RASI;
ASIS_Dist = mean(sqrt(sum(pelvic_w.^2,2))) / 1000; % divide by 1000 to get meters

% These are the magic numbers from lab manual
X_dist = 0.1288 * Leg_L - 0.04856; % no idea what this means
C = Leg_L * 0.115 - 0.0153; % some constant?
B = deg2rad(18); % angle in radians
T = deg2rad(28.4); % another angle

% Hip joint center calcs - this is confusing
LHJC_L = [(C*sin(T) - 0.5*ASIS_Dist); 
         (-X_dist - R_mark)*cos(B) + (C*cos(T)*sin(B)); 
         (-X_dist - R_mark)*sin(B) - (C*cos(T)*cos(B))] * 1000;
         
RHJC_L = [-(C*sin(T) - 0.5*ASIS_Dist); 
         (-X_dist - R_mark)*cos(B) + (C*cos(T)*sin(B)); 
         (-X_dist - R_mark)*sin(B) - (C*cos(T)*cos(B))] * 1000;

% Get pelvis origin
pelvis_org = (RASI + LASI) / 2; % middle of front markers
LHJC = zeros(3, f_max); 
RHJC = zeros(3, f_max);

% Loop thru all frames - rotation matrix stuff I think
for f = 1:f_max
    % Calculate some reference frame? Not 100% sure I understand this
    i = (RASI(f,:).' - pelvis_org(f,:).') / norm(RASI(f,:) - pelvis_org(f,:));
    L2 = pelvis_org(f,:).' - RPSI(f,:).'; 
    L2 = L2 / norm(L2);
    k = cross(i, L2); 
    j = cross(k, i); 
    GRL = [i, j, k]; % rotation matrix I think?
    
    % Calculate hip joint centers
    LHJC(:,f) = GRL * LHJC_L + pelvis_org(f,:).';
    RHJC(:,f) = GRL * RHJC_L + pelvis_org(f,:).';
end

% Get COM position from hip centers (not sure if this is right approach)
COM = (LHJC + RHJC) / 2;
COM_Z = COM(3,:) / 1000; % just z-component in meters

% Get initial height by averaging first second
first_second_idx = find(t <= 1);
initial_height_gold_standard = mean(COM_Z(first_second_idx));

%% Kinematic Method - easier to understand than gold standard
% F = ma so a = F/m
net_force = combined_normal - BW; % subtract bodyweight
acceleration = net_force / mass;

% Integrate accel to get velocity
velocity = cumtrapz(t_force, acceleration);

% Find takeoff and landing - when force is low
foot_off_idx = find(combined_normal < 50, 1, 'first');
foot_land_idx = find(combined_normal < 50, 1, 'last');

% Get max velocity at takeoff
max_velocity = max(velocity(1:foot_off_idx));

% Integrate velocity to get displacement/position
displacement_kin = cumtrapz(t_force, velocity) + initial_height_gold_standard;

%% Projectile Motion - physics formulas
% Time in air
flight_time = t_force(foot_land_idx) - t_force(foot_off_idx);
mid_flight_time = t_force(foot_off_idx) + flight_time/2; % middle point

% Get flight phase times
proj_t = t_force(t_force >= t_force(foot_off_idx) & t_force <= t_force(foot_land_idx));
proj_t_rel = proj_t - t_force(foot_off_idx); % make relative to takeoff

% Use projectile formula: h = v0*t - 0.5*g*t^2 + h0
proj_disp = max_velocity .* proj_t_rel - 0.5 * grav .* proj_t_rel.^2 + initial_height_gold_standard;

%% Energy Method - simplest approach
% PE = KE so mgh = 0.5mv^2
% h = v^2/2g 
energy_height = (max_velocity^2)/(2*grav) + initial_height_gold_standard;

%% Plot the different height calculation methods
figHeights = figure('Name', 'SJ Height Methods - Trial 1');

% Plot 1 - Gold Standard
subplot(2, 2, 1);
plot(t, COM_Z, 'k', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Gold Standard Method');
grid on;

% Plot 2 - Kinematic Method
subplot(2, 2, 2);
plot(t_force, displacement_kin, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Kinematic Method');
grid on;

% Plot 3 - Projectile Motion
subplot(2, 2, 3);
plot(proj_t, proj_disp, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Projectile Motion');
grid on;

% Plot 4 - Energy Method (constant line)
subplot(2, 2, 4);
plot(t_force, ones(size(t_force)) * energy_height, 'm--', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Energy Method');
grid on;

%% Plot force data with jump phases
figForce = figure('Name', 'Force vs Time - SJ Trial 1');
plot(t_force, combined_normal, 'k', 'LineWidth', 1.5);
hold on;

% Add lines for takeoff and landing
xline(t_force(foot_off_idx), 'r--', 'Foot Off', 'LineWidth', 1.5);
xline(t_force(foot_land_idx), 'b--', 'Foot Land', 'LineWidth', 1.5);

xlabel('Time (s)'); 
ylabel('Force (N)');
title('Force-Time Graph - SJ Trial 1');
legend('Vertical GRF', 'Foot Off', 'Foot Land');
grid on;

%% Plot all the parameters in one figure
% Make a figure with accel, vel and disp
figParams = figure('Name', 'Movement Parameters - SJ Trial 1');

% Accel subplot
subplot(2, 2, 1);
plot(t_force, acceleration, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Accel (m/s^2)');
title('Acceleration vs Time');
grid on;

% Velocity subplot
subplot(2, 2, 2);
plot(t_force, velocity, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Velocity (m/s)');
title('Velocity vs Time');
grid on;

% Displacement subplot
subplot(2, 2, 3);
plot(t_force, displacement_kin, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Displacement vs Time');
grid on;

% Combined plot with all three parameters
subplot(2, 2, 4);
plot(t_force, acceleration, 'r', 'LineWidth', 1.5); 
hold on;
plot(t_force, velocity, 'g', 'LineWidth', 1.5);
plot(t_force, displacement_kin, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Values');
title('All Parameters Together');
legend('Accel', 'Vel', 'Disp');
grid on;