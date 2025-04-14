%% Countermovement Jump (CMJ) Trial 1

close all; 
clear; 
clc; 

%% Set constants and measurements

mass = 75.5;          % Mass in kg
gravity = 9.81;       % Gravity
% Leg lengths
L_LL = 0.840;         
R_LL = 0.840;         
% Get average leg length
Leg_L = (L_LL + R_LL) / 2; 
R_mark = 0.012;       % Marker radius

% Import data from Group 3 CMJ Trial 1 file
opts = detectImportOptions('Group 3 S CMJ 02 Trial 1  trajectories.csv', 'NumHeaderLines', 5);
traj = readmatrix('Group 3 S CMJ 02 Trial 1  trajectories.csv', opts);

opts = detectImportOptions('Group 3 S CMJ 02 Trial 1 force.csv', 'NumHeaderLines', 5);
force = readmatrix('Group 3 S CMJ 02 Trial 1 force.csv', opts);
%% 

% Process data - still working on understanding all this
BW1 = mean(force(force(:,2)==0, 5)); % weight
f_max1 = size(traj,1);
Freq1 = 100; % Hz frequency 
t1 = (0:f_max1-1)/Freq1; % time 

% Get force data - not sure if this is the best way but it works
ReSamp1 = find(force(:,2) == 0); 
combined_normal1 = -force(ReSamp1,5); % vertical force (negative because of direction)
t_force1 = (0:length(combined_normal1)-1)/1000; % Convert from 1000Hz to right time scale
%% Davis Hip Model

% Get marker positions
LASI = traj(:,72:74); % left anterior
RASI = traj(:,75:77); % right anterior
LPSI = traj(:,78:80); % left posterior
RPSI = traj(:,81:83); % right posterior

% Gold Standard Method - not sure about all these calculations but they came from the lab manual
pelvic_w1 = LASI - RASI;
ASIS_Dist1 = mean(sqrt(sum(pelvic_w1.^2,2))) / 1000; % need to divide by 1000 to convert to m

% Calculations from lab manual
X_dist1 = 0.1288 * Leg_L - 0.04856;
C1 = Leg_L * 0.115 - 0.0153;
B1 = deg2rad(18); 
T1 = deg2rad(28.4);

% hip joint center calculations 
LHJC_L1 = [(C1*sin(T1) - 0.5*ASIS_Dist1); (-X_dist1 - R_mark)*cos(B1) + (C1*cos(T1)*sin(B1)); (-X_dist1 - R_mark)*sin(B1) - (C1*cos(T1)*cos(B1))] * 1000;
RHJC_L1 = [- (C1*sin(T1) - 0.5*ASIS_Dist1); (-X_dist1 - R_mark)*cos(B1) + (C1*cos(T1)*sin(B1)); (-X_dist1 - R_mark)*sin(B1) - (C1*cos(T1)*cos(B1))] * 1000;

% pelvis calculations - need to check if this is correct?
pelvis_org1 = (RASI + LASI) / 2;
LHJC1 = zeros(3, f_max1); 
RHJC1 = zeros(3, f_max1);

% Loop through frames - this bit is confusing
for f = 1:f_max1
    % get vectors for each frame
    i1 = (RASI(f,:).' - pelvis_org1(f,:).') / norm(RASI(f,:) - pelvis_org1(f,:));
    L21 = pelvis_org1(f,:).' - RPSI(f,:).'; 
    L21 = L21 / norm(L21);
    k1 = cross(i1, L21); 
    j1 = cross(k1, i1); 
    GRL1 = [i1, j1, k1];
    
    % Calculate hip joint centers
    LHJC1(:,f) = GRL1 * LHJC_L1 + pelvis_org1(f,:).';
    RHJC1(:,f) = GRL1 * RHJC_L1 + pelvis_org1(f,:).';
end

% Find COM position
COM1 = (LHJC1 + RHJC1) / 2;
COM_Z1 = COM1(3,:) / 1000; % convert to meters

% Get initial height from first second of data
first_second_idx1 = find(t1 <= 1);
initial_height_gold_standard1 = mean(COM_Z1(first_second_idx1));

% Kinematic Method - this is from class notes
net_force1 = combined_normal1 - BW1; % subtract bodyweight
acceleration1 = net_force1 / mass; % F=ma
velocity1 = cumtrapz(t_force1, acceleration1); % integrate accel to get velocity

% Find foot off and landing times - when force drops below 50N
foot_off_idx1 = find(combined_normal1 < 50, 1, 'first');
foot_land_idx1 = find(combined_normal1 < 50, 1, 'last');

% Find max velocity
max_velocity1 = max(velocity1(1:foot_off_idx1));

% Calculate displacement by integrating velocity
displacement_kin1 = cumtrapz(t_force1, velocity1) + initial_height_gold_standard1;

% Projectile Motion Analysis - this is from biomechanics formulas
flight_time1 = t_force1(foot_land_idx1) - t_force1(foot_off_idx1);
mid_flight_time1 = t_force1(foot_off_idx1) + flight_time1/2; % not actually used but kept it
proj_t1 = t_force1(t_force1 >= t_force1(foot_off_idx1) & t_force1 <= t_force1(foot_land_idx1));
proj_t_rel1 = proj_t1 - t_force1(foot_off_idx1);
proj_disp1 = max_velocity1 .* proj_t_rel1 - 0.5 * gravity .* proj_t_rel1.^2 + initial_height_gold_standard1;

% Energy Method - simple way to find max height
energy_height1 = (max_velocity1^2)/(2*gravity) + initial_height_gold_standard1;

%% Plotting for CMJ Trial 1
% Make first figure with 4 different method comparisons
fig1 = figure('Name', 'Height Methods Comparison - Trial 1');

% Plot Gold Standard
subplot(2, 2, 1);
plot(t1, COM_Z1, 'k', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Gold Standard Method');
grid on;

% Plot Kinematic Method 
subplot(2, 2, 2);
plot(t_force1, displacement_kin1, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Kinematic Method');
grid on;

% Plot Projectile Motion
subplot(2, 2, 3);
plot(proj_t1, proj_disp1, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Projectile Motion');
grid on;

% Plot Energy Method - just a horizontal line
subplot(2, 2, 4);
plot(t_force1, ones(size(t_force1)) * energy_height1, 'm--', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Energy Method');
grid on;

%% Force time graph showing phases
figure('Name', 'Force-Time Graph - CMJ Trial 1');
plot(t_force1, combined_normal1, 'k', 'LineWidth', 1.5);
hold on;

% Add vertical lines to show foot off and landing
xline(t_force1(foot_off_idx1), 'r--', 'Foot Off', 'LineWidth', 1.5);
xline(t_force1(foot_land_idx1), 'b--', 'Foot Land', 'LineWidth', 1.5);

xlabel('Time (s)'); 
ylabel('Force (N)');
title('Force-Time Graph with CMJ Phases - Trial 1');
legend('Vertical GRF', 'Foot Off', 'Foot Land');
grid on;

%% Plot all the motion parameters together
% Make a figure with acceleration, velocity and displacement
fig3 = figure('Name', 'Motion Parameters - CMJ Trial 1');

% Acceleration plot
subplot(2, 2, 1);
plot(t_force1, acceleration1, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Accel (m/s^2)');
title('Acceleration vs Time');
grid on;

% Velocity plot
subplot(2, 2, 2);
plot(t_force1, velocity1, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Velocity (m/s)');
title('Velocity vs Time');
grid on;

% Displacement plot
subplot(2, 2, 3);
plot(t_force1, displacement_kin1, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Height (m)');
title('Displacement vs Time');
grid on;

% Combined plot with all three
subplot(2, 2, 4);
plot(t_force1, acceleration1, 'r', 'LineWidth', 1.5); 
hold on;
plot(t_force1, velocity1, 'g', 'LineWidth', 1.5);
plot(t_force1, displacement_kin1, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('Values');
title('All Parameters vs Time');
legend('Accel', 'Vel', 'Disp');
grid on;