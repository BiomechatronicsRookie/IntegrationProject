%% Integration Project - July 2021 
% Authors: Víctor Antonio Carmona Ortiz - SN 2623277
%          Alejandro Lopez Tellez - SN 2586266
%% Parameter Definition
load('iddatamir.mat')
time = (0:1:length(u)-1)*ts;
Km = 0.3;                           % Motor Constant [N/A]
R = 0.35;                           % Armature Resistance [Ohm]
Motor_d = 0.022;                    % Motor diameter [m]
coil_m = 0.036;                     % Moving Mass of coil [kg]
plate_dim = [0.03 0.03 0.01];       % Mirror support dimensions {length,width,height} [m, m, m]
mirror_thick = 0.0003;              % Mirror thikness [m]
dist_between_motors = 0.07;         % Heart to Heart motors [m]
dist_moth_mirr = 0.03;              % Heart motors to mirror [m]
w_spr_dim = [0.045, 0.0006];        % Wire Spring Dimensions {length,thickness} [m, m]
lf_spr_dim = [0.05, 0.02, 0.0001];  % Leafspring dimesnions {length,width,height} [m, m, m]
dx = 0.3;                           % Distance moment in Y
dy = 0.35;                          % Distance moment in X

%% Plot given data
figure()
subplot(2,1,1)
plot(time,u,'linewidth',1.5)
ylabel('Input []')
title('Input and Output')
subplot(2,1,2)
plot(time,y,'linewidth',1.5)
xlabel('Time [s]')
ylabel('Output []')
legend('x','y')

%% Q1 
K = [0.1389, 0; 0, 0.2731];                                 % Stiffnes matrix 
B0 = [dx^2 * Km/R, -dx^2*Km/R; dy^2*Km/R, dy^2*Km/R];       % Mapping from voltage to torque 
C0 = [0.01, 0; 0, 0.01];                                    % Measurement martix
D = [dx^2*Km^2/R, dx^2*Km^2/R; dy^2*Km^2/R, dy^2*Km^2/R];   % Damping matrix
Mass_supp = 0.0243;                                         % Mass of the support
J_s = 1/12*Mass_supp*(plate_dim(1)^2+plate_dim(3)^2);       % Support Inertia
M = [J_s+coil_m*dy^2*2, 0;0, J_s+coil_m*2*dx^2];                % Inertia Matrix

%ppeneeeeeeeeeeeeeeeeee


