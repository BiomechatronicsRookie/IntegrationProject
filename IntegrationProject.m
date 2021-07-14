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
B0 = [dx*dx*2 * Km/R, -dx*dx*Km/R; dy*dy*Km/R, dy*dy*Km/R];       % Mapping from voltage to torque 
C0 = [0.01, 0; 0, 0.01];                                    % Measurement martix
D0 = [dx*dx*Km^2/R, dx*dx*Km^2/R; dy*dy*Km^2/R, dy*dy*Km^2/R];  % Damping matrix
Mass_supp = 0.0243;                                         % Mass of the support
J_s = 1/12*Mass_supp*(plate_dim(1)^2+plate_dim(3)^2);       % Support Inertia
M = [J_s+coil_m*dy*dy*2, 0;0, J_s+coil_m*2*dx*dx];            % Inertia Matrix

%% STATE SPACE
A = [-M\D0, -M\K;
     1 0 0 0;
     0 1 0 0];
B = [-M\B0; zeros(2)];
C = [0 0 0.01 0;0 0 0 0.01];
D = zeros(2);

Input = [time',u];
sim model.slx;
%% 
figure()
subplot(2,1,1)
plot(time,y(:,1),'linewidth',1.5)
hold on 
plot(time,out.Output(:,1),'linewidth',1.5)
ylabel('Input []')
title('Output Comparison')
subplot(2,1,2)
plot(time,y(:,2),'linewidth',1.5)
hold on 
plot(time,out.Output(:,2),'linewidth',1.5)
xlabel('Time [s]')
ylabel('Output []') 
legend('sensor1','sensor2')


