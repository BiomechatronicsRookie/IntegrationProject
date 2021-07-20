%% Integration Project - July 2021 
% Authors: Víctor Antonio Carmona Ortiz - SN 2623277
%          Alejandro Lopez Tellez - SN 2586266
clear all;close all;clc

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
dx = 0.03;                           % Distance moment in Y
dy = 0.035;                          % Distance moment in X

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
B0 = [dx*Km/R , -dx*Km/R ; dy*Km/R, dy*Km/R];               % Mapping from voltage to torque 
C0 = [0.01, 0; 0, 0.01];                                    % Measurement martix
D0 = [dx*dx*Km^2/R, 0; 0, dy*dy*Km^2/R];                    % Damping matrix
Mass_supp = 0.0243;                                         % Mass of the support
J_s = 1/12*Mass_supp*(plate_dim(1)^2+plate_dim(3)^2);       % Support Inertia
M = [J_s+coil_m*dy*dy*2, 0;0, J_s+coil_m*2*dx*dx];          % Inertia Matrix

%% STATE SPACE
A = [-M\D0, -M\K;
     1 0 0 0;
     0 1 0 0];
B = [M\B0; zeros(2)];
C = [0 0 0.01 0;0 0 0 0.01];
D = zeros(2);

s = tf('s');
tf_sys = C*inv(eye(4)*s - A)*B;

sys = ss(A,B,C,D);
figure()
bodeplot(sys)

Input = [time',u];
sim('model.slx',max(time));
%% 
figure()
subplot(2,1,1)
plot(time,y(:,1),'linewidth',1.5)
hold on 
plot(time,Output(:,1),'linewidth',1.5)
ylabel('Input []')
title('Sensor 1')
subplot(2,1,2)
plot(time,y(:,2),'linewidth',1.5)
hold on 
plot(time,Output(:,2),'linewidth',1.5)
xlabel('Time [s]')
ylabel('Output []') 
title('Sensor 2')
legend('Data','Sim')

%% Decoupling 
[Ty,Tu] = wadyadicd(sys,1,500,1);
sys_dec = Ty*tf(sys)*Tu;
sys_dec = minreal(sys_dec);

figure()
bodeplot(sys_dec)
G1 = sys_dec(1,1);
G2 = sys_dec(2,2);

%% Profiles 
u = 2;
dt = 2e-4;
t = 0:dt:0.3;
a = zeros(size(t));
n = u/0.025;
a(1:126) = t(1:126).*n;
a(126:251) = u-t(1:126).*n;
a(1001:1126) = -t(1:126).*n;
a(1126:1251) = -u + t(1:126).*n;
v = zeros(size(a));
x = zeros(size(a));

for i = 2:length(a)
    v(i) = v(i-1)+a(i-1)*dt;
    x(i) = x(i-1)+v(i-1)*dt;
end

figure()
subplot(3,1,1)
plot(t,x,'linewidth',1.5)
subplot(3,1,2) 
plot(t,v,'linewidth',1.5)
subplot(3,1,3) 
plot(t,a,'linewidth',1.5)


%% Parallel PID controller
meq1 = 1/G1.Numerator{1}(3);        % equivalent mass sys 1
meq2 = 1/G2.Numerator{1}(3);        % equivalent mass sys21
wc = 643.6596;                      % rad/s - crossover frequency for the controler
beta = 1;                           % Tameness factor
alpha = 0.1;                        % lead factor
tz = sqrt((1/alpha))/wc;            % zero place
ti = beta*tz;                       % Integral place
tp = 1/(wc*sqrt(1/alpha));          % pole place
kp1 = meq1 * wc*wc/sqrt(1/alpha);   % gain 1
kp2 = meq2 * wc*wc/sqrt(1/alpha);   % gain 2

s= tf('s');
Controller1 = kp1*(s*tz+1)*(s*ti+1)/(s*ti*(s*tp+1));
Controller2 = kp2*(s*tz+1)*(s*ti+1)/(s*ti*(s*tp+1));

figure()
margin(G1*Controller1)
grid on

figure()
margin(G2*Controller2)
grid on

F1 = G1.Denominator{1}(3)/(G1.Numerator{1}(3));
F2 = G2.Denominator{1}(5)/(G2.Numerator{1}(5));

%% To simulink
Ref1 = [t',x'];
Ref2 = [t',zeros(size(x))'];

sim('Model_Controlled.slx',max(t));
figure()
subplot(3,1,1)
plot(t,x,'LineWidth',1.5)
ylabel('Amplitude [m]')
title('Reference Signal')
subplot(3,1,2)
plot(Controlled1.Time,Controlled1.Data,'LineWidth',1.5)
ylabel('Amplitude [m]')
title('System Output')
subplot(3,1,3)
plot(Controlled1.Time,Controlled1.Data'-x,'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Amplitude [m]')
title('Error')


figure()
subplot(3,1,1)
plot(t,Ref2(:,2),'LineWidth',1.5)
ylabel('Amplitude [m]')
title('Reference Signal')
subplot(3,1,2)
plot(Controlled2.Time,Controlled2.Data,'LineWidth',1.5)
ylabel('Amplitude [m]')
title('System Output')
subplot(3,1,3)
plot(Controlled2.Time,Controlled2.Data-Ref2(:,2),'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Amplitude [m]')
title('Error')


