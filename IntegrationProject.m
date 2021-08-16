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
plot(y(:,1),'linewidth',1.5)
hold on 
plot(Output(:,1),'linewidth',1.5)
ylabel('Input []')
title('Sensor 1')
subplot(2,1,2)
plot(y(:,2),'linewidth',1.5)
hold on 
plot(Output(:,2),'linewidth',1.5)
xlabel('Time [s]')
ylabel('Output []') 
title('Sensor 2')
legend('Data','Sim')

%% Coupling analysis
H = tf(sys);
H = minreal(H);
fq = evalfr(H,40);
[U,S,V] = svd(fq);
H_dec = U'*H*V;
H_dec = minreal(H_dec);
figure()
bodeplot(H_dec);

G1 = H_dec(1,1);
G2 = H_dec(2,2);

%% Decoupling 
% [Ty,Tu] = wadyadicd(sys,1,500,1);
% sys_dec = Ty*tf(sys)*Tu;
% sys_dec = minreal(sys_dec);
% 
% figure()
% bodeplot(sys_dec)
% G1 = sys_dec(1,1);
% G2 = sys_dec(2,2);

%% Profiles 
p = 2;
dt = 2e-4;
t = 0:dt:0.3;
a = zeros(size(t));
n = p/0.025;
a(1:126) = t(1:126).*n;
a(126:251) = p-t(1:126).*n;
a(1001:1126) = -t(1:126).*n;
a(1126:1251) = -p + t(1:126).*n;
v = zeros(size(a));
x = zeros(size(a));
jrk = zeros(size(a));

for i = 2:length(a)
    v(i) = v(i-1)+a(i-1)*dt;
    x(i) = x(i-1)+v(i-1)*dt;
    jrk(i) = (a(i)-a(i-1))/dt;
end

figure()
subplot(4,1,1)
plot(t,x,'linewidth',1.5)
subplot(4,1,2) 
plot(t,v,'linewidth',1.5)
subplot(4,1,3) 
plot(t,a,'linewidth',1.5)
subplot(4,1,4) 
plot(t,jrk,'linewidth',1.5)

%% Parallel PID controller
%From the assignment data
vel_max= 0.05; % m/s
acc_max= 2; %m/s^2
%Spot manipulation is 1/3 of spot diameter= 6 micrometers: 
error_max= 6e-6;

%The motion profile is Third order reference:
% Vel_max= 2*h_m/tm
% Acc_max= 8*h_m/tm^3
%So h_m and t_m are:
h_m= vel_max^2; %Also from the assignment is equal to 2.5 mm
t_m= 2*vel_max;
jerk_max= (32*h_m)/(t_m^3); % jerk_max= error_max in PID controller


%Obtained from Paper 3 and ref from CDMS
    
meq1 = 1/G1.Numerator{1}(3);                        % equivalent mass sys 1
meq2 = 1/G2.Numerator{1}(3);                        % equivalent mass sys21
beta = 1;                                           % Tameness factor
alpha = 0.1;                                        % lead factor
wc = ((jerk_max*2*(1/alpha))/error_max )^(1/3);     % rad/s - crossover frequency for the controler  
tz = sqrt((1/alpha))/wc;                            % zero place
ti = beta*tz;                                       % Integral place
tp = 1/(wc*sqrt(1/alpha));                          % pole place
kp1 = meq1 * wc*wc/sqrt(1/alpha);                   % gain 1
kp2 = meq2 * wc*wc/sqrt(1/alpha);                   % gain 2

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
F2 = G2.Denominator{1}(3)/(G2.Numerator{1}(3));

%% Predictive error
Kj = beta/(alpha*wc^3);
Ka1 = beta*(D0(1,1)/M(1,1))/(alpha*wc^3);
Ka2 = beta*(D0(2,2)/M(2,2))/(alpha*wc^3);
Kv1 = 0;
Kv2 = 0;
err_1 = Kj*jrk + Ka1*a + Kv1*v;
err_2 = Kj*jrk + Ka2*a + Kv2*v;

figure()
subplot(2,1,1)
plot(t,err_1,'linewidth',1.5)
ylabel('Error [m]')
subplot(2,1,2)
plot(t,err_2,'linewidth',1.5)
ylabel('Error [m]')
xlabel('Time [s]')


%% To simulink
Ref1 = [t',x'];
Ref2 = [t',x'];

sim('Model_Controlled.slx',max(t));
figure()
subplot(3,1,1)
plot(t,x,'LineWidth',1.5)
ylabel('Amplitude [m]')
title('Reference Signal')
subplot(3,1,2)
plot(t,Controlled1,'LineWidth',1.5)
ylabel('Amplitude [m]')
title('System Output')
subplot(3,1,3)
plot(t,x'-Controlled1,'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Amplitude [m]')
title('Error')

figure()
subplot(3,1,1)
plot(t,Ref2(:,2),'LineWidth',1.5)
ylabel('Amplitude [m]')
title('Reference Signal')
subplot(3,1,2)
plot(t,Controlled2,'LineWidth',1.5)
ylabel('Amplitude [m]')
title('System Output')
subplot(3,1,3)
plot(t,x'-Controlled2,'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Amplitude [m]')
title('Error')

%% Comparison between simulated and predicted errors
figure()
subplot(2,1,1)
plot(t,x'-Controlled1,'linewidth',1.5)
hold on 
plot(t,err_1,'linewidth',1.5)
ylabel('Error [m]')
title('X plane')
subplot(2,1,2)
plot(t,x'-Controlled2,'linewidth',1.5)
hold on 
plot(t,err_2,'linewidth',1.5)
title('Y Plane')
legend('Simulatd','Predicted')
xlabel('Time [s]')
ylabel('Error [m]')


%% Controler discretisation 
C1_d = c2d(Controller1,ts);
C2_d = c2d(Controller2,ts);

load('Order12.mat');
% Use system identification to identify the plant in discrete time
DataStructure = iddata(y,u,ts);
%  systemIdentification
figure()
bodeplot(tf_sys)
hold on
bodeplot(Order12)

%% Discrete system decoupling
tf_sys_id = tf(Order12);
tf_sys_id = minreal(tf_sys_id);

[Ty,Tu] = wadyadicd(d2c(tf(Order12)),49,33,0);
H_DYAD_dec = Ty*d2c(tf(Order12))*Tu;
H_DYAD_dec = c2d(minreal(H_DYAD_dec),ts);

fq = evalfr(tf_sys_id,49);
[U,S,V] = svd(fq);
H_SVD_dec = U'*H*V;
H_SVD_dec = minreal(H_SVD_dec);
figure()
bodeplot(H_SVD_dec);

G1_d = H_SVD_dec(1,1);
G2_d = H_SVD_dec(2,2);

figure()
bodeplot(H)
hold on
bodeplot(tf_sys_id)
bodeplot(H_SVD_dec)
bodeplot(H_DYAD_dec)
legend('Continous','Discrete','SVD Dec.','Dyad. Dec.')

%% Discrete model
sim('Discrete_Model_Controlled.slx',0.3)
% Results plots
figure()
subplot(3,1,1)
plot(t,Ref2(:,2))
title('Reference Signal')
ylabel('Amplitude [m]')
subplot(3,1,2)
plot(Discrete_Res(:,1))
ylim([-0.01 0.01])
ylabel('Amplitude [m]')
subplot(3,1,3)
plot(Discrete_Res(:,2))
ylim([-0.01 0.01])
ylabel('Amplitude [m]')

% No stability, no robust stability(?)

%% Check characteristic loci
figure()
pzmap(Order12) % Two poles on the right half plane --> unstable on its own



