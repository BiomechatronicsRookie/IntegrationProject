%% Integration Project - July 2021 
% Authors: Víctor Antonio Carmona Ortiz - SN 2623277
%          Alejandro Lopez Tellez - SN 2586266
%% Parameter Definition
%load('iddatamir.mat')
%time = (0:1:length(u)-1)*ts;
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
B0 = [dx^2*Km/R, -dx^2*Km/R; dy^2*Km/R, dy^2*Km/R]; % Mapping from voltage to torque 
C0 = [0.01, 0; 0, 0.01];
C0_test= [0  0.01; 0.01  0  ]
% Measurement martix
D = [dx^2*Km^2/R, 0; 0, dy^2*Km^2/R];   % Damping matrix
Mass_supp = 0.0243;                                         % Mass of the support
J_s = 1/12*Mass_supp*(plate_dim(1)^2+plate_dim(3)^2);       % Support Inertia
M = [J_s+coil_m*dy^2*2, 0;0, J_s+coil_m*2*dx^2];             % Inertia Matrix
M_correct= 0.1*M;

%% Q2 lets do some testing
A= [-inv(M)*D , -inv(M)*K ; eye(2), zeros(2)];
B= [inv(M)*B0; zeros(2)];
C= [zeros(2), C0];

C_test=[C0_test zeros(2)];
A_test2= [ zeros(2) eye(2) ;  -inv(M)*K  -inv(M)*D];
B_test2= [zeros(2) ; inv(M_correct)*B0]; 
A_sy= [-inv(M_correct)*D , -inv(M_correct)*K ; eye(2), zeros(2)];
B_sy= [inv(M_correct)*B0; zeros(2)];

A_test=10*A; 
B_test= 10*B;

s=tf('s');
% tf_test= C_test*inv(eye(4)*s - A_test2)*B_test2;
% tf_sys_sy= C*inv(eye(4)*s - A_sy)*B_sy;
tf_sys= C*inv(eye(4)*s - A)*B;
% tf_sys_test= C*inv(s*eye(4) - A_test)*B_test;
sys=ss(A,B,C,[]);

% figure(1)
% bodeplot(tf_sys_sy)
figure(2)
bodeplot(tf_sys)
% figure(3)
% bodeplot(tf_sys_test)

%% Q3

% two methods
[Ty,Tu]=wadyadicd(sys,1,500,1);


sys_dec = Ty*tf_sys*Tu;

figure()
bodeplot(sys_dec)
% A_decoupled= [Ty, zeros(2); zeros(2), Tu]*A*;
% B_decoupled= Ty*D*Tu;
% C_decoupled= Ty*C*Tu;

%% Q4

%From the assignment data
vel_max= 0.05; % m/s
acc_max= 2; %m/s^2
%Spot manipulation is 1/3 of spot diameter= 6 micrometers: 
error_max= 6*10^-6;

%The motion profile is Third order reference:
% Vel_max= 2*h_m/tm
% Acc_max= 8*h_m/tm^3
%So h_m and t_m are:
h_m= vel_max^2; %Also from the assignment is equal to 2.5 mm
t_m= 2*vel_max;

jerk_max= (32*h_m)/(t_m^3); % jerk_max= error_max in PID controller

%Given parameters for controller 
beta=1; %integral action
alfa=0.1; % phase lead equal, distance between pole and zeros in the lead phase controller

%Obtained from Paper 3 and ref from CDMS
wc=( (jerk_max*2*(1/alfa))/error_max )^(1/3)

%Nominals plants
G1= minreal(sys_dec(1,1));
G2= minreal(sys_dec(2,2));

%From paper 3, the equivalent mass of the plant is equal to the inverse of
%the gain in HF(w=inf)
meq1=1/(G1.Numerator{1}(3));
meq2=1/(G2.Numerator{1}(3));

%Parameters from equation 24 of the assignment:

kp_1=(meq1*wc^2)/(sqrt(1/alfa));
kp_2=(meq2*wc^2)/(sqrt(1/alfa));

t_z=sqrt(1/alfa)/wc; %zero of the lead controller
t_p= 1/(wc*sqrt(1/alfa)); %pole of the lead controller
t_i=beta*t_z; %integrator pole


K_1= kp_1*( (s*t_z +1)*(s*t_i +1) ) / (s*t_i*(s*t_p +1) );
K_2= kp_2*( (s*t_z +1)*(s*t_i +1) ) / (s*t_i*(s*t_p +1) );

L_1=minreal(minreal(sys_dec(1,1))*K_1);
T_1=minreal(feedback(L_1,1));

L_2=minreal(minreal(sys_dec(2,2))*K_2);
T_2=minreal(feedback(L_2,1));

figure()
title('Control Action K1 in G1')
bodeplot(L_1)
hold on
bodeplot(T_1)
legend('OpenLoop(L1)','ClosedLoop(T1)')

figure()
title('Control Action K2 in G2')
bodeplot(L_2)
hold on
bodeplot(T_2)
legend('OpenLoop(L2)','ClosedLoop(T2)')

%% Q7
%Creation of data for SysIdent Toolbox
data2= iddata(y,10*u,ts);
