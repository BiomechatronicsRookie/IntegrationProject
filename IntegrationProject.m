%% Integration Project - July 2021 
% Authors: VÃ­ctor Antonio Carmona Ortiz - SN 2623277
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
B0 = [dx*Km/R, -dx*Km/R; dy*Km/R, dy*Km/R]; % Mapping from voltage to torque 
C0 = [0.025, 0; 0, 0.025];
% C0_test= [0  0.01; 0.01  0  ]
% Measurement martix
D = [dx^2*Km^2/R, 0; 0, dy^2*Km^2/R];   % Damping matrix
Mass_supp = 0.0243;                                         % Mass of the support
J_s = 1/12*Mass_supp*(plate_dim(1)^2+plate_dim(3)^2);       % Support Inertia
M = [J_s+coil_m*dy^2*2, 0;0, J_s+coil_m*2*dx^2];             % Inertia Matrix
% M_correct= 0.1*M;

%% Q2 lets do some testing
A= [-inv(M)*D , -inv(M)*K ; eye(2), zeros(2)];
B= [inv(M)*B0; zeros(2)];
C= [zeros(2), C0];

s=tf('s');
tf_sys= C*inv(eye(4)*s - A)*B;
sys=ss(A,B,C,[]);


figure(2)
bodeplot(tf_sys)
grid on 

w_resonance= sqrt(inv(M)*K)
%% Q3

% two methods
[Ty,Tu]=wadyadicd(sys,1,500,1);


sys_dec = Ty*tf_sys*Tu;


bodeplot(sys_dec)


%% Q4

%From the assignment data
vel_max= 0.05; % m/s
acc_max= 2; %m/s^2
%Spot manipulation is 1/3 of spot diameter= 6 micrometers: 
error_max= (6*10^-6);

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
wc=( (jerk_max*2*(1/alfa))/error_max )^(1/3);

%Estimated plants
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
margin(L_1)
hold on
margin(T_1)
legend('OpenLoop(L1)','ClosedLoop(T1)')

figure()
title('Control Action K2 in G2')
margin(L_2)
hold on
margin(T_2)
legend('OpenLoop(L2)','ClosedLoop(T2)')

%Feedforward assuming perfect stifness means LF only interaction.

F1= G1.Denominator{1}(3)/(G1.Numerator{1}(3));
F2= G2.Denominator{1}(5)/(G2.Numerator{1}(5));

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

Ref1 = [t',x'];

%% Q6
%First lets construct Q filter:

%Lets take a butherworth filter of 2 degree.
s=tf('s');
tau=0.0005; %Double the cross-over frecuency of C(s)
Q=(2*tau*s+1)/(( (s*tau)^3+ 2*(s*tau)^2 + 2*s*tau)+1);

 figure()
 margin(Q) 
% 
% S_observer=feedback(1,Q);
 hold on
 bode(1-Q)
% bode(L_observer)

%Some norms to check stability:
norm(Q,inf) %Butherworth norm 
norm(1-Q,inf)


% G nominal coming from a dominant mass:
G_nominal= ( 1)/ (meq1*s^2);
G_nominal_inv=inv(G_nominal);


QPn=minreal(Q*G_nominal_inv);
Qminus=minreal(inv(1-Q));

Feed=minreal(Qminus*Q*s^2)

%PD controller
%Calculation 
w_controller=sqrt((acc_max*sqrt(1/alfa))/error_max);

kp_1=(meq1*w_controller^2)/(sqrt(1/alfa));
kp_2=(meq2*w_controller^2)/(sqrt(1/alfa));

t_z=sqrt(1/alfa)/w_controller; %zero of the lead controller
t_p= 1/(w_controller*sqrt(1/alfa)); %pole of the lead controller

C_1= kp_1*(s*t_z +1)/(s*t_p +1);


L_controller= G_nominal*C_1;
norm(QPn*Qminus*G1,inf)

T= minreal((L_controller+Q)/(L_controller+1))
norm(T,inf)

% S_controller= feedback(1,L_observer);
% T_controller= feedback(L_observer,1);
% figure()
% margin(Q)
% 
% S_observer=feedback(1,Q);
% hold on
% bode(1-Q)
% bode(L_observer)
% bode(S_controller)
% bode(T_controller)

% legend()
% figure()
% pzmap(Q)


%% Display of Simulink

sim('Model_Controlled_2019b.slx',max(t));
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



%% Q7
%Creation of data for SysIdent Toolbox
dataset1= iddata(y,u,ts);

test_order12=n4sid(dataset1,12);
%% Q9 H2 Mixed Sensitivity function:
dataset1= iddata(y,u,ts);

%% Reduction of System
%Identified System:


ss_reduced_or4= reduce(ss4,4);
% ss_reduced_or7= reduce(ss7,4);
ss_reduced_or6= reduce(ss4,6);

figure(10)
bodeplot(ss4)
hold on
% bodeplot(ss7)
% bodeplot(ss_reduced_or7)
bodeplot(ss_reduced_or4)
legend()

figure(11)
iopzmap(d2c(ss4))
hold on
iopzmap(d2c(ss_reduced_or4))
legend()
%Lets take the case to continous time:
ss_reduced_or4_ct= d2c(ss_reduced_or4);
% ss_reduced_or7_ct= d2c(ss_reduced_or7);

[Tyex9,Tuex9]=wadyadicd(ss_reduced_or4_ct,32.85,48.26,1);
% [Tyex91,Tuex91]=wadyadicd(ss_reduced_or7_ct,38.2,42,1);

sys_order4_decoupled = Tyex9*tf(ss_reduced_or4_ct)*Tuex9
% sys_order4_ss7decoupled = Tyex91*tf(ss_reduced_or7_ct)*Tuex91

sys_order4_decoupled_red= minreal(sys_order4_decoupled)
% sys_order4_ss7decoupled_red= minreal(sys_order4_ss7decoupled,0.0003)

%Verify the correct decoupling:
figure()
bodeplot(sys_order4_decoupled_red)

figure()
iopzmap(sys_order4_decoupled_red)
%% H2 controller:

G1=sys_order4_decoupled_9red(1,1);
G2=sys_order4_decoupled_red(2,2);

%Checking actual performance for the controller:
pole(G1)
zero(G1)

L=G1;
T=minreal((L)/(L+1));
S=minreal((1)/(L+1));

figure(92)
bodeplot(T)
hold on
bodeplot(S)
bodeplot(G1)
legend()

pole(G2)
zero(G2)

L=G2;
T=minreal((L)/(L+1));
S=minreal((1)/(L+1));

figure(93)
bodeplot(T)
hold on
bodeplot(S)
bodeplot(G2)
hold off
legend()



%% Controller G1 code:
G1.u = 'u2';
G1.y = 'y';
s=tf('s')
% W1 = 1/(s+0.00001); %Good one
%Try 1 s+0.001 100 veces mas alejado de s
alfa=10;%How far we are from the integrator
beta=0.01; %Gain added to sensitivy function:
delta= 0.01; %mas alto aumenta T, baja peak en S y t
W1 = (s+0.0001)/((s+alfa)); 
W1.u = 'y2';
W1.y = 'y11';

W2 = tf(beta); 
W2.u = 'u2';
W2.y = 'y12';

S = sumblk('y2 = u1 - y');
 
P = connect(G1,S,W1,W2,{'u1','u2'},{'y11','y12','y2'});
[K,CL,gamma] = h2syn(P,1,1);
gamma
K=minreal(tf(K)) % may have to test with tolerance
L=G1*K;
L_inf=norm(L,inf)
T=minreal((L)/(L+1));
S=minreal((1)/(L+1));
Q=minreal(-K/(L+1));
figure()
bodeplot(T)
hold on
bodeplot(S)
bodeplot(L)
bode(K)
legend()

str = sprintf('%d in 1/(s+alfa) and beta %d for W2(s)', alfa,beta);
title(str)

figure()
nichols(L)
title(str)
figure()
nyquist(L)
title(str)
figure()
step(T)
hold on
legend()
hold off
title(str)
bandwit=bandwidth(T)
T_inf=norm(T,inf)
S_inf=norm(S,inf)
K_inf=norm(tf(K),inf)

figure()
pzmap(K)

B_Controller=bandwidth(K)

figure()
sigma(Q)
%% Controller G2 code:

G2.u = 'u2';
G2.y = 'y';

W21 = 1/(s+0.00001); 
W21.u = 'y2';
W21.y = 'y11';

W22 = tf(0.0000001); 
W22.u = 'u2';
W22.y = 'y12';

S = sumblk('y2 = u1 - y');
 
P1 = connect(G2,S,W21,W22,{'u1','u2'},{'y11','y12','y2'});
[K2,CL,gamma] = h2syn(P1,1,1)

L2=G2*tf(K2);
T2=minreal((L2)/(L2+1))
S2=minreal((1)/(L2+1));
figure()
bodeplot(T2)
hold on
bodeplot(S2)
bodeplot(L2)

legend()
figure()
nichols(L2)

figure()
nyquist(L2)

figure()
step(T2)
hold off
bandwidth(T2)
norm(T,inf);
%%
s=tf('s');
W2=(s+0.0001)/(s+1); % trial and error: c Æ 0.1, r Æ 0
W1=1/(s+0.001); % standard H2 solution needs "+0.001"
V1=G1*(s+1);
G=[W1*V1 , W1*G1;
0 , W2; 
-V1 , -G1]
G=minreal(G);
[K_test,CL,gamma]=h2syn(G,1,1);
K_test=minreal(tf(K),1e-7); % may have to test with tolerance
L=G1*K_test;
T=L/(1+L);
bandwidth(T)

%%
s=tf('s');

disp('New design')
alfa=10; %Aumenta la velocidad del close-loop mas bandwith y mas rapido, pero L es mas grande
V1=((s+alfa)/(s+0.000001)^2);
W2=0.01*((s+0.001))/(s+alfa); % trial and error: c Æ 0.1, r Æ 0
W1=10; % standard H2 solution needs "+0.001"
G=[W1*V1 W1*G1;
0 W2;
-V1 -G1];
G=minreal(G);
[K,CL,gamma]=h2syn(G,1,1);
K=minreal(tf(K),1e-7); % may have to test with tolerance
gamma
L=minreal(G1*K);
L_inf=norm(L,inf)
T=minreal((L)/(L+1));
S=minreal((1)/(L+1));
Q=minreal(-K/(L+1));
figure()
bodeplot(T)
hold on
bodeplot(S)
bodeplot(L)
bodeplot(G1)
bodeplot(K)
legend()

str = sprintf('%d in 1/(s+alfa) and beta %d for W2(s)', alfa,beta);
title(str)

figure()
nichols(L)
title(str)
figure()
nyquist(L)
title(str)
figure()
step(T)
hold on
legend()
hold off
title(str)
bandwit=bandwidth(T)
T_inf=norm(T,inf)
S_inf=norm(S,inf)
K_inf=norm(tf(K),inf)

figure()
pzmap(K)

B_Controller=bandwidth(K)

figure()
sigma(Q)

%%

s=tf('s');

disp('New design')
E=zpk(pole(G1),[],1)
E=tf(E);
M=tf((s^4 + 2.613*s^3+3.414*s^2+2.613*s +1))
V=M/E;
bodeplot(1/V)
%%
alfa=10; %Aumenta la velocidad del close-loop mas bandwith y mas rapido, pero L es mas grande
V1=((s+alfa)/(s+0.000001)^2);
W2=0.01%*((s+0.001))/(s+alfa); % trial and error: c Æ 0.1, r Æ 0
W1=10; % standard H2 solution needs "+0.001"
G=[W1*V W1*G1;
0 W2;
-V -G1];
G=minreal(G);
[K,CL,gamma]=h2syn(G,1,1);
K=minreal(tf(K),1e-7); % may have to test with tolerance
gamma
L=minreal(G1*K);
L_inf=norm(L,inf)
T=minreal((L)/(L+1));
S=minreal((1)/(L+1));
Q=minreal(-K/(L+1));
figure()
bodeplot(T)
hold on
bodeplot(S)
bodeplot(L)
bodeplot(G1)
bodeplot(K)
legend()

str = sprintf('%d in 1/(s+alfa) and beta %d for W2(s)', alfa,beta);
title(str)

figure()
nichols(L)
title(str)
figure()
nyquist(L)
title(str)
figure()
step(T)
hold on
legend()
hold off
title(str)
bandwit=bandwidth(T)
T_inf=norm(T,inf)
S_inf=norm(S,inf)
K_inf=norm(tf(K),inf)

figure()
pzmap(K)

B_Controller=bandwidth(K)

figure()
sigma(Q)
