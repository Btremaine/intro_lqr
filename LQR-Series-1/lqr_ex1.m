% lqr_param.
%
% Example using lqr gains with simple moving mass
clear

% start in continuous time, use MKS units
% states: x1=pos, x2=vel, w1=bias
% calculate ss matricies from parameters, then simulate as ss models.
% assume no external bias or offsets

M1 = 10;   % mass kg
kv = 0.1;  % viscous damping
kf = 2.0;  % force constant

ref = 10;  % reference position
w = 2.0;   % external bias
Nx = [1; 0];

%% c.t. plant model with no bias
% two states: pos, vel 
A = [0  1
     0 -kv/M1];
B = [0 
    kf];
C = [1 0];
D = 0;
sys_plant=ss(A,B,C,D);

%% discrete plant
% two states: pos, vel
Ts= 0.25;
sysd_plant= c2d(sys_plant,Ts);
Ad = sysd_plant.A;
Bd = sysd_plant.B;
Cd = sysd_plant.C;
Dd = sysd_plant.D;

co = ctrb(Ad,Bd);

%% LQR gain for feedback, two states
Qf = [0.1   0
      0     2.0];
Rf = 1.0;
[Kgain] = lqr(sysd_plant,Qf,Rf);

%% LQ time-vartying case
N = 80;         % final time step
S= Qf;          % adjoint boundary condition
Gain = [];
for k= 1:N
    P = Rf + Bd'*S*Bd;
    M = S - S*Bd*inv(P)*Bd'*S;
    K = inv(P)*Bd'*S*Ad;
    Gain = [Gain; K];
    S = Ad'*M*Ad + Qf;
end
Gain = flip(Gain);
figure(1); plot(Gain,'x-'); title('Time-optimal gains')
xlabel('time step'); ylabel('vel & pos gains')



