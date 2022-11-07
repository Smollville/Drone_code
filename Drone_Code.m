% Drone

clear all, close all, clc

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;

%% 1.1.8 2-D Model

clear X Y
    
g = 9.8065; % gravity force [m/s^2]
beta = 1/2*1.225*0.3*0.1; % 1/2*p*Cd*A where p is air density, Cd is drag coefficient, A is frontal Area
v1e = 3;
v1e_x = 1; % speed of the follower drone
v1e_y = sqrt(v1e^2-v1e_x^2);
v0e_x = v1e_x ; % speed of the leader drone
v0e_y = v1e_y;
p10 = 0.5;
p10_x = 0.1;
p10_y = sqrt(p10^2-p10_x^2);
m = 0.1; % mass in kg
we_x = 0; % typ. wind speed [m/s]
we_y = 0;
theta1e = atan((beta/(m*g))*(v0e_x-we_x)^2); % pitch angle in rad
theta1_e = rad2deg(theta1e),
phi1e = atan(-(beta/(m*g))*cos(theta1e)*(v0e_y-we_y)^2);% roll angle in rad
phi1_e = rad2deg(phi1e)
Al = [0 0 -1 0; 0 0 0 -1; 0 0 -2*beta*v1e_x/m 0;0 0 0 -2*beta*v1e_y/m];
Bl = [0 0; 0 0;(g*sec(theta1e)^2) 0;(-g*sec(theta1e)*tan(theta1e)*tan(phi1e)) (-g*sec(phi1e)^2*sec(theta1e))];
El = [0 0 0 1;0 0 0 1; 2*beta*we_x/m 0 0 0; 0 2*beta*we_y/m 0 0];
Cl = eye(4); % Cl = [1 0 0 0; 0 1 0 0] but we wanted to access speed
Dl = 0;

% Discrete-time state-space model approximation

Ts = 0.1; % sampling time [s]

Ad = [1 0 -Ts 0; 0 1 0 -Ts; 0 0 1-2*beta*v1e_x*Ts/m 0;0 0 0 1-2*beta*v1e_y*Ts/m];
Bd = [0 0; 0 0;Ts*(g*sec(theta1e)^2) 0;Ts*(-g*sec(theta1e)*tan(theta1e)*tan(phi1e)) Ts*(-g*sec(phi1e)^2*sec(theta1e))];
Ed = [0 0 0 Ts;0 0 0 Ts; 2*beta*we_x*Ts/m 0 0 0; 0 2*beta*we_y*Ts/m 0 0];
Cd = eye(4);  % Cd = [1 0 0 0; 0 1 0 0] but we wanted to access speed
Dd = 0;

%  Stability, controllability, and observability of the discrete time model.

[eigenvectors,eigenvalues] = eig(Ad); % stability -> eigenvalues <= 1
[similarity_transform,J] = jordan(Ad);
disp('System marginally stable because eigenvectors <= 1');
disp(eigenvalues);
disp('and all jordan blocks associated with eigenvalues = 1 are 1x1 blocks');
disp(J);

Co = ctrb(Ad,Bd); % controlable - rank of controlability matrix equals n = 2
if rank(Co) == size(Ad,1)
   disp('Controlable') 
else
    disp('Not Controlable') 
end    

Ob = obsv(Ad,Cd); % observable - rank of observability matrix equals n = 2
if rank(Ob) == size(Ad,1)
   disp('Observable') 
else
    disp('Not Observable') 
end   

% Simulate system response

% Continuous response
sys_c = ss(Al,Bl,Cl,Dl); 

tmax = 5; % duration of the simulation
t = 0:Ts:tmax;
u = 0.1*ones(length(t),2);
x0 = [0 0 0 0]'; % [p10x p10y v1x v1y]'
y_c = lsim(sys_c,u,t,x0);
y_c(:,1) = y_c(:,1) + p10_x;
y_c(:,2) = y_c(:,2) + p10_y;
y_c(:,3) = y_c(:,3) + v1e_x;
y_c(:,4) = y_c(:,4) + v1e_y;

% Discrete response
Ts = 0.1;

sys_d = ss(Ad,Bd,Cd,Dd,Ts);
y_d = lsim(sys_d,u,t,x0);

figure('Name','Pole Zero Map - 2-D Linearized and Discretized System','NumberTitle','off');
set(gcf,'defaultTextInterpreter','tex');
pzplot(sys_d,'b');
grid on;
axis equal;
% change marker size
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
    set(a(i),'color',sstblue) %change marker size
end
title('Pole Zero Map');

y_d(:,1) = y_d(:,1) + p10_x;
y_d(:,2) = y_d(:,2) + p10_y;
y_d(:,3) = y_d(:,3) + v1e_x;
y_d(:,4) = y_d(:,4) + v1e_y;

% sys_d_z = c2d(sys_c, Ts, 'zoh');
% y_d_z = lsim(sys_d_z,u,t,x0);
% y_d_z(:,1) = y_d_z(:,1) + p10_x;
% y_d_z(:,2) = y_d_z(:,2) + p10_y;
% y_d_z(:,3) = y_d_z(:,3) + v1e_x;
% y_d_z(:,4) = y_d_z(:,4) + v1e_y;
% 
% sys_d_t = c2d(sys_c, Ts, 'tustin');
% y_d_t = lsim(sys_d_t,u,t,x0);
% y_d_t(:,1) = y_d_t(:,1) + p10_x;
% y_d_t(:,2) = y_d_t(:,2) + p10_y;
% y_d_t(:,3) = y_d_t(:,3) + v1e_x;
% y_d_t(:,4) = y_d_t(:,4) + v1e_y;
Ny = length(y_d(:,1));
y_c_total = zeros(Ny,2);
y_d_total = zeros(Ny,2);
for i = 1:size(y_d(:,1))
    y_c_total(i,1) = sqrt(y_c(i,1)^2 + y_c(i,2)^2);
    y_c_total(i,2) = sqrt(y_c(i,3)^2 + y_c(i,4)^2);
    y_d_total(i,1) = sqrt(y_d(i,1)^2 + y_d(i,2)^2);
    y_d_total(i,2) = sqrt(y_d(i,3)^2 + y_d(i,4)^2);
end

figure('Name','Phase Plot - 2-D Linearized and Discretized System','NumberTitle','off');
subplot(2,1,1)
plot(y_d(:,1),y_d(:,3),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$p_{10_x}$$');
ylabel('$$v_{1_x}$$');
legend('trajectory');
title('Phase plot');

subplot(2,1,2)
plot(y_d(:,2),y_d(:,4),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$p_{10_y}$$');
ylabel('$$v_{1_y}$$');
legend('trajectory');
title('Phase plot');

figure('Name','2-D Simulation: 4 states','NumberTitle','off');
subplot(2,1,1);
plot(t,y_c(:,1),'Color',sstgray,'LineStyle','-');
hold on
stairs(t,y_d(:,1),'Color',sstlightgray,'LineStyle','-');
hold on
plot(t,y_c(:,2),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,2),'Color',sstblue,'LineStyle','-');
hold on
stairs(t,y_c_total(:,1),'Color','r','LineStyle','--');
hold on
stairs(t,y_d_total(:,1),'Color','r','LineStyle','-');
% stairs(t,y_d_z(:,1),'Color',sstgray,'LineStyle','-');
% hold on
% stairs(t,y_d_t(:,1),'Color',sstdarkblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$p_{10}$');
title('Position');
legend('$p_{10_x}$ cont.','$p_{10_x}$ disc.','$p_{10_y}$ cont.','$p_{10_y}$ disc.','$p_{10}$ cont.','$p_{10}$ disc.');

subplot(2,1,2);
plot(t,y_c(:,3),'Color',sstgray,'LineStyle','-');
hold on
stairs(t,y_d(:,3),'Color',sstlightgray,'LineStyle','--');
hold on
plot(t,y_c(:,4),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,4),'Color',sstblue,'LineStyle','--');
hold on
stairs(t,y_c_total(:,2),'Color','r','LineStyle','--');
hold on
stairs(t,y_d_total(:,2),'Color','r','LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$v_{1}$');
title('Velocity');
legend('$v_{1_x}$ cont.','$v_{1_x}$ disc.','$v_{1_y}$ cont.','$v_{1_y}$ disc.','$v_1$ cont.','$v_1$ disc.');
hold off;

%% 1.1 Linear 1-D Model 

clear X Y
    
g = 9.8065; % gravity force [m/s^2]
beta = 1/2*1.225*0.3*0.1; % 1/2*p*Cd*A where p is air density, Cd is drag coefficient, A is frontal Area
v1e = 3; % speed of the follower drone
v0e = 3; % speed of the leader drone
m = 0.1; % mass in kg
we = 0; % typ. wind speed [m/s]
theta1e = atan((beta/(m*g))*(v0e-we)^2); % pitch angle in rad


Al = [0 -1; 0 -2*beta*v1e/m];
Bl = [0 g*sec(theta1e)^2]';
El = [0 1; 2*beta*we/m 0];
Cl = eye(2); % Cl = [1 0] but we wanted to access speed
Dl = 0;

% Discrete-time state-space model approximation

Ts = 0.1; % sampling time [s]

Ad = [1 -Ts; 0 1-2*beta*v1e*Ts/m];
Bd = [0 Ts*g*sec(theta1e)^2]';
Ed = [0 Ts; 2*beta*we*Ts/m 0];
Cd = eye(2);  % Cd = [1 0] but we wanted to access speed
Dd = 0;

%  Stability, controllability, and observability of the discrete time model.

[eigenvectors,eigenvalues] = eig(Ad); % stability -> eigenvalues <= 1
[similarity_transform,J] = jordan(Ad);
disp('System marginally stable because eigenvectors <= 1');
disp(eigenvalues);
disp('and all jordan blocks associated with eigenvalues = 1 are 1x1 blocks');
disp(J);

Co = ctrb(Ad,Bd); % controlable - rank of controlability matrix equals n = 2
if rank(Co) == size(Ad,1)
   disp('Controlable') 
else
    disp('Not Controlable') 
end    

Ob = obsv(Ad,Cd); % observable - rank of observability matrix equals n = 2
if rank(Ob) == size(Ad,1)
   disp('Observable') 
else
    disp('Not Observable') 
end   

% Simulate system response

% Continuous response
sys_c = ss(Al,Bl,Cl,Dl); 

tmax = 5; % duration of the simulation
t = 0:Ts:tmax;
u = 0.1*ones(length(t),1);
x0 = [0 0]';
y_c = lsim(sys_c,u,t,x0);
y_c(:,1) = y_c(:,1) + 0.5;
y_c(:,2) = y_c(:,2) + v1e;

% Discrete response
Ts = 0.1;

sys_d = ss(Ad,Bd,Cd,Dd,Ts);
y_d = lsim(sys_d,u,t,x0);

figure('Name','Pole Zero Map - Linearized and Discretized System','NumberTitle','off');
set(gcf,'defaultTextInterpreter','tex');
pzplot(sys_d,'b');
grid on;
axis equal;
% change marker size
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
    set(a(i),'color',sstblue) %change marker size
end
title('Pole Zero Map');

y_d(:,1) = y_d(:,1) + 0.5;
y_d(:,2) = y_d(:,2) + v1e;

sys_d_z = c2d(sys_c, Ts, 'zoh');
y_d_z = lsim(sys_d_z,u,t,x0);
y_d_z(:,1) = y_d_z(:,1) + 0.5;
y_d_z(:,2) = y_d_z(:,2) + v1e;

sys_d_t = c2d(sys_c, Ts, 'tustin');
y_d_t = lsim(sys_d_t,u,t,x0);
y_d_t(:,1) = y_d_t(:,1) + 0.5;
y_d_t(:,2) = y_d_t(:,2) + v1e;

figure('Name','Phase Plot - Linearized and Discretized System','NumberTitle','off');
plot(y_d(:,1),y_d(:,2),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('trajectory');
title('Phase plot');

figure('Name','Simulation of both states for u = 0.1 rad','NumberTitle','off');
subplot(2,1,1);
plot(t,y_c(:,1),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,1),'Color',sstblue,'LineStyle','-');
hold on
stairs(t,y_d_z(:,1),'Color',sstgray,'LineStyle','-');
hold on
stairs(t,y_d_t(:,1),'Color',sstdarkblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$p_{10}$');
title('Position $p_{10}$');
legend('cont.','disc. euler','disc. zoh', 'dic. tustin');

subplot(2,1,2);
plot(t,y_c(:,2),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,2),'Color',sstblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$v_{1}$');
title('Velocity $v_{1}$');
legend('cont.','disc.');
hold off;

% Same but with possibility to simulate some disturbances
nu = size(Bd,2);
nk = tmax/Ts;
TX = 1:nk+1;
x0 = [0 0]';
U = 0.1*ones(nu,nk);
Ud = [0 1]'; % wind and v0 disturbances 
X(:,1) = x0;
X(:,2) = x0;
for k = 2:nk
    X(:,k+1) = Ad*X(:,k) + Bd*U(:,k) + Ed*Ud;
    Y(:,k+1) = Cd*X(:,k+1);
end
Y(1,:) = Y(1,:) + 0.5;
figure('Name','Simulation of both states with disturbances','NumberTitle','off');
subplot(2,1,1);
stairs(t,y_d(:,1),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,Y(1,:),'Color',sstblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$p_{10}$');
title('Position $p_{10}$');
legend('disc.','disc.with E matrix');

v1d = v1e + X(2,:);
subplot(2,1,2);
stairs(t,y_d(:,2),'Color',sstblue,'LineStyle','-');
hold on;
stairs(t,v1d,'-','Color',sstgreen);
grid on;
hold off;
xlabel('t [s]');
ylabel('$$v_{1}$$');
legend('without perturb.','with perturb.');
title('$$v_{1}$$');



%% 1.2 Unconstrained MPC

set(0,'defaultTextInterpreter','latex');

Cl = [1 0];
Cd = [1 0];
clear X Xd Y U dUopt Uopt;

N = 5;
P = 10*eye(1);
Q = 10*eye(1);
R = .1;
xd0 = [0 0]';
nx = size(Bd,1);
nu = size(Bd,2);

% compute matrices
[F,G,Qb,Rb,H,Fd,Gd,Hd,A,B,C] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
Fb = H*F;
Gb = H*G;
Fdb = Hd*Fd;
Gdb = Hd*Gd;

% compute final cost function matrices and control gains
Rt = Gb'*Qb*Gb + Rb,
St = Gb'*Qb,
Ky = Rt^(-1)*St,
K  = Rt^(-1)*St*Fb,

%simulate controlled system:
nk = 200;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = square_dt(nk+N,80); % square reference of amplitude -1 1 on the real system but -1.5 0.5 on simulated linearized system because p10e = 0.5
dist_x1 = 0*0.5*ones(size(ref)).*(Tref>=0);
x0 = [xd0*0 ; Cd*xd0];
U = zeros(nu,nk);
dU = zeros(nu,nk);
Xd(:,1) = xd0;
X(:,1) = x0;
Y(:,1) = Cd*xd0;
Xd(:,2) = xd0;
X(:,2) = x0;
Y(:,2) = Cd*xd0;
for k = 2:nk
    
    % compute initial condition and current reference sequence
    Yb = ref(:,k:k+N)';
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; Cd*Xd(:,k)];
    xk = X(:,k);
    
    % compute optimal incremental control sequence:
    dUopt(:,:,k) = reshape(-(K*xk-Ky*Yb) ,nu,N);

    % set MPC control policy and simulate system:
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    uk = Uopt(:,1,k);
    U(:,k) = uk;
    
    % simulate original system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k) + [dist_x1(:,k);0];
    Y(:,k+1) = Cd*Xd(:,k+1);
        
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; Cd*Xd(:,k+1)];
    
figure('Name','Phase Plot - Unconstrained MPC','NumberTitle','off');
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('trajectory');
title('Phase plot');

figure('Name','State Evolution - Unconstrained MPC','NumberTitle','off');
plot(Tref(1:end-5),ref(1:end-5),'gs-');
grid on;
hold on;
plot(TX(1:end-1),Xd(1,1:end-1),'s-','Color',sstblue);
plot(TX(1:end-1),Xd(2,1:end-1),'d-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('ref','$$x_1$$','$$x_2$$','Location','SouthEast');
title('State evolution');

figure('Name','Control Action - Unconstrained MPC','NumberTitle','off');
plot(1:length(U(1,:,:)),U(1,:,:),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('input');
title('Input');


figure('Name','$$\Delta$$ State evolution - Unconstrained MPC','NumberTitle','off');
plot(TX(1:end-1),X(1,1:end-1),'s-','Color',sstblue);
grid on;
hold on;
plot(TX(1:end-1),X(2,1:end-1),'d-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$\Delta x(t_k)$$');
legend('$$\Delta x_1$$','$$\Delta x_2$$','Location','SouthEast');
title('$$\Delta$$ State evolution');


figure('Name','Simulation - Unconstrained MPC','NumberTitle','off');
subplot(2,1,1)
plot(TX(1:end-1),Xd(1,1:end-1),'s-','Color',sstblue);
grid on;
hold on;
plot(Tref(1:end-5),ref(1,1:end-5),'s-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$y(t_k)$$');
legend('output','reference');
title('Output');
subplot(2,1,2)
plot(1:length(U(1,:,:)),U(1,:,:),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('input');
title('Input');

%% 1.3 Constrained MPC

clear X Xd Xd2 Y U U2 dUopt dUopt2 Uopt Uopt2;

N = 5; %5
P = 10*eye(1); %0.8
Q = 10*eye(1); %0.8
R = 0.1; %0.01
xd0 = [0 0]';
nx = size(Bd,1);
nu = size(Bd,2);

% compute batch matrices
[F,G,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
Fb = H*F;
Gb = H*G;
Fdb = Hd*Fd;
Gdb = Hd*Gd;
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^(-1)*St,
K = Rt^(-1)*St*Fb,

% compute constraints matrices:
u_max = deg2rad(30);
y_max = 5-0.5; %p10max - p10e (to adjust for linearized system)
y_min = 0.2-0.5;
U_max = kron(u_max,ones(N,1));
Y_max = kron(y_max,ones(N+1,1));
Y_min = kron(y_min,ones(N+1,1));
M3 = tril(ones(N*nu));
M4 = ones(N*nu,nu);
Mu = [-M3;M3];
My = [-Gb;Gb];
M = [Mu;My];

%simulate controlled system:
nk = 200;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = square_dt(nk+N,80);% square reference of amplitude -1 1, on the real system this would be -.5 to 1.5 because p10e = 0.5
dist_x1 = 0*.5*ones(size(ref)).*(Tref>=0);
x0 = [xd0*0 ; Cd*xd0];
U = zeros(nu,nk);
U2 = zeros(nu,nk);
Xd(:,1) = xd0;
Xd2(:,1) = xd0;
X(:,1) = x0;
Xd(:,2) = xd0;
Xd2(:,2) = xd0;
X(:,2) = x0;
for k = 2:nk
    
    % compute initial conditions
    Yb = ref(:,k:k+N)';    
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; Cd*Xd(:,k)];
    xk = X(:,k);
    u_1 = U(:,k-1);
    wu = [U_max + M4*u_1;U_max - M4*u_1];
    wy = [-Y_min + Fb*xk;Y_max - Fb*xk];
    w = [wu;wy];
    Dxdk2 = Xd2(:,k)-Xd2(:,k-1);
    xk2 = [ Dxdk2; Cd*Xd2(:,k)];
    
    % compute constrained optimal incremental control sequence and MPC policy
    [dUo,Jo,exitflag,output,lambda] = quadprog(2*Rt,2*St*(Fb*xk-Yb),M,w);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end
    dUopt(:,:,k) = reshape( dUo ,nu,N);
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    U(:,k) = Uopt(:,1,k);
    
    % simulate original system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k) + Ed*[dist_x1(:,k);0];
    
    % compute unconstrained optimal sequence and MPC policy
    dUopt2(:,:,k) = reshape(-(K*xk2-Ky*Yb) ,nu,N);
    Uopt2(:,:,k) = U2(:,k-1) + dUopt2(:,:,k);
    U2(:,k) = Uopt2(:,1,k);
    
    % simulate system 2 (for comparison):
    Xd2(:,k+1) = Ad*Xd2(:,k) + Bd*U2(:,k);
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; Cd*Xd(:,k+1)];
    
figure('Name','Phase Plot - Constrained MPC','NumberTitle','off');
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
grid on;
hold on;
plot(Xd2(1,:),Xd2(2,:),'o-','Color',sstgray);
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('const.','unconst.');
title('Phase plot');

figure('Name','State Evolution - Constrained MPC','NumberTitle','off');
plot(Tref(1,1:end-5),ref(1,1:end-5),'-');
grid on;
hold on;
plot(TX(1:end-1),Xd(1,1:end-1),'s-','Color',sstblue);
plot(TX(1:end-1),Xd(2,1:end-1),'s-','Color',sstgreen);
plot(TX(1:end-1),Xd2(1,1:end-1),'-','Color',sstgray);
plot(TX(1:end-1),Xd2(2,1:end-1),'-','Color',sstlightgray);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('ref','$$x_1$$ const.','$$x_2$$ const.','$$x_1$$ unc.','$$x_2$$ unc.','Location','SouthEast');
title('State evolution');

figure('Name','Simulation - Constrained MPC','NumberTitle','off');
subplot(2,1,1)
plot(TU,U,'-','Color',sstblue);
grid on;
hold on;
plot(TU,U2,'-','Color',sstgray);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('const.','unc.');
title('Input');

subplot(2,1,2)
plot(TX(1:end-1),Xd(1,1:end-1),'-','Color',sstblue);
grid on;
hold on;
plot(Tref(1,1:end-5),ref(1,1:end-5),'-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$y(t_k)$$');
legend('output','reference');
title('Output');





