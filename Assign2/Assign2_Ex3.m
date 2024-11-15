% Spacecraft Guidance and Navigation (2022/2023)
% Assignment #2 - Guidance - Ex 3 Sequential Filers
% Author: Matteo Luciardello Lecardi
% Personal code: 10572870
% Matricola MSc Degree: 977701
% Total expected execution time ~ 1500s 
% Worst run ~ 2350s

%% Load kernels
% cd 'D:\Courses\5\1 Spacecraft guidance and navigation\2022_23\Assignment02'
% addpath addpath Assignment02/simulator_ex3
addpath('simulator_ex3')
cspice_furnsh 'assignment.tm';
clearvars; close all; clc; format long
% opengl hardware
ex3tic = tic;
format compact

%% EX 3.1 & 3.2 REALITY SIMULATION
tic
% build Camera Model
foc = 30; % [mm]
dens = 54; % [pix/mm]
b = 1;
p0 = [1920/2;1200/2]; % center pixel location;
Cframe = [1, 0,  0;
          0, 0, -1;
          0, 1,  0];
R = 10;

Cam.f = foc;
Cam.d = dens;
Cam.p0 = p0;
Cam.b = b;
Cam.Cframe = Cframe;
Cam.R = R;

% mean motion of GEO orbit computation:
mu = 398600;
R = 42241.08;
n = sqrt(mu/R^3);

% Target SGN-1 geometrical parameters
m = 213000; % [kg]
l = 10; % [m]
d = 3;
H = 5;
J = m/12*diag([d^2 + H^2,l^2 + H^2,l^2 + d^2]); % inertia matrix

% Data @t0 for nominal trajectory
r0 = [12,-60,0]'; % [m]
v0 = [1e-4,-2*n*r0(1),-1.2e-3]'; % [m/s]
q0 = [0.674156352338764 0.223585877389611 0.465489474399161 0.528055032413102]'; % [-]
w0 = [-0.001262427155865 0.001204540074343 -0.000039180139156]'; % [rad/s]
t0 = cspice_str2et('2023-04-01T14:55:12.023'); % [s]

% nominal trajectory
t = (0:86400) + t0; % [s] measurements at 1Hz, over a time window of 1 day, from t0
x = [r0;v0;q0;w0];
X = zeros(length(t),13);
X(1,:) = x;
prop_tic = tic;
fprintf('\nLoading propagation...              ')
for k = 1:length(t)-1
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%5.d/86401s\n',k+1)
    x = f(x,t(k),t(k+1),n,J); % continuity is implemented
    X(k+1,:) = x;
end
fprintf('\nPropagation nominal trajectory + quaternions time: %.2fs\n',toc(prop_tic))

ANS.SIM(1).fullPropagation = X;

x = X(:,1:6);
r = X(:,1:3)';
q = X(:,7:10)';
 
C_TI = quat2dcm(q');

%% C_LI

C_LI = [cos(n*1) sin(n*1) 0; -sin(n*1) cos(n*1) 0; 0 0 1];
ANS.ROT(1).C_LI = C_LI;

for k = 1:length(t)
    ANS.SIM(k).Propagation = X(k,:);
    ANS.ROT(k).C_TI = C_TI(:,:,k);
end

%% STM

% State Transition Matrix at each time step of 1s
tt = 1; % [s]
STMrr = [4-3*cos(n*tt) 0 0; 6*(sin(n*tt)-n*tt) 1 0; 0 0 cos(n*tt)];
STMrv = [1/n*sin(n*tt) 2/n*(1-cos(n*tt)) 0; 2/n*(cos(n*tt)-1) 1/n*(4*sin(n*tt)-3*n*tt) 0; 0 0 1/n*sin(n*tt)];
STMvr = [3*n*sin(n*tt) 0 0; 6*n*(cos(n*tt)-1) 0 0; 0 0 -n*sin(n*tt)];
STMvv = [cos(n*tt) 2*sin(n*tt) 0; -2*sin(n*tt) 4*cos(n*tt)-3 0; 0 0 cos(n*tt)];
STM = [STMrr,STMrv;STMvr,STMvv];
ANS.ROT(1).STM = STM;

%% REALITY SIMULATION
fprintf('\nLoading simulation...              ')
tic
for k = 1:length(t)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%5.d/86401s\n',k)
    qn = q(:,k)/norm(q(:,k)); % sueperfluous, it is slighlty improving precision
    meas = meas_sim_pvt(n,r(:,k),qn,0,t(k),Cam);
    ANS.SIM(k).meas = meas;
end
fprintf('\nSimulation time: %.2fs\n',toc)

%% Data collection

id_visible = zeros(length(t),8);
number_visible = zeros(length(t),1);
pixels = zeros(3,8,length(t));

for k = 1:length(t)
    id_visible(k,1:length(ANS.SIM(k).meas.visible)) = ANS.SIM(k).meas.visible;
    id = id_visible(k,1:length(ANS.SIM(k).meas.visible));
    number_visible(k) = length(ANS.SIM(k).meas.visible);
    num = number_visible(k);
    for j = 1:num
        pixels(:,id(j),k) = ANS.SIM(k).meas.y(:,j);
    end
    ANS.SIM(k).pixels = pixels(:,:,k);
end


%% Plot simulated motion

figure(1)
close(1)
figure('Name','1: Ex3.1 Nominal Trajectory','NumberTitle','off')
title('LVLH @SGN-1','Nominal trajectory: chaser OSR relative motion','FontSize',15,'FontWeight','bold')
axis equal padded
hold on 
grid on
box on
xlabel('X [m]','FontSize',15,'FontWeight','bold')
ylabel('Y [m]','FontSize',15,'FontWeight','bold')
zlabel('Z [m]','FontSize',15,'FontWeight','bold')

plot3(0,0,0,'Marker','o','MarkerFaceColor',"#EDB120",'MarkerEdgeColor','k','LineStyle','none','MarkerSize',14)
plot3(r(1,1),r(2,1),r(3,1),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','LineStyle','none','MarkerSize',14)
plot3(r(1,:),r(2,:),r(3,:),'k','LineWidth',1.5)
plot3(r(1,end),r(2,end),r(3,end),'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k','LineStyle','none','MarkerSize',6)

legend('target SGN-1','departure','','arrival','FontSize',15)

%% Plot number of visible vertices

figure(2)
close(2)
figure('Name','2: Ex3.1 Number of visible vertices','NumberTitle','off')
title('Number of visible vertices from simulation','FontSize',15,'FontWeight','bold')
hold on 
grid on
box on
% xlim tight
ylim padded

plot((t-t0)/3600,number_visible,'b','LineWidth',1.5)

xlabel('time [h]','FontSize',15,'FontWeight','bold')
xticks(linspace(0,24,9))
xlim([0 24])

%% Plot Pixels V1

figure(3)
close(3)
figure('Name','3: Ex3.2 Verify visible features 1','NumberTitle','off')
% hold on 
% grid on
T = tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
title(T,'Simulated visible features (uncoupled)','FontSize',22,'FontWeight','bold');
for p = 1:8
    ANS.Visibility(p).Indexes = find(pixels(1,p,:)); % visibility indexes
    vis = ANS.Visibility(p).Indexes;
    nexttile(p)
    hold on
    grid minor
    box on
    xlim([0,24])
    ylim padded
    xTicks = linspace(0,24,9);
    xticks(xTicks)

    plot(0,0,'Color',"#0072BD",'LineWidth',2)
    plot(0,0,'Color',"#77AC30",'LineWidth',2)
    plot(0,0,'Color',"#D95319",'LineWidth',2)
    yline(1920,'Color','b','LineWidth',2,'LineStyle','--')
    yline(1200,'Color',"#77AC30",'LineWidth',2,'LineStyle','--')
    plot((t(vis)-t0)/3600,squeeze(pixels(1,p,vis)),'Marker','.','LineStyle','none','Color',"#0072BD",'MarkerSize',10)
    plot((t(vis)-t0)/3600,squeeze(pixels(2,p,vis)),'Marker','.','LineStyle','none','Color',"#77AC30",'MarkerSize',10)
    plot((t(vis)-t0)/3600,squeeze(pixels(3,p,vis)),'Marker','.','LineStyle','none','Color',"#D95319",'MarkerSize',10)

    title(['Point ',num2str(p)],'FontSize',20,'FontWeight','bold')
    xlabel('time [h]','FontSize',19,'FontWeight','bold')
    ylabel('[pix]','FontSize',19,'FontWeight','bold')
    legend('Hpixel','Vpixel','Baseline','Hsize','Vsize','FontSize',20,'Location','northeastoutside')
end

%% Plot Pixels V2

d = char(916);
figure(4)
close(4)
figure('Name','4: Ex3.2 Verify visible features 2','NumberTitle','off')
T = tiledlayout(2,4,'TileSpacing','tight','Padding','tight');
% title(T,'Simulated visible features (coupled)','FontSize',22,'FontWeight','bold')
% hold on 
% grid on
for p = 1:8
    vis = ANS.Visibility(p).Indexes;
    nexttile(p)
    hold on
    grid minor
    box on
    xlim([0,1920])
    ylim([0,1200])
    daspect([20 20 1])
    title(['Point ',num2str(p)],'FontSize',22,'FontWeight','bold')
    scatter3(squeeze(pixels(1,p,vis)),squeeze(pixels(2,p,vis)),squeeze(pixels(3,p,vis)),10,(t(vis)-t0)/3600,'filled')

    c = colorbar(gca,'Location','east','Ticks',linspace(0,24,9),'Limits',[0 24]);
    c.Label.String = '$\mathbf{t\ [h]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 22;
    c.Label.FontWeight = 'bold';

    xlabel('$\mathbf{p_x\ [pix]}$ ','interpreter','latex','FontSize',22)
    ylabel('$\mathbf{p_y\ [pix]}$','interpreter','latex','FontSize',22)
    zlabel('$\mathbf{d\ [pix]}$','interpreter','latex','FontSize',22)

end
nexttile(2)
plot(-1,-1,'Marker','.','LineStyle','none','LineWidth',0.1,'MarkerSize',1,'Color','none')
legend('','Simulated visible features (coupled)','FontSize',22,'Box','off','Color','none','FontWeight','bold','Location','northoutside')

%% EX 3.3 KALMAN FILTERS
% Estimate of initial relative state & its covariance @t0
x0_hat = [15.792658268071492 -59.044939772661586 3.227106250277039 ...
          -0.053960274403210 -0.053969644762889 -0.089140748762173]'; % [m,m/s]
P0 = diag([10,10,10,0.1,0.1,0.1]); % [m^2,m^2/s, m^2/s^2]

ANS.Filters.EKF.x_hat = zeros(length(t),6);
ANS.Filters.EKF.P = zeros(6,6,length(t));
ANS.Filters.UKF.x_hat = zeros(length(t),6);
ANS.Filters.UKF.P = zeros(6,6,length(t));

% EKF - EXTENDED
x_hat_plus_ekf = x0_hat; 
P_plus_ekf = P0;
ANS.Filters.EKF.x_hat(1,:) = x0_hat';
ANS.Filters.EKF.P(:,:,1) = P0;
ANS.Filters.EKF.P_read(1).P_read = P0;

tic_ekf = tic;
% PREDICTION - CORRECTION X(k-1)+ --> X(k)- --> X(k)+
fprintf('\nLoading Extended Kalman Filter...              ')
for k = 1:length(t)-1
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%5.d/86401s\n',k+1)
    [x_hat_plus_ekf,P_plus_ekf] = EKF(x_hat_plus_ekf,P_plus_ekf,t(k),n,J,Cam,ANS.ROT(k+1),C_LI,STM,ANS.SIM(k+1).meas); 
    ANS.Filters.EKF.x_hat(k+1,:) = x_hat_plus_ekf; 
    ANS.Filters.EKF.P(:,:,k+1) = P_plus_ekf;
    ANS.Filters.EKF.P_read(k+1).P_read = P_plus_ekf;
end
fprintf('\nExecution time EKF: %.2fs\n',toc(tic_ekf))

%% UKF - UNSCENTED
x_hat_plus_ukf = x0_hat; 
P_plus_ukf = P0;
ANS.Filters.UKF.x_hat(1,:) = x0_hat';
ANS.Filters.UKF.P(:,:,1) = P0;
ANS.Filters.UKF.P_read(1).P_read = P0;

tic_ukf = tic;
% PREDICTION - CORRECTION X(k-1)+ --> X(k)- --> X(k)+
fprintf('\nLoading Unscented Kalman Filter...              ')
for k = 1:length(t)-1
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%5.d/86401s\n',k+1)
    [x_hat_plus_ukf,P_plus_ukf] = UKF(x_hat_plus_ukf,P_plus_ukf,t(k),n,J,Cam,ANS.ROT(k+1),C_LI,ANS.SIM(k+1).meas);
    ANS.Filters.UKF.x_hat(k+1,:) = x_hat_plus_ukf;
    ANS.Filters.UKF.P(:,:,k+1) = P_plus_ukf;
    ANS.Filters.UKF.P_read(k+1).P_read = P_plus_ukf;
end
    

fprintf('\nExecution time UKF: %.2fs\n',toc(tic_ukf))

%% EX 3.4 ERRORS

x_ekf = ANS.Filters.EKF.x_hat;
x_ukf = ANS.Filters.UKF.x_hat;
P_ekf = ANS.Filters.EKF.P;
P_ukf = ANS.Filters.UKF.P;

% error definition
e_ekf = x - x_ekf;
e_ukf = x - x_ukf;

STD_ekf = zeros(length(t),6);
STD_ukf = zeros(length(t),6);

% estimated standard deviation 3σ
for k = 1 : length(t)
    STD_ekf(k,:) = 3*sqrt(diag(P_ekf(:,:,k))); 
    STD_ukf(k,:) = 3*sqrt(diag(P_ukf(:,:,k)));
end

ANS.Filters.EKF.Error = e_ekf;
ANS.Filters.UKF.Error = e_ukf;
ANS.Filters.EKF.STD = STD_ekf;
ANS.Filters.UKF.STD = STD_ukf;

%% ekf error only
% x_ekf = ANS.Filters.EKF.x_hat;
% % x_ukf = ANS.Filters.UKF.x_hat;
% P_ekf = ANS.Filters.EKF.P;
% % P_ukf = ANS.Filters.UKF.P;
% 
% e_ekf = x - x_ekf;
% % e_ukf = x - x_ukf;
% 
% STD_ekf = zeros(length(t),6);
% % STD_ukf = zeros(length(t),6);
% 
% % estimated standard deviation 3σ
% for k = 1 : length(t)
%     STD_ekf(k,:) = 3*sqrt(diag(P_ekf(:,:,k))); 
% %     STD_ukf(k,:) = 3*sqrt(diag(P_ukf(:,:,k)));
% end
%% ukf error only
% % x_ekf = ANS.Filters.EKF.x_hat;
% x_ukf = ANS.Filters.UKF.x_hat;
% % P_ekf = ANS.Filters.EKF.P;
% P_ukf = ANS.Filters.UKF.P;
% 
% % e_ekf = x - x_ekf;
% e_ukf = x - x_ukf;
% 
% % STD_ekf = zeros(length(t),6);
% STD_ukf = zeros(length(t),6);
% 
% % estimated standard deviation 3σ
% for k = 1 : length(t)
% %     STD_ekf(k,:) = 3*sqrt(diag(P_ekf(:,:,k))); 
%     STD_ukf(k,:) = 3*sqrt(diag(P_ukf(:,:,k)));
% end

%% EKF PLOT

figure(5)
close(5)
figure('Name','5: Ex 3.4 EKF Error along navigation window','NumberTitle','off')
T = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
title(T,'$\textbf{EKF\ estimation\ errors}$','Interpreter','latex','FontSize',25)

% j = [1 4 2 5 3 6];
str0 = {'{r}_X','{r}_Y','{r}_Z','{v}_X','{v}_Y','{v}_Z'};
str1 = {'[m]','[m]','[m]','[m/s]','[m/s]','[m/s]'};
lims = [0.03 0.15 0.03 6e-6 6e-6 6e-6];

for i = 1:6
% k = j(i);
nexttile(i)
hold on
grid on
box on
title(['\boldmath${',str0{i},'}$'],'Interpreter','latex','FontSize',25)
semilogy((t-t0)/3600,e_ekf(:,i),'b','LineWidth',2)
semilogy((t-t0)/3600,STD_ekf(:,i),'k','LineStyle','--','LineWidth',2)
semilogy((t-t0)/3600,-STD_ekf(:,i),'k','LineStyle','--','LineWidth',2)

ylim([-lims(i) lims(i)])
xlim([0 24])
% ylim padded
xticks(xTicks)
xlabel('$\textbf{time\ [h]}$','Interpreter','latex','FontSize',20)
ylabel(['$\textbf{',str1{i},'}$'],'Interpreter','latex','FontSize',20)
legend('$\textbf{Error}$','\boldmath${3\sigma\ \textbf{boundaries}}$','interpreter','latex','FontSize',22)
end

%% UKF PLOT

figure(6)
close(6)
figure('Name','6: Ex 3.4 UKF Error along navigation window','NumberTitle','off')
T = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
title(T,'$\textbf{UKF\ estimation\ errors}$','Interpreter','latex','FontSize',25)

% j = [1 4 2 5 3 6];
str0 = {'{r}_X','{r}_Y','{r}_Z','{v}_X','{v}_Y','{v}_Z'};
str1 = {'[m]','[m]','[m]','[m/s]','[m/s]','[m/s]'};
lims = [0.03 0.15 0.03 6e-6 6e-6 6e-6];

for i = 1:6
% k = j(i);
nexttile(i)
hold on
grid on
box on
title(['\boldmath${',str0{i},'}$'],'Interpreter','latex','FontSize',25)
semilogy((t-t0)/3600,e_ukf(:,i),'b','LineWidth',2)
semilogy((t-t0)/3600,STD_ukf(:,i),'k','LineStyle','--','LineWidth',2)
semilogy((t-t0)/3600,-STD_ukf(:,i),'k','LineStyle','--','LineWidth',2)

ylim([-lims(i) lims(i)])
xlim([0 24])
% ylim padded
xticks(xTicks)
xlabel('$\textbf{time\ [h]}$','Interpreter','latex','FontSize',20)
ylabel(['$\textbf{',str1{i},'}$'],'Interpreter','latex','FontSize',20)
legend('$\textbf{Error}$','\boldmath${3\sigma\ \textbf{boundaries}}$','interpreter','latex','FontSize',22)
end


fprintf('\n\nEx 3 Total execution time: %.2fs\n\n',toc(ex3tic))

%% Clear Kernel Pool
cspice_kclear
% opengl('save','none') 

%% FUNCTIONS

function dy = propagationfun(t,f,n,J)

% Clohessy-Wiltshire linear equations in LVLH of chaser @target

x = f(1);
z = f(3);
xdot = f(4);
ydot = f(5);
zdot = f(6);
dy = [3*n^2*x + 2*n*ydot; -2*n*xdot; -n^2*z];
dy = [xdot;ydot;zdot;dy];

    if length(f) > 6
     
        q = [f(7) f(8) f(9) f(10)]';    % [-]
        w = [f(11) f(12) f(13)]';       % [rad/s]
 %   Euler's equation of rigid body motion
        dw = J\(cross(-w,J*w));
        om = [-w(1) -w(2) w(3) -w(3) -w(2) w(1)];
        W = zeros(4);
        W(triu(true(4) == 1,1)) = om;
        W = W - W';
%   Quaternions dynamics
        dq = 1/2*W*q;
        dy = [dy;dq;dw];
    end
end

%--------------------------------------------------------------------------

function [out,x,t] = f(f,t0,tf,n,J)
% Propagation function
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);
[t,x] = ode113(@propagationfun,[t0,tf],f,opt,n,J);
out = x(end,:)';
end

%--------------------------------------------------------------------------

function [y,H] = h(x,Cam,C_TI,C_LI) 
% Measurement function (loop-less for performance)
% INPUTs: - x     [3 x 1] : Chaser mean position, Target Centered, LVLH Fixed
%         - CAM   struct with the following fields:
%                 - f       a scalar that represents the focal length [mm]
%                 - d       a scalar that represents the pixel density [pix/mm]
%                 - p0      a 2x1 vector containing the coordinates of the center pixel (u0;v0) [pix]
%                 - b       a scalar representing the stereo camera baseline [m]
%                 - Cframe  a 3x3 director cosine matrix representing the rotation necessary to express a vector in LVLH frame to camera frame [-]
%         - C_TI  [3 x 3] : director cosine matrix, time dependent, Inertial --> Target Body Frame
%         - C_LI  [3 x 3] : director cosine matrix, time dependent, Inertial --> LVLH Frame

% OUTPUTs: - y [3M x 1] : triplet Horizontal, Vertical and Baseline Pixels, for each of the M vertices
%          - H [3M x 6] : dh(x)/dx Jacobian of measurement function 

D = Cam.d;
f = Cam.f;
b = Cam.b;
u0 = Cam.p0(1);
v0 = Cam.p0(2);
C_CL = Cam.Cframe;

% Geometrical features of Target
l = 10; % [m]
d = 3;  % [m]
h = 5;  % [m]

% Retreive Vertices position with a set of rotations
% Target Centered, Target Body Fixed
p = [l/2 d/2 h/2];
P = [1 -1 -1; 1 1 -1; 1 1 1; 1 -1 1; -1 -1 -1; -1 1 -1; -1 1 1; -1 -1 1].*p;
% Target Centered, LVLH fixed --> Concatenation of 2 Rotations :
% Body -> Intertial -> LVLH
rP = C_LI * C_TI' * P';

% X Chaser Centered, Camera Fixed --> Concatenation of Translation + Rotation
% Target Centered -> Chaser Centered
% LVLH -> Camera Fixed
x_hat = [x(1) x(2) x(3)]';
Xcccf = C_CL*(rP - x_hat);

X = Xcccf(1,:);
Y = Xcccf(2,:);
Z = Xcccf(3,:);

% stereoscopic imaging model
y = [u0 - D*f * Y./Z; v0 + D*f*X./Z; b*D*f./Z]; % [3 x M]
y = reshape(y,[numel(y),1]); % [3M x 1]

if nargout > 1
%  Jacobian: concatenation of derivatives
%  dh(x_hat)/dx_hat = dh(x_hat)/dXcccf * dXcccf/dx_hat
    O = zeros(1,length(P));
    dh_xhat__dXcccf = [O -D*f./Z D*f*Y./Z.^2; D*f./Z O -D*f*X./Z.^2; O O -b*D*f./Z.^2]; % [3 x M,M,M]
    dh_xhat__dXcccf = reshape(dh_xhat__dXcccf,[length(dh_xhat__dXcccf),3]);           % [3M x 3]
    dXcccf__dx_hat = zeros(3,6);
    dXcccf__dx_hat(:,1:3) = - C_CL;       % [3 x 6]
    H = dh_xhat__dXcccf * dXcccf__dx_hat; % [3M x 6]
end
end

%--------------------------------------------------------------------------

function [x,P] = EKF(x0,P0,t0,n,J,Cam,ROT,C_LI,STM,meas)
% predicted state X(k-1)+ --> X(k)-
x_pred = f(x0,t0,t0+1,n,J);  % 1s time step
% predicted covariance of errors P(k-1)+ --> P(k)-
P_pred = STM * P0 * STM'; % Linear system has analytical STM for Δt = 1s time step

mv = meas.visible;
M = length(mv); % number of visible/available vertices
if ~isempty(mv) % if correction data is available
% Measurement function
    C_TI = ROT.C_TI;
%     C_LI = ROT.C_LI;
    [y_pred,H_pred] = h(x_pred,Cam,C_TI,C_LI);

% UPDATE
% extract only M visible terms: y_pred --> y_hat / H_pred --> H
    y_hat = zeros(3*M,1);
    H = zeros(3*M,6);
    for j = 1:M
        id = mv(j);
        y_hat(3*j-2:3*j,1) = y_pred(3*id-2:3*id,1); % pick only 3M entries
        H(3*j-2:3*j,:) = H_pred(3*id-2:3*id,:); % pick only 3M rows
    end
    
    R = Cam.R * eye(3*M); % measurement Noise covariance
    % K = P_pred * H'*inv(H * P_pred * H' + R);
    K = P_pred * H'/(H * P_pred * H' + R);
    y = reshape(meas.y,[3*M,1]); % real measurement data

    x = x_pred + K * (y - y_hat);
    P = (eye(6) - K * H) * P_pred;
else % if correction is impossible due to data unavailability
    x = x_pred;
    P = P_pred; 
end
P = (P + P')/2; % Forcing symmetry
end

%--------------------------------------------------------------------------

function [x,P] = UKF(x0,P0,t0,n,J,Cam,ROT,C_LI,meas)
% SIGMA POINTS χ
alpha = 1e-3; % spread of σ points
N = 6;        % size of state
K = 0;
beta = 2;
c = alpha^2*(N+K);

% compute matrix square root of P0

srP = sqrtm(c*P0);

% compose sigma points χ(i)
% loop-less for efficiency
X0 = zeros(2*N+1,N);
X0(1,:) = x0;
X0(2:N+1,:) = (x0 + srP)';
X0(N+2 : 2*N+1,:) = (x0 - srP)';

% PROPAGATION OF SIGMA POINTS χ
X_pred = zeros(2*N+1,N);
for k = 1 : 2*N+1
    X_pred(k,:) = f(X0(k,:),t0,t0+1,n,J);
end

% Compute weights
wm(1,1) = 1-N/c;
wm(2 : 2*N+1,1) = 1/(2*c);
wc = wm;
wc(1,1) = (2-alpha^2+beta)-N/c;

Xi = wm .* X_pred;

x_pred = sum(Xi,1)';

Xi_pred = X_pred' - x_pred;
Pi = zeros(N,N,2*N+1);
for k = 1 : 2*N+1
    Pi(:,:,k) = wc(k) * Xi_pred(:,k)*Xi_pred(:,k)';
end
P_pred = sum(Pi,3);


mv = meas.visible;
M = length(mv); % number of visible/available vertices
if ~isempty(mv) % if correction data is available
% Measurement function
    C_TI = ROT.C_TI;
%     C_LI = ROT.C_LI;
    g_pred = zeros(2*N+1,24);
    for k = 1 : 2*N+1
        g_pred(k,:) = h(X_pred(k,:),Cam,C_TI,C_LI);
    end

%  extract only M visible terms: g_pred --> g
    g = zeros(2*N+1,3*M);
    for k = 1 : 2*N+1
        for j = 1:M
            id = mv(j);
            g(k,3*j-2:3*j) = g_pred(k,3*id-2:3*id); % pick only 3M entries
        end
    end

    gi = wm .* g;
    y_pred = sum(gi,1)';

    gi_pred = g' - y_pred;

    R = Cam.R * eye(3*M); % measurement Noise covariance
    Peei = zeros(3*M,3*M,2*N+1);
    Pxyi = zeros(N,3*M,2*N+1);

    for k = 1 : 2*N+1
        Peei(:,:,k) = wc(k) * gi_pred(:,k)*gi_pred(:,k)';
        Pxyi(:,:,k) = wc(k) * Xi_pred(:,k)*gi_pred(:,k)';
    end

    Pee = sum(Peei,3) + R;
    Pxy = sum(Pxyi,3);

    % K = Pxy * inv(Pee);
    K = Pxy / Pee;
    y = reshape(meas.y,[3*M,1]); % real measurement data
    x = x_pred + K * (y - y_pred);
    P = P_pred - K * Pee * K';
else % if correction is impossible due to data unavailability
    x = x_pred;
    P = P_pred;
end
P = (P + P')/2; % Forcing symmetry
end
