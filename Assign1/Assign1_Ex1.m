% Spacecraft Guidance and Navigation (2022/2023)
% Assignment #1 - Guidance - Ex 1 Periodic Orbits
% Author: Matteo Luciardello Lecardi
% Personal code: 10572870
% Matricola MSc Degree: 977701
% Total expected execution time ~ 80s

%% Load kernels

% cd 'D:\Courses\5\1 Spacecraft guidance and navigation\2022_23\Assignment01'
cspice_furnsh 'assignment.tm';
% opengl hardware

%% Ex 1.1
clearvars; close all; clc; format long g;
figure('Name','1: Earth-Moon Collinear Langrange points and potential curve','NumberTitle','off')

%Sun-Earth
mu = 3.0359*10^-6;
dUdx = @(x) x-(1-mu)*(x+mu)./abs(x+mu).^3-mu*(x+mu-1)./abs(x+mu-1).^3 ;

% 10-digit accuracy
options.TolX = 1e-10;

% Libration Points
eps = 10^-4;
x0_L1 = fzero(dUdx,[-mu+eps 1-mu-eps],options);
x0_L2 = fzero(dUdx,[1-mu+eps 2],options); 
x0_L3 = fzero(dUdx,[-2 -mu-eps],options);

P1 = -mu;  % Earth barycenter
P2 = 1-mu; % Moon barycenter

x1 = -2:eps:-mu-eps ; x2 = -mu+eps:eps:1-mu-eps ; x3 = 1-mu+eps:eps:2 ;

% PLOT
tiledlayout(1,2,"TileSpacing","compact","Padding","tight")

% subplot(1,2,1)
nexttile(1)
title('Earth','FontSize',15,'FontWeight','bold')
grid on
hold on
box on
plot(x1, dUdx(x1), x2, dUdx(x2),'Color','k','LineWidth',2); axis([-2 0.99 -80 80])
xline(P1,'--','LineWidth',1.8,'Color',[0.7 0.7 0.7])
plot(P1,0,'Marker','o','MarkerFaceColor',[0.07 0.62 1],'MarkerEdgeColor','k','MarkerSize',10)
plot(x0_L3,0,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8)
text(x0_L3-0.05,-5.22,'L3','FontSize',14)

legend('\boldmath$\frac{\partial U}{\partial x}$','interpreter','latex','FontSize',20,'Location','southwest')
xlabel('X [-]','FontSize',15,'FontWeight','bold')
ylabel('Y [-]   ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')

% subplot(1,2,2)
nexttile(2)
title('Moon','FontSize',15,'FontWeight','bold')
grid on
hold on
box on

plot(x2, dUdx(x2), x3, dUdx(x3),'Color','k','LineWidth',2) ; axis([0.981 1.019 -0.35 0.35]);
xline(P2,'--','LineWidth',1.8,'Color',[0.7 0.7 0.7])
plot(P2,0,'Marker','o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','k','MarkerSize',9)
plot(x0_L1,0,'o',x0_L2,0,'o','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8)
text(x0_L1-0.0005,-0.022,'L1','FontSize',14)
text(x0_L2-0.0005,-0.022,'L2','FontSize',14)

xlabel('X [-]','FontSize',15,'FontWeight','bold')
legend('\boldmath$\frac{\partial U}{\partial x}$','interpreter','latex','FontSize',20,'Location','southwest')

ylabel('Y [-]   ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')

fprintf('Ex 1: \n\n x0_L1 = %.10f[-] \n\n x0_L2 = %.10f[-] \n\n x0_L3 = %.10f[-] \n\n',x0_L1,x0_L2,x0_L3)
clear eps

%% EX 1.2+1.3

ticHalo = tic;
x0 = [1.008296144180133; 0; 0.001214294450297; 0; 0.010020975499502; 0]; % initial state
z0 = linspace(x0(3),0.0046,10);
ex_t = zeros(size(z0,2),1);


for i = 1:length(z0)
    tic
    x0(3) = z0(i);
    fprintf('\n                                                     ')
    [ex_t(i), x0, orbit] = halo_orbit(x0,mu,i);
    ANS(i).orbit = orbit;
    ex_t(i) = toc;
end
%% 

figure(2)
close figure 2

plotHalo(P2,x0_L2,z0,ANS)

for i = 1:length(z0)
  
  fprintf('\nExecution time orbit n. %d: %.2f seconds',i,ex_t(i))
  
end

fprintf('\n\nEx 1 total execution time halo orbits: %.2f seconds\n',toc(ticHalo))


%% Functions

function [value, isterminal, direction] = ev_1_orbit2d(t,XM,mu)

    value = XM(2);
    isterminal = 1;
    direction = -1;

end

%--------------------------------------------------------------------------

function dXMdt = fxm3d(t,XM,mu)
% features STM propagation

    x = XM(1);
    y = XM(2);
    z = XM(3);
    u = XM(4);
    v = XM(5);
    w = XM(6);
    
    M = reshape(XM(7:42),[6 6])';
    r1 = sqrt((x+mu)^2+y^2+z^2);
    r2 = sqrt((x+mu-1)^2+y^2+z^2);
    
    % r3bp dynamics
    
    dxdt = [u; 
            v; 
            w; 
            2*v+x-(1-mu)*(x+mu)/r1^3-mu*(x+mu-1)/r2^3;
            -2*u+y-(1-mu)*y/r1^3-mu*y/r2^3;
            -(1-mu)*z/r1^3-mu*z/r2^3];
    
    % variational equations    
    
    omegaxx = 1-(1-mu)/r1^3+3*(1-mu)*(x+mu)^2/r1^5-mu/r2^3+3*mu*(x-1+mu)^2/r2^5;
    omegaxy = 3*(1-mu)*(x+mu)*y/r1^5+3*mu*(x+mu-1)*y/r2^5;
    omegayx = omegaxy;
    omegayy = 1-(1-mu)/r1^3+3*(1-mu)*y^2/r1^5-mu/r2^3+3*mu*y^2/r2^5;
    omegazz = -(1-mu)/r1^3-mu/r2^3+3*z^2*(1-mu)/r1^5+3*mu*z^2/r2^5;
    omegazx = 3*mu*z*(mu+x-1)/r2^5+3*z*(1-mu)*(mu+x)/r1^5;
    omegaxz = omegazx;
    omegazy = 3*mu*y*z/r2^5+3*y*z*(1-mu)/r1^5;
    omegayz = omegazy;
    
    A = [      0       0       0  1 0 0
               0       0       0  0 1 0
               0       0       0  0 0 1
         omegaxx omegaxy omegaxz  0 2 0
         omegayx omegayy omegayz -2 0 0
         omegazx omegazy omegazz  0 0 0];
    
    dMdt = A*M;
    
    % packed output
    
    dXMdt = [dxdt; reshape(dMdt,[36 1])];

end 
    
%--------------------------------------------------------------------------    

function phi = RTPB_propagator(x0,t,mu)
    options = odeset('Abstol',1e-12,'Reltol',1e-12);
    [~,phi] = ode113(@odefun3d, [0 t], x0, options, mu);       
end

%--------------------------------------------------------------------------

function dy = odefun3d(t,f,mu)
% Propagation function of CRTBP

    x = f(1);
    y = f(2);
    z = f(3);
    u = f(4);
    v = f(5);
    w = f(6);
    
    r1 = sqrt((x+mu)^2+y^2+z^2);
    r2 = sqrt((x+mu-1)^2+y^2+z^2);
    
    dy = [u; 
          v;  
          w;
          2*v+x-(1-mu)*(x+mu)/r1^3-mu*(x+mu-1)/r2^3;
          -2*u+y-(1-mu)*y/r1^3-mu*y/r2^3;
          -(1-mu)*z/r1^3-mu*z/r2^3];

end

%--------------------------------------------------------------------------

function [ex_t, xx0, phi] = halo_orbit(x0,mu,i)

    tic
    tol = 1e-8;
    Vxe = 1;
    Vze = 1;
    n = 0;
    xx0 = [x0; reshape(eye(6),[36 1])];
    
    N = [1172 1226 1272 1301 1315 1283 1241 1185 1080 1178];
    
    if i == 1
        k = '1st';
    elseif i == 2
        k = '2nd';
    elseif i == 3
        k = '3rd';
    else
        k = sprintf('%2.d°',i);
    end

    while abs(Vxe) >= tol && abs(Vze) >= tol

        fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' ...
            ' %s/10 halo orbit has been corrected %4.d/%d times'],k,n,N(i))


        options = odeset('Abstol',1e-9,'Reltol',1e-7,'Events',@ev_1_orbit2d);
        [~,~,te,xxMMe,~] = ode113(@fxm3d,[0 2*pi],xx0,options,mu);

        % final time variables
        T = te; % half period
        xxf = xxMMe(1:6); % final state
        MMf = reshape(xxMMe(7:42),[6 6]);

        % Variational approach via linear system

        % state vector components at y = 0 (half period)
        x = xxf(1); y = xxf(2); z = xxf(3);
        u = xxf(4); v = xxf(5); w = xxf(6);

        b = [-xxMMe(4); -xxMMe(6)];
        A = [MMf(4,1), MMf(4,5); MMf(6,1), MMf(6,5)];
        LINEAR = A\b;
        xx0(1) = xx0(1) + LINEAR(1);
        xx0(5) = xx0(5) + LINEAR(2);

        Vxe = u;
        Vze = w;

        n = n+1;

    end

% Orbit propagation from corrected initial condition

    phi = RTPB_propagator(xx0(1:6),2*T,mu);
%     phi_e = phi(end,1:6)';
    xx0 = xx0(1:6);
    ex_t = toc;
end

%--------------------------------------------------------------------------

function plotHalo(P2,x0_L2,z0,ANS)

% Plot halo orbit
figure('Name','2: 10 halo orbits around L2','NumberTitle','off')

az = [0 90 0 37.5];
el = [0 0 90 30];
%
DU = 3.84405e5; % [Km] Distance unite Earth-Moon

radii = cspice_bodvrd('MOON','RADII',3);
reM = radii(1)/DU;
[XM,YM,ZM] = ellipsoid(0,0,0,reM,reM,reM,60);
T = tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
title(T,'Earth-Moon rotating frame','FontSize',15,'FontWeight','bold')
str0 = {'Y [-]','Y [-]','Y [-]        ','Y [-]'};
str1 = {'top','top','middle','top'};
str2 = {'on','on','on','off'};
for j = 1:length(az)

%     subplot(2,2,j)
    nexttile(j);
    view(az(j),el(j))

    hold on
    grid on
    axis equal padded

    s(j) = surf(XM+P2,YM,ZM,'EdgeColor',[0.8 0.8 0.8],'FaceColor','none');

    plot3(x0_L2,0,0,'Marker','o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k','LineStyle','none')
    plot3(P2,0,0,'Marker','o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])

    xlabel('X [-]','FontSize',17,'FontWeight','bold')
    ylabel(str0{j},'FontSize',17,'FontWeight','bold','Rotation',0,'VerticalAlignment',str1{j},'HorizontalAlignment','center')
    zlabel('Z [-]      ','FontSize',17,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
    ax = gca;
    ax.XAxis.Exponent = -3;
    pbaspect(ax,[1 1 1])
    box(str2{j})

    for i = 1 : length(z0)
        phi = ANS(i).orbit;
        m = plot3(phi(1,1),phi(1,2),phi(1,3),'o','MarkerSize',4,'MarkerFaceColor',"#4DBEEE",'MarkerEdgeColor','k');
        p = plot3(phi(:,1),phi(:,2),phi(:,3),'Color',"#4DBEEE",'LineWidth',1.4);
        if i == 1
            m.MarkerFaceColor = 'k';
            p.Color = 'k';
        end
    end

end
nexttile(1)

text(phi(1,1)-0.0025,phi(1,2),phi(1,3)+0.002,'\boldmath${z_0 = 0.0046}$','FontSize',20,'Interpreter','latex')
txtp1 = [phi(1,1)-0.0025,phi(1,2),phi(1,3)+0.001];
txtp2 = [phi(1,1),phi(1,2),phi(1,3)];
legend('Moon','L2 point','','corrected initial point','1st Halo orbit','','L2 family','Location','south','FontSize',11)
end

