% Spacecraft Guidance and Navigation (2022/2023)
% Assignment #2 - Guidance - Ex 1 Uncertainty Propagation
% Author: Matteo Luciardello Lecardi
% Personal code: 10572870
% Matricola MSc Degree: 977701
% Total expected execution time ~ 90s

%% 

clc
clear
close all

% cd 'D:\Courses\5\1 Spacecraft guidance and navigation\2022_23\Assignment02\Assignment02'
cspice_furnsh 'assignment.tm';
% opengl hardware
tic2_1 = tic;
tic

%% EX 1.1

% DATA
% ref pericenter epoch t0 [UTC]
t0_utc = '2022-11-11T19:08:49.824';
t0 = cspice_str2et(t0_utc);             % [seconds]
% mean state x0
r0 = [6054.30795817484; -3072.03883303992; -133.115352431876]; % [Km]
v0 = [4.64750094824087; 9.18608475681236; -0.62056520749034];  % [Km/s]
x0 = [r0;v0];
% covariance P0 [km^2, km^2/s, km^2/s^2]
P0 = [+5.6e-3 +3.5e-3 -7.1e-4 0      0      0
      +3.5e-3 +9.7e-3 +7.6e-4 0      0      0
      -7.1e-4 +7.6e-4 +8.1e-4 0      0      0
       0       0       0     +2.8e-7 0      0
       0       0       0      0     +2.7e-7 0
       0       0       0      0      0     +9.6e-8];

ANS.Data.x0 = x0;
ANS.Data.P0 = P0;

%% COMPUTE ORBIT PROPAGATION (Pure Keplerian motion)
figure(1)
figure(2)
close figure 1
close figure 2

mu = cspice_bodvrd('Earth','GM',1);         % [km^3/s^2] GM
h = cross(r0,v0);                           % [km^2/s] specific angular momentum
e = (cross(v0,h)/mu) - (r0/norm(r0));       % eccentricity vector
ne = norm(e);                               % eccentricity
p = dot(h,h)/mu;                            % [Km] semi-latus rectum
a = p/(1 - ne^2);                           % [Km] semimajor axis
T = 2*pi*sqrt(a^3/mu);                      % [s] orbital period
tf = t0 + 4*T;                              % final time
ANS.Data.orbital_mechanics.mu = mu; 
ANS.Data.orbital_mechanics.h = h; 
ANS.Data.orbital_mechanics.e_vec = e; 
ANS.Data.orbital_mechanics.e = ne; 
ANS.Data.orbital_mechanics.p = p; 
ANS.Data.orbital_mechanics.a = a; 
ANS.Data.orbital_mechanics.T = T; 

% CONFIDENCE INTERVAL [%]
s3 = 99.7;
ANS.Data.confidence_interval = s3;
S3 = ['3\sigma: ',num2str(s3),'%'];
% Confidence interval ellipsoid properties
type = 'xy'; % equatorial plane choice for 'error_ellipse_type' fun
               
[~,x,t] = orbit_propagation(x0,t0,t0+T);
r_apo = -a*(1+ne)*e/ne;      % apocenter position
r_peri = a*(1-ne)*e/ne;      % pericenter position
ANS.Data.orbital_mechanics.r_apo = r_apo;
ANS.Data.orbital_mechanics.r_peri = r_peri;

[~,xx,tt,xe,te] = orbit_propagation(x0,t0,tf);

Te = cspice_timout(te','\itYYYY-MON-DD-HR:MN:SC.#### (UTC) ::UTC');
ANS.t_events.T0_UTC = Te(1,9:end); 
for i = 1:4
    ANS.t_events.te(i).apogee = te(2*i);
    ANS.t_events.Te_UTC(i).apogee = Te(2*i,9:end);
    ANS.t_events.te(i).perigee = te(2*i+1);
    ANS.t_events.Te_UTC(i).perigee = Te(2*i+1,9:end);
end

% PLOT orbit

figure('Name','1: First Keplerian orbit','NumberTitle','off')
plot3(x(:,1),x(:,2),x(:,3),'b')
hold on
grid on 
box on
axis equal
plot3(0,0,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
plot3(r_peri(1),r_peri(2),r_peri(3),'LineStyle','none','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(r_apo(1),r_apo(2),r_apo(3),'LineStyle','none','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
title('@ECI: J2000','FontSize',15)
xlabel('X [Km]','FontSize',15,'FontWeight','bold')
ylabel('Y [Km]      ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
[az,el] = view([h(1) h(2) h(3)]);
view([0 el])
legend('','Earth','perigee & apogee','Location','northeast','FontSize',12)

figure('Name','2: Full Keplerian orbit propagation + initial samples','NumberTitle','off')
tiledlayout(1,2,'TileSpacing','compact','Padding','tight')

for k = 1:2
    
    nexttile(k)
%     subplot(1,2,k)
    plot3(0,0,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
    hold on
    grid on  
    box on
    axis equal
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]      ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    view([0 90])
end

nexttile(1)
% subplot(1,2,1)
for i = 2:length(te)
    plot3(xe(i,1),xe(i,2),xe(i,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
end
plot3(xx(:,1),xx(:,2),xx(:,3),'Color','b')
title('@ECI: J2000','FontSize',15)

nexttile(2)
% subplot(1,2,2)
plot3(r0(1),r0(2),r0(3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
H0 = plot3(x(:,1),x(:,2),x(:,3),'Color','b');
title('Initial position',Te(1,:),'Fontsize',15)
ax = gca;
ax.TitleHorizontalAlignment = 'right';
ax.XAxis.Exponent = 4;
ax.YAxis.Exponent = 4;
ax.ZAxis.Exponent = 4;

%% TAKING SAMPLES (multivariate normal random numbers)

rng('default') % For reproducibility
n = 100;       % n. of samples
R0n = mvnrnd(x0',P0,n);
ANS.MC.InitialState_Samples = R0n;

% Plot samples

figure(2)
for i = 1:2
    nexttile(i)
%     subplot(1,2,i)
    H1 = plot3(R0n(:,1),R0n(:,2),R0n(:,3),'+','Color',[0.3010 0.7450 0.9330]); % Samples distribution around mean Initial Mean State

end

nexttile(1)
% subplot(1,2,1)
legend('Earth','perigee & apogee','','','','','','','','orbit','Initial samples')

tol = 4e-1;
nexttile(2)
% subplot(1,2,2)
xlim([x0(1)-tol x0(1)+tol]);
ylim([x0(2)-tol x0(2)+tol]);
zlim([x0(3)-tol x0(3)+tol]);

% Confidence interval ellipse/ellipsoid
% toggle ellipse 'xy' / ellipsoid ('-oid')
% other entries: 'yz', 'zx', 'all', 'all2d'

% type_in = 'xy';
type_in = '-oid';

if strcmp(type_in,'xy')
    H = error_ellipse_type(type_in,P0(1:3,1:3),x0(1:3),s3/100);
    H.Color = [0.4940 0.1840 0.5560];
elseif strcmp(type_in,'-oid')
    H = error_ellipse_type(type_in,P0(1:3,1:3),x0(1:3),s3/100);
    H.EdgeColor = [0.4940 0.1840 0.5560];
    H.FaceColor = 'none';
    H.LineWidth = 0.1;
    H0.LineWidth = 1.2;
    H1.LineWidth = 1.2;
    H1.MarkerSize = 10;
end

% plot(r0(1),r0(2),'Color',[0.4940 0.1840 0.5560])
legend('','initial position = perigee','','Multivariate normal random samples',S3,'Location','southeast','FontSize',12)

%% LINCOV

% MEAN
% Propagated initial mean will remain mean of final distribution

figure('Name','3: Keplerian propagation at each revolution + Apogee and Perigee Samples','NumberTitle','off')

xxx = x0';
t00 = t0;
n_revolution = {'1st','2nd','3rd','4th'};
tiledlayout(3,4,'TileSpacing','compact','Padding','tight')
for i = 1:4
    
    nexttile(i)
%     subplot(3,4,i)
    tff = t00 + T;
    [~,xxx,ttt,xxe,tte] = orbit_propagation(xxx(end,:),t00,tff);

    ANS.Propagation(i).rv_mean = xxx;

    t00 = tff;
    plot3(xxx(:,1),xxx(:,2),xxx(:,3),'b')
    hold on
    grid on   
    box on
    switch i
        case 1
              ANS.Propagation(1).rv_apogee = xe(2,:);
              ANS.Propagation(1).rv_perigee = xe(3,:);
        otherwise
              ANS.Propagation(i).rv_apogee = xe(2*i,:);
              ANS.Propagation(i).rv_perigee = xe(2*i+1,:);
    end
    plot3(0,0,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({[n_revolution{i},' revolution: LINCOV']},'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    axis equal
    view([0 90])

end

% COVARIANCE

% Compute STM for Covariance propagation at apogees and perigees
% Numerical approach

PHI_t = zeros(6,6,length(te)-1);
P = zeros(6,6,length(te)-1);

for i = 1:length(te)-1
    
    % STM Φ(t0,t) 
    PHI_t(:,:,i) = STM(x0,t0,te(i+1));
    
    % Covariance P
    P(:,:,i) = PHI_t(:,:,i) * P0 * PHI_t(:,:,i)';

end
for i = 1:4
    
    ANS.LINCOV.STM(i).apogee = PHI_t(:,:,2*i-1);
    ANS.LINCOV.STM(i).perigee = PHI_t(:,:,2*i);
    ANS.LINCOV.P(i).apogee = P(:,:,2*i-1);
    ANS.LINCOV.P(i).perigee = P(:,:,2*i);

end

%% Evaluation of LINCOV performance

% Take samples

Rn = zeros(n,6,length(te)-1);
for i = 1:length(te)-1
    
    Rn(:,:,i) = mvnrnd(xe(i+1,:),P(:,:,i),n);

end

for i = 1:4
    
    ANS.LINCOV.Samples(i).apogee = Rn(:,:,2*i-1);
    ANS.LINCOV.Samples(i).perigee = Rn(:,:,2*i);

end

% Plot and compare SAMPLES

figure(3)
lgnd = {'','','LINCOV samples'};
for i = 1:4
    
    nexttile(i)
%     subplot(3,4,i)
    plot3(Rn(:,1,2*i-1),Rn(:,2,2*i-1),Rn(:,3,2*i-1),'+','Color',[0.3010 0.7450 0.9330]) % Samples distribution around mean APOGEE
    plot3(Rn(:,1,2*i),Rn(:,2,2*i),Rn(:,3,2*i),'+','Color',[0.3010 0.7450 0.9330])       % Samples distribution around mean PERIGEE
    plot3(ANS.Propagation(i).rv_apogee(1,1),ANS.Propagation(i).rv_apogee(1,2),ANS.Propagation(i).rv_apogee(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    plot3(ANS.Propagation(i).rv_perigee(1,1),ANS.Propagation(i).rv_perigee(1,2),ANS.Propagation(i).rv_perigee(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    legend(lgnd,'Location','northeast','FontSize',12)

end


figure('Name','4: APOGEE Samples','NumberTitle','off')
tiledlayout(3,4,'TileSpacing','compact','Padding','tight')

figure('Name','5: PERIGEE Samples','NumberTitle','off')
tiledlayout(3,4,'TileSpacing','compact','Padding','tight')

tol = [1.2e3 9.2e3 3.5e2];
lgndapo = {'','LINCOV samples','LIN. mean apogee',S3};
lgndperi = {'','LINCOV samples','LIN. mean perigee',S3};

for i = 1:4

figure(4)
    
    nexttile(i)
%     subplot(3,4,i)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on    
    box on
    plot3(Rn(:,1,2*i-1),Rn(:,2,2*i-1),Rn(:,3,2*i-1),'+','Color',[0.3010 0.7450 0.9330]) % Samples distribution around mean APOGEE
    plot3(ANS.Propagation(i).rv_apogee(1,1),ANS.Propagation(i).rv_apogee(1,2),ANS.Propagation(i).rv_apogee(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    H = error_ellipse_type(type,P(1:3,1:3,2*i-1),xe(2*i,1:3),s3/100);
    H.Color = [0.4940 0.1840 0.5560];
    view([0 90])
    xlim([ANS.Propagation(i).rv_apogee(1,1)-tol(3) ANS.Propagation(i).rv_apogee(1,1)+tol(3)])
    ylim([ANS.Propagation(i).rv_apogee(1,2)-tol(1) ANS.Propagation(i).rv_apogee(1,2)+tol(1)])
    zlim([ANS.Propagation(i).rv_apogee(1,3)-tol(1) ANS.Propagation(i).rv_apogee(1,3)+tol(1)])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')    
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['APOGEE ',n_revolution{i},' revolution: LINCOV']},Te(2*i,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    legend(lgndapo,'Location','southeast','FontSize',12)

figure(5) 

    nexttile(i)
%     subplot(3,4,i)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on    
    box on
    plot3(Rn(:,1,2*i),Rn(:,2,2*i),Rn(:,3,2*i),'+','Color',[0.3010 0.7450 0.9330])       % Samples distribution around mean PERIGEE
    plot3(ANS.Propagation(i).rv_perigee(1,1),ANS.Propagation(i).rv_perigee(1,2),ANS.Propagation(i).rv_perigee(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    H = error_ellipse_type(type,P(1:3,1:3,2*i),xe(2*i+1,1:3),s3/100);
    H.Color = [0.4940 0.1840 0.5560];
    view([0 90])
    xlim([ANS.Propagation(i).rv_perigee(1,1)-tol(2) ANS.Propagation(i).rv_perigee(1,1)+tol(2)])
    ylim([ANS.Propagation(i).rv_perigee(1,2)-tol(2) ANS.Propagation(i).rv_perigee(1,2)+tol(2)])
    zlim([ANS.Propagation(i).rv_perigee(1,3)-tol(2) ANS.Propagation(i).rv_perigee(1,3)+tol(2)])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')    
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['PERIGEE ',n_revolution{i},' revolution: LINCOV']},Te(2*i+1,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    ax.ZAxis.Exponent = 4;
    legend(lgndperi,'Location','southeast','FontSize',12)

end

% Square root of covariance trace 
% & Eigenvalues of covariance

for i = 1:4

    ANS.LINCOV.srcov_trace(i).rr_apo = sqrt(trace(ANS.LINCOV.P(i).apogee(1:3,1:3)));
    ANS.LINCOV.srcov_trace(i).rr_peri = sqrt(trace(ANS.LINCOV.P(i).perigee(1:3,1:3)));
    ANS.LINCOV.srcov_trace(i).vv_apo = sqrt(trace(ANS.LINCOV.P(i).apogee(4:6,4:6)));
    ANS.LINCOV.srcov_trace(i).vv_peri = sqrt(trace(ANS.LINCOV.P(i).perigee(4:6,4:6)));

    apo_revolutions(i,:) = eig(ANS.LINCOV.P(i).apogee);
    peri_revolutions(i,:) = eig(ANS.LINCOV.P(i).perigee);

end

for j = 1:6

    ANS.LINCOV.eigenvalues.apogee(j).first = apo_revolutions(1,j);
    ANS.LINCOV.eigenvalues.apogee(j).second = apo_revolutions(2,j);
    ANS.LINCOV.eigenvalues.apogee(j).third = apo_revolutions(3,j);
    ANS.LINCOV.eigenvalues.apogee(j).fourth = apo_revolutions(4,j);

    ANS.LINCOV.eigenvalues.perigee(j).first = peri_revolutions(1,j);
    ANS.LINCOV.eigenvalues.perigee(j).second = peri_revolutions(2,j);
    ANS.LINCOV.eigenvalues.perigee(j).third = peri_revolutions(3,j);
    ANS.LINCOV.eigenvalues.perigee(j).fourth = peri_revolutions(4,j);

end

%%
fprintf('\nExecution time LINCOV: %.2fs\n',toc)
tic

%% UNSCENTED TRANSFORM

% Initial PDF (Probability Density Function)

alpha = 1e-3; % spread of σ points
N = 6;        % size of state
K = 0;
beta = 2;
c = alpha^2*(N+K);

% compute matrix square root of P0

srP = sqrtm(c*P0);

% compose σ points χ(i)

X = zeros(13,6);
X(1,:) = x0;

for i = 1:N

    X(i+1,:) = x0 + srP(:,i);
    X(i+N+1,:) = x0 - srP(:,i);

end

ANS.UT.X_SIGMA_points = X;

% compute weights for mean and covariance

Wm = zeros(2*N+1,1);
Wc = Wm;
Wm(1) = 1-N/c;
Wc(1) = (2-alpha^2+beta)-N/c;
ANS.UT.Weights(1).Wm = Wm(1);
ANS.UT.Weights(1).Wc = Wc(1);

for i = 1:2*N
    Wm(i+1) = 1/(2*c);
    Wc(i+1) = Wm(i+1);

    ANS.UT.Weights(i+1).Wm = Wm(i+1);
    ANS.UT.Weights(i+1).Wc = Wc(i+1);

end

% propagate σ points: γ = g(χ)

g = zeros(2*N+1,N,8);
for j = 1:8
    for i = 1:2*N+1
        g(i,:,j) = orbit_propagation(X(i,:),t0,te(j+1));
    end
end
for i = 1:4
        ANS.UT.G_SIGMA_points(i).apogee = g(:,:,2*i-1);
        ANS.UT.G_SIGMA_points(i).perigee = g(:,:,2*i);
end

% compute weighted sample mean

y = zeros(N,8);
Y = zeros(2*N+1,N);
for j = 1:8
    for i = 1:2*N+1
        Y(i,:,j) = Wm(i)*g(i,:,j);
    end
    y(:,j) = sum(Y(:,:,j),1);
end
for i = 1:N
    
    ANS.UT.y_weighted_MEAN.apogee(i).first = y(i,1);
    ANS.UT.y_weighted_MEAN.apogee(i).second = y(i,3);
    ANS.UT.y_weighted_MEAN.apogee(i).third = y(i,5);
    ANS.UT.y_weighted_MEAN.apogee(i).fourth = y(i,7);

    ANS.UT.y_weighted_MEAN.perigee(i).first = y(i,2);
    ANS.UT.y_weighted_MEAN.perigee(i).second = y(i,4);
    ANS.UT.y_weighted_MEAN.perigee(i).third = y(i,6);
    ANS.UT.y_weighted_MEAN.perigee(i).fourth = y(i,8);

end

% compute weighted sample covariance

PY = zeros(N,N,2*N+1,8);
Py = zeros(N,N,8);
for j = 1:8
    for i = 1:2*N+1
        PY(:,:,i,j) = Wc(i)*(g(i,:,j)' - y(:,j))*(g(i,:,j)' - y(:,j))';
    end
    Py(:,:,j) = sum(PY(:,:,:,j),3);
end

for i = 1:4
    
    ANS.UT.Py_weighted_COV(i).apogee = Py(:,:,2*i-1);
    ANS.UT.Py_weighted_COV(i).perigee = Py(:,:,2*i);

end

%% Evaluation of UT performance

% Take SAMPLES
clear Rn
Rn = zeros(n,N,8);
for i = 1:8
    
    Rn(:,:,i) = mvnrnd(y(:,i),(Py(:,:,i)+Py(:,:,i)')/2,n);

end

for i = 1:4
    
    ANS.UT.Samples(i).apogee = Rn(:,:,2*i-1);
    ANS.UT.Samples(i).perigee = Rn(:,:,2*i);

end
% Plot and compare SAMPLES
lgnd = {'','UT samples'};
lgndapo = {'','UT samples','UT mean apogee',S3};
lgndperi = {'','UT samples','UT mean perigee',S3};

for i = 1:4

figure(3)

    nexttile(i+4)
%     subplot(3,4,i+4)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    axis equal
    hold on
    grid on    
    box on
    plot3(Rn(:,1,2*i-1),Rn(:,2,2*i-1),Rn(:,3,2*i-1),'+','Color',[0.3010 0.7450 0.9330]) % Samples distribution around mean APOGEE
    plot3(Rn(:,1,2*i),Rn(:,2,2*i),Rn(:,3,2*i),'+','Color',[0.3010 0.7450 0.9330])       % Samples distribution around mean PERIGEE
    plot3(y(1,2*i-1),y(2,2*i-1),y(3,2*i-1),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    plot3(y(1,2*i),y(2,2*i),y(3,2*i),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    plot3(0,0,0,'Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
    view([0 90])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({[n_revolution{i},' revolution: UT']},'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    legend(lgnd,'Location','northeast','FontSize',12)

figure(4)

    nexttile(i+4)
%     subplot(3,4,i+4)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on
    box on
    plot3(Rn(:,1,2*i-1),Rn(:,2,2*i-1),Rn(:,3,2*i-1),'+','Color',[0.3010 0.7450 0.9330]) % Samples distribution around mean APOGEE
    plot3(y(1,2*i-1),y(2,2*i-1),y(3,2*i-1),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    H = error_ellipse_type(type,Py(1:3,1:3,2*i-1),y(1:3,2*i-1),s3/100);
    H.Color = [0.4940 0.1840 0.5560];
    view([0 90])
    xlim([y(1,2*i-1)-tol(3) y(1,2*i-1)+tol(3)])
    ylim([y(2,2*i-1)-tol(1) y(2,2*i-1)+tol(1)])
    zlim([y(3,2*i-1)-tol(1) y(3,2*i-1)+tol(1)])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['APOGEE ',n_revolution{i},' revolution: UT']},Te(2*i,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    legend(lgndapo,'Location','southeast','FontSize',12)

figure(5) 

    nexttile(i+4)
%     subplot(3,4,i+4)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on
    box on
    plot3(Rn(:,1,2*i),Rn(:,2,2*i),Rn(:,3,2*i),'+','Color',[0.3010 0.7450 0.9330])       % Samples distribution around mean PERIGEE
    plot3(y(1,2*i),y(2,2*i),y(3,2*i),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    H = error_ellipse_type(type,Py(1:3,1:3,2*i),y(1:3,2*i),s3/100);
    H.Color = [0.4940 0.1840 0.5560];
    view([0 90])
    xlim([y(1,2*i)-tol(2) y(1,2*i)+tol(2)])
    ylim([y(2,2*i)-tol(2) y(2,2*i)+tol(2)])
    zlim([y(3,2*i)-tol(2) y(3,2*i)+tol(2)])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['PERIGEE ',n_revolution{i},' revolution: UT']},Te(2*i+1,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    ax.ZAxis.Exponent = 4;
    legend(lgndperi,'Location','southeast','FontSize',12)

end

% Square root of covariance trace 
% & Eigenvalues of covariance

for i = 1:4

    ANS.UT.srcov_trace(i).rr_apo = sqrt(trace(ANS.UT.Py_weighted_COV(i).apogee(1:3,1:3)));
    ANS.UT.srcov_trace(i).rr_peri = sqrt(trace(ANS.UT.Py_weighted_COV(i).perigee(1:3,1:3)));
    ANS.UT.srcov_trace(i).vv_apo = sqrt(trace(ANS.UT.Py_weighted_COV(i).apogee(4:6,4:6)));
    ANS.UT.srcov_trace(i).vv_peri = sqrt(trace(ANS.UT.Py_weighted_COV(i).perigee(4:6,4:6)));

    apo_revolutions(i,:) = eig(ANS.UT.Py_weighted_COV(i).apogee);
    peri_revolutions(i,:) = eig(ANS.UT.Py_weighted_COV(i).perigee);

end

for j = 1:6

    ANS.UT.eigenvalues.apogee(j).first = apo_revolutions(1,j);
    ANS.UT.eigenvalues.apogee(j).second = apo_revolutions(2,j);
    ANS.UT.eigenvalues.apogee(j).third = apo_revolutions(3,j);
    ANS.UT.eigenvalues.apogee(j).fourth = apo_revolutions(4,j);

    ANS.UT.eigenvalues.perigee(j).first = peri_revolutions(1,j);
    ANS.UT.eigenvalues.perigee(j).second = peri_revolutions(2,j);
    ANS.UT.eigenvalues.perigee(j).third = peri_revolutions(3,j);
    ANS.UT.eigenvalues.perigee(j).fourth = peri_revolutions(4,j);

end


%%
fprintf('Execution time UT: %.2fs\n',toc)
tic

%% EX 1.2 MONTE CARLO

% Retrieve 100 samples from initial state and propagate them at each epoch

% propagate initial state samples: γ = g(R0)

g_mc = zeros(n,N,8);
for j = 1:8
    for i = 1:n
        g_mc(i,:,j) = orbit_propagation(R0n(i,:)',t0,te(j+1));
    end
end
for i = 1:4
        ANS.MC.SamplesMap(i).apogee = g_mc(:,:,2*i-1);
        ANS.MC.SamplesMap(i).perigee = g_mc(:,:,2*i);
end

% Estimate mean and covariance from mapped samples at each epoch

% Mean
y_mc = zeros(6,8);
for j = 1:8
    y_mc(:,j) = mean(g_mc(:,:,j),1);
end

for i = 1:N
    
    ANS.MC.y_mc_MEAN.apogee(i).first = y_mc(i,1);
    ANS.MC.y_mc_MEAN.apogee(i).second = y_mc(i,3);
    ANS.MC.y_mc_MEAN.apogee(i).third = y_mc(i,5);
    ANS.MC.y_mc_MEAN.apogee(i).fourth = y_mc(i,7);

    ANS.MC.y_mc_MEAN.perigee(i).first = y_mc(i,2);
    ANS.MC.y_mc_MEAN.perigee(i).second = y_mc(i,4);
    ANS.MC.y_mc_MEAN.perigee(i).third = y_mc(i,6);
    ANS.MC.y_mc_MEAN.perigee(i).fourth = y_mc(i,8);

end

% Covariance

Py_mc = zeros(6,6,8);

for j = 1:8
   Py_mc(:,:,j) = cov(g_mc(:,:,j));     
end
for i = 1:4
    ANS.MC.Py_mc_COV(i).apogee = Py_mc(:,:,2*i-1);
    ANS.MC.Py_mc_COV(i).perigee = Py_mc(:,:,2*i);
end
%% 

% Plot and compare SAMPLES
lgnd = {'','MC samples'};
lgndapo = {'','MC samples','MC mean apogee',S3};
lgndperi = {'','MC samples','MC mean perigee',S3};

for i = 1:4

figure(3)

    nexttile(i+8)
%     subplot(3,4,i+8)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    axis equal
    hold on
    grid on
    box on
    plot3(g_mc(:,1,2*i-1),g_mc(:,2,2*i-1),g_mc(:,3,2*i-1),'+','Color',[0.3010 0.7450 0.9330]) % Samples distribution around mean APOGEE
    plot3(g_mc(:,1,2*i),g_mc(:,2,2*i),g_mc(:,3,2*i),'+','Color',[0.3010 0.7450 0.9330])       % Samples distribution around mean PERIGEE
    plot3(y_mc(1,2*i-1),y(2,2*i-1),y_mc(3,2*i-1),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    plot3(y_mc(1,2*i),y_mc(2,2*i),y_mc(3,2*i),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    plot3(0,0,0,'Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
    view([0 90])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({[n_revolution{i},' revolution: MC']},'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    legend(lgnd,'Location','northeast','FontSize',12)

figure(4)
    
    nexttile(i+8)
%     subplot(3,4,i+8)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on
    box on
    plot3(g_mc(:,1,2*i-1),g_mc(:,2,2*i-1),g_mc(:,3,2*i-1),'+','Color',[0.3010 0.7450 0.9330]) % Samples distribution around mean APOGEE
    plot3(y_mc(1,2*i-1),y_mc(2,2*i-1),y_mc(3,2*i-1),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    H = error_ellipse_type(type,Py_mc(1:3,1:3,2*i-1),y_mc(1:3,2*i-1),s3/100);
    H.Color = [0.4940 0.1840 0.5560];
    view([0 90])
    xlim([y_mc(1,2*i-1)-tol(3) y_mc(1,2*i-1)+tol(3)])
    ylim([y_mc(2,2*i-1)-tol(1) y_mc(2,2*i-1)+tol(1)])
    zlim([y_mc(3,2*i-1)-tol(1) y_mc(3,2*i-1)+tol(1)])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['APOGEE ',n_revolution{i},' revolution: MC']},Te(2*i,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    legend(lgndapo,'Location','southeast','FontSize',12)

figure(5) 
    
    nexttile(i+8)
%     subplot(3,4,i+8)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on
    box on
    plot3(g_mc(:,1,2*i),g_mc(:,2,2*i),g_mc(:,3,2*i),'+','Color',[0.3010 0.7450 0.9330])       % Samples distribution around mean PERIGEE
    plot3(y_mc(1,2*i),y_mc(2,2*i),y_mc(3,2*i),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    H = error_ellipse_type(type,Py_mc(1:3,1:3,2*i),y_mc(1:3,2*i),s3/100);
    H.Color = [0.4940 0.1840 0.5560];
    view([0 90])
    xlim([y_mc(1,2*i)-tol(2) y_mc(1,2*i)+tol(2)])
    ylim([y_mc(2,2*i)-tol(2) y_mc(2,2*i)+tol(2)])
    zlim([y_mc(3,2*i)-tol(2) y_mc(3,2*i)+tol(2)])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['PERIGEE ',n_revolution{i},' revolution: MC']},Te(2*i+1,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    ax.ZAxis.Exponent = 4;
    legend(lgndperi,'Location','southeast','FontSize',12)

end

% Square root of covariance trace 
% & Eigenvalues of covariance

for i = 1:4

    ANS.MC.srcov_trace(i).rr_apo = sqrt(trace(ANS.MC.Py_mc_COV(i).apogee(1:3,1:3)));
    ANS.MC.srcov_trace(i).rr_peri = sqrt(trace(ANS.MC.Py_mc_COV(i).perigee(1:3,1:3)));
    ANS.MC.srcov_trace(i).vv_apo = sqrt(trace(ANS.MC.Py_mc_COV(i).apogee(4:6,4:6)));
    ANS.MC.srcov_trace(i).vv_peri = sqrt(trace(ANS.MC.Py_mc_COV(i).perigee(4:6,4:6)));

    apo_revolutions(i,:) = eig(ANS.MC.Py_mc_COV(i).apogee);
    peri_revolutions(i,:) = eig(ANS.MC.Py_mc_COV(i).perigee);

end

for j = 1:6

    ANS.MC.eigenvalues.apogee(j).first = apo_revolutions(1,j);
    ANS.MC.eigenvalues.apogee(j).second = apo_revolutions(2,j);
    ANS.MC.eigenvalues.apogee(j).third = apo_revolutions(3,j);
    ANS.MC.eigenvalues.apogee(j).fourth = apo_revolutions(4,j);

    ANS.MC.eigenvalues.perigee(j).first = peri_revolutions(1,j);
    ANS.MC.eigenvalues.perigee(j).second = peri_revolutions(2,j);
    ANS.MC.eigenvalues.perigee(j).third = peri_revolutions(3,j);
    ANS.MC.eigenvalues.perigee(j).fourth = peri_revolutions(4,j);

end


%%
fprintf('Execution time MC: %.2fs\n\n',toc)

%% COMPREHENSIVE PLOT

% Apogee

figure(6)
close(6)
figure('Name','6: Comprehensive plot with MC samples vs Lincov & UT covariance ellipses (apocenters)','NumberTitle','off')
tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

lgndMCapo = {'Reference orbit','MC samples','MC mean apogee',['LINCOV ',S3],['UT ',S3],['MC ',S3]};
for i = 1:4
    nexttile(i)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on
    box on
    
%  plot3(ANS.Propagation(i).rv_apogee(1,1),ANS.Propagation(i).rv_apogee(1,2),ANS.Propagation(i).rv_apogee(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
 plot3(g_mc(:,1,2*i-1),g_mc(:,2,2*i-1),g_mc(:,3,2*i-1),'+','Color',[0.3010 0.7450 0.9330],'LineWidth',1,'MarkerSize',10) % Samples distribution around mean APOGEE
 plot3(y_mc(1,2*i-1),y_mc(2,2*i-1),y_mc(3,2*i-1),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])


    % LINCOV
    H = error_ellipse_type(type,P(1:3,1:3,2*i-1),xe(2*i,1:3),s3/100);
    H.Color = [0.4667 0.6745 0.1882];
    H.LineStyle = '--';
    H.LineWidth = 1.5;
    % UT
    H = error_ellipse_type(type,Py(1:3,1:3,2*i-1),y(1:3,2*i-1),s3/100);
    H.Color = [3 2 1]/3;
    H.LineStyle = '-.';
    H.LineWidth = 1.5;

    % MC
    H = error_ellipse_type(type,Py_mc(1:3,1:3,2*i-1),y_mc(1:3,2*i-1),s3/100);
    H.Color = [0.8784 0.2275 0.2275];
    H.LineStyle = '-';
    H.LineWidth = 1.5;

    view([0 90])
    xlim([y_mc(1,2*i-1)-tol(3) y_mc(1,2*i-1)+tol(3)])
    ylim([y_mc(2,2*i-1)-tol(1) y_mc(2,2*i-1)+tol(1)])
    zlim([y_mc(3,2*i-1)-tol(1) y_mc(3,2*i-1)+tol(1)])    
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')    
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['APOGEE ',n_revolution{i},' revolution']},Te(2*i,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    ax.ZAxis.Exponent = 4;
    legend(lgndMCapo,'Location','southeast','FontSize',12)
end

% Perigee

figure(7)
close(7)
figure('Name','7: Comprehensive plot with MC samples vs Lincov & UT covariance ellipses (pericenters)','NumberTitle','off')
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')
lgndMCperi = {'Reference orbit','MC samples','MC mean perigee',['LINCOV ',S3],['UT ',S3],['MC ',S3]};

for i = 1:4
    nexttile(i)
    plot3(ANS.Propagation(i).rv_mean(:,1),ANS.Propagation(i).rv_mean(:,2),ANS.Propagation(i).rv_mean(:,3),'b')
    hold on
    grid on
    box on
    axis equal
    
%  plot3(ANS.Propagation(i).rv_apogee(1,1),ANS.Propagation(i).rv_apogee(1,2),ANS.Propagation(i).rv_apogee(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(g_mc(:,1,2*i),g_mc(:,2,2*i),g_mc(:,3,2*i),'+','Color',[0.3010 0.7450 0.9330],'LineWidth',1,'MarkerSize',10)       % Samples distribution around mean PERIGEE 
plot3(y_mc(1,2*i),y_mc(2,2*i),y_mc(3,2*i),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])

    % LINCOV
    H = error_ellipse_type(type,P(1:3,1:3,2*i),xe(2*i+1,1:3),s3/100);
    H.Color = [0.4667 0.6745 0.1882];
    H.LineStyle = '--';
    H.LineWidth = 1.5;
    % UT
    H = error_ellipse_type(type,Py(1:3,1:3,2*i),y(1:3,2*i),s3/100);
    H.Color = [3 2 1]/3;
    H.LineStyle = '-.';
    H.LineWidth = 1.5;
    % MC
    H = error_ellipse_type(type,Py_mc(1:3,1:3,2*i),y_mc(1:3,2*i),s3/100);
    H.Color = [0.8784 0.2275 0.2275];
    H.LineStyle = '-';
    H.LineWidth = 1.5;

    view([0 90])
    xlim([y_mc(1,2*i)-tol(2) y_mc(1,2*i)+tol(2)])
    ylim([y_mc(2,2*i)-tol(2) y_mc(2,2*i)+tol(2)])
    zlim([y_mc(3,2*i)-tol(2) y_mc(3,2*i)+tol(2)])
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')    
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    title({['PERIGEE ',n_revolution{i},' revolution']},Te(2*i+1,:),'Fontsize',15)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    ax.ZAxis.Exponent = 4;
    legend(lgndMCperi,'Location','southeast','FontSize',12)
end

%% 
fprintf('Ex 1 execution time: %.1fs\n\n',toc(tic2_1))

%% Clear SPICE kernel pool
cspice_kclear();

%% Functions

function dxdt = keplerian(t,x)

% Keplerian unperturbed motion @ Earth

    mu = cspice_bodvrd('Earth','GM',1);   % [km^3/s^2] GM
    r_ECI = [x(1);x(2);x(3)];
    rdot = [x(4);x(5);x(6)];
     
    norm_r_ECI = norm(r_ECI);
    vdot = - mu*r_ECI/(norm_r_ECI^3);

    dxdt = [rdot ; vdot];

end

%---------------------------------------------------------------

function [value,isterminal,direction]= eventfun(t,x)

r = [x(1);x(2);x(3)];
v = [x(4);x(5);x(6)];
value = dot(r,v);
isterminal = 0;
direction = 0;

end

%---------------------------------------------------------------

function [phi,x,tt,xe,te] = orbit_propagation(x0,t0,t)
% outputs:  
%         phi = [N x 1] final column state x at time t (N = dim of state)
%           x = [I x N] row state propagation along orbit, (I = n. of states)
%          tt = [I x 1] time corresponding to each row of x
%          xe = [9 x N] row states at each epoch corresponding to eventfun
%          te = [9 x 1] time corresponding to each row of xe
    
    if nargout >= 4
        options = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@eventfun);
        [tt,x,te,xe] = ode113(@keplerian,[t0 t],x0,options);
    else
        options = odeset('reltol', 1e-12, 'abstol', 1e-12);
        [tt,x] = ode113(@keplerian,[t0 t],x0,options);
    end
    phi = x(end,:)';

end

%---------------------------------------------------------------

function PHI_t = STM(x0,t0,t)
    % NUMERICAL APPROACH
    % Initialize perturbations and STM
    dx = zeros(length(x0));
    PHI_t = zeros(length(x0));
    % Reference solution
    phi_0 = orbit_propagation(x0,t0,t);
    % Computation of STM by columns
    for i = 1 : length(x0)
        dx(i,i) = sqrt(eps)*max(1, abs(x0(i)));     % Perturbation
        phi_x = orbit_propagation(x0+dx(:,i),t0,t); % Perturbed solution
        PHI_t(:,i) = (phi_x - phi_0) / dx(i,i);     % Forward differences
    end

end

%---------------------------------------------------------------
%---------------------------------------------------------------

function h = error_ellipse_type(type,varargin)
%  Copyright (c) 2004, AJ Johnson
%  All rights reserved.    
%
%    ERROR_ELLIPSE - plot an error ellipse, or ellipsoid, defining confidence region
%    
%    ERROR_ELLIPSE(type,...) - 'string' for 3d case, choose which ellipse to get among 
%    'xy','yz','zx','-oid','all','all2d';
%    leave it [] for 2d case
%    
%    ERROR_ELLIPSE(...,C22) - Given a 2x2 covariance matrix, plot the
%    associated error ellipse, at the origin. It returns a graphics handle
%    of the ellipse that was drawn.
%
%    ERROR_ELLIPSE(...,C33) - Given a 3x3 covariance matrix, plot the
%    associated error ellipsoid, at the origin, as well as its projections
%    onto the three axes. Returns a vector of 4 graphics handles, for the
%    three ellipses (in the X-Y, Y-Z, and Z-X planes, respectively) and for
%    the ellipsoid.
%
%    ERROR_ELLIPSE(...,C,MU) - Plot the ellipse, or ellipsoid, centered at MU,
%    a vector whose length should match that of C (which is 2x2 or 3x3).
%
%    ERROR_ELLIPSE(...,'Property1',Value1,'Name2',Value2,...) sets the
%    values of specified properties, including:
%      'C' - Alternate method of specifying the covariance matrix
%      'mu' - Alternate method of specifying the ellipse (-oid) center
%      'conf' - A value betwen 0 and 1 specifying the confidence interval.
%        the default is 0.5 which is the 50% error ellipse.
%      'scale' - Allow the plot the be scaled to difference units.
%      'style' - A plotting style used to format ellipses.
%      'clip' - specifies a clipping radius. Portions of the ellipse, -oid,
%        outside the radius will not be shown.
%
%    NOTES: C must be positive definite for this function to work properly.

default_properties = struct(...
  'C', [], ... % The covaraince matrix (required)
  'mu', [], ... % Center of ellipse (optional)
  'conf', 0.5, ... % Percent confidence/100
  'scale', 1, ... % Scale factor, e.g. 1e-3 to plot m as km
  'style', '', ...  % Plot style
  'clip', inf); % Clipping radius

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.C = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.mu = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.conf = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.scale = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & ~ischar(varargin{1})
  error('Invalid parameter/value pair arguments.') 
end

prop = getopt(default_properties, varargin{:});
C = prop.C;

if isempty(prop.mu)
  mu = zeros(length(C),1);
else
  mu = prop.mu;
end

conf = prop.conf;
scale = prop.scale;
style = prop.style;

if conf <= 0 | conf >= 1
  error('conf parameter must be in range 0 to 1, exclusive')
end

[r,c] = size(C);
if r ~= c | (r ~= 2 & r ~= 3)
  error(['Don''t know what to do with ',num2str(r),'x',num2str(c),' matrix'])
end

x0=mu(1);
y0=mu(2);

% Compute quantile for the desired percentile
k = sqrt(qchisq(conf,r)); % r is the number of dimensions (degrees of freedom)

hold_state = get(gca,'nextplot');

if r==3 & c==3
  z0=mu(3);
  
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end

  % C is 3x3; extract the 2x2 matricies, and plot the associated error
  % ellipses. They are drawn in space, around the ellipsoid; it may be
  % preferable to draw them on the axes.
  Cxy = C(1:2,1:2);
  Cyz = C(2:3,2:3);
  Czx = C([3 1],[3 1]);

  if strcmp(type,'xy') || strcmp(type,'all')||strcmp(type,'all2d')
  [x,y,z] = getpoints(Cxy,prop.clip);
  h1=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  end 
  if strcmp(type,'yz') || strcmp(type,'all')||strcmp(type,'all2d')
  [y,z,x] = getpoints(Cyz,prop.clip);
  h2=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  end 
  if strcmp(type,'zx') || strcmp(type,'all')||strcmp(type,'all2d')
  [z,x,y] = getpoints(Czx,prop.clip);
  h3=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  end
  if  strcmp(type,'-oid')||strcmp(type,'all')
  [eigvec,eigval] = eig(C);

  [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
  XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
  
  X(:) = scale*(k*XYZ(:,1)+x0);
  Y(:) = scale*(k*XYZ(:,2)+y0);
  Z(:) = scale*(k*XYZ(:,3)+z0);
  h4=surf(X,Y,Z);
  colormap gray
  alpha(0.3)
  camlight
  end
  if nargout
      if strcmp(type,'xy')
        h = h1;
      elseif strcmp(type,'yz')
        h = h2;
      elseif strcmp(type,'zx')
        h = h3;
      elseif strcmp(type,'-oid')
        h = h4;
      elseif strcmp(type,'all')
        h = [h1 h2 h3 h4];
      elseif strcmp(type,'all2d')
        h = [h1 h2 h3];
      end
  end
elseif r==2 & c==2
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end

  [x,y,z] = getpoints(C,prop.clip);
  h1=plot(scale*(x0+k*x),scale*(y0+k*y),prop.style);
  set(h1,'zdata',z+1)
  if nargout
    h=h1;
  end
else
  error('C (covaraince matrix) must be specified as a 2x2 or 3x3 matrix)')
end
%axis equal

set(gca,'nextplot',hold_state);
end
%---------------------------------------------------------------
% getpoints - Generate x and y points that define an ellipse, given a 2x2
%   covariance matrix, C. z, if requested, is all zeros with same shape as
%   x and y.
function [x,y,z] = getpoints(C,clipping_radius)

n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle

[eigvec,eigval] = eig(C); % Compute eigenvector and eigenvalue
xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x = xy(:,1);
y = xy(:,2);
z = zeros(size(x));

% Clip data to a bounding radius
if nargin >= 2
  r = sqrt(sum(xy.^2,2)); % Euclidian distance (distance from center)
  x(r > clipping_radius) = nan;
  y(r > clipping_radius) = nan;
  z(r > clipping_radius) = nan;
end
end
%---------------------------------------------------------------
function x=qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end

s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;

for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end
end
%---------------------------------------------------------------
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));

if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end
end
%---------------------------------------------------------------
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
end
%---------------------------------------------------------------
function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
  arg = varargin{ii};
  if isempty(TargetField)
    if ~ischar(arg)
      error('Propery names must be character strings');
    end
    f = find(strcmp(prop_names, arg));
    if length(f) == 0
      error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
    end
    TargetField = arg;
  else
    % properties.(TargetField) = arg; % Ver 6.5 and later only
    properties = setfield(properties, TargetField, arg); % Ver 6.1 friendly
    TargetField = '';
  end
end
if ~isempty(TargetField)
  error('Property names and values must be specified in pairs.');
end
end