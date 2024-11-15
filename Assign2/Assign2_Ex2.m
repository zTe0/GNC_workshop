% Spacecraft Guidance and Navigation (2022/2023)
% Assignment #2 - Guidance - Ex 2 Batch Filers
% Author: Matteo Luciardello Lecardi
% Personal code: 10572870
% Matricola MSc Degree: 977701
% Total expected execution time ~ 500s (about 8.3 min, ~ 50s for ex 2.1 & 2.2)
% worst try ~ 1000s

%% 

clc
clear
close all

% cd 'D:\Courses\5\1 Spacecraft guidance and navigation\2022_23\Assignment02\Assignment02'
cspice_furnsh 'assignment.tm';
addpath('sgp4')
% opengl hardware
% warning off
ticex_2 = tic;

%% EX 2.1 Compute visibility windows
tic
% Data
% mean state x0
r0 = [6054.30795817484; -3072.03883303992; -133.115352431876]; % [Km]
v0 = [4.64750094824087; 9.18608475681236; -0.62056520749034];  % [Km/s]
x0 = [r0;v0];

% Orbital elements
mu = cspice_bodvrd('Earth','GM',1);         % [km^3/s^2] GM
h = cross(r0,v0);                           % [km^2/s] specific angular momentum
e = (cross(v0,h)/mu) - (r0/norm(r0));       % eccentricity vector
ne = norm(e);                               % eccentricity
p = dot(h,h)/mu;                            % [Km] semi-latus rectum
a = p/(1 - ne^2);                           % [Km] semimajor axis
T1 = 2*pi*sqrt(a^3/mu);                     % [s] orbital period

% Compute time interval for visibility windows
t0_utc = '2022-11-12T04:30:00.000'; 
tf_utc = '2022-11-14T16:30:00.000';

t0 = cspice_str2et(t0_utc); %[s]
tf = cspice_str2et(tf_utc); %[s]

norbits = ceil((tf-t0)/T1);   % number of orbits around Earth
t_grid = t0:60:tf;      % [s] one point per minute grid

TFORMAT = 'YYYY-MON-DD-HR:MN:SC.### ::UTC';
T1_grid = cspice_timout(t_grid + 1e-4 ,TFORMAT); % 1e-4 sec offset for table readibility and indexing

for i = 1:length(t_grid) 
    ANS.Data.Time_grid(i).UTC_grid = T1_grid(i,:);                   
    ANS.Data.Time_grid(i).ephemeris_grid = t_grid(:,i) + 1e-4; % [s]
end

ANS.Data.x0 = x0;
ANS.Data.t0 = t0_utc;
ANS.Data.tf = tf_utc;

% Orbit propagation
x_grid = zeros(length(t_grid),6);
Te = zeros(length(t_grid),1);
x_grid(1,:) = x0';
xx0 = x0;
tt0 = t0;
for i = 1:length(t_grid)-1
    [x_grid(i+1,:),~,~,~,te] = orbit_propagation(xx0,tt0,t_grid(i+1));
    if ~isempty(te) 
        Te(i,:) = te;
    end
    xx0 = x_grid(i+1,:);
    tt0 = t_grid(i+1);
end
ANS.Propagation.ECI = x_grid;

% retrieve t_ grid indexes at the beginning of each revolution
TTe = find(Te);
TTe(5) = length(t_grid);
n_revolutions = length(TTe)-1;
for i = 1:length(TTe)-1
ANS.Data.Revolutions(i).UTC = T1_grid(TTe(i),:);
ANS.Data.Revolutions(i).ephemeris = t_grid(TTe(i)) + 1e-4;
ANS.Data.Revolutions(i).index = TTe(i);
end

% Compute orbit state frame transformation from ECI to ECEF

ROT6_ECI2ECEF = zeros(6,6,length(t_grid));
state_ECEF = zeros(length(t_grid),6);
for i = 1:length(x_grid) 
    ROT6_ECI2ECEF(:,:,i) = cspice_sxform('J2000','ITRF93',t_grid(i)); % rotation matrix ECI to ECEF convertion
%     ROT6_ECI2ECEF(:,:,i) = cspice_sxform('J2000','IAU_EARTH',t_grid(i)); % rotation matrix ECI to ECEF convertion
    state_ECEF(i,:) = ROT6_ECI2ECEF(:,:,i) * x_grid(i,:)';
end

% Extract position

ROT_ECI2ECEF = ROT6_ECI2ECEF(1:3,1:3,:);
pos_ECEF = state_ECEF(:,1:3);

ANS.Propagation.ECEF = pos_ECEF;

% Get Earth radii (equatorial and polar) 
% (requires pck00010.tpc kernel)
radii = cspice_bodvrd('EARTH','RADII',3);
re = radii(1); rp = radii(3);
% Compute flattening
flat = (re - rp)/re;
% flat = 0;

% Define stations coordinates
% KOUROU   
lat(1,1) =   5.25144;  % [deg] 
lon(1,1) = -52.80466;  % [deg]
alt(1,1) = -14.67;     % [m] 
% PERTH
lat(2,1) = -31.80252;  % [deg] 
lon(2,1) = 115.88516;  % [deg]
alt(2,1) =  22.16;     % [m]

lat_rad = zeros(2,1);
lon_rad = zeros(2,1);
alt_km = zeros(2,1);
pos_station_ECEF = zeros(2,3);
ECEF2TOPO = zeros(3,3,2);

for i = 1:2
    % Convert to radians and km
    lat_rad(i,1) = deg2rad(lat(i)); 
    lon_rad(i,1) = deg2rad(lon(i));
    alt_km(i,1) = alt(i)/1000;
    % Compute station pos wrt Earth center (ECEF)
    pos_station_ECEF(i,:) = cspice_pgrrec('EARTH',lon_rad(i),lat_rad(i),alt_km(i),re,flat);
    % Compute ECEF2TOPO rotation matrix
    ECEF2TOPO(:,:,i) = cspice_eul2m(lat_rad(i)-pi, pi-lon_rad(i),pi/2,2,1,2);
end

% Compute orbit position frame transformation from ECEF to TOPO for both stations

pos_TOPO = zeros(length(t_grid),3,2);
for i = 1:length(t_grid)
    for j = 1:2     
        pos_TOPO(i,:,j) = ECEF2TOPO(:,:,j)*(pos_ECEF(i,:)' - pos_station_ECEF(j,:)');
    end
end

ANS.Propagation.Topocentric.KOUROU = pos_TOPO(:,:,1);
ANS.Propagation.Topocentric.PERTH = pos_TOPO(:,:,2);

% Azimuth and elevation

Az = zeros(length(t_grid),2);
El = zeros(length(t_grid),2);
stationName = {'KOUROU','PERTH'};

for i = 1:length(t_grid)
    for j = 1:2
        Az(i,j) = atan2d(pos_TOPO(i,2,j),pos_TOPO(i,1,j));
        El(i,j) = asind(pos_TOPO(i,3,j)/norm(pos_TOPO(i,:,j)));
%         [Az_ap(i,j),El_ap(i,j)] = antenna_pointing(stationName{j},t_grid(i),[x_grid(i,1:3)';x_grid(i,4:6)']);
%         [~, Az_ap(i,j),  El_ap(i,j)] = cspice_reclat(pos_TOPO(i,:,j)');
    end
end

% Visibility Windows

% Compute indexes of visibility, accounting for minimum elevation
i_visibility = false(length(t_grid),2);
%     i_visibility(:,i) = El(:,i) > 0;
i_visibility(:,1) = El(:,1) >= 10; % [deg] minimum elevation for KOUROU station
i_visibility(:,2) = El(:,2) >= 5;  % [deg] minimum elevation for PERTH station

% Compute visible azimuth & elevation
aaz = zeros(length(t_grid),2);
eel = zeros(length(t_grid),2);
for i = 1:2
    aaz(i_visibility(:,i),i) = Az(i_visibility(:,i),i);
    eel(i_visibility(:,i),i) = El(i_visibility(:,i),i);
end
for i = 1:length(t_grid)
    ANS.Visibility.Azimuth(i).KOUROU = aaz(i,1);
    ANS.Visibility.Azimuth(i).PERTH = aaz(i,2);
    ANS.Visibility.Elevation(i).KOUROU = eel(i,1);
    ANS.Visibility.Elevation(i).PERTH = eel(i,2);
end

% Compute visibility time grid
T = zeros(length(t_grid),2);
for i = 1:2
    T(i_visibility(:,i),i) = t_grid(i_visibility(:,i));
end
for i = 1:length(t_grid)
    ANS.Visibility.Time(i).KOUROU = T(i,1);
    ANS.Visibility.Time(i).PERTH = T(i,2);
end

% Compute visibility position in ECI & ECEF
vis_pos_ECI = zeros(length(t_grid),3,2);
vis_pos_ECEF = zeros(length(t_grid),3,2);
for i = 1:2
    vis_pos_ECI(i_visibility(:,i),:,i) = x_grid(i_visibility(:,i),1:3);
    vis_pos_ECEF(i_visibility(:,i),:,i) = pos_ECEF(i_visibility(:,i),:);
end
for i = 1:length(t_grid)
    ANS.Visibility.ECI(i).KOUROU = vis_pos_ECI(i,:,1);
    ANS.Visibility.ECI(i).PERTH = vis_pos_ECI(i,:,2);
    ANS.Visibility.ECEF(i).KOUROU = vis_pos_ECEF(i,:,1);
    ANS.Visibility.ECEF(i).PERTH = vis_pos_ECEF(i,:,2);
end

% Visibility table

% Manually 
% retrieve indexes of Start & End times from spy(i_visibility) plot in figure(3)

% F_start = [16 913 1142 1420 2300 2948 3427 3547];
% F_end = [668 1063 1150 2271 2857 3392 3457 3601];
% StationName = [stationName,'KOUROU',stationName,stationName,'PERTH'];

% Automatic process 
% to retrieve indexes of Start & End times of visibility windows and their associated Station

% [F_start,F_end,StationName,ANS.Visibility.Windows] = automatic_visibility_table(i_visibility,stationName,t_grid);
[F_start,F_end,StationName] = automatic_visibility_table(i_visibility,stationName);

% Data collection

T_start = cell(8,1);
T_end = cell(8,1);
for j = 1:length(F_start)
    T_start{j} = T1_grid(F_start(j),:);
    T_end{j} = T1_grid(F_end(j),:);
end

ANS.Visibility.Windows = struct;
for i = 1:length(F_start)
    ANS.Visibility.Windows(i).Pass = i;
    ANS.Visibility.Windows(i).StationName = StationName{i};
    ANS.Visibility.Windows(i).StartTime_UTC = T_start{i};
    ANS.Visibility.Windows(i).EndTime_UTC = T_end{i};
end
Windows = struct2table(ANS.Visibility.Windows);
for i = 1:width(Windows)
    Windows.(i) = categorical(Windows.(i));
end
Windows.Properties.VariableNames = {'Pass#','Station Name','Start Time(UTC)','End Time(UTC)'}
ANS.Visibility.Windows = Windows;

%% PLOTS

figure(1)
close figure 1
figure(2)
close figure 2
figure(3)
close figure 3
figure(4)
close figure 4
figure(5)
close figure 5
figure(6)
close figure 6


figure('Name','1: Orbit propagation','NumberTitle','off')
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
for i = 1:4
    nexttile(i)
%     subplot(2,2,i)
    grid on
    hold on
    box on
%     view([90 0])
    axis equal tight
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    ax = gca;
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;   
    ax.ZAxis.Exponent = 4;   
    switch i
        case {1,2}
            plot3(0,0,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
            view([0 90])
        case {3}
            plot3(0,0,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.9290 0.6940 0.1250])
            view([90 0])
        case {4}
            plot3(0,0,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.72 0.27 1])
            view([90 0])
    end
end

nexttile(1)
% subplot(2,2,1)
plot3(x_grid(:,1),x_grid(:,2),x_grid(:,3),'k')
plot3(r0(1),r0(2),r0(3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(x_grid(end,1),x_grid(end,2),x_grid(end,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
title('@ECI J2000','FontSize',15,'FontWeight','bold')
legend('Earth','orbit','r_{0}','r_{end}','Location','northeastoutside','FontSize',15)

nexttile(2)
% subplot(2,2,2)
plot3(pos_ECEF(:,1),pos_ECEF(:,2),pos_ECEF(:,3),'k')
plot3(pos_ECEF(1,1),pos_ECEF(1,2),pos_ECEF(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(pos_ECEF(end,1),pos_ECEF(end,2),pos_ECEF(end,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
plot3(pos_station_ECEF(1,1),pos_station_ECEF(1,2),pos_station_ECEF(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.9290 0.6940 0.1250])
plot3(pos_station_ECEF(2,1),pos_station_ECEF(2,2),pos_station_ECEF(2,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.72 0.27 1])
title('@ECEF ITRF93','FontSize',15,'FontWeight','bold')
legend('Earth','orbit','r_{0}','r_{end}','KOUROU','PERTH','Location','northeastoutside','FontSize',15)

nexttile(3)
% subplot(2,2,3)
plot3(pos_TOPO(:,1,1),pos_TOPO(:,2,1),pos_TOPO(:,3,1),'k')
plot3(pos_TOPO(1,1,1),pos_TOPO(1,2,1),pos_TOPO(1,3,1),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(pos_TOPO(end,1,1),pos_TOPO(end,2,1),pos_TOPO(end,3,1),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
title('@KOUROU TOPO','FontSize',15,'FontWeight','bold')
legend('KOUROU','orbit','r_{0}','r_{end}','Location','northeastoutside','FontSize',15)

nexttile(4)
% subplot(2,2,4)
plot3(pos_TOPO(:,1,2),pos_TOPO(:,2,2),pos_TOPO(:,3,2),'k')
plot3(pos_TOPO(1,1,2),pos_TOPO(1,2,2),pos_TOPO(1,3,2),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(pos_TOPO(end,1,2),pos_TOPO(end,2,2),pos_TOPO(end,3,2),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
title('@PERTH TOPO','FontSize',15,'FontWeight','bold')
legend('PERTH','orbit','r_{0}','r_{end}','Location','northeastoutside','FontSize',15)


figure('Name','2: Indexes of pericenter along time grid','NumberTitle','off')
spy(Te,'*')
pbaspect auto
grid on
box on
xticks(1)
yticks(TTe)
ylim([0 3601])
xlabel({['pericenter count = ',num2str(nnz(Te))]})
ylabel('i_{th} element of t_{grid}')
title('@ECI orbit reference','FontSize',15,'FontWeight','bold')
legend('pericenter indexes','FontSize',15)


figure('Name','3: Indexes for Visibility Windows','NumberTitle','off')
hold on
grid on
box on
xticks([1 2])
title('Retrieve indexes of T_{start} & T_{end} for each pass')
yyaxis left
spy([i_visibility(:,1) zeros(length(t_grid),1)],'*')
ylabel('i_{th} element of t_{grid} for KOUROU')
yticks([F_start(1),F_end(1),F_start(2),F_end(2),F_start(3),F_end(3),F_start(4),F_end(4),F_start(5),F_end(5),F_start(6),F_end(6),F_start(7),F_end(7),F_start(8),F_end(8)])
ax = gca;
ax.YTickLabel = {F_start(1),F_end(1),'','',F_start(3),F_end(3),F_start(4),F_end(4),'','',F_start(6),F_end(6),'','','',''};
yyaxis right
grid on
box on
spy([zeros(length(t_grid),1) i_visibility(:,2)],'r*')
ylabel('i_{th} element of t_{grid} for PERTH')
yticks([F_start(2),F_end(2),F_start(5),F_end(5),F_start(7),F_end(7),F_start(8),F_end(8)])
pbaspect auto
xlabel('')
legend(stationName)
ax.GridColor = [.15 .15 .15];


figure('Name','4: ECI orbit w/ visibility windows per each revolution','NumberTitle','off')
tiledlayout(2,2,"TileSpacing","compact","Padding","tight")
n_revolution = {'1st','2nd','3rd','4th'};
color = {'#0072BD','#A2142F'};
for j = 1:n_revolutions           % j_{th} revolution
    
    nexttile(j)
%     subplot(2,2,j)
    [X,Y,Z] = ellipsoid(0,0,0,re,re,rp,30);
    h = surf(X,Y,Z);
    h.EdgeColor = [0.07 0.62 1];
    h.FaceColor = 'none';
    % h.FaceColor = [0.07 0.62 1];
    % h.FaceAlpha = 0.5;
    % camlight
    hold on
    grid on
    box on
    axis equal
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    plot3(x_grid(TTe(j):TTe(j+1),1),x_grid(TTe(j):TTe(j+1),2),x_grid(TTe(j):TTe(j+1),3),'k')
    for i = 1:2
        plot(0,0,'Color',color{i},'LineWidth',2)
        plot3(vis_pos_ECI(TTe(j):TTe(j+1),1,i),vis_pos_ECI(TTe(j):TTe(j+1),2,i),vis_pos_ECI(TTe(j):TTe(j+1),3,i),'Color',color{i},'Marker','.','LineStyle','none')
    end
    view([0 90])
    plot(0,0,'Color',[0.07 0.62 1],'Marker','*','Linestyle','none')
    plot3(x_grid(1,1),x_grid(1,2),x_grid(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    title({'@ECI J2000 ',[n_revolution{j},' revolution']},'FontSize',15,'FontWeight','bold')
    ax = gca;
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    ax.ZAxis.Exponent = 4;

    switch j
        case {1}
            legend('Earth','s/c orbit','KOUROU visibility','','PERTH visibility','','','r_{0}','Location','northeast','FontSize',15)
        case {2}
            legend('Earth','s/c orbit','KOUROU visibility','','','','','pericenter','Location','northeast','FontSize',15)
        case {3}
            legend('Earth','s/c orbit','KOUROU visibility','','PERTH visibility','','','pericenter','Location','northeast','FontSize',15)
        case {4}
            plot3(x_grid(end,1),x_grid(end,2),x_grid(end,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
            legend('Earth','s/c orbit','','','PERTH visibility','','','pericenter','r_{end}','Location','northeast','FontSize',15)
    end
end

rev = 'next revolution';
figure('Name','5: ECEF orbit w/ position of stations on earth surface & visibility windows','NumberTitle','off')
[X,Y,Z] = ellipsoid(0,0,0,re,re,rp,30);
h = surf(X,Y,Z);
h.EdgeColor = [0.07 0.62 1];
h.FaceColor = 'none';
% h.FaceColor = [0.07 0.62 1];
% h.FaceAlpha = 0.5;
% camlight
hold on
grid on
box on
axis equal
xlabel('X [Km]','FontSize',15,'FontWeight','bold')
ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
view([0 90])
plot3(pos_ECEF(:,1),pos_ECEF(:,2),pos_ECEF(:,3),'k')
for i = 1:2
    plot(0,0,'Color',color{i},'LineWidth',2)
    plot3(vis_pos_ECEF(:,1,i),vis_pos_ECEF(:,2,i),vis_pos_ECEF(:,3,i),'Color',color{i},'Marker','.','LineStyle','none')
end
plot(0,0,'Color',[0.07 0.62 1],'Marker','*','Linestyle','none')
plot3(pos_ECEF(1,1),pos_ECEF(1,2),pos_ECEF(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(pos_ECEF(TTe(2:4),1),pos_ECEF(TTe(2:4),2),pos_ECEF(TTe(2:4),3),'LineStyle','none','Marker','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot3(pos_ECEF(end,1),pos_ECEF(end,2),pos_ECEF(end,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
plot3(pos_station_ECEF(1,1),pos_station_ECEF(1,2),pos_station_ECEF(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.9290 0.6940 0.1250])
plot3(pos_station_ECEF(2,1),pos_station_ECEF(2,2),pos_station_ECEF(2,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.72 0.27 1])
title('@ECEF ITRF93','FontSize',15,'FontWeight','bold')
legend('Earth','s/c orbit','KOUROU visibility','','PERTH visibility','','','r_{0}',rev,'r_{end}','KOUROU','PERTH','Location','northwest','FontSize',15)
ax = gca;
ax.XAxis.Exponent = 4;
ax.YAxis.Exponent = 4;
ax.ZAxis.Exponent = 4;


figure('Name','6: Time windows, Azimuth & Elevation','NumberTitle','off')
tiledlayout(4,1,"TileSpacing","compact","Padding","compact")
for i = 1:4
    nexttile(i)
%     subplot(4,1,i)
    hold on
    grid on
    box on
    xlim tight
end
for i = 1:2
    nexttile(1)
%     subplot(4,1,1)
        plot((t_grid-t0)/cspice_spd(),i_visibility(:,i),'Color',color{i},'LineWidth',2)
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[Logical]','FontSize',15,'FontWeight','bold')
        title('Visibility','FontSize',15,'FontWeight','bold')
        yticks([0,1])
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        ylim([-0.3 1.3])
    nexttile(2)
%     subplot(4,1,2)
        plot((t_grid-t0)/cspice_spd(),aaz(:,i),'Color',color{i},'LineWidth',2)
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[deg]','FontSize',15,'FontWeight','bold')
        title('Azimuth','FontSize',15,'FontWeight','bold')
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        yTicks = linspace(-180,180,5);
        ylim([yTicks(1)-1e1 yTicks(end)+1e1])
        yticks(yTicks);
    nexttile(3)
%     subplot(4,1,3)
        plot((t_grid-t0)/cspice_spd(),eel(:,i),'Color',color{i},'LineWidth',2)
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[deg]','FontSize',15,'FontWeight','bold')
        title('Elevation','FontSize',15,'FontWeight','bold')
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        yTicks = linspace(0,90,4);
        ylim([yTicks(1) yTicks(end)])
        yticks(yTicks)
    nexttile(4)   
%     subplot(4,1,4)
        hold on
        grid on
        box on
        plot(aaz(i_visibility(:,i),i),eel(i_visibility(:,i),i),'*','Color',color{i})
        xlabel('Azimuth [deg]','FontSize',15,'FontWeight','bold')
        ylabel('Elevation [deg]','FontSize',15,'FontWeight','bold')
        title('Polar','FontSize',15,'FontWeight','bold')
        xTicks = linspace(-180,180,13);
%         xlim([xTicks(1) xTicks(end)])
        xticks(xTicks); 
        xlim([-180 180])
        yTicks = linspace(0,90,4);
        ylim([yTicks(1) yTicks(end)])
        yticks(yTicks)
end
for i = 1:4
    nexttile(i)
%     subplot(4,1,i)
    switch i
        case 4
        otherwise
            xline((t_grid(TTe(2:4))-t0)/cspice_spd(),'--','LineWidth',1.8,'Color',[0.7 0.7 0.7])
    end
    switch i
        case {1,4}
            legend('KOUROU','PERTH',rev,'Location','east','FontSize',15,'FontWeight','bold')
        otherwise
            legend('KOUROU','PERTH',rev,'Location','northeast','FontSize',15,'FontWeight','bold')
    end
end

%
fprintf('\nExecution time ex2.1: %.2fs\n',toc)
tic

%% EX 2.2 Simulate measurements
%% Ex 2.2.a Compute position and expected measurements

% 3LE
% ARIANE 5 R/B   
l1 = '1 87654U 22110B   22316.00967942  .00000002  00000-0  32024-3 0  9990';
l2 = '2 87654   3.6309 137.1541 8138191 196.1632  96.6141  1.26411866   834';

% initialize the satrec structure
satrec = twoline2rv(l1,l2,'u','e','a',72);
ANS.TLE.satrec = satrec;

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
t_sat_utc = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
clear min
t_sat = cspice_str2et(t_sat_utc);
ANS.TLE.t_sat = t_sat;

t_since = t_grid - t_sat; % [s] -> consider SGP4 needs [min]
ANS.TLE.t_since = t_since';

% Evaluate the TLE along t_grid
% TEME
rteme = zeros(length(t_grid),3);
vteme = zeros(length(t_grid),3);
ateme = zeros(length(t_grid),3);
for i = 1:length(t_grid)
t_since = t_grid(i) - t_sat; % [s] -> consider SGP4 needs [min]
[satrec,rteme(i,:),vteme(i,:)] = sgp4(satrec,t_since/60); 

end
ANS.TLE.rteme = rteme;
ANS.TLE.vteme = vteme;

% % Centuries from TDT 2000 January 1 00:00:00.000
% ttt = cspice_unitim(t_grid(i),'ET','TDT')/cspice_jyear()/100;

% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

% Retrieve celestial pole offsets dPsi & dEps along t_grid (3 Days)
% from EARTH ORIENTATION PARAMETERS DATA at TLE reference epoch

T_grid = cspice_timout(t_grid,TFORMAT); % w/o offset for correct day indexing when choosing dPsi & dEps
dPsi = [-.113564 -.113334 -.112965] * arcsec2rad; % [rad]
dEps = [-.007029 -.006981 -.006952] * arcsec2rad; % [rad]

% Manually retrieve indexes of beginning/end of each day

% i1 = 1:1171;        % indexes day1
% i2 = 1172:2611;     % indexes day2
% i3 = 2612:3601;     % indexes day3
% n_days = 3;
% 
% I = zeros(max(diff([0 i1(end) i2(end) i3(end)])),n_days);
% for i = 1:3
%     I(1:length(i1),1) = i1;
%     I(1:length(i2),2) = i2;
%     I(1:length(i3),3) = i3;
% end

% Automatic process
[I,Day,n_days] = automatic_day_indexing(T_grid);

% Data collection
for i = 1:n_days 
    ANS.Data.Days(i).UTC = T1_grid(I(1,i),:);
    ANS.Data.Days(i).ephemeris = t_grid(I(1,i)) + 1e-4;
    ANS.Data.Days(i).index = I(1,i);
end

% ECI
reci = zeros(length(t_grid),3);
veci = zeros(length(t_grid),3);
for j = 1:n_days
    for i = nonzeros(I(:,j))'
        % Centuries from TDT 2000 January 1 00:00:00.000
        ttt = cspice_unitim(t_grid(i),'ET','TDT')/cspice_jyear()/100;
        [reci(i,:),veci(i,:)] = teme2eci(rteme(i,:)',vteme(i,:)',ateme(i,:)',ttt,dPsi(j),dEps(j));
    end
end

ANS.TLE.Propagation.ECI = reci;

% ECEF
recef = zeros(length(t_grid),3);
for i = 1:length(t_grid)
    recef(i,:) = ROT_ECI2ECEF(:,:,i) * reci(i,:)';
end
ANS.TLE.Propagation.ECEF = recef;

% TOPO
rtopo = zeros(length(t_grid),3,2);
for i = 1:length(t_grid)
    for j = 1:2
       rtopo(i,:,j) = ECEF2TOPO(:,:,j)*(recef(i,:)' - pos_station_ECEF(j,:)');
    end
end
ANS.TLE.Propagation.Topocentric.KOUROU = rtopo(:,:,1);
ANS.TLE.Propagation.Topocentric.PERTH = rtopo(:,:,2);

% Measurements
% One - way Range
RHO = zeros(length(t_grid),2);
AZ = zeros(length(t_grid),2);
EL = zeros(length(t_grid),2);
for i = 1:length(t_grid)
    for j = 1:2
        RHO(i,j) = sqrt((recef(i,:) - pos_station_ECEF(j,:))*(recef(i,:) - pos_station_ECEF(j,:))');
        AZ(i,j) = atan2d(rtopo(i,2,j),rtopo(i,1,j));
        EL(i,j) = asind(rtopo(i,3,j)/norm(rtopo(i,:,j)));
%         [AZ(i,j),EL(i,j),RHO] = antenna_pointing(stationName{j},t_grid(i),[x_grid(i,1:3)';x_grid(i,4:6)']);
%         [RHO(i,j), AZ(i,j),  EZ(i,j)] = cspice_reclat(pos_TOPO(i,:,j)');
    end
end

% Elevation visibility within time windows
eL = zeros(length(t_grid),2);
for i = 1:2
    eL(i_visibility(:,i),i) = EL(i_visibility(:,i),i);
end
% Visibility dependency on measured elevation lower boundary
real_visibility = false(length(t_grid),2);
real_visibility(:,1) = eL(:,1) > 10;
real_visibility(:,2) = eL(:,2) > 5;

rho = zeros(length(t_grid),2);
az = zeros(length(t_grid),2);
el = zeros(length(t_grid),2);
for i = 1:2
    rho(real_visibility(:,i),i) = RHO(real_visibility(:,i),i);
    az(real_visibility(:,i),i) = AZ(real_visibility(:,i),i);
    el(real_visibility(:,i),i) = EL(real_visibility(:,i),i);
end

% Compute visibility time grid
T = zeros(length(t_grid),2);
for i = 1:2
    T(real_visibility(:,i),i) = t_grid(real_visibility(:,i));
end
for i = 1:length(t_grid)
    ANS.TLE.Visibility.Time(i).KOUROU = T(i,1);
    ANS.TLE.Visibility.Time(i).PERTH = T(i,2);
end

% Compute visibility measured position in ECI & ECEF
Reci = zeros(length(t_grid),3,2);
Recef = zeros(length(t_grid),3,2);
for i = 1:2
    Reci(real_visibility(:,i),:,i) = reci(real_visibility(:,i),1:3);
    Recef(real_visibility(:,i),:,i) = recef(real_visibility(:,i),:);
end
for i = 1:length(t_grid)
    ANS.TLE.Visibility.ECI(i).KOUROU = Reci(i,:,1);
    ANS.TLE.Visibility.ECI(i).PERTH = Reci(i,:,2);
    ANS.TLE.Visibility.ECEF(i).KOUROU = Recef(i,:,1);
    ANS.TLE.Visibility.ECEF(i).PERTH = Recef(i,:,2);
end

%% Plot
 
figure(7)
close figure 7
figure(8)
close figure 8
figure(9)
close figure 9


figure('Name','7: ECI measured orbit w/ visibility windows per each revolution','NumberTitle','off')
% n_revolution = {'1st','2nd','3rd','4th'};
% color = {'#0072BD','#A2142F'};
tiledlayout(2,2,"TileSpacing","compact","Padding","compact")
for j = 1:n_revolutions           % j_{th} revolution
    nexttile(j)
%     subplot(2,2,j)
    [X,Y,Z] = ellipsoid(0,0,0,re,re,rp,30);
    h = surf(X,Y,Z);
    h.EdgeColor = [0.07 0.62 1];
    h.FaceColor = 'none';
    % h.FaceColor = [0.07 0.62 1];
    % h.FaceAlpha = 0.5;
    % camlight
    hold on
    grid on
    box on
    axis equal
    xlabel('X [Km]','FontSize',15,'FontWeight','bold')
    ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
    zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
    plot3(reci(TTe(j):TTe(j+1),1),reci(TTe(j):TTe(j+1),2),reci(TTe(j):TTe(j+1),3),'k')
    for i = 1:2
        plot(0,0,'Color',color{i},'LineWidth',2)
        plot3(Reci(TTe(j):TTe(j+1),1,i),Reci(TTe(j):TTe(j+1),2,i),Reci(TTe(j):TTe(j+1),3,i),'Color',color{i},'Marker','.','LineStyle','none')
    end
    view([0 90])
    plot(0,0,'Color',[0.07 0.62 1],'Marker','*','Linestyle','none')
    plot3(reci(1,1),reci(1,2),reci(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
    title({'@ECI J2000 ',[n_revolution{j},' revolution']},'FontSize',15,'FontWeight','bold')
    ax = gca;
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    ax.ZAxis.Exponent = 4;
    switch j
        case {1}
            legend('Earth','s/c orbit','KOUROU visibility','','PERTH visibility','','','r_{0}','Location','northeast','FontSize',15)
        case {2}
            legend('Earth','s/c orbit','KOUROU visibility','','','','','r_{0}','Location','northeast','FontSize',15)
        case {3}
            legend('Earth','s/c orbit','KOUROU visibility','','PERTH visibility','','','r_{0}','Location','northeast','FontSize',15)
        case {4}
            plot3(reci(end,1),reci(end,2),reci(end,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
            legend('Earth','s/c orbit','','','PERTH visibility','','','r_{0}','r_{end}','Location','northeast','FontSize',15)
    end
end


figure('Name','8: ECEF measure orbit w/ visibility windows','NumberTitle','off')
[X,Y,Z] = ellipsoid(0,0,0,re,re,rp,30);
h = surf(X,Y,Z);
h.EdgeColor = [0.07 0.62 1];
h.FaceColor = 'none';
% h.FaceColor = [0.07 0.62 1];
% h.FaceAlpha = 0.5;
% camlight
hold on
grid on
box on
axis equal
xlabel('X [Km]','FontSize',15,'FontWeight','bold')
ylabel('Y [Km]','FontSize',15,'FontWeight','bold')
zlabel('Z [Km]','FontSize',15,'FontWeight','bold')
view([0 90])
plot3(recef(:,1),recef(:,2),recef(:,3),'k')
for i = 1:2
    plot(0,0,'Color',color{i},'LineWidth',2)
    plot3(Recef(:,1,i),Recef(:,2,i),Recef(:,3,i),'Color',color{i},'Marker','.','LineStyle','none')
end
plot(0,0,'Color',[0.07 0.62 1],'Marker','*','Linestyle','none')
plot3(recef(1,1),recef(1,2),recef(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot3(recef((TTe(2:4)),1),recef((TTe(2:4)),2),recef((TTe(2:4)),3),'LineStyle','none','Marker','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot3(recef(end,1),recef(end,2),recef(end,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
plot3(pos_station_ECEF(1,1),pos_station_ECEF(1,2),pos_station_ECEF(1,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.9290 0.6940 0.1250])
plot3(pos_station_ECEF(2,1),pos_station_ECEF(2,2),pos_station_ECEF(2,3),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.72 0.27 1])
title('@ECEF ITRF93','FontSize',15,'FontWeight','bold')
legend('Earth','s/c orbit','KOUROU visibility','','PERTH visibility','','','r_{0}',rev,'r_{end}','KOUROU','PERTH','Location','southeast','FontSize',15)
ax = gca;
ax.XAxis.Exponent = 4;
ax.YAxis.Exponent = 4;
ax.ZAxis.Exponent = 4;


f9 = figure('Name','9: Simulated Range (one-way), Azimuth & Elevation','NumberTitle','off');
tiledlayout(6,1,"TileSpacing","compact","Padding","compact")
for i = 1:6
    nexttile(i)
%     subplot(6,1,i)
    hold on
    grid on
    box on
    xlim tight
end
for i = 1:2
    nexttile(1)
%     subplot(6,1,1)
        plot((t_grid-t0)/cspice_spd(),i_visibility(:,i),'Color',color{i},'LineWidth',2);
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[Logical]','FontSize',15,'FontWeight','bold')
        title('Simulation windows','FontSize',15,'FontWeight','bold')
        yticks([0,1])
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        ylim([-0.3 1.3])
    
    nexttile(2)    
%     subplot(6,1,2)
        plot((t_grid-t0)/cspice_spd(),real_visibility(:,i),'Color',color{i},'LineWidth',2)
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[Logical]','FontSize',15,'FontWeight','bold')
        title('Actual Visibility','FontSize',15,'FontWeight','bold')
        yticks([0,1])
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        ylim([-0.3 1.3])

    nexttile(3)    
%     subplot(6,1,3)
%         plot(t_grid/cspice_spd(),RHO(:,i),'b')
        plot((t_grid-t0)/cspice_spd(),rho(:,i),'Color',color{i},'LineWidth',2)
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[Km]','FontSize',15,'FontWeight','bold')
        title('Range measurement','FontSize',15,'FontWeight','bold')
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
    
    nexttile(4)
%     subplot(6,1,4)
%         plot(t_grid/cspice_spd(),Az(:,i),'b')
        plot((t_grid-t0)/cspice_spd(),az(:,i),'Color',color{i},'LineWidth',2)
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[deg]','FontSize',15,'FontWeight','bold')
        title('Azimuth measurement','FontSize',15,'FontWeight','bold')
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        yTicks = linspace(-180,180,5);
        ylim([yTicks(1)-1e1 yTicks(end)+1e1])
        yticks(yTicks);
    
    nexttile(5)
%     subplot(6,1,5)
%         plot(t_grid/cspice_spd(),El(:,i),'b')
        plot((t_grid-t0)/cspice_spd(),el(:,i),'Color',color{i},'LineWidth',2)
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[deg]','FontSize',15,'FontWeight','bold')
        title('Elevation measurement','FontSize',15,'FontWeight','bold')
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        yTicks = linspace(0,90,4);
        ylim([yTicks(1) yTicks(end)])
        yticks(yTicks)
    
    nexttile(6)   
%     subplot(6,1,6)
        hold on
        grid on
        box on
%         plot(Az(:,i),El(:,i),'*','Color','b')
        plot(az(real_visibility(:,i),i),el(real_visibility(:,i),i),'*','Color',color{i})
        xlabel('Azimuth [deg]','FontSize',15,'FontWeight','bold')
        ylabel('Elevation [deg]','FontSize',15,'FontWeight','bold')
        title('Polar','FontSize',15,'FontWeight','bold')
        xTicks = linspace(-180,180,13);
        xlim([xTicks(1) xTicks(end)])
        xticks(xTicks); 
        yTicks = linspace(0,90,4);
        ylim([yTicks(1) yTicks(end)])
        yticks(yTicks)
end
for i = 1:6
    nexttile(i)
%     subplot(6,1,i)
    switch i
        case 6
        otherwise
            xline((t_grid(TTe(2:4))-t0)/cspice_spd(),'--','LineWidth',1.8,'Color',[0.7 0.7 0.7])
    end
    switch i
        case {1,2,3,6}
            legend('KOUROU','PERTH',rev,'Location','east','FontSize',15)
        case 4
            legend('KOUROU','PERTH',rev,'Location','southeast','FontSize',15)
        case 5
            legend('KOUROU','PERTH',rev,'Location','northeast','FontSize',15)
    end
end


%
fprintf('Execution time ex2.2.a: %.2fs\n',toc)
tic

%% Ex. 2.2.b ADD NOISE

% Measurement noise data
sigma_Az_El = 100;  % [mdeg]
sigma_Az_El = sigma_Az_El*1e-3; % [deg]
sigma_range = 0.01; % [km]

% Diagonal noise matrix R
SIGMA = eye(3).*[sigma_Az_El^2,sigma_Az_El^2,sigma_range^2];
% SIGMA = diag([sigma_Az_El^2,sigma_Az_El^2,sigma_range^2]);
% SIGMA = [sigma_Az_El^2,sigma_Az_El^2,sigma_range^2];
MU = zeros(length(t_grid),3,2);
MMU = zeros(length(t_grid),3,2);
R = zeros(length(t_grid),3,2);
RR = zeros(length(t_grid),3,2);
NOISE = zeros(length(t_grid),3,2);

rng('default') % For reproducibility
for j = 1:2
    MMU(:,:,j) = [az(:,j),el(:,j),rho(:,j)];
    R(:,:,j) = mvnrnd(MMU(:,:,j),SIGMA);
    RR(i_visibility(:,j),:,j) = R(i_visibility(:,j),:,j);
    NOISE(:,:,j) = RR(:,:,j) - MMU(:,:,j);
end
for i = 1:length(t_grid)
    ANS.Noise(i).Mean = MMU(i,:);
    ANS.Noise(i).StandardDeviation = RR(i,:);
    ANS.Noise(i).Error = NOISE(i,:);
end

% PLOT noise
% colorN = [0.4940 0.1840 0.5560]; % violet
% colorN = "#EDB120";              % yellow
colorN = "#77AC30";              % dark green           

% Expand the noise so to visualize it better
% Naz = 8e1;
Nel = 3e1;
Nrho = 2e5;
Naz = Nel;

%%
mksz = 3.5;
f9.Name = '9: Simulated Range (one-way), Azimuth & Elevation + random noise over expected time windows';
for i = 1:2 
        nexttile(3)
%         subplot(6,1,3)
            plot((t_grid(i_visibility(:,i))-t0)/cspice_spd(),rho(i_visibility(:,i),i) + Nrho*NOISE(i_visibility(:,i),3,i),'Marker','.','Color',colorN,'MarkerSize',mksz,'LineStyle','none')
            plot(0,0,'Color',colorN,'LineWidth',2)
    
        nexttile(4)
%         subplot(6,1,4)
           plot((t_grid(i_visibility(:,i))-t0)/cspice_spd(),az(i_visibility(:,i),i) + Naz*NOISE(i_visibility(:,i),1,i),'Marker','.','Color',colorN,'MarkerSize',mksz,'LineStyle','none')
           plot(0,0,'Color',colorN,'LineWidth',2)

        
        nexttile(5)
%         subplot(6,1,5)
            plot((t_grid(i_visibility(:,i))-t0)/cspice_spd(),el(i_visibility(:,i),i) + Nel*NOISE(i_visibility(:,i),2,i),'Marker','.','Color',colorN,'MarkerSize',mksz,'LineStyle','none')
            plot(0,0,'Color',colorN,'LineWidth',2)
            ylim([-10 yTicks(end)])

     
       nexttile(6)
%      subplot(6,1,6)
            plot(az(real_visibility(:,i),i) + Naz*NOISE(real_visibility(:,i),1,i),el(real_visibility(:,i),i) + Nel*NOISE(real_visibility(:,i),2,i),'Marker','.','Color',colorN,'MarkerSize',mksz,'LineStyle','none')
            plot(0,-100,'Marker','.','MarkerSize',50,'Color',colorN,'LineStyle','none')
end
for i = 1:6
    nexttile(i)
%     subplot(6,1,i)
    switch i
        case {1,2}
            legend('KOUROU','PERTH',rev,'Location','east','FontSize',15)
        case 3
%             legend({'KOUROU','PERTH',['random' newline 'noise']},'Location','east')
%             legend({'KOUROU','PERTH',rev,'','','',['RND NOISE x' num2str(Nrho,'%.e')]},'Location','east','FontSize',15)
            legend({'','','','','','',['RND NOISE x' num2str(Nrho,'%.e')]},'Location','east','FontSize',15)
        
        case 4
%             legend({'KOUROU','PERTH',['random' newline 'noise']},'Location','southeast')
%             legend({'KOUROU','PERTH',rev,'','','',['RND NOISE x',num2str(Naz,'% .e')]},'Location','southeast','FontSize',15)
            legend({'','','','','','',['RND NOISE x',num2str(Naz,'% .e')]},'Location','southeast','FontSize',15)
        case 5
%             legend({'KOUROU','PERTH',['random' newline 'noise']},'Location','northeast')
%             legend({'KOUROU','PERTH',rev,'','','',['RND NOISE x',num2str(Nel,'% .e')]},'Location','southeast','FontSize',15)
            legend({'','','','','','',['RND NOISE x',num2str(Nel,'% .e')]},'Location','southeast','FontSize',15)
        case 6
%             legend({'KOUROU','PERTH',['random' newline 'noise']},'Location','east')
%             legend({'KOUROU','PERTH','',['RND NOISE x',num2str(Nel,'% .e')]},'Location','east','FontSize',15)
            legend({'','','',['RND NOISE x',num2str(Nel,'% .e')]},'Location','east','FontSize',15)
    end
end

%
fprintf('Execution time ex2.2.b: %.2fs\n',toc)
tic

%% EX 2.3 Solve the navigation problem
%% EX 2.3.a PERTH ONLY for optimization
tic
x_0 = [reci(1,:),veci(1,:)]';

lb = [];
ub = [];
opt = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Diagnostics','off','Display',...
                                                                                            ...'iter-detailed');
                                                                                            'off');

[x_A,resnorm,residual,exitflag(1),~,~,jac] = lsqnonlin(@cost,x_0,lb,ub,opt,'a',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO,RR);
Resnorm(1) = resnorm;

% Covariance

Jac = full(jac);
P_ls(:,:,1) = resnorm / (length(residual)-length(x_0)) .* inv(Jac'*Jac);

%
fprintf('Execution time ex2.3.a: %.2fs\n',toc)
tic

%% EX 2.3.b Both stations for optimization
tic
opt = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Diagnostics','off','Display',...
                                                                                            ...'iter-detailed');
                                                                                            'off');
% opt.FunctionTolerance = 1e-8;
% opt.StepTolerance = 1e-8;
[x_B,resnorm,residual,exitflag(2),~,~,jac] = lsqnonlin(@cost,x_0,lb,ub,opt,'b',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO,RR);
Resnorm(2) = resnorm;

% Covariance

Jac = full(jac);
P_ls(:,:,2) = resnorm / (length(residual)-length(x_0)) .* inv(Jac'*Jac);

%
fprintf('Execution time ex2.3.b: %.2fs\n',toc)
tic

%% EX 2.3.c Both stations + J2 perturbed motion
tic
opt = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Diagnostics','off','Display',...
                                                                                            ...'iter-detailed');
                                                                                            'off');
opt.FunctionTolerance = 1e-8;
opt.StepTolerance = 1e-8;
[x_C,resnorm,residual,exitflag(3),~,~,jac] = lsqnonlin(@cost,x_0,lb,ub,opt,'c',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO,RR);
Resnorm(3) = resnorm;
% Covariance

Jac = full(jac);
P_ls(:,:,3) = resnorm / (length(residual)-length(x_0)) .* inv(Jac'*Jac);

%
fprintf('Execution time ex2.3.c: %.2fs\n',toc)
tic

%% Data collection
tic
% Square root of covariance trace 
% & Eigenvalues of covariance
P_ls_acell = cell(6);
P_ls_bcell = P_ls_acell;
P_ls_ccell = P_ls_acell;
ANS.FIT.covmatrix.a = struct;
ANS.FIT.covmatrix.b = struct;
ANS.FIT.covmatrix.c = struct;

for i = 1:6
    for k = 1:6
        P_ls_acell(i,k) = {num2str(P_ls(i,k,1),'%.3e')};
        P_ls_bcell(i,k) = {num2str(P_ls(i,k,2),'%.3e')};
        P_ls_ccell(i,k) = {num2str(P_ls(i,k,3),'%.3e')};
    end
end
for i = 1:6
    ANS.FIT.covmatrix.a(i).one = P_ls_acell(i,1);
    ANS.FIT.covmatrix.a(i).two = P_ls_acell(i,2);
    ANS.FIT.covmatrix.a(i).three = P_ls_acell(i,3);
    ANS.FIT.covmatrix.a(i).four = P_ls_acell(i,4);
    ANS.FIT.covmatrix.a(i).five = P_ls_acell(i,5);
    ANS.FIT.covmatrix.a(i).six = P_ls_acell(i,6);

    ANS.FIT.covmatrix.b(i).one = P_ls_bcell(i,1);
    ANS.FIT.covmatrix.b(i).two = P_ls_bcell(i,2);
    ANS.FIT.covmatrix.b(i).three = P_ls_bcell(i,3);
    ANS.FIT.covmatrix.b(i).four = P_ls_bcell(i,4);
    ANS.FIT.covmatrix.b(i).five = P_ls_bcell(i,5);
    ANS.FIT.covmatrix.b(i).six = P_ls_bcell(i,6);

    ANS.FIT.covmatrix.c(i).one = P_ls_ccell(i,1);
    ANS.FIT.covmatrix.c(i).two = P_ls_ccell(i,2);
    ANS.FIT.covmatrix.c(i).three = P_ls_ccell(i,3);
    ANS.FIT.covmatrix.c(i).four = P_ls_ccell(i,4);
    ANS.FIT.covmatrix.c(i).five = P_ls_ccell(i,5);
    ANS.FIT.covmatrix.c(i).six = P_ls_ccell(i,6);
end
cmA = struct2table(ANS.FIT.covmatrix.a);
cmB = struct2table(ANS.FIT.covmatrix.b);
cmC = struct2table(ANS.FIT.covmatrix.c);

for i = 1:width(cmA)
    cmA.(i) = categorical(cmA.(i));
    cmB.(i) = categorical(cmB.(i));
    cmC.(i) = categorical(cmC.(i));
end
ANS.FIT.covmatrix.a = cmA;
ANS.FIT.covmatrix.b = cmB;
ANS.FIT.covmatrix.c = cmC;

for i = 1:3
    ANS.FIT.sqr_cov_tr(i).rr = num2str(trace(P_ls(1:3,1:3,i)),'%.6e');
    ANS.FIT.sqr_cov_tr(i).vv = num2str(trace(P_ls(4:6,4:6,i)),'%.6e');
end
eigA = eig(P_ls(:,:,1));
eigB = eig(P_ls(:,:,2));
eigC = eig(P_ls(:,:,3));

for j = 1:6
    ANS.FIT.x0(j).a = num2str(x_A(j),'%.6e');
    ANS.FIT.x0(j).b = num2str(x_B(j),'%.6e');
    ANS.FIT.x0(j).c = num2str(x_C(j),'%.6e');
    ANS.FIT.eig(j).a = num2str(eigA(j),'%.6e');
    ANS.FIT.eig(j).b = num2str(eigB(j),'%.6e');
    ANS.FIT.eig(j).c = num2str(eigC(j),'%.6e');
end

%% TEST SOLUTIONS

% First guess for A & B
[~,az0,el0,rho0] = cost(x_0,'b',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO);
% A solution
[~,azA,elA,rhoA] = cost(x_A,'b',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO);
% B solution
[~,azB,elB,rhoB] = cost(x_B,'b',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO);
% First guess for C
[~,az0C,el0C,rho0C] = cost(x_0,'c',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO);
% C solution
[~,azC,elC,rhoC] = cost(x_C,'c',t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO);

%% PLOT FITTING of solutions

figure(10)
close figure 10
figure('Name','10: Fit A, B & C solutions to simulated measurements','NumberTitle','off')
tiledlayout(3,1,"TileSpacing","compact","Padding","compact")

color_firstguess = "#EDB120";
color_A = "#77AC30";
color_B = "#4DBEEE";
color_firstguessC = "#7E2F8E";
color_C = "#D95319";

for i = 1:3
    nexttile(i)
%     subplot(3,1,i)
    hold on
    grid on
    box on
    xlim tight
end
for j = 1:2
    nexttile(1)
%     subplot(3,1,1)
        plot((t_grid(i_visibility(:,j))-t0)/cspice_spd(),rho(i_visibility(:,j),j) + NOISE(i_visibility(:,j),3,j),'Marker','*','Color',color{j},'MarkerSize',10,'LineStyle','none','LineWidth',2)
        plot(-1e3,0,'Marker','*','Color',color{2},'LineStyle','none','MarkerSize',10,'LineStyle','none','LineWidth',2)
        xline((t_grid(TTe(2:4))-t0)/cspice_spd(),'--','LineWidth',1.8,'Color',[0.7 0.7 0.7])
        plot((t_grid-t0)/cspice_spd,rho0(:,j),'Color',color_firstguess,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,rhoA(:,j),'Color',color_A,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,rhoB(:,j),'Color',color_B,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,rho0C(:,j),'Color',color_firstguessC,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,rhoC(:,j),'Color',color_C,'LineWidth',2)
        title('Range','FontSize',15,'FontWeight','bold')
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[Km]','FontSize',15,'FontWeight','bold')
        xlim([(t_grid(1)-t0)/cspice_spd() (t_grid(end)-t0)/cspice_spd()])
%---------------------------------------------------------------------------        
        % Example in report
        xlim([0.733960852	0.737583259])
        ylim([59123.98278	59251.59261])
%---------------------------------------------------------------------------    
    nexttile(2)
%     subplot(3,1,2)
        plot((t_grid(i_visibility(:,j))-t0)/cspice_spd(),az(i_visibility(:,j),j) + NOISE(i_visibility(:,j),1,j),'Marker','*','Color',color{j},'MarkerSize',10,'LineStyle','none','LineWidth',2) 
        plot(0,-1e3,'Marker','*','Color',color{2},'LineStyle','none','MarkerSize',10,'LineStyle','none','LineWidth',2)
        xline((t_grid(TTe(2:4))-t0)/cspice_spd(),'--','LineWidth',1.8,'Color',[0.7 0.7 0.7])
        plot((t_grid-t0)/cspice_spd,az0(:,j),'Color',color_firstguess,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,azA(:,j),'Color',color_A,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,azB(:,j),'Color',color_B,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,az0C(:,j),'Color',color_firstguessC,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,azC(:,j),'Color',color_C,'LineWidth',2)
        title('Azimuth','FontSize',15,'FontWeight','bold')
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[deg]','FontSize',15,'FontWeight','bold')
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        yTicks = linspace(-180,180,5);
        ylim([yTicks(1)-1e1 yTicks(end)+1e1])
        yticks(yTicks)
%---------------------------------------------------------------------------        
        % Example in report
        xlim([0.733960852	0.737583259])
        ylim([-16.13570612	-13.58179509])
        ax = gca;
        ax.YTickMode = 'auto';
%---------------------------------------------------------------------------    
    nexttile(3)    
%     subplot(3,1,3)
        plot((t_grid(i_visibility(:,j))-t0)/cspice_spd(),el(i_visibility(:,j),j) + NOISE(i_visibility(:,j),2,j),'Marker','*','Color',color{j},'MarkerSize',10,'LineStyle','none','LineWidth',2)
        plot(0,-1e3,'Marker','*','Color',color{2},'LineStyle','none','MarkerSize',10,'LineStyle','none','LineWidth',2)
        xline((t_grid(TTe(2:4))-t0)/cspice_spd(),'--','LineWidth',1.8,'Color',[0.7 0.7 0.7])
        plot((t_grid-t0)/cspice_spd,el0(:,j),'Color',color_firstguess,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,elA(:,j),'Color',color_A,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,elB(:,j),'Color',color_B,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,el0C(:,j),'Color',color_firstguessC,'LineWidth',2)
        plot((t_grid-t0)/cspice_spd,elC(:,j),'Color',color_C,'LineWidth',2)
        title('Elevation','FontSize',15,'FontWeight','bold')
        xlabel('t [days]','FontSize',15,'FontWeight','bold')
        ylabel('[deg]','FontSize',15,'FontWeight','bold')
%         xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
        yTicks = linspace(0,90,4);
        ylim([-10 yTicks(end)])
        yticks(yTicks)
%---------------------------------------------------------------------------
        % Example in report
        xlim([0.733960852	0.737583259])
        ylim([52.93138372	53.36384609])
        ax = gca;
        ax.YTickMode = 'auto';
%---------------------------------------------------------------------------
end
for i = 1:3
    nexttile(i)
%     subplot(3,1,i)
        legend('KOUROU simulation','PERTH simulation',rev,'','','first guess','A solution','B solution','J2 first guess','C solution','Location','east','FontSize',15)
end

%
fprintf('Execution time ex2.3 test: %.2fs\n',toc)

fprintf('Ex 2 total execution time: %.2fs\n',toc(ticex_2))

%% Clear Kernel Pool
cspice_kclear();


%% Functions

function [value,isterminal,direction]= eventfun(t,x)

r = [x(1);x(2);x(3)];
v = [x(4);x(5);x(6)];
value = dot(r,v);
isterminal = 0;
direction = 1;

end

%---------------------------------------------------------------

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

function dxdt = J2pert_keplerian(t,x)

Re = cspice_bodvrd('399','RADII',3);  
Re = Re(1);                           % [km] Mean Earth radius
mu = cspice_bodvrd('Earth','GM',1);   % [km^3/s^2] GM

% Keplerian Perturbed motion (J2)

J2 = 0.0010826269;

% SPICE Position method

pos_ECI = [x(1);x(2);x(3)];
ROT_ECI2ECEF = cspice_pxform('J2000','ITRF93',t); % rotation matrix ECI to ECEF convertion
pos_ECEF = ROT_ECI2ECEF*pos_ECI;
r = pos_ECEF;
aJ2 = 3/2*mu*J2*(r/norm(r)^3)*(Re/norm(r))^2.*(5*(r(3)/norm(r))^2-[1;1;3]);
ROT_ECEF2ECI = cspice_pxform('ITRF93','J2000',t);
aJ2_ECI = ROT_ECEF2ECI*aJ2;
norm_r_ECI = norm(pos_ECI);
dxdt = [x(4);x(5);x(6); - mu*pos_ECI/(norm_r_ECI^3) + aJ2_ECI];

end

%---------------------------------------------------------------

function [phi,x,tt,xe,te] = orbit_propagation(x0,t0,t,J2)
% outputs:  
%         phi = [N x 1] final column state x at time t (N = dim of state)
%           x = [I x N] row state propagation along orbit, (I = n. of states)
%          tt = [I x 1] time corresponding to each row of x
%          xe = [9 x N] row states at each epoch corresponding to eventfun
%          te = [9 x 1] time corresponding to each row of xe
% Special Input: 
%          J2 = add char vector 'J2' to select propagation with J2 acceleration
    
    if nargin == 3
        J2 = [];
    end
    if nargout >= 4 && (~strcmp(J2,'J2') || nargin == 3)
        options = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@eventfun);
        [tt,x,te,xe] = ode113(@keplerian,[t0 t],x0,options);
    elseif nargout < 4 && (~strcmp(J2,'J2') || nargin == 3)
        options = odeset('reltol', 1e-12, 'abstol', 1e-12);
        [tt,x] = ode113(@keplerian,[t0 t],x0,options);
    elseif nargout >= 4 && strcmp(J2,'J2') && nargin == 4
        options = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@eventfun);
        [tt,x,te,xe] = ode113(@J2pert_keplerian,[t0 t],x0,options);
    elseif nargout < 4 && strcmp(J2,'J2') && nargin == 4
        options = odeset('reltol', 1e-12, 'abstol', 1e-12);
        [tt,x] = ode113(@J2pert_keplerian,[t0 t],x0,options);
    end
    phi = x(end,:)';
end

%---------------------------------------------------------------

function [F_start,F_end,StationName,Windows] = automatic_visibility_table(i_visibility,stationName,t_grid)
% INPUT: 
%       i_visibility [Nxn] = N number of elements of t_grid, n number of stations, logical grid of visibility.
%       stationName [1xn] = cell array with names of each station.
%       t_grid [Nx1] = [s] time grid ephimeris time, required for output Windows only.
% OUTPUT:
%       F_start [n*kx1] = k counter of visibility windows per each station, vis. win. starting indexes
%       F_end [n*kx1] = visibility windows ending indexes.
%       StationName [n*kx1] = stations associated with each time window, in order
%       TABLE = visibility windows table, displayed in command window and
%               stored in 'ANS.Visibility.Windows' structure, needs t_grid input.

n = size(i_visibility,2); % number of stations
    for j = 1:n % station number
        F2 = 1;
        k = 1;  % visibility window counter for each station
        while 1
            I = i_visibility(F2:end,j);
            F1 = min(find(I)) + F2 - 1;
            if isempty(F1)
                break
            end
            f_start(k,j) = F1;
            I = i_visibility(F1:end,j);
            F2 = min(find(~I)) + F1 - 1;
            if isempty(F2)
                f_end(k,j) = length(i_visibility);
                break
            end
            f_end(k,j) = F2 - 1;
            k = k + 1;
        end
    end
    F_start = sort(nonzeros(f_start));
    F_end = sort(nonzeros(f_end));

    c_start = zeros(length(F_start),1);
    for i = 1:length(F_start)
        [~,c_start(i)] = find(F_start(i) == f_start);
    end
    StationName = cell(length(F_start),1);
    for i = 1:length(F_start)
        StationName{i} = stationName(c_start(i)); 
    end

    
    if nargout == 4

        T_start = cell(8,1);
        T_end = cell(8,1);
        TFORMAT = 'YYYY-MON-DD-HR:MN:SC.### ::UTC';
        
            for j = 1:length(F_start)
                T_start{j} = cspice_timout(t_grid(F_start(j)),TFORMAT);
                T_end{j} = cspice_timout(t_grid(F_end(j)),TFORMAT);
            end
        
        ANS.Visibility.Windows = struct;
        for i = 1:length(F_start)
            ANS.Visibility.Windows(i).Pass = i;
            ANS.Visibility.Windows(i).StationName = StationName{i};
            ANS.Visibility.Windows(i).StartTime_UTC = T_start{i};
            ANS.Visibility.Windows(i).EndTime_UTC = T_end{i};
        end
        Windows = struct2table(ANS.Visibility.Windows);
        Windows.Properties.VariableNames = {'Pass#','Station Name','Start Time(UTC)','End Time(UTC)'}
%         ANS.Visibility.Windows = Windows;
    end
end

%---------------------------------------------------------------

function [I,Day,n_days] = automatic_day_indexing(t_grid)
% INPUT: both are accepted 
%       t_grid [1xN] = [s] ephemeris time
%       T_grid [Nx24] = char array, N number of elements, UTC time grid,
%                       day must be a number 'DD' in position 10:11.
% OUTPUT:
%       I [nxm] = n max number of indexes in any of the days, m number of days,
%                 contains indexes of each epoch, divided per each column day.
%       Day [mx24] = char array, UTC time for beginning of orbit and next days.
%       n_days [1x1] = scalar, number of different days 
    
    if isnumeric(t_grid)
        TFORMAT = 'YYYY-MON-DD-HR:MN:SC.### ::UTC';
        T_grid = cspice_timout(t_grid,TFORMAT);
    elseif ischar(t_grid)
        T_grid = t_grid;
    end
    days = str2num(T_grid(:,10:11));
    ith = find(diff(days));
    n_days = length(ith)+1;
    Day = zeros(1,n_days);
    Day(1) = days(1);
    Day(end) = days(end);
    i1 = 1:ith(1);
    iend = ith(end)+1:length(T_grid);
    I = zeros(max(diff([0 ith' length(T_grid)])),n_days);
    I(1:length(i1),1) = i1';
    I(1:length(iend),end) = iend';
    for i = 1:n_days-2
        Day(1,i+1) = days(ith(i)+1);
        I(:,i+1) = ith(i)+1:ith(i+1);
    end
end

%---------------------------------------------------------------

function [residual,az,el,rho] = cost(x_0,sol,t_grid,i_visibility,pos_station_ECEF,ROT_ECI2ECEF,ECEF2TOPO,RR)

% Choose sol: char array, select between 'a','b','c'

if strcmp(sol,'a')
    J = 2;
    J2 = [];
elseif strcmp(sol,'b')
    J = 1:2;
    J2 = [];
elseif strcmp(sol,'c')
    J = 1:2;
    J2 = 'J2';
end

% cost function for non linear least-squares problems
    residual = [];
% ECI PROPAGATION with pure keplerian motion
    x_propagation_eci = zeros(6,length(t_grid),2);
    for j = J
        t0 = t_grid(1);
        x0 = x_0;
        for i = find(i_visibility(:,j))'
            x_propagation_eci(:,i,j) = orbit_propagation(x0,t0,t_grid(i),J2);
            x0 = x_propagation_eci(:,i,j);
            t0 = t_grid(i);
        end

    end
% ECEF
    recef = zeros(3,length(t_grid));
    for j = J
        for i = find(i_visibility(:,j))'
            recef(:,i,j) = ROT_ECI2ECEF(:,:,i) * x_propagation_eci(1:3,i,j);
        end
    end
% TOPO
    rtopo = zeros(3,length(t_grid),2);
    for j = J
        for i = find(i_visibility(:,j))'
            rtopo(:,i,j) = ECEF2TOPO(:,:,j)*(recef(:,i,j) - pos_station_ECEF(j,:)');
        end
    end

    
% PREDICTED MEASUREMENTS

    RHO = zeros(length(t_grid),2);
    AZ = zeros(length(t_grid),2);
    EL = zeros(length(t_grid),2);
    for j = J
        for i = find(i_visibility(:,j))'
            RHO(i,j) = sqrt((recef(:,i,j) - pos_station_ECEF(j,:)')'*(recef(:,i,j) - pos_station_ECEF(j,:)'));
            AZ(i,j) = atan2d(rtopo(2,i,j),rtopo(1,i,j));
            EL(i,j) = asind(rtopo(3,i,j)/norm(rtopo(:,i,j)));
        end
    end
% Visibility dependency on measured elevation lower boundary
    real_visibility = false(length(t_grid),2);
    real_visibility(:,1) = EL(:,1) > 10;
    real_visibility(:,2) = EL(:,2) > 5;

% Compute visible range, azimuth & elevation
    rho = zeros(length(t_grid),2);
    az = zeros(length(t_grid),2);
    el = zeros(length(t_grid),2);
    for j = J
        rho(real_visibility(:,j),j) = RHO(real_visibility(:,j),j);
        az(real_visibility(:,j),j) = AZ(real_visibility(:,j),j);
        el(real_visibility(:,j),j) = EL(real_visibility(:,j),j);
    end

% COMPUTE RESIDUAL
    if nargin == 8
        W_m = diag(1./[100*1e-3 100*1e-3 1e-2]);
% Obtain indexes over time grid,so to build one vector for each measurement
        k = [];
        aZ = [];
        eL = [];
        rhO = [];
        sim = [];
        for j = J
            k = sort([k;find(i_visibility(:,j))]);
        end
        for i = k'
            if isempty(nonzeros(az(i,:)))
                aZ = [aZ;0];
                eL = [eL;0];
                rhO = [rhO;0];
            else
                aZ = [aZ;nonzeros(az(i,:))];
                eL = [eL;nonzeros(el(i,:))];
                rhO = [rhO;nonzeros(rho(i,:))];
            end
            if isempty(nonzeros(RR(i,:,:)))
                sim = [sim;0,0,0];
            else
                sim = [sim;nonzeros(RR(i,:,:))'];
            end
        end

        for i = 1 : length(aZ)
            diff_meas_weigthed = W_m*([aZ(i);eL(i);rhO(i)] - sim(i,:)');
            residual = [residual;diff_meas_weigthed];
        end
    end
end
