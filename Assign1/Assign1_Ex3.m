% Spacecraft Guidance and Navigation (2022/2023)
% Assignment #1 - Guidance - Ex 3 Continuous guidance
% Author: Matteo Luciardello Lecardi
% Personal code: 10572870
% Matricola MSc Degree: 977701
% Total expected execution time: variable, on average 20s ~ 60s for Nth = [4,3]

%% Load kernels

% cd 'D:\Courses\5\1 Spacecraft guidance and navigation\2022_23\Assignment01'
cspice_furnsh 'assignment.tm';
clearvars; close all; clc; format long
% opengl hardware
ex3tic = tic;
% warning off

%% CONVERSION of UNITS
% execute orbit_propagation with [KM,S,KG] --> [AU,YEARS,KG]
% instead of using [AU,YEARS,KG] for input {λ0,tf} only


%% EX 3.1
ticsolve = tic;
Nth = [4;3]; % Number of active thrusters (max 4 per execution)
% Nth = [4;3;2;1]; % examine all faulties' cases
% Nth = [6;8;10];
% Nth = 50;
m0 = 1500;  % [kg] s/c mass at t0

T0 = '2022-08-03-12:45:20.000 UTC';
t0 = cspice_str2et(T0); % [s]
x0Earth = cspice_spkezr('Earth',t0,'ECLIPJ2000','none','Sun'); % [km][Km/s]

%     x0Earth(1:3) = cspice_convrt(x0Earth','KM','AU')';
    km2au = cspice_convrt(1,'KM','AU');
    sec2year = cspice_convrt(1,'SECONDS','YEARS');

    x0Earth(1:3) = x0Earth(1:3)*km2au;
    x0Earth(4:6) = x0Earth(4:6)*km2au/sec2year;
    t0 = t0*sec2year;

% xMars = cspice_spkezr('Mars',tf,'ECLIPJ2000','none','Sun');

x0 = [x0Earth;m0;t0]; % initial boundary conditions + t0



%% FSOLVE

opt = optimoptions('fsolve', ...
               ... refinements for best results, using F [8x1] w/ hemiltonian
                   ...'FunctionTolerance',9e-44, ...
                   ...'OptimalityTolerance',8e-10, ...
                   ...
               ... refinements for F [7x1] w/o hemiltonian
               ... for 'levenberg-marquardt' algorithm, OptTol = 1e-4*FunTol
               ... 'FunctionTolerance',1e-7, ... refinement which causes more exflag 4 
               ... good fix if exitflag = 4 --> cannot reach F-O-Opt 1e-10
               ... add if doing Nth 2 and 1
                  'FunctionTolerance',1e-5, ... 
               ... good fix if exitflag = 4 --> reduces the step faster than F-O-Opt
                   'StepTolerance',1e-50, ...
                   ...
                   ...'FiniteDifferenceType','central', ...
                   ...'FiniteDifferenceStepSize',1e-10, ...
                   ...
                   ...'MaxIterations',1e6, ...
                   'MaxFunctionEvaluations',2e3, ...
                   'Display', ...
                   'off');
                   ...'iter-detailed');


e = 0;
disc = 0;
for k = 1:length(Nth)
    n = 0;

ticNth = tic;
    while 1
        tic
        n = n + 1;

% F [7x1] --> w/o hemiltonian (Hf-HM), nor Hf
          x = rand(7,1); % {λ0} --> [AU,Kg]

%   x = a + (b-a)*rand, search x between interval (a,b)
%   search Δt within interval (t1,t2) months

          if Nth(k) == 4 
              t1 = 7;  % [months]
              t2 = 10; % [months]
          elseif Nth(k) == 3
              t1 = 12;  
              t2 = 15;
          elseif Nth(k) == 2
              t1 = 21;
              t2 = 24;
          elseif Nth(k) == 1
              t1 = 34;
              t2 = 37;
          else
              t1 = 0;
              t2 = 7;
          end
          x(8) = t1/12 + (t2/12 - t1/12)*rand; % {Δtf} --> [YEARS]


          %             x(8) = abs(x(8));
%         end
%         best first guess, good for both configurations with 4 & 3 trhusters
%         FunctionTol = 9e-44;
%         OptimalityTol = 8e-10;
%         F [8x1] --> hemiltonian Hf - HM
%         x = [-0.12193514258764015512070955082891 0.08786591446944001793095679886392 0.00608782895565807957022341767583 -0.02496348110866202299473748382752 -0.01166378555229282548699742960707 -0.00033133353382922621933187601329 0.00637121163311578277838043149472 0.76825865619307276421778851727140]';
           
%         F [7x1] --> w/o hemiltonian Hf - HM
          % x = [0.83719859616159586490624633370317 0.75296087798745203745198750766576 0.36452816913868701931278337724507 0.59964360141248329583163467759732 0.34273249154068041733012250915635 0.43647975954081552618646355767851 0.25822074929657723441067673775251 0.97940044300902995288993224676233]';


        D = char(916);
        L = char(955);

        inguess = sprintf(['\nNumber of Thrusters: %d\n\nFIRST GUESS:\n'...  
                 '%sr: %.10f %.10f %.10f; \n%sv: %.10f %.10f %.10f; \n' ...
                 '%sm: %.10f; %st: %.10f[years]\n' ...
                 'x = [%.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f]'';\n\n'],...
                 Nth(k),L,x(1:3),L,x(4:6),L,x(7),D,x(8),x);
        
        fprintf('%s',inguess)
%         tic
        [y,F,exitflag,output] = fsolve(@boundariesfun,x,opt,Nth(k),x0);
        y(7) = abs(y(7));

        foo = output.firstorderopt;
        Tf = cspice_convrt(t0+y(8),'YEARS','SECONDS');
        Tf = cspice_timout(Tf,'YYYY-MM-DD-HR:MM:SC.### UTC');
        outguess = sprintf(['INITIAL CONDITION:\n'...  
                 '%sr: %.10f %.10f %.10f; \n%sv: %.10f %.10f %.10f; \n' ...
                 '%sm: %.10f; %st: %.10f [years]\nArrival Date: %s\n' ...
                 'y = [%.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f]'';\n\n'],...
                 L,y(1:3),L,y(4:6),L,y(7),D,y(8),Tf,y);
        
        
%         fprintf('%s',inguess)
        fprintf('exit flag: %d\n\n',exitflag)
        
        fprintf('%s\n',outguess)

        
        % errors % [AU,YEARS]
        Er = abs(F(1:3)); % position error
        Ev = abs(F(4:6)); % velocity error
        er = norm(Er);
        ev = norm(Ev);
        r = norm(F)^2; 
        au2km = cspice_convrt(1,'AU','KM');
        year2sec = cspice_convrt(1,'YEARS','SECONDS');
        error = sprintf(['errors:\n' ...
                 'er: %.10e[Km] %.10e[Km] %.10e[Km]\n', ...
                 'ev: %.10e[Km/s] %.10e[Km/s] %.10e[Km/s]\n' ...
                 '||er||: %.10e[Km]; ||ev||: %.10e[Km/s]\n||F||^2: %.10e\n' ...
                 'First Order Optimality: %.10e\n'],Er*au2km,Ev*au2km/year2sec,er*au2km,ev*au2km/year2sec,r,foo);
        fprintf('%s',error)
%         fprintf(['\nExecution time input combination number: %d  --> ' ...
%                  '%.2fs\nElapsed so far %.2fs\n\n'],n,toc,toc(ticsolve))
        
        ANS(k).initialCondition.firstguess_x = x;
        ANS(k).initialCondition.y = y;
        ANS(k).errors.Er = Er;
        ANS(k).errors.Er_km = Er*au2km;
        ANS(k).errors.er = er;
        ANS(k).errors.Ev = Ev;
        year2sec = cspice_convrt(1,'YEARS','SECONDS');
        ANS(k).errors.Ev_kms = Ev*au2km/year2sec;
        ANS(k).errors.FirstOderOpt = foo;
        ANS(k).inguess = inguess;
        ANS(k).outguess = outguess;
        ANS(k).error = error;

        fprintf('Execution time N. %d Thrusters try n. %d = %.2fs\n',Nth(k),n,toc)

        if exitflag == 1
            fprintf('\nExecution time N. %d Thrusters = %.2fs\n\n\n',Nth(k),toc(ticNth))
            break
        elseif exitflag == 0 || exitflag == -2 || exitflag == -3
            e = e + 1; % number of mistakes
        elseif exitflag == 4 || exitflag == 3 || exitflag == 2
            disc = disc + 1; % number of discarded good guesses
        end
    end
end
fprintf('\nExecution time fsolve (%d bad guess/es, %d good discarded) = %.2fs\n',e,disc,toc(ticsolve))

%% ERROR TABLE
errortab = struct;
errorTab = table;
prec = '%.10e';
for k = 1:length(Nth)
    errortab(k).Nth = Nth(k);
    errortab(k).FOO = num2str(ANS(k).errors.FirstOderOpt,'%.3e');
    errortab(k).dt = num2str(ANS(k).initialCondition.y(8)*12,'%.2f');
    errortab(k).r = num2str(norm(ANS(k).errors.Er_km),prec);
    errortab(k).v = num2str(norm(ANS(k).errors.Ev_kms),prec);
    errortab(k).rx = num2str(ANS(k).errors.Er_km(1),prec);
    errortab(k).ry = num2str(ANS(k).errors.Er_km(2),prec);
    errortab(k).rz = num2str(ANS(k).errors.Er_km(3),prec);
    errortab(k).vx = num2str(ANS(k).errors.Ev_kms(1),prec);
    errortab(k).vy = num2str(ANS(k).errors.Ev_kms(2),prec);
    errortab(k).vz = num2str(ANS(k).errors.Ev_kms(3),prec);
end
errorTab = struct2table(errortab);
for i = 1 : width(errorTab)
    if height(errorTab) == 1
        break
    end
    errorTab.(i) = categorical(errorTab.(i));
end

errorTab.Properties.VariableNames = {'N. Thrusters','First Order Optimality', ...
    [D,'t[months]'],['||',D,'r||[km]'],['||',D,'v||[km/s]'],[D,'r_x[km]'], ...
    [D,'r_y[km]'],[D,'r_z[km]'],[D,'v_x[km/s]'],[D,'v_y[km/s]'],[D,'v_z[km/s]']}

ANS(1).errortab = errorTab;
ANS = orderfields(ANS, [6,1:5]);

%% Propagation
%--------------------------------------------------------------------
% % EXAMPLEs OF POSSIBLE INITIAL CONDITIONS
% y = zeros(8,4);
% Nth = [4;3;2;1];
% % 4 thrusters
% y(:,1) = [-17.39653860539522511885479616466910 12.53448590369614557005206734174863 0.86353989597166680436401975384797 -3.55887968025170220442987556452863 -1.65818669525369322137464678235119 -0.04601682335200094758898003988179 0.00205368265645299811344992590989 0.77161174258527831515408479390317]';
% 
% % 3 thrusters
% y(:,2) = [-3.34115510846615615392352083290461 2.76490762969343739641203683277126 0.04568214904854188407323789533621 -0.69363416071156824838084276052541 -0.36542856256628292177524031103530 0.00533118315384050345945299298478 0.00034638942553806865640164680542 1.20090864058272961401030443084892]';
% 
% % 2 thrusters
% % y(:,3) = [-20.68526015504369652830973791424185 24.79591121589701785410397860687226 0.04205196500014871802131111167000 -4.42170823728118378426188428420573 -3.57868332823426849031989149807487 0.00295765750642380548934839978870 0.00319912517502370964683766629832 1.88012633052037791969723912188783]';
% % y(:,3) = [4.18064550407858703806596167851239 -3.69229936214240650471651861153077 0.26197659802934941675189861598483 0.88256105639763104964146123165847 0.54895083501713048956816010104376 0.00904705591089080204958428055306 0.00189886945509386826681785809257 1.94710317879797045748091477435082]';
% 
% % 1 thruster
% % y(:,4) = [1.36558913569362871953671856317669 -0.64130208653953280073523046667106 0.28628248900086128436726085055852 0.13442113316601098094160704476963 -0.05426820654042494029845045133698 0.01093973110993226756371132779577 0.00087735128623492124582544660072 2.99278221104349206882488942937925]';
% % y(:,4) = [7.37532307797099306867494306061417 -3.46356744880173028633407739107497 1.54616479640205684020770604547579 0.72598650625345739761939967138460 -0.29309368802039076351562130184902 0.05908369451444801107031778997225 0.00473843048376602037075056017557 2.99278221104116548545448495133314]';
% 
% for k = 1:length(Nth) 
%     ANS(k).initialCondition.y = y(:,k);
% end
%--------------------------------------------------------------------


km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');

for k = 1:length(Nth)
    y = ANS(k).initialCondition.y;
    tf = cspice_convrt(t0+y(8),'YEARS','SECONDS');
    xM = cspice_spkezr('Mars',tf,'ECLIPJ2000','none','Sun');

    xM(1:3) = xM(1:3)*km2au;
% xM(4:6) = xM(4:6)*km2au/sec2year;


    [~,out,t_out] = orbit_propagator([x0(1:7);y(1:7)],t0,t0+y(8),Nth(k));
    
    % costates
    Lr = out(:,8:10);
    Lv = out(:,11:13);
    Lm = out(:,14);

    % norm of costates
    lr = zeros(length(t_out),1);
    lv = zeros(length(t_out),1);
    for j = 1:length(t_out)
        lr(j) = norm(Lr(j,:));
        lv(j) = norm(Lv(j,:));
    end

    % primer vector

    alpha = - Lv ./ lv;
    % alpha = zeros(length(t_out),3);
    % for j = 1:height(out)
    %     alpha(j,:) = - Lv(j,:)/norm(Lv(j,:)); % primer vector
    % end
    ANS(k).out = out;
    ANS(k).t_out = t_out;
    ANS(k).Lr = Lr;
    ANS(k).lr = lr;
    ANS(k).Lv = Lv;
    ANS(k).lv = lv;
    ANS(k).Lm = Lm;
    ANS(k).alpha = alpha;
    ANS(k).xM.xM = xM;
end
ANS = orderfields(ANS, [1:2,15,3:14]);

%% PLOT Orbit

color_orbit = ["b","k",'r','g'];
color = ["#0072BD","#A2142F",'r','g'];
% Environment
figure(1)
close(1)
f = figure('Name','1: Ex3 Heliocentric Earth-Mars transfer','NumberTitle','off');
hold on
grid on
box on
axis equal tight
zlim([-1,1])
xlabel('X [AU]','FontSize',15)
ylabel('Y [AU]      ','FontSize',15,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
zlabel('Z [AU]','FontSize',15)
title('ECLIPJ2000 @Sun','FontSize',15)

plot3(0,0,0,'Marker','o','MarkerFaceColor','#EDB120','MarkerEdgeColor','k','MarkerSize',12,'LineStyle','none')
% RS = cspice_bodvrd('Sun','RADII',3)*km2au;
% [X,Y,Z] = ellipsoid(0,0,0,RS(1),RS(1),RS(1),50);
% surf(X,Y,Z,'EdgeColor','#EDB120','FaceColor','none')

plot3(x0Earth(1),x0Earth(2),x0Earth(3),'Marker','o','MarkerFaceColor',"#0072BD",'MarkerEdgeColor','k','MarkerSize',10,'LineStyle','none')


STR_ = {'\textbf{Sun}','$\mathbf{Earth\ @t_0}$','$\mathbf{Mars\ @t_f}$'};
STR = STR_;    
str = cell(1,1);
thr = cell(1,1);

for k = 1:length(Nth)
    xM = ANS(k).xM.xM;
    out = ANS(k).out;
    thr{k} = num2str(Nth(k));
    
    plot3(xM(1),xM(2),xM(3),'Marker','o','MarkerFaceColor',"#D95319",'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
    plot3(out(:,1),out(:,2),out(:,3),'Color',color_orbit(k),'LineWidth',2)
%     plot3(out(end,1),out(end,2),out(end,3),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4,'LineStyle','none')
    str(2*k-1:2*k,:) = {['\boldmath$\underline{r}\textrm{\ :\textsf{\ ',thr{k},'}\ thrusters}$'],''};
    STR = {STR{:},str{2*k-1:2*k}};

    legend(STR,'Location','south', ...
        'Interpreter','latex',...
        'FontSize',17);
end

%% Prime Vector
figure(1)
f.Name = sprintf('1: Ex3 Heliocentric Earth-Mars transfer & primer vector %s',char(945));
STR = STR(1:end-1);
str = cell(1,1);
title('\textbf{Primer Vector / Thrust Angles} \boldmath${\underline{\hat{\alpha}}^*}$','\textbf{ECLIPJ2000 @Sun}',...
      'FontSize',17,'Interpreter','latex')
for k = 1:length(Nth)
    out = ANS(k).out;
    alpha = ANS(k).alpha;
    quiver3(out(:,1),out(:,2),out(:,3),alpha(:,1),alpha(:,2),alpha(:,3),'Color',color(k),'LineWidth',1);
    
    str(k,:) = {['\boldmath${\underline{\hat{\alpha}}^*}\textrm{:\textsf{\ ',thr{k},'}\ thrusters}$']};
    STR = {STR{:},str{k}};

    legend(STR,'Location','south', ...
           'Interpreter','latex',...
           'FontSize',17);
end

%% Vectorial Costates

figure(2)
close(2)
L = char(955);
figure('Name',sprintf('2: costate %sr',L),'NumberTitle','off')
hold on
grid on
box on
axis equal tight
zlim([-1,1])
xlabel('X [AU]','FontSize',15)
ylabel('Y [AU]      ','FontSize',15,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
zlabel('Z [AU]','FontSize',15)
title('\boldmath${\underline{\lambda}_r}$','\textbf{ECLIPJ2000 @Sun}',...
      'FontSize',17,'Interpreter','latex')


plot3(0,0,0,'Marker','o','MarkerFaceColor','#EDB120','MarkerEdgeColor','k','MarkerSize',12,'LineStyle','none')
plot3(x0Earth(1),x0Earth(2),x0Earth(3),'Marker','o','MarkerFaceColor',"#0072BD",'MarkerEdgeColor','k','MarkerSize',10,'LineStyle','none')

for k = 1:length(Nth)
    xM = ANS(k).xM.xM;
    plot3(xM(1),xM(2),xM(3),'Marker','o','MarkerFaceColor',"#D95319",'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
end
for k = 1:length(Nth)
    out = ANS(k).out;
    plot3(out(:,1),out(:,2),out(:,3),'Color',color_orbit(k),'LineWidth',2)
end
for k = 1:length(Nth)
    out = ANS(k).out;
    Lr = ANS(k).Lr;
    quiver3(out(:,1),out(:,2),out(:,3),Lr(:,1),Lr(:,2),Lr(:,3),'Color',color(k),'LineWidth',1);
end


STR_ = {'\textbf{Sun}','$\mathbf{Earth\ @t_0}$','$\mathbf{Mars\ @t_f}$'};
for k = 1:length(Nth)
    str0(k,:) = {''};
    str1(k,:) = {['\boldmath$\underline{r}\textrm{\ :\textsf{\ ',thr{k},'}\ thrusters}$']};
    str2(k,:) = {['\boldmath${\underline{\lambda}_r}\textrm{:\textsf{\ ',thr{k},'}\ thrusters}$']};
    str3(k,:) = {['\boldmath${\underline{\lambda}_v}\textrm{:\textsf{\ ',thr{k},'}\ thrusters}$']};
end
STR = {STR_{:},str0{1:end-1},str1{:},str2{:}};
legend(STR,'Location','south', ...
       'Interpreter','latex',...
       'FontSize',17);


figure(3)
close(3)
figure('Name',sprintf('3: costate %sv',L),'NumberTitle','off')
hold on
grid on
box on
axis equal tight
zlim([-1,1])
xlabel('X [AU]','FontSize',15)
ylabel('Y [AU]      ','FontSize',15,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
zlabel('Z [AU]','FontSize',15)
title('\boldmath${\underline{\lambda}_v}$','\textbf{ECLIPJ2000 @Sun}', ...
      'FontSize',17,'Interpreter','latex')


plot3(0,0,0,'Marker','o','MarkerFaceColor','#EDB120','MarkerEdgeColor','k','MarkerSize',12,'LineStyle','none')
plot3(x0Earth(1),x0Earth(2),x0Earth(3),'Marker','o','MarkerFaceColor',"#0072BD",'MarkerEdgeColor','k','MarkerSize',10,'LineStyle','none')


for k = 1:length(Nth)
    xM = ANS(k).xM.xM;
    plot3(xM(1),xM(2),xM(3),'Marker','o','MarkerFaceColor',"#D95319",'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
end
for k = 1:length(Nth)
    out = ANS(k).out;
    plot3(out(:,1),out(:,2),out(:,3),'Color',color_orbit(k),'LineWidth',2)
end
for k = 1:length(Nth)
    out = ANS(k).out;
    Lv = ANS(k).Lv;
    quiver3(out(:,1),out(:,2),out(:,3),Lv(:,1),Lv(:,2),Lv(:,3),'Color',color(k),'LineWidth',1);
end

STR = {STR_{:},str0{1:end-1},str1{:},str3{:}};
legend(STR,'Location','south', ...
       'Interpreter','latex',...
       'FontSize',17);

%% Plot Mass / Switching function & scalar costates
for k = 1:length(Nth)
    out = ANS(k).out;
    m = out(:,7)/m0; % [-] mass ratio
    Lv = ANS(k).Lv;
    Lm = ANS(k).Lm;
    St = zeros(length(m),1);
    for j = 1:length(out)
        St(j) = switchingfun(m(j),Lv(j,:),Lm(j,:));
    end
    ANS(k).m = m;
    ANS(k).St = St;
end

figure(4)
close(4)
figure('Name','4: Thrust ratio, Mass ratio, Switching function & scalar costates','NumberTitle','off')

t = tiledlayout(6,1,'TileSpacing','tight','Padding','compact');
title(t,'Time-Optimal quantities & costates','FontSize',15,'FontWeight','bold')
xlabel(t,'time[years]','FontSize',15,'FontWeight','bold')
for i = 1:6
    nexttile(i)
    grid on
    hold on
    box on
    xlim tight
    ylim padded
end

    Tr = max(Nth/4);
for k = 1:length(Nth)
    t_out = ANS(k).t_out;
    nexttile(1)
    tr = ones(length(t_out),1)/4*Nth(k);
    plot(t_out-t0,tr,'Color',color(k),'LineWidth',1.5);
%     plot(t_out-t0,0.75*ones(length(t_out),1),'Color',color(2),'LineWidth',1.5);
    title('Thrust ratio (ref. 4 thrusters)','FontSize',14)
    ylim([0 Tr+0.2])
    yticks([0 1])

    nexttile(2)
    m = ANS(k).m;
    plot(t_out-t0,m,'Color',color(k),'LineWidth',1.5);
    title('Mass ratio (ref. m_0)','FontSize',14)
    nexttile(3)
    St = ANS(k).St;
    plot(t_out-t0,St,'Color',color(k),'LineWidth',1.5)
    title('Switching function S_t','FontSize',14)
    nexttile(4)
    lr = ANS(k).lr;
    plot(t_out-t0,lr,'Color',color(k),'LineWidth',1.5);
    title('\boldmath$|\!|{\underline{\lambda}_r}|\!|$','Interpreter','latex','FontSize',17)
    nexttile(5)
    lv = ANS(k).lv;
    plot(t_out-t0,lv,'Color',color(k),'LineWidth',1.5);
    title('\boldmath$|\!|{\underline{\lambda}_v}|\!|$','Interpreter','latex','FontSize',17)
    nexttile(6)
    Lm = ANS(k).Lm;
    plot(t_out-t0,Lm,'Color',color(k),'LineWidth',1.5);
    title('\boldmath${\lambda_m}$','Interpreter','latex','FontSize',17)
    str(k,:) = {[thr{k},' thrusters']};
end

for i = 1:6
    nexttile(i)
%     legend('4 thrusters','3 thrusters','Location','east','FontSize',14)
    legend(str{:},'Location','east','FontSize',14)
end


fprintf('\nEx 3 execution time = %.2fs\n',toc(ticsolve))


%% Clear Kernel Pool
cspice_kclear();


%% Functions

function [mu,ustar,Tmax,Isp,g0] = constants % [AU,YEARS,KG]
mu = cspice_bodvrd('Sun','GM',1); % [km^3/s^2]
ustar = 1; % thrust throttle factor for time-optimal problem
Tmax = 150*1e-6; % [kN --> kg*Km/s^2] Max Thrust
Isp = 3000; % [s] % Specific impulse
g0 = 9.80665*1e-3; % [km/s^2] gravitational acceleration at sea level

km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');

mu = mu * km2au^3/sec2year^2; % [AU^3/years^2]
Tmax = Tmax * km2au/sec2year^2; % [kg*AU/years^2]
Isp = Isp * sec2year; % [years]
g0 = g0 * km2au/sec2year^2; % [AU/years^2]

end

%----------------------------------------------------------------------

function dy = TPBVP(t,x,Nth) %  2 point boundary value problem
% works with input units triplet[AU,YEARS,KG]

% variables
r = x(1:3);
v = x(4:6);
m = x(7);
Lr = x(8:10);
Lv = x(11:13);
Lm = x(14);

[mu,ustar,Tmax,Isp,g0] = constants;
% [AU,YEARS,KG]

% dynamics
rdot = v;
vdot = -mu*r/norm(r)^3 - ustar*Nth*Tmax/m * Lv/norm(Lv);
mdot = -ustar*Nth*Tmax/(Isp*g0);
Lrdot = -3*mu/norm(r)^5*dot(r,Lv)*r + mu/norm(r)^3*Lv;
Lvdot = -Lr;
Lmdot = -ustar*norm(Lv)*Nth*Tmax/m^2;

dy = [rdot;vdot;mdot;Lrdot;Lvdot;Lmdot];

end

%----------------------------------------------------------------------

function [y,out,t_out] = orbit_propagator(x0,t0,tf,Nth)
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);
[t_out,out] = ode113(@TPBVP,[t0 tf],x0',opt,Nth);
y = out(end,:)';
end

%----------------------------------------------------------------------

function F = boundariesfun(x,Nth,x_start) 
% input units triplet [AU,Year,Kg] for performance
% {λ0,tf} zero-finding variables 
% x_start [Km,s]


Lm0 = abs(x(7)); % [AU,kg] λm guess > 0
Dtf = abs(x(8)); % [Years] Δt guess > 0


Lx = x(1:6);

y0 = [x_start(1:7);Lx;Lm0];
t0 = x_start(8);
tf = Dtf + t0;
Y = orbit_propagator(y0,t0,tf,Nth);
rf = Y(1:3);
vf = Y(4:6);
mf = Y(7);
Lrf = Y(8:10);
Lvf = Y(11:13);
Lmf = Y(14);

tf = cspice_convrt(tf,'YEARS','SECONDS');

xM = cspice_spkezr('Mars',tf,'ECLIPJ2000','none','Sun'); % [km][Km/s]

km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');

xM(1:3) = xM(1:3)*km2au;
xM(4:6) = xM(4:6)*km2au/sec2year;
tf = tf * sec2year;

rM = xM(1:3);
vM = xM(4:6);

[mu,ustar,Tmax,Isp,g0] = constants; % [AU,YEARS,KG]
aM = - mu*rM/norm(rM)^3;
Stf = switchingfun(mf,Lvf,Lmf);
% Hemiltonian
Hf = 1 + dot(Lrf,vf) - mu/norm(rf)^3*dot(rf,Lvf) + Nth*Tmax/(Isp*g0)*ustar*Stf;
HM = dot(Lrf,vM) + dot(Lvf,aM);

% F = [rf - rM; vf - vM; Lmf; Hf - HM];
% F = [rf - rM; vf - vM; Lmf; Hf];
F = [rf - rM; vf - vM; Lmf];

end

%----------------------------------------------------------------------

function St = switchingfun(m,Lv,Lm) % Switching function
[~,~,~,Isp,g0] = constants;
St = - norm(Lv)*Isp*g0/m - Lm;
end

