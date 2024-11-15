% Spacecraft Guidance and Navigation (2022/2023)
% Assignment #1 - Guidance - Ex 2 Impulsive guidance (shooting methods)
% Author: Matteo Luciardello Lecardi
% Personal code: 10572870
% Matricola MSc Degree: 977701

% Total expected execution time ~ 650 -> 800s

%% Load kernels

% cd 'D:\Courses\5\1 Spacecraft guidance and navigation\2022_23\Assignment01'
% cspice_furnsh 'assignment.tm';
clearvars; close all; clc; format long
clear memoizer        % initialize memoize()
% opengl hardware
ex2tic = tic;

%% CUSTOMIZATION
% Choose MEMOIZE() ON/OFF

% Simple Shooting
memo1 = {1,'ON'};   
% memo1 = {0,'OFF'};

% Multiple Shooting
% memo2 = {1,'ON'};
memo2 = {0,'OFF'};

% Best combination is ON/OFF

% Choose number of mesh points for Multiple Shooting
% NN = 4;
% NN = [4 23];
NN = [4 23 30];

%% EX 2.1 FIRST GUESS
fprintf('\nEx 2.1: Initializing First Guess...')
tic
a = 1.5*pi; % [rad] angle defined on the Earth circular parking orbit
b = 1.41;   % [sqrt(km)/s] initial-to-circular velocity ratio
d = 7;      % [days] transfer duration
ti = 0;     % [days] initial time
tf = ti + d;% [days] arrival time

mi = 1.21506683*10^-2;   % Earth-Moon mass parameter 
Re = 6378;               % [km] Mean Earth/s radius
Rm = 1738;               % [km] Mean Moon's radius
hi = 167;                % [km] Altitude of departure orbit 
hf = 100;                % [km] Altitude of arrival orbit
DU = 3.84405*10^5;       % [km] Distance Unit Earth-Moon
oms = -9.25195985*10^-1; % [s^-1] Scaled angular velocity of the Sun
ms = 3.28900541*10^5;    % Scaled mass of the Sun
r = 3.88811143*10^2;     % Scaled Sun-(Earth+Moon) distance
omem = 2.66186135*10^-6; % [s^-1] Scaled Earth-Moon angular velocity
TU = 4.34811305;         % [days] Time Unit
T = 23*TU;               % [days] Maximum prescribed transfer time 

ri = (Re+hi)/(DU);    % [km] initial circular orbit scaled radius
rf = (Rm+hf)/(DU);    % [km] final circular orbit scaled radius

% INITIAL GUESS SOLUTION

v0 = b*sqrt((1-mi)/ri);

x0 = double(ri*cos(sym(a))-mi);      % r0*cos(a)-mi
y0 = ri*sin(a);
x0dot = -(v0-ri)*sin(a);
y0dot = 0;               % (v0-r0)*cos(a)

X0 = [x0; y0; x0dot; y0dot]; % initial transfer state (alpha, beta)
fprintf('\n x0: %.10f[-]\n y0: %.10f[-]\n u0: %.10f[-/s]\n v0: %.10f[-/s]\n',X0)
% PLOT of transfer state in the rotating frame
P1 = -mi;  % Earth barycenter
P2 = 1-mi; % Moon barycenter

[XF,out,t_ex] = orbit_propagation(X0,ti,tf);

Dvi = sqrt((X0(3)-X0(2))^2+(X0(4)+X0(1)+mi)^2)-sqrt((1-mi)/ri); % cost for initial maneuver
Dvf = sqrt((XF(3)-XF(2))^2+(XF(4)+XF(1)+mi-1)^2)-sqrt((mi)/rf); % cost for arrival maneuver
Dv = abs(Dvi) + abs(Dvf);                                       % total cost of transfer
fprintf('\n %sV = %.10f[-/s]',char(916),Dv)

% PLOT

plot1(mi,P1,P2,X0,XF,ri,rf,out,t_ex)

fprintf('\n Ex 2.1: First Guess execution time %.3fs\n\n',toc)


%% EX 2.2 SIMPLE SHOOTING
% a) Simple shooting w/o derivatives
% b) derivatives with STM estimation (variational equations apporach)

fprintf('\nEx 2.2: Initializing Simple Shooting...\n')
figure(2)
close figure 2
ANS.SimpleShooting = [];
ANS.VariableSimpleShooting = [];

ex2_2 = tic;

x0 = [X0;ti;tf];

options = optimoptions('fmincon', 'Algorithm', ...
                       'active-set', ...
                       ...'interior-point', ...
                       ...'UseParallel',true, ...  % Suggested by Oshima for multiple shooting
                       ...
                       ... PROBLEM TYPE
                       ...'ConstraintTolerance',1e-10, ... % Suggested by Oshima for multiple shooting, 10 digits accuracy
                       ...'SpecifyObjectiveGradient',true, ...
                       ...'SpecifyConstraintGradient',true, ...
                       ...'CheckGradients',true, ...    % won't actually work
                       ...
                       ... FOR PERFORMANCE
                       'FiniteDifferenceType','central', ... % accuracy
                       'FiniteDifferenceStepSize',1e-10, ...
                       ...
                       ... ADDITIONAL PROPERTIES
                       'MaxFunctionEvaluations',1e6, ...
                       'MaxIterations',1e6, ...
                       ...'StepTolerance',eps, ...     % 16 digits accuracy input
                       ...'FunctionTolerance',eps, ... % 16 digits accuracy output
                       'Display', ...
                       ...'iter-detailed');
                       'off');


% Choose MEMOIZE() ON/OFF
memo = memo1;
fprintf('MEMOIZE: %s\n\n',memo{2})

% Different precision executions
N1 = {'double' '16digits'};
N2 = {'best accuracy' 'double' sprintf('16dig.%sCons',char(8711)) '16digits'};
str = {'Manually computed matrix multiplications' ...
       ' --> good floating-point accuracy...' ...
       'Derivatives by auto matrix moltiplications' ...
       ' --> poor floating-point accuracy...' ...
       ' --> improved floating-point accuracy...'};
str = {[sprintf('\b\b')] [sprintf('\b%s',str{5})] [] [];[str{1},newline,str{2}] [str{3},newline,str{4}] [str{3},newline,str{5}] [str{3},newline,str{5}]};
       
% Tolerances
ConstTol1 = [1e-10 1e-10];
% ConstTol2 = [1e-10 2e-4 1e-4 1e-4]; % for max accuracy at "best accuracy"
% case --> as in A.2e table
ConstTol2 = [1e-7 2e-4 1e-4 1e-4]; % for speed & relatively high accuracy 

ConstTol = {ConstTol1 ConstTol2};

% Gradients selector
Gsel = {false true};

% Choose problem of Simple Shooting
a = N1; % w/o derivatives
b = N2; % w/ gradients of objective and constraints functions
N = {a b};

problem = {'2.2.a) W/O Derivatives' '2.2.b) W/ Gradients'};
problem1 = [0 1];
% fprintf('\n Initialization Simple Shooting\n')

A = [];
Aeq = [];
b = [];
beq = [];
lb = [-Inf,-Inf,-Inf,-Inf,0,0];
ub = [Inf,Inf,Inf,Inf,2*pi/abs(oms),2*pi/abs(oms)+T];
nonlcon = @cons;

for k = 1:length(N)
    fprintf(' %s...\n\n',problem{k})
    for i = 1:length(N{k})
        clearAllMemoizedCaches
        clear memoizer
        options.SpecifyObjectiveGradient = Gsel{k};
        options.SpecifyConstraintGradient = Gsel{k};
        options.ConstraintTolerance = ConstTol{k}(i);
        
        fprintf('Precision: %s...\n %s',N{k}{i},str{k,i})

        tic

        [x0ab,f,exflag,output] = fmincon(@shootobjfun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options,mi,ri,rf,memo{1},N{k}{i});
        
        fprintf('\n Execution time: %.3fs\n\n',toc)

        ANS.SimpleShooting(i+2*k-2).Problem = problem1(k);
        ANS.SimpleShooting(i+2*k-2).Precision = N{k}{i};
        ANS.SimpleShooting(i+2*k-2).Iter = output.iterations;
        ANS.SimpleShooting(i+2*k-2).t = num2str(toc,'%.2f');
        ANS.SimpleShooting(i+2*k-2).Dv = num2str(f,'%.10f');
        ANS.SimpleShooting(i+2*k-2).FirstOrderOpt = num2str(output.firstorderopt,'%.10e');
        ANS.SimpleShooting(i+2*k-2).ConstrViol = num2str(output.constrviolation,'%.10e');
        ANS.VariableSimpleShooting(i+2*k-2).x0ab = x0ab;
    end
end

SimpleShooting = struct2table(ANS.SimpleShooting);
for i = 1:width(SimpleShooting)
    SimpleShooting.(i) = categorical(SimpleShooting.(i));
end
SimpleShooting.Properties.VariableNames = {char(8711) 'Precision' 'Iterations' 'Exe. Time[s]' [char(916) 'V[-/s]'] 'First Order Optimality' 'Constraints Violation'}
ANS.SimpleShooting = SimpleShooting;

x0a = ANS.VariableSimpleShooting(2).x0ab;
x0b = ANS.VariableSimpleShooting(3).x0ab;

disp(['Choose precision from "Constraints Violation" lowest value:' newline ' a) 16digits' newline ' b) best accuracy'])

% PLOT 


[XFa,out_a,t_ex_a] = orbit_propagation(x0a(1:4),x0a(5),x0a(6));
[XFb,out_b,t_ex_b] = orbit_propagation(x0b(1:4),x0b(5),x0b(6));

 
plot2(mi,P1,P2,x0a(1:4),x0b(1:4),XFa,XFb,ri,rf,out_a,out_b,t_ex_a,t_ex_b)
%% 


fprintf('\n Ex 2.2: Simple Shooting execution time: %.1fs\n',toc(ex2_2))

%% EX 2.3 Multiple Shooting

ANS.MultipleShooting = [];
ANS.VariableMultipleShooting = [];

fprintf('\n\nEx 2.3: Initializing Multiple Shooting...')

ex2_3tic = tic;

% Choose number of mesh points
% NN = 30;
% NN = [4 23 30];
% NN = [4 10 23];

% Choose MEMOIZE() ON/OFF
memo = memo2;

fprintf('\nMEMOIZE: %s\n\n',memo{2})
for k = 1 : length(NN)
clearAllMemoizedCaches
clear memoizer
N = NN(k);
tic
% Setup initial conditions of NLP variables vector
tj = zeros(1,N);
Xj = zeros(4,N);
xj0 = zeros(1,4*N+2);
tj(1) = ti;
Xj(:,1) = X0;
xj0(1:4) = X0;
for j = 2 : N
    tj(j) = ti + (j - 1)*(tf - ti)/(N - 1);
    Xj(:,j) = orbit_propagation(Xj(:,j-1),tj(j-1),tj(j));
    xj0(4*j-3:4*j) = Xj(:,j);
end
xj0(end-1:end) = [ti,tf];


% FMINCON

% Setup Parallel pool
% pp = gcp;
% if pp.Connected ~= 1
%    parpool;
% end

options = optimoptions('fmincon', 'Algorithm', ...
                       'active-set', ...
                       ...'UseParallel',true, ...  % Suggested by Oshima for multiple shooting
                       ...
                       ... PROBLEM TYPE
                       'ConstraintTolerance',1e-10, ... % Suggested by Oshima for multiple shooting, 10 digits accuracy
                       'SpecifyObjectiveGradient',true, ...
                       'SpecifyConstraintGradient',true, ...
                       ...
                       ... FOR PERFORMANCE
                       'FiniteDifferenceType','central', ... % accuracy, suggested by Topputo 2013
                       'FiniteDifferenceStepSize',1e-10, ...
                       ...
                       ... ADDITIONAL PROPERTIES
                       'MaxFunctionEvaluations',1e6, ...
                       'MaxIterations',1e6, ...
                       ...'StepTolerance',1e-10, ...     % 10 digits accuracy input
                       ...'FunctionTolerance',1e-10, ... % 10 digits accuracy output
                       'Display', ...
                       ...'iter-detailed');
                       'off');

A = [];
Aeq = [];
b = [];
beq = [];
lb = [ones(1,4*N)*-Inf,0,0];
ub = [ones(1,4*N)*Inf,2*pi/abs(oms),2*pi/abs(oms)+T];
nonlcon = @cons;

fprintf('Computing solution with N = %d points...\n',N)
[xjF,f,exflag,output] = fmincon(@shootobjfun,xj0,A,b,Aeq,beq,lb,ub,nonlcon,options,mi,ri,rf,memo{1});
% output;

 fprintf(' %sV = %.10f[-/s]\n Execution time: %.3fs\n\n',char(916),f,toc)
 
% Table
    ANS.MultipleShooting(k).N = N;
    ANS.MultipleShooting(k).Grad = 1;
    ANS.MultipleShooting(k).Iter = output.iterations;
    ANS.MultipleShooting(k).t = num2str(toc,'%.2f');
    ANS.MultipleShooting(k).Dv = num2str(f,'%.10f');
    ANS.MultipleShooting(k).FirstOrderOpt = num2str(output.firstorderopt,'%.10e');
    ANS.MultipleShooting(k).ConstrViol = num2str(output.constrviolation,'%.10e');
    ANS.VariableMultipleShooting(k).xj = xjF;

end

fprintf('\n Ex 2.3: Multiple Shooting execution time: %.1fs\n',toc(ex2_3tic))

MultipleShooting = struct2table(ANS.MultipleShooting);
 
for i = 1 : width(MultipleShooting)
    if height(MultipleShooting) == 1
        break
    end
    MultipleShooting.(i) = categorical(MultipleShooting.(i));
end
MultipleShooting.Properties.VariableNames = {'N' char(8711) 'Iterations' 'Exe. Time[s]' [char(916) 'V[-/s]'] 'First Order Optimality' 'Constraints Violation'}
ANS.MultipleShooting = MultipleShooting;

%% 

%PLOT
figure(3), figure(4), figure(5)
close figure 3, close figure 4, close figure 5
for k = 1 : length(NN)
    xj = ANS.VariableMultipleShooting(k).xj; 
    plot3(xj,P1,P2,ri,rf,mi,k)
end
 
fprintf('\nConfiguration MEMOIZE: %s/%s\n',memo1{2},memo2{2})
fprintf('\nEx 2 total execution time: %.1fs\n',toc(ex2tic))

%% 
% Close parallel pool
%  delete(gcp('nocreate'))

%%
% If code slows down, it may be caused by cache saturation.
% It's good practice to clear it after each use
clearAllMemoizedCaches
clear memoizer

%% Functions

function [xunit,yunit] = circles(x,y,r)
% hold on
th = 0:pi/50:2*pi;
xunit = zeros(1,length(th));
yunit = zeros(1,length(th));
   
    for i = 1:length(x)
        xunit = r(i) * cos(th) + x(i);
        yunit = r(i) * sin(th) + y(i);
        if nargout == 0
            hold on
            plot(xunit, yunit,'Color',[0.6 0.6 0.6],'LineWidth',0.3);
        end
    end
end

%--------------------------------------------------------------------------

function dy = pcr4bp(t,f)

    mi = 1.21506683*10^-2;   
    oms = -9.25195985*10^-1;
    ms = 3.28900541*10^5;    
    r = 3.88811143*10^2;
    
%   ...  om3(x,y) = 1/2*(x^2+y^2)+(1-mi)/r1(x,y)+mi/r2(x,y)+1/2*mi*(1-mi);
%   ...  om4(x,y,t_) = (1/2*(x^2+y^2)+(1-mi)/(((x+mi)^2+y^2)^(1/2))+mi/(((x+mi-1)^2+y^2)^(1/2))+1/2*mi*(1-mi))+ms/(((x-r*cos(oms*t_))^2+(y-r*sin(oms*t_))^2)^(1/2))-(ms/r^2)*(x*cos(oms*t_)+y*sin(oms*t_));

    x = f(1);
    y = f(2);
    xdot = f(3);
    ydot = f(4);
     
    % 4-body effective potential partial derivatives wrt dx and dy
    OM4x = x - (ms*(x - r*cos(oms*t)))/((x - r*cos(oms*t))^2 + ...
           (y - r*sin(oms*t))^2)^(3/2) - (mi*(mi + x - 1))/((mi + x - 1)^2 + ...
           y^2)^(3/2) - (ms*cos(oms*t))/r^2 + ((mi + x)*(mi - 1))/((mi + x)^2 + y^2)^(3/2);
    OM4y = y - (ms*(y - r*sin(oms*t)))/((x - r*cos(oms*t))^2 + ...
           (y - r*sin(oms*t))^2)^(3/2) - (mi*y)/((mi + x - 1)^2 + y^2)^(3/2) - ...
           (ms*sin(oms*t))/r^2 + (y*(mi - 1))/((mi + x)^2 + y^2)^(3/2);
     
    dy = [       xdot;
                 ydot;
          2*ydot+OM4x;
         -2*xdot+OM4y];

    if length(f) > 4  % Build STMdot
    
        OM4xx = 1 - ms/((x - r*cos(oms*t))^2+(y - r*sin(oms*t))^2)^(3/2) + ...
                (3*ms*(x - r*cos(oms*t))^2)/((x - r*cos(oms*t))^2 + ...
                (y - r*sin(oms*t))^2)^(5/2) - (1 - mi)/((x + mi)^2 + y^2)^(3/2) + ...
                (3*(1 - mi)*(x + mi)^2)/((x + mi)^2 + y^2)^(5/2) - ...
                mi/((x + mi - 1)^2 + y^2)^(3/2) + (3*mi*(x + mi - 1)^2)/((x + mi - 1)^2 + y^2)^(5/2);
        OM4xy = (3*ms*(x - r*cos(oms*t))*(y - r*sin(oms*t)))/((y - r*sin(oms*t))^2 + ...
                (x - r*cos(oms*t))^2)^(5/2) + (3*(1 - mi)*(x + mi)*y)/(y^2+(x + mi)^2)^(5/2) + ...
                (3*mi*(x + mi - 1)*y)/(y^2 + (x + mi - 1)^2)^(5/2);
        OM4yx = OM4xy;
        OM4yy = 1 - ms/((y - r*sin(oms*t))^2 + (x-r*cos(oms*t))^2)^(3/2) + ...
                (3*ms*(y - r*sin(oms*t))^2)/((y - r*sin(oms*t))^2 + ...
                (x - r*cos(oms*t))^2)^(5/2) - (1 - mi)/(y^2 + (x + mi)^2)^(3/2)+ ...
                (3*(1 - mi)*y^2)/(y^2 + (x + mi)^2)^(5/2) - mi/(y^2 + ...
                (x + mi - 1)^2)^(3/2) + (3*mi*y^2)/(y^2 + (x + mi - 1)^2)^(5/2);
    
        A = [0     0      1 0;
             0     0      0 1;
             OM4xx OM4xy  0 2;
             OM4yx OM4yy -2 0];
    
        dPHI = A * reshape(f(5:20),[4 4]);  % variational equations
        dy = [dy;reshape(dPHI,[16 1])];   

    end
end

%--------------------------------------------------------------------------

function [phi,x,tt,xe,te] = orbit_propagation(x0,t0,t)
% outputs:  
%         phi = [N x 1] final column state x at time t (N = dim of state)
%           x = [I x N] row state propagation along orbit, (I = n. of states)
%          tt = [I x 1] time corresponding to each row of x
    
options = odeset('reltol', 1e-12, 'abstol', 1e-12); % suggested in Oshima, Topputo, Yanao CMDA 2019
[tt,x] = ode113(@pcr4bp,[t0 t],x0,options);
phi = x(end,:)';

end

%--------------------------------------------------------------------------

function [xf,PHI] = STM(xi,ti,tf)
    
% VARIATIONAL APPROACH

    xx0 = [xi;reshape(eye(4),[16 1])];
    
    phi = orbit_propagation(xx0,ti,tf);
    xf = phi(1:4);
    PHI = reshape(phi(5:20),[4 4]);

end

% function [phi_0,PHI] = STMnum(xi,ti,tf)
%     
% % NUMERICAL APPROACH
% % Initialize perturbations and STM
%     dx = zeros(length(xi));
%     PHI_t = zeros(length(xi));
% % Reference solution
%     phi_0 = orbit_propagation(xi,ti,tf);
% % Computation of STM by columns
%     for i = 1 : length(xi)
%         dx(i,i) = sqrt(eps)*max(1, abs(xi(i)));      % Perturbation
%         phi_x = orbit_propagation(xi+dx(:,i),ti,tf); % Perturbed solution
%         PHI_t(:,i) = (phi_x - phi_0) / dx(i,i);      % Forward differences
%     end
% end

%--------------------------------------------------------------------------

function [x2,PHI,f1,f2] = memoSTMfun(x)

% INPUT: 
%       x [6] = (X,t1,t2) variables of simple shooting problem 
% OUTPUT: 
%       {vars} [4] = variables to be shared within gradient functions of 
%                    Objective and Constraints functions
% Performance: MEMOIZE() can be used so to cache combinations of input / output 
%              of this function, therefore the propagation of STM(t1,t2) 
%              for both gradient functions can be avoided

    x1 = x(1:4);
    t1 = x(5);
    t2 = x(6);

    if nargout == 1
        x2 = orbit_propagation(x1,t1,t2);
    else
        [x2,PHI] = STM(x1,t1,t2);    
    
        xdot1 = pcr4bp(t1,x1);
        f1 = xdot1;
        xdot2 = pcr4bp(t2,x2);
        f2 = xdot2;
    end
end

%--------------------------------------------------------------------------

function  [MEMOSTM,STATS,x2,PHI,f1,f2] = memoSTM(x)
    
% MEMOIZE() improves performance
% Performance: MEMOIZE() can be used so to cache combinations of input / output 
%              of this function, therefore the propagation of STM(t1,t2) 
%              for both gradient functions can be avoided

    MEMOSTM = memoize(@memoSTMfun); 
    MEMOSTM.CacheSize = 1e5;
    if nargout <= 3
        x2 = MEMOSTM(x);
    else
        [x2,PHI,f1,f2] = MEMOSTM(x);
    end
    STATS = MEMOSTM.stats;
end

%--------------------------------------------------------------------------

function [dv,gradf] = shootobjfun(x,mi,ri,rf,memo,n)
    
    x1 = x(1:4);
    t1 = x(end-1);
    t2 = x(end);
    
    N = (length(x)-2)/4;

    if N == 1           % Simple Shooting case
        if memo == 1
            [~,s,x2] = memoSTM([x1;t1;t2]); % call to MEMOIZE() for efficiency
            assignin("base","memoize_stats",s); % Retreive memoize info within base workspace
        elseif memo == 0
            x2 = memoSTMfun([x1;t1;t2]); % test w/o MEMOIZE()
        end
    elseif N > 1        % Multiple Shooting case
        x2 = x(end-5 : end-2);  % --> is not time dependent
    end

    if nargin <= 5
        n = [];
    end
    if strcmp('16digits',n)
        x1 = sym(x1);
        x2 = sym(x2);
    end

    Dvi = sqrt((x1(3)-x1(2))^2+(x1(4)+x1(1)+mi)^2)-sqrt((1-mi)/ri); % [scalar] cost for initial maneuver
    Dvf = sqrt((x2(3)-x2(2))^2+(x2(4)+x2(1)+mi-1)^2)-sqrt(mi/rf);   % [scalar] cost for arrival maneuver
    dv = double(abs(Dvi)+abs(Dvf));

% Gradient of the objective function:
    
    if nargout > 1

        gradf = grad_objfun(x,mi,n,N,memo);
    
    end

end

%--------------------------------------------------------------------------

function gof = grad_objfun(x,mi,n,N,memo)
    
    x1 = x(1:4);
    t1 = x(end-1);
    t2 = x(end);
    
    
    if N == 1
        if memo == 1
            [~,s,x2,PHI,f1,f2] = memoSTM([x1;t1;t2]); % call to MEMOIZE() for efficiency
            assignin("base","memoize_stats",s)
        elseif memo == 0
            [x2,PHI,f1,f2] = memoSTMfun([x1;t1;t2]);
        end
    elseif N > 1        % Multiple Shooting case
        x2 = x(end-5 : end-2);  % --> final variable x_N, which is not time dependent
    end
    
    % LEC 08 + OSHIMA (2019)

    if strcmp('16digits',n)
        x1 = sym(x1);
        x2 = sym(x2);
        mi = sym(mi);
        PHI = sym(PHI);  % IMPLEMENTED ONLY FOR SIMPLE SHOOTING --> don't use n for multiple shooting
        f1 = sym(f1);
        f2 = sym(f2);
    end
    
    P1 = [x1(4)+x1(1)+mi,  x1(2)-x1(3),  x1(3)-x1(2),  x1(4)+x1(1)+mi]/...       % dJ/dx1
          sqrt((x1(3)-x1(2))^2+(x1(4)+x1(1)+mi)^2);
    PN = [x2(4)+x2(1)+mi-1,  x2(2)-x2(3),  x2(3)-x2(2),  x2(4)+x2(1)+mi-1]/...   % dJ/dx2
          sqrt((x2(3)-x2(2))^2+(x2(4)+x2(1)+mi-1)^2);
    

    if N == 1
        P2 = PN;
        dJdt1 = P2*(-PHI*f1);   % dJ/dx2 * dx2/dt1
        dJdt2 = P2*f2;          % dJ/dx2 * dx2/dt2
        gof = double([P1 dJdt1 dJdt2]');
    elseif N > 1
        gof = [P1 zeros(1,4*(N-2)) PN zeros(1,2)]';
    end
end

%--------------------------------------------------------------------------

function [c,ceq,DC,DCeq] = cons(x,mi,ri,rf,memo,n)
    x1 = x(1:4);
    t1 = x(end-1);
    t2 = x(end);
    N = (length(x)-2)/4;

%   Initialization
    if N == 1 % Simple Shooting case
        X1 = [];
        tj = [];
        if memo == 1
        [~,s,x2] = memoSTM([x1;t1;t2]); % call to MEMOIZE() for efficiency
        elseif memo == 0
            x2 = memoSTMfun([x1;t1;t2]);
        end
    elseif N > 1 % Multiple Shooting case
        xj = x(1:end-2);
        tj = zeros(1,N);
        X1 = zeros(4,N);
        X2 = zeros(4,N); % --> phi(:,j)(X1(:,j),t(j),t(j+1))
        x2 = x(end-5 : end-2); % --> final variable x_N
        tj(1) = t1;
        X1(:,1) = xj(1:4); 
        
        for j = 2 : N
            tj(j) = t1 + (j - 1)*(t2 - t1)/(N - 1);
            X1(:,j) = xj(4*j-3:4*j);
            if memo == 1
                [~,s,X2(:,j)] = memoSTM([X1(:,j-1);tj(j-1);tj(j)]); %  First column will remain empty
            elseif memo == 0
                X2(:,j) = memoSTMfun([X1(:,j-1);tj(j-1);tj(j)]);
            end
        end
    end
    if memo == 1
        assignin("base","memoize_stats",s); % Retreive memoize info within base workspace
    end
    if nargin <= 5
       n = [];
    end
    

% Inequality constraints
    Tau = t1 - t2;
    if  N == 1
        c = Tau;
    elseif N > 1
        Etaj = zeros(1,2*N);
        Re = 6378/3.84405e5;
        Rm = 1738/3.84405e5;
        for j = 1 : N
           ETAJ =  [Re^2 - (X1(1,j) + mi)^2 - X1(2,j)^2, Rm^2 - (X1(1,j) + mi - 1)^2 - X1(2,j)^2];
           Etaj(1,2*j-1:2*j) = ETAJ;
        end
        c = [Etaj,Tau];
    end

% Equality constraints
    
    Psi1 = [(x1(1)+mi)^2+x1(2)^2-ri^2; (x1(1)+mi)*(x1(3)-x1(2))+x1(2)*(x1(4)+x1(1)+mi)];    
    PsiN = [(x2(1)+mi-1)^2+x2(2)^2-rf^2; (x2(1)+mi-1)*(x2(3)-x2(2))+x2(2)*(x2(4)+x2(1)+mi-1)];
    if N == 1
        Psi2 = PsiN;
        ceq = [Psi1; Psi2];
    elseif N > 1
        zj = zeros(1,4*(N-1));
        for j = 1 : N - 1
            ZJ = X2(:,j+1) - X1(:,j+1);   % psi(x(j),t(j),t(j+1)) - x(j+1)
            zj(1,4*j-3:4*j) = ZJ'; 
        end
        ceq = [zj,Psi1',PsiN']; 

    end

% Gradient of the constraints:
    
    if nargout > 2
    
        [DCeq, DC] = grad_con(x,mi,n,N,X1,tj,memo);
            
    end
end

%--------------------------------------------------------------------------

function [gc_eq, gc_ineq] = grad_con(x,mi,n,N,X1,tj,memo)
  
    x1 = x(1:4);
    t1 = x(end-1);
    t2 = x(end);


   if N == 1
       if memo == 1
           [~,s,x2,PHI,f1,f2] = memoSTM([x1;t1;t2]); % call to MEMOIZE() for efficiency
       elseif memo == 0
           [x2,PHI,f1,f2] = memoSTMfun([x1;t1;t2]); % test w/o MEMOIZE()
       end
   elseif N > 1
       x2 = x(end-5 : end-2); % --> final variable x_N
       gc_eq = zeros(4*(N-1)+4,4*N+2);
       for j = 1 : N - 1
           if memo == 1
               [~,s,~,PHI,f1,f2] = memoSTM([X1(:,j);tj(j);tj(j+1)]); % call to MEMOIZE() for efficiency (expected)
           elseif memo == 0
               [~,PHI,f1,f2] = memoSTMfun([X1(:,j);tj(j);tj(j+1)]);  % test w/o MEMOIZE()
           end
           Q1j = -(N-j)/(N-1)*PHI*f1 + (N-j-1)/(N-1)*f2; % dZeta(j)/dt1  [4x1]
           QNj = -(j-1)/(N-1)*PHI*f1 + j/(N-1)*f2;       % dZeta(j)/dtN  [4x1]
           
           gc_eq(4*j-3:4*j,4*j-3:4*j) = PHI;             % PHI(t(j),t(j+1)) [4x4]
           gc_eq(4*j-3:4*j,4*j+1:4*j+4) = -eye(4);       % [4x4]
           gc_eq(4*j-3:4*j,4*N+1) = Q1j;
           gc_eq(4*j-3:4*j,4*N+2) = QNj;

       end
   end
       
   if memo == 1
       assignin("base","memoize_stats",s); % Retreive memoize info within base workspace
   end

    R1 = [2*(x1(1)+mi)  2*x1(2)          0     0;    % dPsi1/dx1
                 x1(3)    x1(4)   x1(1)+mi x1(2)];

    RN = [2*(x2(1)+mi-1)  2*x2(2)           0     0; % dPsi2/dx2
                   x2(3)    x2(4)  x2(1)+mi-1 x2(2)];  
 
if N > 1
    
    gc_eq(end-3:end-2,1:4) = R1;
    gc_eq(end-1:end,end-5:end-2) = RN;
    gc_eq = gc_eq';

elseif N == 1
    R2 = RN;
% TESTED BEST ACCURACY
% MIX of compact and element by element construction, 
% avoiding matrix multiplications such as [2x2]*[2x1] & [2x4]*[4x4]
% which seem to be problematic for fmincon

    if strcmp('best accuracy',n)
    
%   dPsi2/dx1 = dPsi2/dx2 * dx2/dx1 = (R2*PHI)'
    
        dPsi2dx1 = [2*(x2(1) + mi - 1)*PHI(1,1) + 2*x2(2)*PHI(2,1),  PHI(1,1)*(x2(3) - x2(2)) + (x2(1) + mi - 1)*(PHI(3,1) - PHI(2,1)) + PHI(2,1)*(x2(4) + x2(1) + mi - 1) + x2(2)*(PHI(4,1) + PHI(1,1));
                    2*(x2(1) + mi - 1)*PHI(1,2) + 2*x2(2)*PHI(2,2),  PHI(1,2)*(x2(3) - x2(2)) + (x2(1) + mi - 1)*(PHI(3,2) - PHI(2,2)) + PHI(2,2)*(x2(4) + x2(1) + mi - 1) + x2(2)*(PHI(4,2) + PHI(1,2));
                    2*(x2(1) + mi - 1)*PHI(1,3) + 2*x2(2)*PHI(2,3),  PHI(1,3)*(x2(3) - x2(2)) + (x2(1) + mi - 1)*(PHI(3,3) - PHI(2,3)) + PHI(2,3)*(x2(4) + x2(1) + mi - 1) + x2(2)*(PHI(4,3) + PHI(1,3));
                    2*(x2(1) + mi - 1)*PHI(1,4) + 2*x2(2)*PHI(2,4),  PHI(1,4)*(x2(3) - x2(2)) + (x2(1) + mi - 1)*(PHI(3,4) - PHI(2,4)) + PHI(2,4)*(x2(4) + x2(1) + mi - 1) + x2(2)*(PHI(4,4) + PHI(1,4))];
   
        dphdt1 = -PHI*f1;
        dPsi1dx1 = R1';
%   dPsi2/dt1 = dPsi2/dx2 * dx2/dt1 = (R2 * (-PHI*f1))'
        dPsi2dt1 = [2*(x2(1) + mi - 1)*dphdt1(1) + 2*x2(2)*dphdt1(2),   dphdt1(1)*x2(3) + dphdt1(2)*x2(4) + dphdt1(3)*(x2(1) + mi - 1) + dphdt1(4)*x2(2)];

        dPsi1dt1 = [2*(x1(1) + mi)*x1(3) + 2*x1(2)*x1(4),  x1(3)^2 + x1(4)^2 + f1(3)*(x1(1) + mi) + f1(4)*x1(2)];
        dPsi1dt2 = [0, 0];
        dPsi2dt2 = [2*(x2(1) + mi-1)*x2(3) + 2*x2(2)*x2(4),  x2(3)^2 + x2(4)^2 + f2(3)*(x2(1) + mi - 1) + f2(4)*x2(2)];

        gc_eq = [dPsi1dx1 dPsi2dx1; dPsi1dt1 dPsi2dt1; dPsi1dt2 dPsi2dt2];

    
  

    elseif strcmp('double',n) || strcmp('16digits',n) || strcmp(sprintf('16dig.%sCons',char(8711)),n)
            if strcmp('16digits',n) || strcmp(sprintf('16dig.%sCons',char(8711)),n)
                x1 = sym(x1);
                x2 = sym(x2);
                mi = sym(mi);
                PHI = sym(PHI);
                f1 = sym(f1);
                f2 = sym(f2);
            end

        dPsi1dx1_1 = R1';
        dPsi2dx1_1 = (R2*PHI)'; %  dPsi2/dx1 = dPsi2/dx2 * dx2/dx1
        dPsi1dt1_1 = [2*(x1(1) + mi)*x1(3) + 2*x1(2)*x1(4),  x1(3)^2 + x1(4)^2 + f1(3)*(x1(1) + mi) + f1(4)*x1(2)];
        dPsi2dt1_1 = (R2 * (-PHI*f1))';  %  dPsi2/dt1 = dPsi2/dx2 * dx2/dt1
        dPsi1dt2_1 = [0, 0];
        dPsi2dt2_1 = [2*(x2(1) + mi-1)*x2(3) + 2*x2(2)*x2(4),  x2(3)^2 + x2(4)^2 + f2(3)*(x2(1) + mi - 1) + f2(4)*x2(2)];

        gc_eq = double([dPsi1dx1_1 dPsi2dx1_1; dPsi1dt1_1 dPsi2dt1_1; dPsi1dt2_1 dPsi2dt2_1]);
    
    end
end
    % Inequality constraints gradient
    St = [1 -1];   % dTau/dt1  dTau/dt2
    if N == 1
        gc_ineq = [zeros(1,4) St]';
    elseif N > 1
        gc_ineq = zeros(2*N+1,4*N+2);
        for j = 1 : N
            Sj = [-2*(X1(1,j) + mi) -2*X1(2,j) 0 0; -2*(X1(1,j) + mi - 1) -2*X1(2,j) 0 0];   % dEta(j)/dx(j)
            gc_ineq(2*j-1:2*j,4*j-3:4*j) = Sj; 
        end
        gc_ineq(end,end-1:end) = St;
        gc_ineq = gc_ineq';
    end
   
end

%--------------------------------------------------------------------------

function plot1(mi,P1,P2,X0,XF,ri,rf,out,t_ex)

figure(1)
close figure 1

figure('Name','1: Ex2.1 First guess solution pcr4bp','NumberTitle','off')
tiledlayout(1,2,"TileSpacing","compact","Padding","tight")
nexttile(1)
% subplot(1,2,1)
title('Earth-Moon rotating frame','FontSize',15,'FontWeight','bold')
hold on
grid on
box on
axis equal padded
plot(P1,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
plot(P2,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])
plot(X0(1),X0(2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(XF(1),XF(2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot(out(:,1),out(:,2),'LineWidth',1.6,'Color','k')

circles([P1 P2],[0 0],[ri rf]);

text(P1-0.12,+0.24,'P1','FontSize',14)
text(P2-0.1,+0.22,'P2','FontSize',14)
text(X0(1)-0.06,X0(2)-0.25,'X_i','FontSize',14)
text(XF(1)-0.06,XF(2)-0.2,'X_f','FontSize',14)



legend('Earth','Moon','Initial transfer position','Final transfer position','S/C first guess transfer','Location','southeast','FontSize',15)

xlabel('X [-]','FontSize',15,'FontWeight','bold')
ylabel('Y [-]       ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
xticks([-2.5 -2 -1.5 -1 -0.5 P1 0.5 P2 1.5 2])


% PLOT of transfer state in Earth-centered (P1) inertial frame

x1 = (out(:,1)+mi).*cos(t_ex)-out(:,2).*sin(t_ex);
y1 = (out(:,1)+mi).*sin(t_ex)+out(:,2).*cos(t_ex);
x1dot = (out(:,3)-out(:,2)).*cos(t_ex)-(out(:,4)+out(:,1)+mi).*sin(t_ex);
y1dot = (out(:,3)-out(:,2)).*sin(t_ex)+(out(:,4)+out(:,1)+mi).*cos(t_ex);
X1 = [x1,y1,x1dot,y1dot];


nexttile(2)
% subplot(1,2,2)
title('ECIrf','FontSize',15,'FontWeight','bold')
hold on
grid on
box on
axis equal padded

plot((P1+mi)*cos(t_ex),(P1+mi)*sin(t_ex),'Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1],'LineStyle','none')
plot((P2+mi)*cos(t_ex(1)),(P2+mi)*sin(t_ex(1)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
plot((P2+mi)*cos(t_ex(end)),(P2+mi)*sin(t_ex(end)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])
plot(X1(1,1),X1(1,2),'Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(X1(end,1),X1(end,2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot((P2+mi)*cos(t_ex),(P2+mi)*sin(t_ex),':','LineWidth',2,'Color',[0.8 0.8 0.8])
plot(X1(:,1),X1(:,2),'LineWidth',1.6,'Color','k')


circles([(P1+mi)*cos(t_ex(1)),(P2+mi)*cos(t_ex(end))],[(P1+mi)*sin(t_ex(1)),(P2+mi)*sin(t_ex(end))],[ri, rf]);

text((P1+mi)*cos(t_ex(1))-0.06,(P1+mi)*sin(t_ex(1))+0.15,'P1','FontSize',14)
text((P2+mi)*cos(t_ex(1))-0.06,(P2+mi)*sin(t_ex(1))+0.1,'P2_i','FontSize',14)
text((P2+mi)*cos(t_ex(end))-0.06,(P2+mi)*sin(t_ex(end))+0.1,'P2_f','FontSize',14)
text(X1(1,1)-0.03,X1(1,1)-0.15,'X_i','FontSize',14)
%text(X1(1,1)*cos(t_ex)-0.125,(P2+mi)*t_ex./t_ex+0.07,'Moon','FontSize',14)
text(X1(end,1)-0.04,X1(end,2)-0.15,'X_f','FontSize',14)

xlabel('X [-]','FontSize',15,'FontWeight','bold')
ylabel('Y [-]       ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')

legend('Earth','Initial Moon position','Final Moon position','','','Moon orbit','S/C first guess transfer','Location','southeast','FontSize',15)
end

%--------------------------------------------------------------------------

function plot2(mi,P1,P2,X0a,X0b,Xfa,Xfb,ri,rf,out_a,out_b,t_ex_a,t_ex_b)

figure(2)
close figure 2

figure ('Name','2: Ex2.2 Simple shooting: a) w/o derivatives b) w/ objective fun & constraints fun gradients','NumberTitle','off')
tiledlayout(2,3,"TileSpacing","compact","Padding","tight")
nexttile(1,[2,1])
% subplot(2,3,[1,4])
title('Earth-Moon rotating frame','FontSize',15,'FontWeight','bold')
hold on
grid on
box on
axis equal padded

plot(P1,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
plot(P2,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])
plot(X0a(1),X0a(2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(X0b(1),X0b(2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(Xfa(1),Xfa(2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot(Xfb(1),Xfb(2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot(out_a(:,1),out_a(:,2),'LineWidth',1.6,'Color','b')
plot(out_b(:,1),out_b(:,2),'LineWidth',1.6,'Color',"#D95319")

circles([P1 P2],[0 0],[ri rf]);

text(P1-0.1,+0.22,'P1','FontSize',14)
text(P2-0.1,+0.22,'P2','FontSize',14)
text(X0b(1)-0.08,X0b(2)-0.3,'X_i','FontSize',14)
text(Xfb(1)-0.08,Xfb(2)-0.3,'X_f','FontSize',14)

legend({'Earth','Moon','Initial transfer position','','Final transfer position','',['a) S/C Simple Shooting w/o ',char(8711)],['b) S/C Simple Shooting w/ ',char(8711)]},'Location','southeast','FontSize',15)

xlabel('X [-]','FontSize',15,'FontWeight','bold')
ylabel('Y [-]       ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')

nexttile(2,[2,1])
% subplot(2,3,[2,5])
title('ECIrf','FontSize',15,'FontWeight','bold')
hold on
grid on
box on
axis equal padded

x1a = (out_a(:,1)+mi).*cos(t_ex_a)-out_a(:,2).*sin(t_ex_a);
y1a = (out_a(:,1)+mi).*sin(t_ex_a)+out_a(:,2).*cos(t_ex_a);
x1dota = (out_a(:,3)-out_a(:,2)).*cos(t_ex_a)-(out_a(:,4)+out_a(:,1)+mi).*sin(t_ex_a);
y1dota = (out_a(:,3)-out_a(:,2)).*sin(t_ex_a)+(out_a(:,4)+out_a(:,1)+mi).*cos(t_ex_a);
X1a = [x1a,y1a,x1dota,y1dota];
     
x1b = (out_b(:,1)+mi).*cos(t_ex_b)-out_b(:,2).*sin(t_ex_b);
y1b = (out_b(:,1)+mi).*sin(t_ex_b)+out_b(:,2).*cos(t_ex_b);
x1dotb = (out_b(:,3)-out_b(:,2)).*cos(t_ex_b)-(out_b(:,4)+out_b(:,1)+mi).*sin(t_ex_b);
y1dotb = (out_b(:,3)-out_b(:,2)).*sin(t_ex_b)+(out_b(:,4)+out_b(:,1)+mi).*cos(t_ex_b);
X1b = [x1b,y1b,x1dotb,y1dotb];
    

plot((P1+mi)*cos(t_ex_b),(P1+mi)*sin(t_ex_b),'Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1],'LineStyle','none')
plot((P2+mi)*cos(t_ex_b(1)),(P2+mi)*sin(t_ex_b(1)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])

plot((P2+mi)*cos(t_ex_a(end)),(P2+mi)*sin(t_ex_a(end)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','#ffdd8f')
plot((P2+mi)*cos(t_ex_b(end)),(P2+mi)*sin(t_ex_b(end)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])

plot(X1a(1,1),X1a(1,2),'Marker','o','MarkerSize',7,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','r')
plot(X1b(1,1),X1b(1,2),'Marker','o','MarkerSize',7,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','r')

plot(X1a(end,1),X1a(end,2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot(X1b(end,1),X1b(end,2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot((P2+mi)*cos(t_ex_b),(P2+mi)*sin(t_ex_b),':','LineWidth',2,'Color',[0.8 0.8 0.8])
plot(X1a(:,1),X1a(:,2),'LineWidth',2,'Color','b')
plot(X1b(:,1),X1b(:,2),'LineWidth',2,'Color',"#D95319")

circles([(P1+mi)*cos(t_ex_b(1)),(P2+mi)*cos(t_ex_a(end)),(P2+mi)*cos(t_ex_b(end))],[(P1+mi)*sin(t_ex_b(1)),(P2+mi)*sin(t_ex_a(end)),(P2+mi)*sin(t_ex_b(end))],[ri,rf,rf]);

text((P1+mi)*cos(t_ex_b(1))-0.08,(P1+mi)*sin(t_ex_b(1))+0.15,'P1','FontSize',14)
text((P2+mi)*cos(t_ex_b(1))-0.06,(P2+mi)*sin(t_ex_b(1))+0.1,'P2_i','FontSize',14)

text((P2+mi)*cos(t_ex_a(end))-0.25,(P2+mi)*sin(t_ex_a(end))+0.14,'P2_f','FontSize',14)
text((P2+mi)*cos(t_ex_b(end))-0.04,(P2+mi)*sin(t_ex_b(end))+0.14,'P2_f','FontSize',14)
text(X1b(1,1)-0.035,X1b(1,1)-0.15,'X_i','FontSize',14)
%text(X1(1,1)*cos(t_ex)-0.125,(P2+mi)*t_ex./t_ex+0.07,'Moon','FontSize',14)

text(X1a(end,1)-0.14,X1a(end,2)-0.2,'X_f','FontSize',14)
text(X1b(end,1)-0.08,X1b(end,2)-0.2,'X_f','FontSize',14)

xlabel('X [-]','FontSize',15,'FontWeight','bold')
ylabel('Y [-]       ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')

legend({'Earth','Initial Moon position','a) Final Moon position','b) Final Moon position','','','','','Moon orbit',['a) S/C Simple Shooting w/o ',char(8711)],['b) S/C Simple Shooting w/ ',char(8711)]},'Location','southeast','FontSize',15)

nexttile(3)
% subplot(2,3,3)

title('ECIrf zoom in on arrival','FontSize',15,'FontWeight','bold')
hold on
grid on
box on
axis equal    

plot((P2+mi)*cos(t_ex_b),(P2+mi)*sin(t_ex_b),':','LineWidth',2,'Color',[0.8 0.8 0.8])

circles([(P2+mi)*cos(t_ex_a(end)),(P2+mi)*cos(t_ex_b(end))],[(P2+mi)*sin(t_ex_a(end)),(P2+mi)*sin(t_ex_b(end))],[rf,rf]);

plot((P2+mi)*cos(t_ex_a(end)),(P2+mi)*sin(t_ex_a(end)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','#ffdd8f')
plot((P2+mi)*cos(t_ex_b(end)),(P2+mi)*sin(t_ex_b(end)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])


plot(X1a(end,1),X1a(end,2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot(X1b(end,1),X1b(end,2),'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840])
plot(X1a(:,1),X1a(:,2),'LineWidth',2,'Color','b')
plot(X1b(:,1),X1b(:,2),'LineWidth',2,'Color',"#D95319")

xlim([0.06-0.06  0.06+0.06])
ylim([0.91 1.03])
xlabel('X [-]','FontSize',15,'FontWeight','bold')
ylabel('Y [-]       ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')

legend({'Moon orbit','Arrival orbit','','a) Final Moon position','b) Final Moon position','Final transfer position','',['a) S/C Simple Shooting w/o ',char(8711)],['b) S/C Simple Shooting w/ ',char(8711)]},'Location','southeast','FontSize',15)

nexttile(6)
% subplot(2,3,6)

title('ECIrf zoom in on departure','FontSize',15,'FontWeight','bold')
hold on
grid on
box on
axis equal

Re = 6378/3.84405e5;
[xc,yc,zc] = ellipsoid((P1+mi)*cos(t_ex_b(1)),(P1+mi)*sin(t_ex_b(1)),0,Re,Re,0,100);
surf(xc,yc,zc,'FaceColor', ...
               ...[0.07 0.62  1])
               'none', ...
               'EdgeColor', ...
               [0.07 0.62  1], ...
               'LineWidth',0.6)
plot((P1+mi)*cos(t_ex_b),(P1+mi)*sin(t_ex_b),'Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
circles((P1+mi)*cos(t_ex_b(1)),(P1+mi)*sin(t_ex_b(1)),ri);

plot(X1a(1,1),X1a(1,2),'Marker','o','MarkerSize',7,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','b')
plot(X1b(1,1),X1b(1,2),'Marker','o','MarkerSize',7,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',"#D95319")

plot(X1a(:,1),X1a(:,2),'LineWidth',2,'Color','b')
plot(X1b(:,1),X1b(:,2),'LineWidth',2,'Color',"#D95319")

xlim([-9.5e-3  9.5e-3])
ylim([-18e-3 1e-3])
xlabel('X [-]','FontSize',15,'FontWeight','bold')
ylabel('Y [-]       ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')

legend({'Earth','','Departure orbit','a) Departure transfer position','b) Departure transfer position',['a) S/C Simple Shooting w/o ',char(8711)],['b) S/C Simple Shooting w/ ',char(8711)]},'Location','east','FontSize',15)

end

%--------------------------------------------------------------------------

function plot3(xjF,P1,P2,ri,rf,mi,n)

N = (length(xjF)-2)/4;
tjopt = zeros(1,N);
    for j = 1 : N
        tjopt(j) = xjF(end-1) + (j - 1)*(xjF(end) - xjF(end-1))/(N - 1);
    end

figure('Name',[num2str(2+n),': Ex.2.3 Multiple Shooting N = ',num2str(N)],'NumberTitle','off')
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
    for i = 1:2
        nexttile(i)
%         subplot(1,2,i)
        hold on
        grid on
        box on
        axis equal padded
        xlabel('X [-]','FontSize',15,'FontWeight','bold')
        ylabel('Y [-]       ','FontSize',15,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
    end
nexttile(1)
% subplot(1,2,1)
title('Earth-Moon rotating frame','FontSize',15,'FontWeight','bold')
circles([P1 P2],[0 0],[ri rf]);
plot(P1,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
plot(P2,0,'LineStyle','none','Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])

nexttile(2)
% subplot(1,2,2)
title('ECIrf','FontSize',15,'FontWeight','bold')
circles([(P1+mi)*cos(xjF(end-1)),(P2+mi)*cos(xjF(end))],[(P1+mi)*sin(xjF(end-1)),(P2+mi)*sin(xjF(end))],[ri, rf]);
plot((P1+mi)*cos(0),(P1+mi)*sin(0),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.07 0.62 1])
plot((P2+mi)*cos(xjF(end-1)),(P2+mi)*sin(xjF(end-1)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
plot((P2+mi)*cos(xjF(end)),(P2+mi)*sin(xjF(end)),'Marker','o','MarkerSize',12,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])

for j = 1 : N - 1
    xj = xjF(4*j-3:4*j);
    [XjF,out,t_ex] = orbit_propagation(xj,tjopt(j),tjopt(j+1));
    nexttile(1)
%     subplot(1,2,1)

    plot(xj(1),xj(2),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6,'LineStyle','none')
    plot(XjF(1),XjF(2),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6,'LineStyle','none')
    plot(out(:,1),out(:,2),'k','LineWidth',1.6)
    nexttile(2)
%     subplot(1,2,2)
 
    % PLOT of transfer state in Earth-centered (P1) inertial frame

    x1 = (out(:,1)+mi).*cos(t_ex)-out(:,2).*sin(t_ex);
    y1 = (out(:,1)+mi).*sin(t_ex)+out(:,2).*cos(t_ex);
    x1dot = (out(:,3)-out(:,2)).*cos(t_ex)-(out(:,4)+out(:,1)+mi).*sin(t_ex);
    y1dot = (out(:,3)-out(:,2)).*sin(t_ex)+(out(:,4)+out(:,1)+mi).*cos(t_ex);
    X1 = [x1,y1,x1dot,y1dot]; % state vector

    plot(X1(1,1),X1(1,2),'LineStyle','none','Marker','o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','r')
    plot(X1(end,1),X1(end,2),'LineStyle','none','Marker','o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','r')
    plot(X1(:,1),X1(:,2),'LineWidth',1.6,'Color','k')
    plot((P2+mi)*cos(t_ex),(P2+mi)*sin(t_ex),':','LineWidth',2,'Color',[0.8 0.8 0.8])
    
end
nexttile(1)
% subplot(1,2,1)
legend('','','Earth','Moon',sprintf('N = %d mesh points',N),'','S/C Multiple Shooting transfer','Location','southeast','FontSize',15)
nexttile(2)
% subplot(1,2,2)
legend('','','Earth','Initial Moon position','Final Moon position',sprintf('N = %d mesh points',N),'','S/C Multiple Shooting transfer','Moon orbit','Location','southeast','FontSize',15)

end

