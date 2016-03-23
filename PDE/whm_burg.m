%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Solving 1-D Burgers equation evolution problem by
%                  Wavelet Haar Method (WHM)
%
%               u_t + u u_x = \nu u_{xx}, u=u(x,t), x\in[-1,1], t>=0,
%               IC: u(x,0)      = -\sin(\pi x),
%               BC: u(+/-1,t)   = 0. 
%
%              coded by Oleg Kravchenko, BMSTU, 2013.07.01
%              UPD1: 2016.03.23
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] Ulo Lepik, "Numerical solution of evolution equations by the Haar
%       wavelet method,"Applied Mathematics and Computation 185 (2007) 
%       695–704  
%     link: http://www.sciencedirect.com/science/article/pii/S0096300306009052
%     link: https://www.researchgate.net/profile/Uelo_Lepik/publication/...
%               251142138_Haar_wavelet_method_for_solving_integral_equations_and_evolution_equations/...
%               links/5591507b08aed6ec4bf829ed.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Evolution problem, Wavelet Haar Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc; %#ok<*CLALL>

addpath 'wbas'

%% Wavelet Haar method

% step 1
% collocation points
J = 7;                          % level of resolution
N = 2^(J + 1);                  % N = 2M number of wavelets
j = 1:N;                        % space grid counter
x = (j - 0.5) / N;              % space grid
tin  = 0.00;                    % start time
tfin = 0.30;                    % final time
dt = 1e-4;                      % time step
ts = tin:dt:tfin;               % time grid
nu = 1 / (400*pi);              % nu parameter

%% Initial vectors
H = zeros(N,N);
P = H;
Q = P;

for i = 1:N    
    H(i,:) = h(x, i, J);
    P(i,:) = p(x, i, J);
    Q(i,:) = q(x, i, J);
    disp(['H,P,Q matrixes formation: ' num2str(i) ' / ' num2str(N)])
end;

u0      = sin(2*pi*x);
d1u0    = 2*pi*cos(2*pi*x);
d2u0    = - 4*pi^2*sin(2*pi*x);
u       = u0;
d1u     = d1u0;
d2u     = d2u0;

as = zeros(1,N);        % constant on a sub-interval [t_s, t_s+1]
E  = ones(N,1);         % unit vector
qt = zeros(N,1) + 0.5;
for i = 2:N
    tmp = qtilde(i,J);	% qtilde vector
    qt(i) = 0.25 / tmp^2;   
end

GMatrix = zeros(length(ts),N,N);
for nt = 1:length(ts)
    for nl = 1:length(x)
        GMatrix(nt,:,nl) = Q(:,nl) - qt(:) * x(nl);
    end; % for nl
    disp(['GMatrix assemble: ' num2str(nt) ' / ' num2str(length(ts))])
end; % for nt

GM = zeros(N,N);
bvect = zeros(1,N);
UMatrix = zeros(length(ts),N);
%% Computation (Main Loop)
for nt = 1:length(ts)
    % GM and bvect computation
    for nl = 1:length(x)
        bvect(nl) = - 0.5 * u(nl) * d1u(nl) + nu * d2u(nl);
    end; % for nl
    GM(:,:) = GMatrix(nt,:,:);
    
    % as vector computation
    as = bvect / GM;
    
    % Next time step
    % u, d1u, d2u recomputation
    for nl = 1:length(x)
        d2u(nl) = dt * as * H(:,nl) + nu * d2u(nl);
        d1u(nl) = dt * as * (P(:,nl) - E(nl) * qt(:)) + d1u(nl);
        u(nl)     = dt * as * (Q(:,nl) - x(nl) * qt(:)) + u(nl);
    end
    
    % BC
%     u(1) = u0(1);
%     u(N) = u0(N);
    
    UMatrix(nt,:) = u;
    
    disp(['time step: ' num2str(nt) ' / ' num2str(length(ts))])
end % for nt

%% Plots
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

% fig:01
figure('color','w')
plot(x,u0,x,u,'r.-',[x(1) x(end)],[u0(1) u0(N)],'ko','MarkerFaceColor','g')
xlabel('$x$'); ylabel('$u(x,t)$')
title('Solution of Burgers equation by WHM')
text(0.65,0.5,['$t = ' num2str(tfin) '$'],'FontSize',12)