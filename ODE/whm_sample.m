%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Solution of 1D nonlinear IVP test problem by
%                  Wavelet Haar Method (WHM)
%
%               y''(x) = - 2 y(x) y'(x) , y(0) = 0, y'(0) = - 1,
%               exact solution: y(x) = - tan(x)
%           
%
%              coded by Oleg Kravchenko, BMSTU, 2012.07.20
%              UPD1: 2016.03.21
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] Siraj-ul-Islam, Imran Aziz, Bozidar Sarler, "The numerical solution 
%       of second-order boundary-value problems by collocation method with 
%       the Haar wavelets,"Mathematical and Computer Modelling, Vol. 52, 
%       No. 9-10, 1577-1590, 2012.
%     link: http://www.sciencedirect.com/science/article/pii/S0895717710003006
% [2] Sahoo, Bishnupriya, "A study on solution of differential equations 
%       using Haar wavelet collocation method, MSc thesis, 2012.
%     link: http://ethesis.nitrkl.ac.in/3276/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Nonlinear IVP, ode, WHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc %#ok<CLALL>

addpath 'wbas' 'models'

%% Initial Parameters
flag = 0;                   % 0 - don't save data files, 1 - save data files
%% Wavelet Haar method

% step 1
% collocation points
J = 3;                              % level of decomposition
N = 2^(J + 1); % N = 2M             % number of basis functions
j = 1:N;                            % index of grid points
x = (j - 0.5) ./ N;                 % grid points

% step 2
% initial values
alpha1 = 0;                         % initial value of a function
beta1  = - 1;                       % initial value of the first derivative
a1     = beta1 - alpha1;

% step 3
% Newton solver
W = zeros(N,N);
f = zeros([N 1]);
a = zeros([N 1]);

epsilon = 1e-4;
r = ones([N 1]);

iter_ind = 0;
tic
while max(r) > epsilon   
    for j = 1:N
        % f(x) computation 
        % H, P1, P2 computation
        H = 0;               
        P1 = 0;
        P2 = 0;
        for i = 1:N
            H  = H  + a(i) * haar(x(j), i, J);
            P1 = P1 + a(i) * p1(x(j), i, J);
            P2 = P2 + a(i) * p2(x(j), i, J);            
        end;
        
        f(j) = 2 * (alpha1 + beta1 * x(j) + P2) * ...
            (beta1 + P1) + H;  
        
        
        % W(x) matrix computation        
        for k = 1:N
            W(j,k) = 2 * p2(x(j),k,J) * (beta1 + P1) + ...
                2 * (alpha1 + beta1 * x(j) + P2) * p1(x(j),k,J) + haar(x(j),k,J);            
        end; % for k
    end; % for j
    
    a_new = W \ (W*a - f);      % linear system solution
    r = abs(a_new - a);         % residual 
    disp(['iteration: ' num2str(iter_ind) ' error Newton: ' num2str(max(r))])   
    
    % Update variables
    a = a_new;
    iter_ind = iter_ind + 1;
end; % while
toc

% Reconstruct approximate solution
y = zeros(N,1);
for j = 1:N    
    S = 0;
    for i = 1:N
        S = S + a(i) * p2(x(j),i,J);
    end
    y(j) = alpha1 + x(j) * beta1 + S;
end; % for

%% Exact solution
yexact = - tan(x);
% critical point pi/2 ~= 1.5708
x_zero1 = 0.5 * pi; 

%% Runge - Kutta method
[x, y1] = ode113('model0', x, [alpha1 beta1]);

%% Plot graphics
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

oft = 0.01;

% fig:01
figure('color','w')
plot(x,yexact,'g',x,y,'rs',x,y1(:,1),'b.')
xlabel('$x$'); ylabel('$y$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Exact','WHM', 'RGK')
axis([-oft 1+oft min(yexact)-oft max(yexact)+oft])

% Absolute errors
rRGK = abs(y1(:,1) - yexact');
rWHM = abs(y - yexact');
rRW = abs(y - y1(:,1));

% fig:02
figure('color','w')
plot(x,rRGK,'b.-',x,rWHM,'r.-',x,rRW,'ms-')
xlabel('$x$'); ylabel('Absolute Error');
title('Absolute Error: $\max|y_{numeric} - y_{analytic}|$')
legend('RGK','WHM','Between RGK and WHM',...
    'Location','northoutside','Orientation','horizontal')
axis([-oft 1+oft min([rRGK; rWHM; rRW])-oft max([rRGK; rWHM; rRW])+oft])

%% Disp Errors
disp(['error RGK: ' num2str(max(rRGK)) ' error WHM: ' num2str(max(rWHM)) ...
    ' error RW: ' num2str(max(rRW))])

%% Save data
if flag == 1    
    cd 'dat'
    
    table0 = [x yexact' y y1(:,1)];
    fid = fopen('table0.txt','w');
    fprintf(fid, '%6.2f %6.2f %6.2f %6.2f\n', table0');
    
    fclose(fid);
    disp('Saved.')
    
    cd '..'
end