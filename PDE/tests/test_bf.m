% Plot basis functions
clear all; close all; clc; %#ok<CLALL>

addpath '../wbas'

nx = 2^11;					% number of grid points
x = 0:1/(nx-1):1;			% spatial grid
J = 2;						% resolution level

figure('color','w')
N = 2^J;
offset = 1e-1;
for i=1:N
    [y,m,k] = h(x,i,J);
    subplot(N,1,i), 
    plot(x,y,'linewidth',2)
    axis([min(x)-offset max(x)+offset min(y)-offset max(y)+offset])
    xlabel('$x$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['$\psi_{(' num2str(m) ',' num2str(k) ')}$'],'interpreter','latex')
end

