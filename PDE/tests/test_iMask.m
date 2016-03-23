% test iMat&iMask
% Assemble IMat&IMask matrix
clear all; close all; clc; %#ok<CLALL>

addpath '../wbas'

i = 10;			% number of basis function
J = 4;			% resolution level

if i == 1
    m = 0;
    k = 0;
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end; % for i
    end; % for j

    [m, k] = find(IMask == i);
    m = 2^(m - 1);
    k = k - 1;
end;

% Display IMask Matrix
disp(IMask)

% Display IMask Matrix
disp(IMat)