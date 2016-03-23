function [dydt] = rhs2(t, y, par)
    %
    % function returns the values of the rhs's
    % for the system:
    %    dy1/dt = y2
    %    dy2/dt = -y1
    %
    dydt=[y(2); -y(1)];
