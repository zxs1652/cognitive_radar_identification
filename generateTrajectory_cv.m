function [pos, vel] = generateTrajectory_cv(height, nEpochs,Ts)
%GENERATETRAJECTORY This function generates random trajectory following
%constant velocity model
%   x0 -- initial starting position
%   nEpochs --  time steps desired
%   tr -- turning rate
%   Ts -- sampling period
% Define larger rectangle [xMin, yMin, width, height]
largerRect = [-3000, -3000, 5000, 5000];

% Define smaller rectangle [xMin, yMin, width, height]
smallerRect = [-1500, -1500, 3000, 3000];

dis = 0;
while min(dis) < 3500
    % random select start point and end point
    start_xy = randomPointInRectExclude(largerRect, smallerRect);
    end_xy = randomPointInRectExclude(largerRect, smallerRect);
    
    
    % generate path
    path_x = linspace(round(start_xy(1)), round(end_xy(1)), nEpochs)';
    path_y = linspace(round(start_xy(2)), round(end_xy(2)), nEpochs)';
    path_z = height * ones(nEpochs, 1);
    pos = [path_x, path_y, path_z];
    dis = norm(start_xy - end_xy);
end


vel_x = (round(end_xy(1)) - round(start_xy(1))) / ((nEpochs-1) * Ts) * ones(nEpochs, 1);
vel_y = (round(end_xy(2)) - round(start_xy(2))) / ((nEpochs-1) * Ts) * ones(nEpochs, 1);
vel_z = zeros(nEpochs, 1);
vel = [vel_x, vel_y, vel_z];
end