function [pos, vel] = generateTrajectory_ctr(pos0, nEpochs, tr, Ts)
%GENERATETRAJECTORY Summary of this function goes here
%   x0 -- initial starting position
%   nEpochs --  time steps desired
%   tr -- turning rate
%   Ts -- sampling period

sigma_x0 = 500;
sigma_v0 = 20;
sigma_pos = 0.1;
sigma_vel = 0.01;
sigma_h = 0.001;
sigma_vh = 0.001;
sigma_vh0 = 0.01;
sigma_tr = 1e-3*pi;

pos = zeros(nEpochs, 3);
vel = zeros(nEpochs, 3);
pos(1, 1:2) = pos0(1:2) + sigma_x0*randn;
pos(1, 3) = pos0(3);
vel(1, 1:2) = sigma_v0*randn(1, 2); 
vel(1, 3) = sigma_vh0*randn; 

for iEpoch = 2: nEpochs
    OTs = tr * Ts;
    B = [1, sin(OTs)/tr   , 0, -(1-cos(OTs))/tr;
         0, cos(OTs)               , 0, -sin(OTs);
         0, (1-cos(OTs))/tr, 1, sin(OTs)/tr;
         0, sin(OTs)               , 0, cos(OTs)];
    s = B * [pos(iEpoch-1, 1); vel(iEpoch-1, 1);...
            pos(iEpoch-1, 2); vel(iEpoch-1, 2)];
    pos(iEpoch, 1) = s(1) + sigma_pos*randn;
    vel(iEpoch, 1) = s(2) + sigma_vel*randn;
    pos(iEpoch, 2) = s(3) + sigma_pos*randn;
    vel(iEpoch, 2) = s(4) + sigma_vel*randn;
    pos(iEpoch, 3) = pos(iEpoch-1, 3) + vel(iEpoch-1, 3)*Ts + sigma_h*randn;
    vel(iEpoch, 3) = vel(iEpoch-1, 3) + sigma_vh*randn;
    tr = tr + sigma_tr*randn;
end
