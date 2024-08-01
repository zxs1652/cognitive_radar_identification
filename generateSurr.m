% This script is used to generate surrogate data for null hypothesis
% -- non-cognitive data -- %
clear
close all
clc

%% Parameter Setting
rng(3)
% ====== radar characteristic ====== %
nRadar = 8;
totalPower = 1e6;
minPower = 0.02*totalPower;
maxPower = 0.6*totalPower;
initPowers = totalPower/nRadar * ones(nRadar, 1);
ftx = 10e6; % sampling frequency for trasmission
Ttx = 1/ftx; % sampling period for transmission
frx = 10e6; % sampling frequency for reception
Trx = 1/frx; % sampling period for reception
c = 299792458; % speed of light
fc = 1e9; % carrier frequency
lambda = c./fc; % wave length
% g = 0.1 * ones(nRadar, 1); % antenna gain
% h = 1 * ones(nRadar, 1); % RCS
% A = 10 * ones(nRadar, 1); % reception area
% loss = 1.4e-13; % other path loss

radius = 2e3; % radar alligned in a circle shape
angle = 2*pi/nRadar.*(0:1:nRadar-1)';
radarPosX = radius*cos(angle);
radarPosY = radius*sin(angle);
radarPos = [radarPosX,radarPosY,zeros(nRadar, 1)];

%% Base Signal Generation
% === radar waveform === %
prf = 1e3; % Pulse repetition frequency
nSampleTx = 1/prf*ftx;
nSampleRx = nSampleTx;
t = 1/ftx: 1/ftx: 1/prf;
s0 = zeros(nRadar, length(t));
f0 = 0; % linear FM sweep frequency
f1 = 50e3;
for iRadar = 1: nRadar  
    s0(iRadar,:) = chirp(t,f0,t(end),f1,'linear',0,'complex');
    f1 = f1 + 50e3;
end
% === received signal === %
% generate original sampled code
ii = (1:1:nSampleRx);
sampleInd0 = 1 + mod(floor((ii-1)/(frx/ftx)), nSampleTx);
code0 = s0(:,sampleInd0);

radarInfo = zeros(nRadar + 4, 3);
radarInfo(1:nRadar,:) = radarPos;
radarInfo(1:nRadar, 4) = lambda;
radarInfo(nRadar + 1, 1) = nSampleRx;
radarInfo(nRadar + 2, 1) = nSampleTx;
radarInfo(nRadar + 3, 1) = frx;
radarInfo(nRadar + 4, 1) = ftx;
radarInfo(nRadar + 5, 1) = fc;

%% generate data
nDataset = 1;
pr = cell(nDataset, 1);
rho = cell(nDataset, 1);
% hwb = waitbar(0,'1','Name','Generating Surrogate Data...');

for iDataset = 1:nDataset
    % ====== display ====== %
    iDataset
    % waitbar(iDataset/nDataset,hwb,...
    %     sprintf('%d out of  %d', iDataset, nDataset))
    % ====== trajrctory parameter ====== %
    ft = 1; % tracking frequency
    T = 1/ft;
    % general trajectory
    x0 = -3000; y0 = 3000;
    xBound = 3000; yBound = 3000;
    xInterval = 100; yInterval = 100;
    height = 1000;
    i = 1;
    tarPosX = zeros(3000, 1);
    tarPosY = zeros(3000, 1);
    tarVelX = zeros(3000, 1);
    tarVelY = zeros(3000, 1);
    tarPosX(1,:) = x0; tarPosY(1,:) = y0; 
    tarVelX(1,:) = 0; tarVelY(1,:) = 0; 
    right = 1; 
    while tarPosX(i) < xBound || tarPosY(i) > -yBound
        if right == 1 && tarPosX(i) < xBound
            tarVelX(i,:) = xInterval + 10*randn; tarVelY(i,:) = 0;
            tarPosX(i+1,:) = tarPosX(i,:) + tarVelX(i,:)*T;
            tarPosY(i+1,:) = tarPosY(i,:) + tarVelY(i,:)*T;
        end
        if right == 1 && tarPosX(i) >= xBound
            right = 0;
            tarVelX(i,:) = 0; tarVelY(i,:) = -yInterval + 10*randn;
            tarPosX(i+1,:) = tarPosX(i,:) + tarVelX(i,:)*T;
            tarPosY(i+1,:) = tarPosY(i,:) + tarVelY(i,:)*T;
            i = i + 1;
            continue;
        end
        if right == 0 && tarPosX(i) > -xBound
            tarVelX(i,:) = -xInterval + 10*randn; tarVelY(i,:) = 0; 
            tarPosX(i+1,:) = tarPosX(i,:) + tarVelX(i,:)*T;
            tarPosY(i+1,:) = tarPosY(i,:) + tarVelY(i,:)*T;

        end
        if right == 0 && tarPosX(i) <= -xBound
            right = 1;
            tarVelX(i,:) = 0; tarVelY(i,:) = -yInterval + 10*randn;
            tarPosX(i+1,:) = tarPosX(i,:) + tarVelX(i,:)*T;
            tarPosY(i+1,:) = tarPosY(i,:) + tarVelY(i,:)*T;
            i = i + 1;
            continue;
        end
        i = i + 1;
    end
    tarPos = [tarPosX,tarPosY,height*ones(size(tarPosX))];
    tarVelX(i,:) = 0; tarVelY(i,:) = -yInterval + 2*randn;
    tarVel = [tarVelX,tarVelY,zeros(size(tarVelX))];
    nTrack = i;
    
    % specific trajectory
    % duration = 300;
    % x0 = -2000; y0 = 2000;
    % tarVelX = 25;
    % tarVelY = -10;
    % % tarVelX = 0;
    % % tarVelY = 0;
    % tMove = (0:T:duration-T)';
    % nTrack = length(tMove);
    % tarPosX = x0 + tarVelX * tMove + 50*randn(300,1);
    % tarPosY = y0 + tarVelY * tMove + 50*randn(300,1);
    % height = 1000;
    % tarPos = [tarPosX,tarPosY,height*ones(length(tMove),1)];
    % tarVel = repmat([tarVelX,tarVelY,0], nTrack, 1);
    % ====== visualize trajrctory ====== %
    % figure(3);% visualize the trajectory
    % plot3(radarPos(:,1),radarPos(:,2),radarPos(:,3),...
    %     'o','Color','b','MarkerSize',10,...
    %     'MarkerFaceColor','#D9FFFF')
    % radar_number = {'1','2','3','4','5','6','7','8'};
    % text(radarPos(:,1),radarPos(:,2),radarPos(:,3),...
    %     radar_number,FontSize=16)
    % hold on
    % plot3(tarPos(:,1), tarPos(:,2), tarPos(:,3),...
    %     '-o', 'MarkerSize', 5)
    % grid on
        
    %% received signal without noise   
    % === power allocation policy1 (evenly distributed) === %
    p1 = zeros(nRadar, nTrack);
    distance = zeros(nRadar, nTrack);
    for iTrack = 1:nTrack
        prop1 = 1 / nRadar * ones(nRadar, 1);
        p1(:, iTrack) = totalPower*prop1;
        georange = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
        distance(:,iTrack) = georange;
    end
    
    %% power estimation 
    %% On Target Side 
    % pRx1 = zeros(nRadar, nTrack);
    pRxEst1 = zeros(nRadar, nTrack);
    for iTrack = 1:nTrack
        % clarify the DOA
        tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        [az, el] = generateDOA(radarPos, tarPos(iTrack,:));
        doa = [az, el]';
        %% policy 1
        % generate received signal
        [rTarget1, ~, ~, alphaTarget1] = ...
            generateRxSignalTarget(p1(:,iTrack), radarInfo,tarInfo,code0);
        % pRx1(:,iTrack) = alphaTarget1.^2;
        % create the sensor array
        array = phased.URA('Size',[11 11],'ElementSpacing',[lambda/2 lambda/2],'ArrayNormal','z');
        % sterring vector
        steeringvec1 = phased.SteeringVector('SensorArray', array);
        % collected signal
        xArray1 = collectPlaneWave(array,rTarget1.',doa,fc);
        X1 = xArray1.';
        % beamforming
        nAnt1 = size(xArray1,2);
        w1 = zeros(nAnt1, nRadar);
        yBeam1 = zeros(nRadar, nSampleRx);
        Hsv1 = zeros(nAnt1, nRadar);
        for iRadar = 1: nRadar
            Hsv1(:, iRadar) = steeringvec1(fc, doa(:, iRadar));
        end
        % correlation matrix
        Rxx1 = X1*X1'/nSampleRx + 1e5*eye(size(Hsv1, 1));
        for iRadar = 1: nRadar
            % identify constraints
            P = zeros(1, nRadar);
            P(iRadar) = 1;
            % beamformer
            Hbeam = Hsv1;
            % compute the beam weights
            w1(:, iRadar) = Rxx1^(-1)*Hbeam*(Hbeam'*Rxx1^(-1)*Hbeam)^(-1)*P';
            % beamer output
            yBeam1(iRadar, :) = w1(:, iRadar)' * X1;
        end
        pRxEst1(:, iTrack) = sum(abs(yBeam1).^2, 2)/nSampleRx;
    end
    
    pr{iDataset} = pRxEst1;
    rho{iDataset} = distance;
end
save('surr_pr.mat','pr')
save('surr_rho.mat', 'rho')