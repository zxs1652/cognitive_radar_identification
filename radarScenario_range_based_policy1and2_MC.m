% implement radar scenario simulation
% the radar is tracking a self-defined trajectory in 3 dimensions
clear
close all
clc

%% Parameter Setting
% rng(3)
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

% ====== tracking parameter ====== %
ft = 1; % tracking frequency
T = 1/ft;
% duration = 300;
% x0 = -2000; y0 = 2000;
% tarVelX = 25;
% tarVelY = -10;
% % tarVelX = 0;
% % tarVelY = 0;
% tMove = (0:T:duration-T)';
% nTrack = length(tMove);
% tarPosX = x0 + tarVelX * tMove;
% tarPosY = y0 + tarVelY * tMove;
% height = 1000;
% tarPos = [tarPosX,tarPosY,height*ones(length(tMove),1)];
% tarVel = repmat([tarVelX,tarVelY,0], nTrack, 1);
pos0 = [0, 0, 1000];
tr = 1e-2*pi;
nTrack = 300;

% pos0 = tarPos(1, :)' + 10*randn(3,1);
% vel0 = tarVel(1, :)' + 2*randn(3,1);
% x0 = [pos0;vel0];
% P0 = blkdiag(10^2*eye(3), 2^2*eye(3));
% F = [1, 0, 0, T, 0, 0;
%      0, 1, 0, 0, T, 0;
%      0, 0, 1, 0, 0, T;
%      0, 0, 0, 1, 0, 0;
%      0, 0, 0, 0, 1, 0;
%      0, 0, 0, 0, 0, 1];
% % Q = blkdiag(1^2*eye(3), 1^2*eye(3));
% Q = zeros(6,6);
% sigmaTau0 = 1e-7;
% R0 = sigmaTau0^2*eye(nRadar);


%% Monte-Carlo Experiments for collect data
nMC = 400;
% memory for test data
prx1 = cell(nMC, 1);
prx2 = cell(nMC, 1);
rho = cell(nMC, 1);
target_pos = cell(nMC, 1);
parfor iMC = 1:nMC
    iMC
    %% Signal Generation
    % === generate trajectory === %
    [tarPos, tarVel] = generateTrajectory(pos0, nTrack, tr, T);
    % ====== trajectory visualization ====== %
    % h3 = figure(3);% visualize the trajectory
    % plot3(radarPos(:,1),radarPos(:,2),radarPos(:,3),...
    %     'o','Color','b','MarkerSize',10,...
    %     'MarkerFaceColor','#D9FFFF')
    % radar_number = {'1','2','3','4','5','6','7','8'};
    % text(radarPos(:,1),radarPos(:,2),radarPos(:,3),...
    %     radar_number,FontSize=16)
    % hold on
    % scatter3(tarPos(:,1), tarPos(:,2), tarPos(:,3))
    % grid on
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
    
    % figure(1) % visualize the spectrum
    % % t = unigrid(0,1/waveform.SampleRate,1/waveform.PRF,'[)');
    % plot(t, real(s0(1,:))) % visualize the waveform
    % figure(2)
    % pspectrum(s0(1,:),ftx,'spectrogram', 'FrequencyLimits',[-500e3, 500e3],...
    %   'OverlapPercent',99,'Leakage',0.85) 
    
    %% received signal without noise
    radarInfo = zeros(nRadar + 4, 3);
    radarInfo(1:nRadar,:) = radarPos;
    radarInfo(1:nRadar, 4) = lambda;
    radarInfo(nRadar + 1, 1) = nSampleRx;
    radarInfo(nRadar + 2, 1) = nSampleTx;
    radarInfo(nRadar + 3, 1) = frx;
    radarInfo(nRadar + 4, 1) = ftx;
    radarInfo(nRadar + 5, 1) = fc;
    sigmaW = 1;
    geoRange = zeros(nRadar, nTrack);
    tau = zeros(nRadar, nTrack);
    for iTrack = 1:nTrack
        tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        geoRange(:,iTrack) = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
        tau(:,iTrack) = 2*geoRange(:,iTrack)./c;
    end
    
    % === power allocation policy1 (evenly distributed) === %
    rRadar_no_noise1 = cell(nTrack, 1); % received signal of radar
    p1 = zeros(nRadar, nTrack);
    alpha1 = zeros(nRadar,nTrack);
    sigmaTau1 = zeros(nRadar,nTrack);
    distance = zeros(nRadar, nTrack);
    for iTrack = 1:nTrack
        % tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        prop1 = 1 / nRadar * ones(nRadar, 1);
        p1(:, iTrack) = totalPower*prop1;
        % [rRadar_no_noise1{iTrack}, ~, ~, alpha1(:, iTrack)] = ...
        %     generateRxSignalRadar(p1(:, iTrack), radarInfo,tarInfo,code0);
        georange = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
        distance(:,iTrack) = georange;
        alpha1(:, iTrack) = sqrt(p1(:, iTrack).* (1./(georange)).^4);
        sigmaTau1(:,iTrack) = sqrt(sigmaW^2./(8*pi^2*fc^2*alpha1(:, iTrack).^2));
    end
    
    % === power allocation policy2 (proportional to the distance^4) === %
    rRadar_no_noise2 = cell(nTrack, 1); % received signal of radar
    p2 = zeros(nRadar, nTrack);
    alpha2 = zeros(nRadar,nTrack);
    sigmaTau2 = zeros(nRadar,nTrack);
    for iTrack = 1:nTrack
        % tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        prop2 = (1./geoRange(:,iTrack))./sum(1./geoRange(:,iTrack));
        p2(:, iTrack) = totalPower*prop2;
        % [rRadar_no_noise2{iTrack}, ~, ~, alpha2(:, iTrack)] = ...
        %     generateRxSignalRadar(p2(:, iTrack), radarInfo,tarInfo,code0);
        georange = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
        alpha2(:, iTrack) = sqrt(p2(:, iTrack).* (1./(georange)).^4);
        sigmaTau2(:,iTrack) = sqrt(sigmaW^2./(8*pi^2*fc^2*alpha2(:, iTrack).^2));
    end
    %% On Target Side 
    pRx1 = zeros(nRadar, nTrack);
    pRxEst1 = zeros(nRadar, nTrack);
    pRx2 = zeros(nRadar, nTrack);
    pRxEst2 = zeros(nRadar, nTrack);
    for iTrack = 1:nTrack
        
        % clarify the DOA
        tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        [az, el] = generateDOA(radarPos, tarPos(iTrack,:));
        doa = [az, el]';
        %% policy 1
        % generate received signal
        [rTarget1, ~, ~, alphaTarget1] = ...
            generateRxSignalTarget(p1(:,iTrack), radarInfo,tarInfo,code0);
        pRx1(:,iTrack) = alphaTarget1.^2;
        % create the sensor array
        array = phased.URA('Size',[11 11],'ElementSpacing',[lambda/2 lambda/2],'ArrayNormal','z');
        % sterring vector
        steeringvec1 = phased.SteeringVector('SensorArray', array);
        % collected signal
        [xArray1, sv] = collectPlaneWave(array,rTarget1.',doa,fc);
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


        %% policy 2
         % generate received signal
        [rTarget2, ~, ~, alphaTarget2] = ...
            generateRxSignalTarget(p2(:,iTrack), radarInfo,tarInfo,code0);
        pRx2(:,iTrack) = alphaTarget2.^2;
        % create the sensor array
        array = phased.URA('Size',[11 11],'ElementSpacing',[lambda/2 lambda/2],'ArrayNormal','z');
        % sterring vector
        steeringvec2 = phased.SteeringVector('SensorArray', array);
        % collected signal
        [xArray2, sv] = collectPlaneWave(array,rTarget2.',doa,fc);
        X2 = xArray2.';
        % beamforming
        nAnt2 = size(xArray2,2);
        w2 = zeros(nAnt2, nRadar);
        yBeam2 = zeros(nRadar, nSampleRx);
        Hsv2 = zeros(nAnt2, nRadar);
        for iRadar = 1: nRadar
            Hsv2(:, iRadar) = steeringvec2(fc, doa(:, iRadar));
        end
        % correlation matrix
        Rxx2 = X2*X2'/nSampleRx + 1e5*eye(size(Hsv2, 1));
        for iRadar = 1: nRadar
            % identify constraints
            P = zeros(1, nRadar);
            P(iRadar) = 1;
            % beamformer
            Hbeam = Hsv2;
            % compute the beam weights
            w2(:, iRadar) = Rxx2^(-1)*Hbeam*(Hbeam'*Rxx2^(-1)*Hbeam)^(-1)*P';
            % beamer output
            yBeam2(iRadar, :) = w2(:, iRadar)' * X2;
        end
        pRxEst2(:, iTrack) = sum(abs(yBeam2).^2, 2)/nSampleRx;

        % h8 = figure(8); % plot beamformer output
        % t = (1/frx:1/frx:nSampleRx/frx);
        % subplot(1,2,1);
        % plot(t,real(yBeam0(1 ,:))); axis tight
        % title('Output of MVDR Beamformer for Antenna Array');
        % xlabel('Time (s)');ylabel('Real(y)');
        % subplot(1,2,2);
        % plot(t,real(r0Target(1 ,:))); axis tight
        % title('Original Transmitted Signal to Target');
        % xlabel('Time (s)');ylabel('Real(y)');
        % 
        % h9 = figure(9); % Plot array response with beamformer weighting
        % azPlot = (-180: 0.5: 180);
        % elPlot = (-90: 0.5: 0);
        % powerPlot = zeros(length(azPlot), length(elPlot));
        % for i = 1: length(azPlot)
        %     for j = 1: length(elPlot)
        %         svPlot = steeringvec0(fc, [azPlot(i); elPlot(j)]);
        %         powerPlot(i, j) = abs(w0(:,1)' * svPlot)^2; 
        %     end
        % end
        % [elPlot, azPlot] = meshgrid(elPlot, azPlot);
        % mesh(elPlot, azPlot, -10*log10(powerPlot))
        % hold on 
        % for i = 1:length(doa(1,:))
        %     plot3([doa(2,i),doa(2,i)], [doa(1,i),doa(1,i)], [100,-5],...
        %         'LineWidth',2,'Color',"#A2142F")
        % end
        % title('Beamformer Negative Attenuation at DOAs');
        % xlabel('Elevation');ylabel('Azimuth');zlabel('Power(dB)')
        % legend('Negative Attenuation', 'True DOAs')

    end
    % record test data
    prx1{iMC} = pRxEst1;
    prx2{iMC} = pRxEst2;
    rho{iMC} = distance;
    target_pos{iMC} = tarPos;
end
%% save test data
save("dataset/prx1.mat", "prx1")
save("dataset/prx2.mat", "prx2")
save("dataset/rho.mat", "rho")
save("dataset/target_pos.mat", "target_pos")
%% Data Analysis
% % ====== load data ====== %
% load('dataset\noncognitive\prx.mat')
% load("dataset\noncognitive\rho.mat")
% load("dataset\cognitive\policy2\prx2.mat")
% load('dataset/surrogate/mutualInfo_surr_3200_distanceRatio.mat')
% % ====== sequence mutual information ======%  
% iDataset = 1;
% pr1_seq = prx{iDataset};
% pr2_seq = prx2{iDataset};
% rho_seq = rho{iDataset};
% % estimate marginal pdf
% pts1 = 0:0.005:1;
% pts2 = 0:0.005:1;
% % normalize power
% prx1EstNormalized = pr1_seq .* (rho_seq.^2);
% prx1EstNormalized = prx1EstNormalized / max(prx1EstNormalized(:));
% prx2EstNormalized = pr2_seq .* (rho_seq.^2);
% prx2EstNormalized = prx2EstNormalized / max(prx2EstNormalized(:));
% rhoNormalized = rho_seq.^4;
% rhoTotal = sum(rhoNormalized, 1);
% rhoNormalized = rhoNormalized ./ repmat(rhoTotal, nRadar, 1);
% 
% mutualInfo1_seq = zeros(nRadar, nTrack-10);
% mutualInfo2_seq = zeros(nRadar, nTrack-10);
% decision1 = zeros(nRadar, nTrack-10);
% decision2 = zeros(nRadar, nTrack-10);
% p_fa1 = zeros(nRadar, nTrack-10);
% p_fa2 = zeros(nRadar, nTrack-10);
% for iTrack = 1:nTrack - 10
%     iTrack
%     I_prx1_rho = zeros(nRadar,1);
%     I_prx2_rho = zeros(nRadar,1);
%     parfor iRadar = 1:nRadar
%         I_prx1_rho(iRadar) = mutualInfo(prx1EstNormalized(iRadar,1:10 + iTrack),...
%             rhoNormalized(iRadar,1:10 + iTrack), pts1, pts2);
%         I_prx2_rho(iRadar) = mutualInfo(prx2EstNormalized(iRadar,1:10 + iTrack),...
%             rhoNormalized(iRadar,1:10 + iTrack), pts1, pts2);
%     end
%     mutualInfo1_seq(:,iTrack) = I_prx1_rho;
%     mutualInfo2_seq(:,iTrack) = I_prx2_rho;
%     parfor iRadar = 1:nRadar
%         [decision1(iRadar, iTrack), prob1] = kstest2(mutualInfo_surr(iRadar,:), mutualInfo1_seq(iRadar, 1:iTrack));
%         p_fa1(iRadar, iTrack) = prob1;
%         [decision2(iRadar, iTrack), prob2] = kstest2(mutualInfo_surr(iRadar,:), mutualInfo2_seq(iRadar, 1:iTrack));
%         p_fa2(iRadar, iTrack) = prob2;
%     end
% end
% 
% % ===== Mutual Information ===== %
% h11 = figure(11);
% iPlot = 1; 
% plot(1:1:nTrack-10, mutualInfo1_seq)
% legend('Radar 1', 'Radar 2', 'Radar 3', 'Radar 4',...
%     'Radar 5', 'Radar 6', 'Radar 7', 'Radar 8')
% h12 = figure(12);
% plot(1:1:nTrack-10, p_fa)
% legend('Radar 1', 'Radar 2', 'Radar 3', 'Radar 4',...
%     'Radar 5', 'Radar 6', 'Radar 7', 'Radar 8')



%% function definition
function [x] = const_vel_transition_function(xPrev, T)
    F = [1, 0, 0, T, 0, 0;
        0, 1, 0, 0, T, 0;
        0, 0, 1, 0, 0, T;
        0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, 1, 0;
        0, 0, 0, 0, 0, 1];
    x = F*xPrev;
end

function [y] = delay_measurement_function(x, nRadar, radarPos)
    c = 299792458;
    pos = x(1:3);
    vel = x(4:6);
    georange = vecnorm(repmat(pos',nRadar,1) - radarPos, 2, 2);
    range = 2*georange;
    y = range./c;

end

function [H_gradient] = gradient_H(tarPos,radarPos)
    c = 299792458;
    nRadar = size(radarPos,1);
    diff = tarPos - radarPos;
    ranges = vecnorm(diff,2,2);
    H_gradient = [2/c * diff./repmat(ranges, 1, 3), zeros(nRadar, 3)];

end


%% test
% for i = 1:6
%     figure(10)
%     plot(real(code0(i,:)));
%     hold on
% 
% end

% SNR
% for iRadar = 1:nRadar
%     a1 = r2Radar{1}(iRadar,:)*r2Radar{1}(iRadar,:)'/nSampleRx;
% 
%     noise = sqrt(1/2)*randn(1,nSampleRx) + 1i*sqrt(1/2)*randn(1,nSampleRx);
%     r22 = r2Target(iRadar,:) + noise;
%     a2 = r22*r22'/nSampleRx;
% 
%     snrRadar(iRadar,1) = 20*log10(a1);
%     snrTarget(iRadar,1) = 20*log10(a2);
% end
% snrRadar'
% snrTarget'

%

