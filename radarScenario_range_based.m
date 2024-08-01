% implement radar scenario simulation
% the radar is tracking a self-defined trajectory in 3 dimensions
clear
close all
clc

%% Parameter Setting
rng(1)
% ====== radar characteristic ====== %
nRadar = 6;
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
duration = 400;
x0 = -2000; y0 = 2000;
tarVelX = 25;
tarVelY = -10;
% tarVelX = 0;
% tarVelY = 0;
tMove = (0:T:duration-T)';
nTrack = length(tMove);
tarPosX = x0 + tarVelX * tMove;
tarPosY = y0 + tarVelY * tMove;
height = 1000;
tarPos = [tarPosX,tarPosY,height*ones(length(tMove),1)];
tarVel = repmat([tarVelX,tarVelY,0], nTrack, 1);

pos0 = tarPos(1, :)' + 10*randn(3,1);
vel0 = tarVel(1, :)' + 2*randn(3,1);
x0 = [pos0;vel0];
P0 = blkdiag(10^2*eye(3), 2^2*eye(3));
F = [1, 0, 0, T, 0, 0;
     0, 1, 0, 0, T, 0;
     0, 0, 1, 0, 0, T;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1];
% Q = blkdiag(1^2*eye(3), 1^2*eye(3));
Q = zeros(6,6);
sigmaTau0 = 1e-7;
R0 = sigmaTau0^2*eye(nRadar);

%% Signal Generation
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
for iTrack = 1:nTrack
    % tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
    prop1 = 1 / nRadar * ones(nRadar, 1);
    p1(:, iTrack) = totalPower*prop1;
    % [rRadar_no_noise1{iTrack}, ~, ~, alpha1(:, iTrack)] = ...
    %     generateRxSignalRadar(p1(:, iTrack), radarInfo,tarInfo,code0);
    georange = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
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
%% Monte-Carlo Experiments for RMSE
nMC = 100;
posEst_MC1 = zeros(3, nTrack+1, nMC);
tauEst_MC1 = zeros(nRadar, nTrack, nMC);
posEst_MC2 = zeros(3, nTrack+1, nMC);
tauEst_MC2 = zeros(nRadar, nTrack, nMC);
posEst_MC3 = zeros(3, nTrack+1, nMC);
% tauEst_MC3 = zeros(nRadar, nTrack, nMC);

Jd_MC1 = zeros(6,6,nTrack,nMC);
Jd_MC2 = zeros(6,6,nTrack,nMC);
Jd_MC3 = zeros(6,6,nTrack,nMC);
parfor iMC = 1:nMC
    iMC
    %% On Radar Side 
    % add noise to receiced signal
    % === power allocation policy1 === %
    rRadar1 = cell(nTrack, 1);
    for iTrack = 1:nTrack
        % noise = sigmaW*(sqrt(1/2)*randn(nRadar,nSampleRx) + 1i*sqrt(1/2)*randn(nRadar,nSampleRx));
        % rRadar1{iTrack} = rRadar_no_noise1{iTrack} + noise;
    end
    % === power allocation policy1 === %
    rRadar2 = cell(nTrack, 1);
    for iTrack = 1:nTrack
        % noise = sigmaW*(sqrt(1/2)*randn(nRadar,nSampleRx) + 1i*sqrt(1/2)*randn(nRadar,nSampleRx));
        % rRadar2{iTrack} = rRadar_no_noise2{iTrack} + noise;
    end
    
    %% 2-Step Tracking Implementation
    % abondoned
    % generate measurement directly
    ii = (1:1:nSampleRx);
    sampleInd = 1 + mod(floor((ii-1)/(frx/ftx)), nSampleTx);

    tauEst1 = zeros(nRadar, nTrack);
    xEst1 = zeros(6, nTrack+1);
    xEst1(:,1) = x0;
    tauEst2 = zeros(nRadar, nTrack);
    xEst2 = zeros(6, nTrack+1);
    xEst2(:,1) = x0;

    cvfunc = @const_vel_transition_function;
    hfunc = @delay_measurement_function;
    ckf1 = trackingCKF(cvfunc, hfunc, x0, 'ProcessNoise', Q, 'MeasurementNoise', R0, 'StateCovariance', P0);
    ckf2 = trackingCKF(cvfunc, hfunc, x0, 'ProcessNoise', Q, 'MeasurementNoise', R0, 'StateCovariance', P0);
    
    for iTrack = 1:nTrack % loop for tracking time step
        % waitbar(iTrack/nTrack, hwb, ['Process tracking on the radar side: ', int2str(iTrack), ' of ', int2str(nTrack), ' completed']);

        %% Tracking of policy 1 and 2
        % === Observable Estimation === %
        for iRadar = 1:nRadar % loop for estimation of each radar
            tauEst1(iRadar, iTrack) = tau(iRadar, iTrack)  + sigmaTau1(iRadar, iTrack)*randn;
            tauEst2(iRadar, iTrack) = tau(iRadar, iTrack)  + sigmaTau2(iRadar, iTrack)*randn;
        end
    end
    % === State Estimation (CKF) === % 
    xEst1(:, 2:end) = ckfTimeSequence(ckf1, tauEst1, sigmaTau1, T, radarPos);
    xEst2(:, 2:end) = ckfTimeSequence(ckf2, tauEst2, sigmaTau2, T, radarPos);

    % record MC result
    posEst_MC1(:,:,iMC) = xEst1(1:3,:);
    tauEst_MC1(:,:,iMC) = tauEst1;
    posEst_MC2(:,:,iMC) = xEst2(1:3,:);
    tauEst_MC2(:,:,iMC) = tauEst2;

    % === PCRB computation === %
    for iTrack = 1:nTrack
        dH = gradient_H(tarPos(iTrack,:),radarPos);
        R1 = diag(sigmaTau1(:, iTrack)).^2;
        Jd_MC1(:,:,iTrack,iMC) = dH.' * inv(R1) * dH;
        R2 = diag(sigmaTau2(:, iTrack)).^2;
        Jd_MC2(:,:,iTrack,iMC) = dH.' * inv(R2) * dH;
    end
    
    %% Tracking of policy 3
    xEst3 = zeros(6, nTrack+1);
    xEst3(:,1) = x0;
    ckf3 = trackingCKF(cvfunc, hfunc, x0, 'ProcessNoise', Q, 'MeasurementNoise', R0, 'StateCovariance', P0);
    p3 = totalPower / nRadar * ones(nRadar, 1);
    J3 = zeros(6,6,nTrack+1);
    J3(:,:,1) = P0^(-1);
    J3_condi = zeros(6,6,nTrack+1);
    J3_condi(:,:,1) = P0^(-1);
    J3_test = zeros(6,6,nTrack+1);
    J3_test(:,:,1) = P0^(-1);
    for iTrack = 1:nTrack
        % === Observable Estimation === %
        tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        % [~, ~, ~, alpha3] = ...
        %     generateRxSignalRadar(p3, radarInfo,tarInfo,s);
        georange = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
        alpha3 = sqrt(p3.* (1./(georange)).^4);
        sigmaTau3 = sqrt(sigmaW^2./(8*pi^2*fc^2*alpha3.^2));
        tauEst3 = zeros(nRadar, 1);
        for iRadar = 1:nRadar % loop for estimation of each radar
            tauEst3(iRadar) = tau(iRadar, iTrack)  + sigmaTau3(iRadar)*randn;
        end
        % === State Estimation (CKF) === % 
        R3 = diag(sigmaTau3).^2;
        ckf3.MeasurementNoise = R3;
        y = tauEst3;
        [xPred, PPred] = predict(ckf3, T);
        [xCorr, PCorr] = correct(ckf3, y, nRadar, radarPos);
        xEst3(:, iTrack+1) = xCorr;
        % === PCRB computation === %
        dH = gradient_H(tarPos(iTrack,:),radarPos);
        Jd_MC3(:,:,iTrack,iMC) = dH.' * inv(R3) * dH;
        J3_test(:,:,iTrack+1) = (Q + F*J3_test(:,:,iTrack)^(-1)*F')^(-1) + Jd_MC3(:,:,iTrack,iMC);
        % === CPCRB computation === %
        dH_condi = gradient_H(xEst3(1:3, iTrack+1)',radarPos);
        Jd = dH_condi.' * inv(R3) * dH_condi;
        J3_condi(:,:,iTrack+1) = (Q + F*J3_condi(:,:,iTrack)^(-1)*F')^(-1) + Jd;
        % === power allocation optimization (policy 3) === %
        p = optimvar('p', nRadar, 1,...
            'LowerBound', 0, 'UpperBound', totalPower);
        objFun = @objFunPCRB;
        Jp = (Q + F*J3_test(:,:,iTrack)^(-1)*F')^(-1);%J3_condi(:,:,iTrack)
        [cost, pSum] = fcn2optimexpr(objFun, p, ...
            tarPos(iTrack,:), radarPos, sigmaW, fc, Jp, ...%xEst3(1:3, iTrack+1)'
            'OutputSize',[1,1],'ReuseEvaluation',true);
        prob = optimproblem;
        prob.Objective = cost;
        prob.Constraints.totalPower = pSum == totalPower;
        initialP.p = p2(:, iTrack);
        opts = optimoptions("fmincon", Display="notify-detailed", MaxFunctionEvaluations=1e4);
        [sol,fval] = solve(prob,initialP, Options=opts);
        p3 = sol.p;
    end
    % record MC result
    posEst_MC3(:,:,iMC) = xEst3(1:3,:);
    
    %% On Target Side 
    % pEst0 = zeros(nRadar, nTrack);
    % pEst2 = zeros(nRadar, nTrack);
    % pRx0 = zeros(nRadar, nTrack);
    % pRx2 = zeros(nRadar, nTrack);
    % estToTrue2 = zeros(nRadar, nTrack);
    % hwb = waitbar(0,' estimate observables...');
    % for iTrack = 1:nTrack
    %     %% DOA estimation -- MUSIC (Target side)
    %     waitbar(iTrack/nTrack, hwb, ['Power estimation tracking on the target side: ', int2str(iTrack), ' of ', int2str(nTrack), ' completed']);
    %     % clarify the DOA
    %     tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
    %     [az, el] = generateDOA(radarPos, tarPos(iTrack,:));
    %     doa = [az, el]';
    %     % create the sensor array
    %     array = phased.URA('Size',[11 11],'ElementSpacing',[lambda/2 lambda/2],'ArrayNormal','z');
    %     % array.Element.FrequencyRange = [50.0e6 500.0e6];
    %     % generate the sensor array recieved signal
    %     % policy 0
    %     for iRadar = 1:nRadar
    %         s0(iRadar,:) = sqrt(p0(iRadar, iTrack))*code0(iRadar, :);
    %     end
    %     [r0Target, ~, ~, alpha] = generateRxSignalTarget(radarInfo,tarInfo,s0);
    %     pRx0(:,iTrack) = alpha.*p0(:,iTrack);
    %     % xArray = collectPlaneWave(array,r2Target.',doa,fc);
    %     [xArray0, sv] = collectPlaneWave(array,r0Target.',doa,fc);
    %     sigma = 100;
    %     noise = sigma*(randn(size(xArray0))+1i*randn(size(xArray0)));
    %     xArrayNoisy0 = xArray0 + noise;
    % 
    %     % implement 2D MUSIC
    %     estimator = phased.MUSICEstimator2D('SensorArray',array,...
    %         'OperatingFrequency',fc,...
    %         'NumSignalsSource','Property',...
    %         'DOAOutputPort',true,'NumSignals',nRadar,...
    %         'AzimuthScanAngles',-180:.5:180,...
    %         'ElevationScanAngles',-90:.5:0);
    %     [~,doas0] = estimator(xArrayNoisy0);
    %     % align with true doa
    %     for m = 1:nRadar
    %         for n = 1:nRadar
    %             if abs(doa(1, m) - doas0(1, n))<1 &&...
    %                     abs(doa(2, m) - doas0(2, n))<1
    %                 estToTrue0(m, iTrack) = n;
    %             end
    %         end
    %     end
    % 
    %     % policy 2
    %     for iRadar = 1:nRadar
    %         s2(iRadar,:) = sqrt(p2(iRadar, iTrack))*code0(iRadar, :);
    %     end
    %     [r2Target, ~, ~, alpha] = generateRxSignalTarget(radarInfo,tarInfo,s2);
    %     pRx2(:,iTrack) = alpha.*p2(:,iTrack);
    %     % xArray = collectPlaneWave(array,r2Target.',doa,fc);
    %     [xArray2, sv] = collectPlaneWave(array,r2Target.',doa,fc);
    %     sigma = 1;
    %     noise = sigma*(randn(size(xArray2))+1i*randn(size(xArray2)));
    %     xArrayNoisy2 = xArray2 + noise;
    % 
    %     % implement 2D MUSIC
    %     estimator = phased.MUSICEstimator2D('SensorArray',array,...
    %         'OperatingFrequency',fc,...
    %         'NumSignalsSource','Property',...
    %         'DOAOutputPort',true,'NumSignals',nRadar,...
    %         'AzimuthScanAngles',-180:.5:180,...
    %         'ElevationScanAngles',-90:.5:0);
    %     [~,doas2] = estimator(xArrayNoisy2);
    %     % align with true doa
    %     for m = 1:nRadar
    %         for n = 1:nRadar
    %             if abs(doa(1, m) - doas0(1, n))<1 &&...
    %                     abs(doa(2, m) - doas0(2, n))<1
    %                 estToTrue2(m, iTrack) = n;
    %             end
    %         end
    %     end
    % 
    %     % plot pseudo-spectrum
    %     % figure(6)
    %     % plotSpectrum(estimator);
    %     % hold on
    %     % for i = 1:length(doa(1,:))
    %     %     plot3([doa(1,i),doa(1,i)], [doa(2,i),doa(2,i)], [-22,10],...
    %     %         'LineWidth',1.5,'Color',"#A2142F")
    %     % end
    %     % hold off
    %     % exportgraphics(gcf,'figs/spectrum_music.gif','Append',true);
    %     % 
    %     % figure(7);% visualize the trajectory
    %     % plot3(radarPos(:,1),radarPos(:,2),radarPos(:,3),...
    %     %     'o','Color','b','MarkerSize',10,...
    %     %     'MarkerFaceColor','#D9FFFF')
    %     % grid on
    %     % ylim([-1500,1500])
    %     % xlim([-1500,1500])
    %     % xlabel('x')
    %     % ylabel('y')
    %     % hold on
    %     % scatter3(tarPos(iTrack,1), tarPos(iTrack,2), tarPos(iTrack,3),'*','LineWidth',1)
    %     % hold off
    %     % legend("Radar", "Target")
    %     % exportgraphics(gcf,'figs/target_pos.gif','Append',true);
    % 
    %     %% Beamforming (Target side)
    %     % policy 0
    %     X0 = xArray0.';
    %     XNoisy0 = xArrayNoisy0.';
    % 
    %     nAnt0 = size(xArray0,2);
    %     % Rxx = noise'*noise/nSampleRx;
    %     w0 = zeros(nAnt0, nRadar);
    %     yBeam0 = zeros(nRadar, nSampleRx);
    % 
    %     Hsv0 = zeros(nAnt0, nRadar);
    %     steeringvec0 = phased.SteeringVector('SensorArray', array);
    %     doasAlign0 = doas0(:,estToTrue0(:, iTrack));
    %     for iRadar = 1: nRadar
    %         % Hsv0(:, iRadar) = steeringvec0(fc, doa(:, iRadar));
    %         Hsv0(:, iRadar) = steeringvec0(fc, doa(:, iRadar));
    %     end
    %     % theoretical correlation matrix
    %     % Rxx = Hsv*(r2Target*r2Target')*Hsv'/nSampleRx + sigma*eye(size(Hsv, 1)) + 1e5*eye(size(Hsv, 1));
    %     % Rxx = X*X'/nSampleRx + sigma*eye(size(Hsv, 1)) + 1e5*eye(size(Hsv, 1));
    %     % estimate correlation matrix
    %     Rxx0 = XNoisy0*XNoisy0'/nSampleRx + 1e5*eye(size(Hsv0, 1));
    % 
    %     for iRadar = 1: nRadar
    %         % identify constraints
    %         % P = zeros(1, nRadar + 4);
    %         % P(iRadar) = 1;
    %         % P(nRadar+1:end) = 1;
    %         P = zeros(1, nRadar);
    %         P(iRadar) = 1;
    %         % compute steering vector candidate
    %         svCan = zeros(nAnt0, 4);
    %         angD = [0.05, 0.05, -0.05, -0.05;
    %                 0.05, -0.05, 0.05, -0.05]*2;
    %         for k = 1: 4
    %             svCan(:, k) = steeringvec0(fc, doas2(:, iRadar) + angD(:, k));
    %         end
    %         % Hbeam = [Hsv, svCan];
    %         Hbeam = Hsv0;
    %         % compute the beam weights
    %         w0(:, iRadar) = Rxx0^(-1)*Hbeam*(Hbeam'*Rxx0^(-1)*Hbeam)^(-1)*P';
    %         % w0(:, iRadar) = Hbeam*(Hbeam'*Hbeam)^(-1)*P';
    %         % beamer output
    %         yBeam0(iRadar, :) = w0(:, iRadar)' * X0;
    %     end
    %     pEst0(:, iTrack) = sum(abs(yBeam0).^2, 2)/nSampleRx;
    % 
    %     % policy 2
    %     X2 = xArray2.';
    %     XNoisy2 = xArrayNoisy2.';
    % 
    %     nAnt2 = size(xArray2,2);
    %     % Rxx = noise'*noise/nSampleRx;
    %     w2 = zeros(nAnt2, nRadar);
    %     yBeam2 = zeros(nRadar, nSampleRx);
    % 
    %     Hsv2 = zeros(nAnt2, nRadar);
    %     steeringvec2 = phased.SteeringVector('SensorArray', array);
    %     doasAlign2 = doas0(:,estToTrue2(:, iTrack));
    %     for iRadar = 1: nRadar
    %          % Hsv2(:, iRadar) = steeringvec2(fc, doa(:, iRadar));
    %         Hsv2(:, iRadar) = steeringvec2(fc, doa(:, iRadar));
    %     end
    %     % theoretical correlation matrix
    %     % Rxx = Hsv*(r2Target*r2Target')*Hsv'/nSampleRx + sigma*eye(size(Hsv, 1)) + 1e5*eye(size(Hsv, 1));
    %     % Rxx = X*X'/nSampleRx + sigma*eye(size(Hsv, 1)) + 1e5*eye(size(Hsv, 1));
    %     % estimate correlation matrix
    %     Rxx2 = XNoisy2*XNoisy2'/nSampleRx + 1e5*eye(size(Hsv2, 1));
    % 
    %     for iRadar = 1: nRadar
    %         % identify constraints
    %         % P = zeros(1, nRadar + 4);
    %         % P(iRadar) = 1;
    %         % P(nRadar+1:end) = 1;
    %         P = zeros(1, nRadar);
    %         P(iRadar) = 1;
    %         % compute steering vector candidate
    %         svCan = zeros(nAnt2, 4);
    %         angD = [0.05, 0.05, -0.05, -0.05;
    %                 0.05, -0.05, 0.05, -0.05]*2;
    %         for k = 1: 4
    %             svCan(:, k) = steeringvec2(fc, doas0(:, iRadar) + angD(:, k));
    %         end
    %         % Hbeam = [Hsv, svCan];
    %         Hbeam = Hsv2;
    %         % compute the beam weights
    %         w2(:, iRadar) = Rxx2^(-1)*Hbeam*(Hbeam'*Rxx2^(-1)*Hbeam)^(-1)*P';
    %         % w2(:, iRadar) = Hbeam*(Hbeam'*Hbeam)^(-1)*P';
    %         % beamer output
    %         yBeam2(iRadar, :) = w2(:, iRadar)' * X2;
    %     end
    %     pEst2(:, iTrack) = sum(abs(yBeam2).^2, 2)/nSampleRx;
    % % 
    % %     % h8 = figure(8); % plot beamformer output
    % %     % t = (1/frx:1/frx:nSampleRx/frx);
    % %     % subplot(1,2,1);
    % %     % plot(t,real(yBeam0(1 ,:))); axis tight
    % %     % title('Output of MVDR Beamformer for Antenna Array');
    % %     % xlabel('Time (s)');ylabel('Real(y)');
    % %     % subplot(1,2,2);
    % %     % plot(t,real(r0Target(1 ,:))); axis tight
    % %     % title('Original Transmitted Signal to Target');
    % %     % xlabel('Time (s)');ylabel('Real(y)');
    % %     % 
    % %     % h9 = figure(9); % Plot array response with beamformer weighting
    % %     % azPlot = (-180: 0.5: 180);
    % %     % elPlot = (-90: 0.5: 0);
    % %     % powerPlot = zeros(length(azPlot), length(elPlot));
    % %     % for i = 1: length(azPlot)
    % %     %     for j = 1: length(elPlot)
    % %     %         svPlot = steeringvec0(fc, [azPlot(i); elPlot(j)]);
    % %     %         powerPlot(i, j) = abs(w0(:,1)' * svPlot)^2; 
    % %     %     end
    % %     % end
    % %     % [elPlot, azPlot] = meshgrid(elPlot, azPlot);
    % %     % mesh(elPlot, azPlot, -10*log10(powerPlot))
    % %     % hold on 
    % %     % for i = 1:length(doa(1,:))
    % %     %     plot3([doa(2,i),doa(2,i)], [doa(1,i),doa(1,i)], [100,-5],...
    % %     %         'LineWidth',2,'Color',"#A2142F")
    % %     % end
    % %     % title('Beamformer Negative Attenuation at DOAs');
    % %     % xlabel('Elevation');ylabel('Azimuth');zlabel('Power(dB)')
    % %     % legend('Negative Attenuation', 'True DOAs')
    % % 
    % end
    % close(hwb)
    % 
    % %% KL Divergence between power and distance
    % % estimate marginal pdf
    % pts1 = 0:0.005:1;
    % pts2 = 0:0.005:1;
    % p0EstNormalized = pEst0 .* (rho.^2);
    % % p0EstNormalized = p0;
    % p0EstNormalized = p0EstNormalized / max(p0EstNormalized(:));%./ repmat(max(p0EstNormalized, [], 2), 1, nTrack);
    % p2EstNormalized = pEst2 .* (rho.^2);
    % % p2EstNormalized = p2;
    % p2EstNormalized = p2EstNormalized / max(p2EstNormalized(:));%./ repmat(max(p2EstNormalized, [], 2), 1, nTrack);
    % rhoNormalized = rho.^4;
    % rhoNormalized = rhoNormalized / max(rhoNormalized(:));%./ repmat(max(rho, [], 2), 1, nTrack);
    % I_p0_rho = zeros(nRadar,1);
    % I_p2_rho = zeros(nRadar,1);
    % for iRadar = 1:nRadar
    %     I_p0_rho(iRadar) = mutualInfo(p0EstNormalized(iRadar,:), rhoNormalized(iRadar,:), pts1, pts2);
    %     I_p2_rho(iRadar) = mutualInfo(p2EstNormalized(iRadar,:), rhoNormalized(iRadar,:), pts1, pts2);
    % end
    % 
    % % ===== Mutual Information ===== %
    % h13 = figure(13);
    % barPlot = [I_p0_rho, I_p2_rho];
    % bar(barPlot)
    % legend("Non-Cognitive", "Cognitive")
    % xlabel("Radar")
    % ylabel("Mutual Information")
    % %% KL Divergence between power and absolute coordinate
    % % estimate marginal pdf
    % % pts1 = -2:0.01:4;
    % % pts2 = -2:0.01:4;
    % % p0EstNormalized = pEst0 .* (rho.^2);
    % % p0Interval = max(p0EstNormalized, [], 2) - min(p0EstNormalized, [], 2);
    % % p0EstNormalized = (p0EstNormalized - repmat(min(p0EstNormalized, [], 2), 1, nTrack)) ./ repmat(p0Interval, 1, nTrack);
    % % p2EstNormalized = pEst2 .* (rho.^2);
    % % p2Interval = max(p2EstNormalized, [], 2) - min(p2EstNormalized, [], 2);
    % % p2EstNormalized = (p2EstNormalized - repmat(min(p2EstNormalized, [], 2), 1, nTrack)) ./ repmat(p2Interval, 1, nTrack);
    % % tarPosXNormalized = tarPosX' / max(tarPosX);
    % % tarPosYNormalized = tarPosY' / max(tarPosY);
    % % tarPosXY = tarPosX.^2 + tarPosY.^2;
    % % tarPosXYNormalized = tarPosXY' / max(tarPosXY);
    % % 
    % % I_p0_x = zeros(nRadar,1);
    % % I_p2_x = zeros(nRadar,1);
    % % I_p0_y = zeros(nRadar,1);
    % % I_p2_y = zeros(nRadar,1);
    % % I_p0_xy = zeros(nRadar,1);
    % % I_p2_xy = zeros(nRadar,1);
    % % for iRadar = 1:nRadar
    % %     I_p0_x(iRadar) = mutualInfo(p0EstNormalized(iRadar,:), tarPosXNormalized, pts1, pts2);
    % %     I_p2_x(iRadar) = mutualInfo(p2EstNormalized(iRadar,:), tarPosXNormalized, pts1, pts2);
    % %     I_p0_y(iRadar) = mutualInfo(p0EstNormalized(iRadar,:), tarPosYNormalized, pts1, pts2);
    % %     I_p2_y(iRadar) = mutualInfo(p2EstNormalized(iRadar,:), tarPosYNormalized, pts1, pts2);
    % %     I_p0_xy(iRadar) = mutualInfo(p0EstNormalized(iRadar,:), tarPosXYNormalized, pts1, pts2);
    % %     I_p2_xy(iRadar) = mutualInfo(p2EstNormalized(iRadar,:), tarPosXYNormalized, pts1, pts2);
    % % end
end
% figure(1), hold on % plot trajectory and tracking result
% plot3(xEst(1,:),xEst(2,:),xEst(3,:), "-*")
%% Data Analysis
xEst_30_3 = load('figs/MC100/posEst30_policy3_3.mat');
Jd_30_3 = load('figs/MC100/Jd30_policy3_3.mat');
posEst_MC3 = xEst_30_3.posEst_MC3;
Jd_MC3 = Jd_30_3.Jd_MC3;
nMC3 = 30;
% ====== trajectory ====== %
meanPosEst1 = mean(posEst_MC1, 3);
meanPosEst2 = mean(posEst_MC2, 3);
meanPosEst3 = mean(posEst_MC3, 3);
h3 = figure(3);% visualize the trajectory
plot3(radarPos(:,1),radarPos(:,2),radarPos(:,3),...
    'o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF')
radar_number = {'1','2','3','4','5','6'};
text(radarPos(:,1),radarPos(:,2),radarPos(:,3),...
    radar_number,FontSize=16)
hold on
scatter3(tarPos(:,1), tarPos(:,2), tarPos(:,3))
grid on
% plot trajectory and tracking result
plot3(meanPosEst1(1,:),meanPosEst1(2,:),meanPosEst1(3,:), "-*")
plot3(meanPosEst2(1,:),meanPosEst2(2,:),meanPosEst2(3,:), "-*")
plot3(meanPosEst3(1,:),meanPosEst3(2,:),meanPosEst3(3,:), "-*")
zlim([-100, 1.2*height]) 
legend('radar','target', 'trajectory (policy1)',...
    'trajectory (policy2)', 'trajectory (policy3)')
title("Trajectory")
% saveas(h3,'figs/trajectory.fig')

% ====== RMSE and Bound vs Time ====== % 
RMSEtime1 = zeros(1,nTrack+1);
RMSEtime1(1) = sqrt(trace(P0(1:3,1:3)));
for iTrack = 1:nTrack
    squareError = 0;
    for iMC = 1:nMC
        error = posEst_MC1(:,iTrack+1,iMC) - tarPos(iTrack,:)';
        squareError = squareError + error'*error;
    end
    RMSEtime1(iTrack) = sqrt(squareError/nMC); 
end

RMSEtime2 = zeros(1,nTrack+1);
RMSEtime2(1) = sqrt(trace(P0(1:3,1:3)));
for iTrack = 1:nTrack
    squareError = 0;
    for iMC = 1:nMC
        error = posEst_MC2(:,iTrack+1,iMC) - tarPos(iTrack,:)';
        squareError = squareError + error'*error;
    end
    RMSEtime2(iTrack) = sqrt(squareError/nMC); 
end

RMSEtime3 = zeros(1,nTrack+1);
RMSEtime3(1) = sqrt(trace(P0(1:3,1:3)));
for iTrack = 1:nTrack
    squareError = 0;
    for iMC = 1:nMC3
        error = posEst_MC3(:,iTrack+1,iMC) - tarPos(iTrack,:)';
        squareError = squareError + error'*error;
    end
    RMSEtime3(iTrack) = sqrt(squareError/nMC3); 
end

mean_Jd1 = mean(Jd_MC1,4);
J1 = zeros(6,6,nTrack+1);
J1(:,:,1) = P0^(-1);
bound1 = zeros(6,6,nTrack);
bound1(:,:,1) = P0;
bound_pos1 = zeros(nTrack+1,1);
for iTrack = 1: nTrack
    J1(:,:,iTrack+1) = (Q + F*J1(:,:,iTrack)^(-1)*F')^(-1) + ...
        mean_Jd1(:,:,iTrack);
    bound1(:,:,iTrack+1) = J1(:,:,iTrack+1)^(-1);
end
for iTrack = 1: nTrack+1
    bound_pos1(iTrack) = trace(bound1(1:3,1:3,iTrack));
end

mean_Jd2 = mean(Jd_MC2,4);
J2 = zeros(6,6,nTrack+1);
J2(:,:,1) = P0^(-1);
bound2 = zeros(6,6,nTrack);
bound2(:,:,1) = P0;
bound_pos2 = zeros(nTrack+1,1);
for iTrack = 1: nTrack
    J2(:,:,iTrack+1) = (Q + F*J2(:,:,iTrack)^(-1)*F')^(-1) + ...
        mean_Jd2(:,:,iTrack);
    bound2(:,:,iTrack+1) = J2(:,:,iTrack+1)^(-1);
end
for iTrack = 1: nTrack+1
    bound_pos2(iTrack) = trace(bound2(1:3,1:3,iTrack));
end

mean_Jd3 = mean(Jd_MC3,4);
J3 = zeros(6,6,nTrack+1);
J3(:,:,1) = P0^(-1);
bound3 = zeros(6,6,nTrack);
bound3(:,:,1) = P0;
bound_pos3 = zeros(nTrack+1,1);
for iTrack = 1: nTrack
    J3(:,:,iTrack+1) = (Q + F*J3(:,:,iTrack)^(-1)*F')^(-1) + ...
        mean_Jd3(:,:,iTrack);
    bound3(:,:,iTrack+1) = J3(:,:,iTrack+1)^(-1);
end
for iTrack = 1: nTrack+1
    bound_pos3(iTrack) = trace(bound3(1:3,1:3,iTrack));
end

% bound3_test = zeros(6,6,nTrack);
% bound3_test(:,:,1) = P0;
% bound_pos3_test = zeros(nTrack+1,1);
% for iTrack = 1: nTrack
%     bound3_test(:,:,iTrack+1) = J3_test(:,:,iTrack+1)^(-1);
% end
% for iTrack = 1: nTrack+1
%     bound_pos3_test(iTrack) = trace(bound3_test(1:3,1:3,iTrack));
% end

h4 = figure(4);
semilogy(1:nTrack+1, RMSEtime1, "-*", 'LineWidth', 1, 'Color', "#d10808")
hold on
semilogy(1:nTrack+1, RMSEtime2, "-*", 'LineWidth', 1, 'Color', '#4cb7ed')
semilogy(1:nTrack+1, RMSEtime3, "-*", 'LineWidth', 1, 'Color', '#00ff00')
semilogy(1:nTrack+1, sqrt(bound_pos1), ':', 'LineWidth', 3, 'Color', '#ff6929')
semilogy(1:nTrack+1, sqrt(bound_pos2), ':', 'LineWidth', 3, 'Color', '#0072bd')
semilogy(1:nTrack+1, sqrt(bound_pos3), ':', 'LineWidth', 3, 'Color', '#169c02')
% semilogy(1:nTrack+1, sqrt(bound_pos3_test), 'LineWidth', 3)
hold on 
title("RMSE vs Time")
legend("RMSE (policy1)", "RMSE (policy2)", "RMSE (policy3)",...
    "Bound (policy1)", "Bound (policy2)", "Bound (policy3)")
% saveas(h4,'figs/RMSEvsTime.fig')

% ====== RMSE CDF ======%
h5 = figure(5);
cdfplot(RMSEtime1)
hold on
cdfplot(RMSEtime2)
cdfplot(RMSEtime3)
legend("RMSE (policy1)", "RMSE (policy2)", "RMSE (policy3)")
title("Emperical CDF")
% saveas(h5,'figs/cdfplot.fig')

% % ===== Power Estimation ===== %
% h10 = figure(10);
% subplot(1, 2, 1)
% pEstAlign0 = zeros(size(pEst0));
% for iTrack = 1:nTrack
%     tmp = pEst0(:,iTrack);
%     pEstAlign0(:,iTrack) = tmp(estToTrue0(:,iTrack));
% end
% plot((1:1:nTrack), pRx0(1,:), "LineWidth", 5), hold on % true received power
% plot((1:1:nTrack), pEst0(1,:), "-o", "LineWidth", 2, "MarkerSize", 3) % estimation
% legend("received power", "power estimation")
% 
% subplot(1, 2, 2)
% pEstAlign = zeros(size(pEst2));
% for iTrack = 1:nTrack
%     tmp = pEst2(:,iTrack);
%     pEstAlign(:,iTrack) = tmp(estToTrue2(:,iTrack));
% end
% plot((1:1:nTrack), pRx2(1,:), "LineWidth", 5), hold on % true received power
% % plot((1:1:nTrack), pEstAlign(1,:), "-o", "LineWidth", 2, "MarkerSize", 3) % estimation
% % plot((1:1:nTrack), p2(1,:), "-o", "LineWidth", 2, "MarkerSize", 3) % transmitted power
% plot((1:1:nTrack), pEst2(1,:), "-o", "LineWidth", 2, "MarkerSize", 3) % estimation
% legend("received power", "power estimation")
% % saveas(h10,'figs/power_estimation.fig')
% 
% % ===== Data Correlation ===== %
% h11 = figure(11); % policy 0
% for iRadar = 1:nRadar
%     subplot(6, 1, iRadar)
%     yyaxis left
%     % pTx0 = pEst0 .* (rho.^2);
%     plot((1:1:nTrack), p0EstNormalized(iRadar,:), "LineWidth", 2)
%     ylim([0, 1])
%     xlabel('time instance')
%     ylabel('transimitted power estimation')
%     yyaxis right
%     plot((1:1:nTrack), rho(iRadar,:), "LineWidth", 2)
%     ylabel('true range')
% end
% % saveas(h11,'figs/data_correlation_policy0.fig')
% pm0 = mean(p0EstNormalized, 2);%mean(pEst0, 2);
% rm = mean(rhoNormalized, 2);%mean(rho, 2);
% corrCoe0 = zeros(nRadar, 1);
% for iRadar = 1:nRadar
%     corrCoe0(iRadar) = sum((p0EstNormalized(iRadar, :) - pm0(iRadar)).*(rhoNormalized(iRadar, :) - rm(iRadar))) / ...
%         sqrt(sum((p0EstNormalized(iRadar, :) - pm0(iRadar)).^2) * sum((rhoNormalized(iRadar, :) - rm(iRadar)).^2));
% end
% 
% h12 = figure(12); % policy 2
% for iRadar = 1:nRadar
%     subplot(6, 1, iRadar)
%     yyaxis left
%     % pTx2 = pEst2 .* (rho.^2);
%     plot((1:1:nTrack), p2EstNormalized(iRadar,:), "LineWidth", 2)
%     ylim([0, 1])
%     xlabel('time instance')
%     ylabel('transimitted power estimation')
%     yyaxis right
%     plot((1:1:nTrack), rho(iRadar,:), "LineWidth", 2)
%     ylabel('true range')
% end
% % saveas(h12,'figs/data_correlation_policy2.fig')
% pm2 = mean(p2EstNormalized, 2);%mean(pEst2, 2);
% rm = mean(rhoNormalized, 2);%mean(rho, 2);
% corrCoe2 = zeros(nRadar, 1);
% for iRadar = 1:nRadar
%     corrCoe2(iRadar) = sum((p2EstNormalized(iRadar, :) - pm2(iRadar)).*(rhoNormalized(iRadar, :) - rm(iRadar))) / ...
%         sqrt(sum((p2EstNormalized(iRadar, :) - pm2(iRadar)).^2) * sum((rhoNormalized(iRadar, :) - rm(iRadar)).^2));
% end
% 
% % pm0 = mean(pEst0, 2);
% % rm = mean(rho, 2);
% % corrCoe0 = zeros(nRadar, 1);
% % for iRadar = 1:nRadar
% %     corrCoe0(iRadar) = sum((pEst0(iRadar, :) - pm0(iRadar)).*(rho(iRadar, :) - rm(iRadar))) / ...
% %         sqrt(sum((pEst0(iRadar, :) - pm0(iRadar)).^2) * sum((rho(iRadar, :) - rm(iRadar)).^2));
% % end
% % ===== Mutual Information ===== %
% h13 = figure(13);
% barPlot = [I_p0_rho, I_p2_rho];
% bar(barPlot)
% legend("Non-Cognitive", "Cognitive")
% xlabel("Radar")
% ylabel("Mutual Information")
% % saveas(h13,'figs/MI_comparison.fig')


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

