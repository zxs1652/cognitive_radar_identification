% implement radar scenario simulation
% the radar is tracking a self-defined trajectory in 3 dimensions
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


%% Monte-Carlo Experiments for collect data
nMC = 1;
posEst_MC3 = zeros(3, nTrack+1, nMC);
% tauEst_MC3 = zeros(nRadar, nTrack, nMC);
Jd_MC3 = zeros(6,6,nTrack,nMC);
% memory for test data
prx = cell(nMC, 1);
rho = cell(nMC, 1);
target_pos = cell(nMC, 1);
for iMC = 1:nMC
    iMC
    %% Signal Generation
    % === generate trajectory === %
    [tarPos, tarVel] = generateTrajectory(pos0, nTrack, tr, T);
    % % ====== trajectory ====== %
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
    p2 = zeros(nRadar, nTrack);
    for iTrack = 1:nTrack
        tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        geoRange(:,iTrack) = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
        tau(:,iTrack) = 2*geoRange(:,iTrack)./c;
        prop2 = (1./geoRange(:,iTrack).^4)./sum(1./geoRange(:,iTrack).^4);
        p2(:, iTrack) = totalPower*prop2;
    end
    %% On Radar Side (Tracking of policy 3)
    xEst3 = zeros(6, nTrack+1);
    pos0 = tarPos(1, :)' + 10*randn(3,1);
    vel0 = tarVel(1, :)' + 2*randn(3,1);
    x0 = [pos0;vel0];
    xEst3(:,1) = x0;
    cvfunc = @const_vel_transition_function;
    hfunc = @delay_measurement_function;
    ckf3 = trackingCKF(cvfunc, hfunc, x0, 'ProcessNoise', Q, 'MeasurementNoise', R0, 'StateCovariance', P0);
    p3 = zeros(nRadar, nTrack);
    p3(:,1) = totalPower / nRadar * ones(nRadar, 1);
    % J3 = zeros(6,6,nTrack+1);
    % J3(:,:,1) = P0^(-1);
    J3_condi = zeros(6,6,nTrack+1);
    J3_condi(:,:,1) = P0^(-1);
    J3_test = zeros(6,6,nTrack+1);
    J3_test(:,:,1) = P0^(-1);

    distance = zeros(nRadar, nTrack);
    % rRadar3 = cell(nTrack, 1);
    for iTrack = 1:nTrack
        % === Observable Estimation === %
        % generate signal
        tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        % [rRadar_no_noise3, ~, ~, alpha3] = ...
        %     generateRxSignalRadar(p3, radarInfo,tarInfo,code0);
        % noise = sigmaW*(sqrt(1/2)*randn(nRadar,nSampleRx) + 1i*sqrt(1/2)*randn(nRadar,nSampleRx));
        % rRadar3{iTrack} = rRadar_no_noise3 + noise;
        georange = vecnorm(repmat(tarPos(iTrack,:),nRadar,1) - radarPos, 2, 2);
        distance(:,iTrack) = georange;
        alpha3 = sqrt(p3(:,iTrack).* (1./(georange)).^4);

        % generate noisy measurement
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
        p3(:,iTrack+1) = sol.p;
    end
    % record MC result
    % posEst_MC3(:,:,iMC) = xEst3(1:3,:);
    %% On Target Side 
    % rTarget3 = cell(nTrack, 1);
    pRx3 = zeros(nRadar, nTrack);
    pRxEst3 = zeros(nRadar, nTrack);
    for iTrack = 1:nTrack
        % clarify the DOA
        tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
        [az, el] = generateDOA(radarPos, tarPos(iTrack,:));
        doa = [az, el]';
        % generate received signal
        [rTarget3, ~, ~, alphaTarget3] = ...
            generateRxSignalTarget(p3(:,iTrack), radarInfo,tarInfo,code0);
        pRx3(:,iTrack) = alphaTarget3.^2;
        % create the sensor array
        array = phased.URA('Size',[11 11],'ElementSpacing',[lambda/2 lambda/2],'ArrayNormal','z');
        % sterring vector
        steeringvec3 = phased.SteeringVector('SensorArray', array);
        % collected signal
        [xArray3, sv] = collectPlaneWave(array,rTarget3.',doa,fc);
        X3 = xArray3.';
        %% Beamforming 
        nAnt3 = size(xArray3,2);
        w3 = zeros(nAnt3, nRadar);
        yBeam3 = zeros(nRadar, nSampleRx);
        Hsv3 = zeros(nAnt3, nRadar);
        for iRadar = 1: nRadar
            Hsv3(:, iRadar) = steeringvec3(fc, doa(:, iRadar));
        end
        % correlation matrix
        Rxx3 = X3*X3'/nSampleRx + 1e5*eye(size(Hsv3, 1));
        for iRadar = 1: nRadar
            % identify constraints
            P = zeros(1, nRadar);
            P(iRadar) = 1;
            % beamformer
            Hbeam = Hsv3;
            % compute the beam weights
            w3(:, iRadar) = Rxx3^(-1)*Hbeam*(Hbeam'*Rxx3^(-1)*Hbeam)^(-1)*P';
            % beamer output
            yBeam3(iRadar, :) = w3(:, iRadar)' * X3;
        end
        pRxEst3(:, iTrack) = sum(abs(yBeam3).^2, 2)/nSampleRx;


    end
    % record test data
    prx{iMC} = pRxEst3;
    rho{iMC} = distance;
    target_pos{iMC} = tarPos;
end
%% save test data
save("dataset/prx.mat", "prx")
save("dataset/rho.mat", "rho")
save("dataset/target_pos.mat", "target_pos")
% figure(1), hold on % plot trajectory and tracking result
% plot3(xEst(1,:),xEst(2,:),xEst(3,:), "-*")
%% Data Analysis
% % ====== trajectory ====== %
% meanPosEst3 = mean(posEst_MC3, 3);
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
% % plot trajectory and tracking result
% plot3(meanPosEst3(1,:),meanPosEst3(2,:),meanPosEst3(3,:), "-*")
% zlim([-100, 1.2*height]) 
% % legend('radar','target', 'trajectory (policy1)',...
% %     'trajectory (policy2)', 'trajectory (policy3)')
% legend('radar','target', 'trajectory (policy3)')
% title("Trajectory")
% % saveas(h3,'figs/trajectory.fig')
% 
% % ====== RMSE and Bound vs Time ====== % 
% RMSEtime3 = zeros(1,nTrack+1);
% RMSEtime3(1) = sqrt(trace(P0(1:3,1:3)));
% for iTrack = 1:nTrack
%     squareError = 0;
%     for iMC = 1:nMC
%         error = posEst_MC3(:,iTrack+1,iMC) - tarPos(iTrack,:)';
%         squareError = squareError + error'*error;
%     end
%     RMSEtime3(iTrack) = sqrt(squareError/nMC); 
% end
% 
% mean_Jd3 = mean(Jd_MC3(:,:,:,1:nMC),4);
% J3 = zeros(6,6,nTrack+1);
% J3(:,:,1) = P0^(-1);
% bound3 = zeros(6,6,nTrack);
% bound3(:,:,1) = P0;
% bound_pos3 = zeros(nTrack+1,1);
% for iTrack = 1: nTrack
%     J3(:,:,iTrack+1) = (Q + F*J3(:,:,iTrack)^(-1)*F')^(-1) + ...
%         mean_Jd3(:,:,iTrack);
%     bound3(:,:,iTrack+1) = J3(:,:,iTrack+1)^(-1);
% end
% for iTrack = 1: nTrack+1
%     bound_pos3(iTrack) = trace(bound3(1:3,1:3,iTrack));
% end
% 
% bound3_test = zeros(6,6,nTrack);
% bound3_test(:,:,1) = P0;
% bound_pos3_test = zeros(nTrack+1,1);
% for iTrack = 1: nTrack
%     bound3_test(:,:,iTrack+1) = J3_test(:,:,iTrack+1)^(-1);
% end
% for iTrack = 1: nTrack+1
%     bound_pos3_test(iTrack) = trace(bound3_test(1:3,1:3,iTrack));
% end
% 
% h4 = figure(4);
% semilogy(1:nTrack+1, RMSEtime3, "-*", 'LineWidth', 1)
% hold on
% semilogy(1:nTrack+1, sqrt(bound_pos3), 'LineWidth', 3)
% semilogy(1:nTrack+1, sqrt(bound_pos3_test), 'LineWidth', 3)
% title("RMSE vs Time")
% legend("RMSE (policy3)","Bound (policy3)", "test")
% % saveas(h4,'figs/RMSEvsTime.fig')
% 
% % ====== RMSE CDF ======%
% h5 = figure(5);
% cdfplot(RMSEtime3)
% legend("RMSE (policy3)")
% title("Emperical CDF")
% % saveas(h5,'figs/cdfplot.fig')
% 
% % ===== Power Estimation ===== %
% h10 = figure(10);
% iPlot = 1;
% plot((1:1:nTrack), pRx3(iPlot,:), "LineWidth", 5), hold on % true received power
% plot((1:1:nTrack), pRxEst3(iPlot,:), "-o", "LineWidth", 2, "MarkerSize", 3) % estimation
% legend("received power", "power estimation")
% % saveas(h10,'figs/power_estimation.fig')
% 


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
% xEst_30_1 = load('posEst30_policy3_1.mat');
% Jd_30_1 = load('Jd30_policy3_1.mat');
% % xEst_30_2 = load('posEst30_policy3_2.mat');
% % Jd_30_2 = load('Jd30_policy3_2.mat');
% xEst_30_3 = load('posEst30_policy3_3.mat');
% Jd_30_3 = load('Jd30_policy3_3.mat');
% % 
% posEst_MC3 = cat(3,xEst_30_1.posEst_MC3, xEst_30_3.posEst_MC3);
% Jd_MC3 = cat(4, Jd_30_1.Jd_MC3, Jd_30_3.Jd_MC3);
% nMC = 60;
% 
% posEst_MC3 = xEst_30_1.posEst_MC3;
% Jd_MC3 = Jd_30_1.Jd_MC3;
% nMC = 30;
% 
% posEst_MC3 = xEst_30_2.posEst_MC3;
% Jd_MC3 = Jd_30_2.Jd_MC3;
% nMC = 30;


