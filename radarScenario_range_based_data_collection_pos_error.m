% implement radar scenario simulation
% the radar is tracking a self-defined trajectory in 3 dimensions
clear
close all
clc

%% Parameter Setting
rng(3)
% ====== radar characteristic ====== %
radarConfig;

% === Radar position uncertainy settings === %
sigma_candidates = (0:5:150);
nCandidates = length(sigma_candidates);
for iCandidate = 1 :nCandidates
    % === Monte Carlo settings === %
    nMC = 1;
    prx1 = cell(nMC, 1);
    prx2 = cell(nMC, 1);
    prx3 = cell(nMC, 1);
    rho = cell(nMC, 1);
    target_pos = cell(nMC, 1);
    radar_pos = cell(nMC, 1);
    sigma_pos_error = sigma_candidates(iCandidate);
    for iMC = 1:nMC
        %% display
        fprintf("Monte Carlo %d of %d...\n", iMC, nMC)
        %% generate random trajectory
        [tarPos, tarVel] = generateTrajectory_cv(hgt, nTrack, T);
        
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
        Q = zeros(6,6);
        sigmaTau0 = 1e-7;
        R0 = sigmaTau0^2*eye(nRadar);
        
        %% compute power of transmitted signal on radar side
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
        for iTrack = 1:nTrack
            % tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
            prop1 = 1 / nRadar * ones(nRadar, 1);
            p1(:, iTrack) = totalPower*prop1;
            % [rRadar_no_noise1{iTrack}, ~, ~, alpha1(:, iTrack)] = ...
            %     generateRxSignalRadar(p1(:, iTrack), radarInfo,tarInfo,code0);
        end
        
        % === power allocation policy2 (proportional to the distance^2) === %
        rRadar_no_noise2 = cell(nTrack, 1); % received signal of radar
        p2 = zeros(nRadar, nTrack);
        alpha2 = zeros(nRadar,nTrack);
        for iTrack = 1:nTrack
            % tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
            prop2 = (1./geoRange(:,iTrack).^2)./sum(1./geoRange(:,iTrack).^2);
            p2(:, iTrack) = totalPower*prop2;
            % [rRadar_no_noise2{iTrack}, ~, ~, alpha2(:, iTrack)] = ...
            %     generateRxSignalRadar(p2(:, iTrack), radarInfo,tarInfo,code0);
        end
    
        % === power allocation policy3 (need tracking implementation) === %
        cvfunc = @const_vel_transition_function;
        hfunc = @delay_measurement_function;
        ckf3 = trackingCKF(cvfunc, hfunc, x0, 'ProcessNoise', Q, 'MeasurementNoise', R0, 'StateCovariance', P0);
        p3 = zeros(nRadar, nTrack);
        p3(:,1) = (1./geoRange(:,1).^4)./sum(1./geoRange(:,1).^4);
        J3 = P0^(-1);
        for iTrack = 1:nTrack
            % === Observable Estimation === %
            % generate signal
            tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
            alpha3 = sqrt(p3(:,iTrack).* (1./(geoRange(:,iTrack))).^4);
    
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
            % === PCRB computation === %
            dH = gradient_H(tarPos(iTrack,:),radarPos);
            Jd3 = dH.' * inv(R3) * dH;
            J3 = (Q + F*J3^(-1)*F')^(-1) + Jd3;
            % === power allocation optimization (policy 3) === %
            p = optimvar('p', nRadar, 1,...
                'LowerBound', minPower, 'UpperBound', totalPower);
            objFun = @objFunPCRB;
            Jp = (Q + F*J3^(-1)*F')^(-1);%J3_condi(:,:,iTrack)
            [cost, pSum] = fcn2optimexpr(objFun, p, ...
                tarPos(iTrack,:), radarPos, sigmaW, fc, Jp, ...%xEst3(1:3, iTrack+1)'
                'OutputSize',[1,1],'ReuseEvaluation',true);
            prob = optimproblem;
            prob.Objective = cost;
            prob.Constraints.totalPower = pSum == totalPower;
            initial_sol = [];
            initial_prop = (1./geoRange(:,iTrack).^4)./sum(1./geoRange(:,iTrack).^4);
            initial_sol.p = totalPower*initial_prop;
            opts = optimoptions("fmincon", Display="off", MaxFunctionEvaluations=1e4);
            [sol,fval] = solve(prob,initial_sol, Options=opts);
            % [sol, fval] = fmincon(prob,initialP, Options=opts);
            p3(:,iTrack+1) = sol.p;
        end
    
    %% Collect data from target side
        %% momery allocation
        pRx1 = zeros(nRadar, nTrack);
        pRxEst1 = zeros(nRadar, nTrack);
        pRx2 = zeros(nRadar, nTrack);
        pRxEst2 = zeros(nRadar, nTrack);
        pRx3 = zeros(nRadar, nTrack);
        pRxEst3 = zeros(nRadar, nTrack);
        %% generate received signal and estimate power on target side
        radarPos_given = radarPos + sigma_pos_error * randn(nRadar, 3);
        for iTrack = 1:nTrack
            %% clarify the DOA
            tarInfo = [tarPos(iTrack,:),tarVel(iTrack,:)];
            [az, el] = generateDOA(radarPos, tarPos(iTrack,:));
            doa = [az, el]';
    
            [az_given, el_given] = generateDOA(radarPos_given, tarPos(iTrack,:));
            doa_given = [az_given, el_given]';
            
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
            xArray1 = collectPlaneWave(array,rTarget1.',doa,fc);
            X1 = xArray1.';
            % beamforming
            nAnt1 = size(xArray1,2);
            w1 = zeros(nAnt1, nRadar);
            yBeam1 = zeros(nRadar, nSampleRx);
            Hsv1 = zeros(nAnt1, nRadar);
            for iRadar = 1: nRadar
                Hsv1(:, iRadar) = steeringvec1(fc, doa_given(:, iRadar));
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
            xArray2 = collectPlaneWave(array,rTarget2.',doa,fc);
            X2 = xArray2.';
            % beamforming
            nAnt2 = size(xArray2,2);
            w2 = zeros(nAnt2, nRadar);
            yBeam2 = zeros(nRadar, nSampleRx);
            Hsv2 = zeros(nAnt2, nRadar);
            for iRadar = 1: nRadar
                Hsv2(:, iRadar) = steeringvec2(fc, doa_given(:, iRadar));
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
    
            %% policy 3
             % generate received signal
            [rTarget3, ~, ~, alphaTarget3] = ...
                generateRxSignalTarget(p3(:,iTrack), radarInfo,tarInfo,code0);
            pRx3(:,iTrack) = alphaTarget3.^2;
            % create the sensor array
            array = phased.URA('Size',[11 11],'ElementSpacing',[lambda/2 lambda/2],'ArrayNormal','z');
            % sterring vector
            steeringvec3 = phased.SteeringVector('SensorArray', array);
            % collected signal
            xArray3 = collectPlaneWave(array,rTarget3.',doa,fc);
            X3 = xArray3.';
            % beamforming
            nAnt3 = size(xArray3,2);
            w3 = zeros(nAnt3, nRadar);
            yBeam3 = zeros(nRadar, nSampleRx);
            Hsv3 = zeros(nAnt3, nRadar);
            for iRadar = 1: nRadar
                Hsv3(:, iRadar) = steeringvec3(fc, doa_given(:, iRadar));
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
        % record data
        prx1{iMC} = pRxEst1;
        prx2{iMC} = pRxEst2;
        prx3{iMC} = pRxEst3;
        rho{iMC} = geoRange;
        target_pos{iMC} = tarPos;
        radar_pos{iMC} = radarPos_given;
    end
    %% save test data
    % save_path_prx1 = sprintf("dataset/uncertain_radar_pos/prx1_%d_rng3_sigma_%d.mat", nMC, sigma_pos_error);
    % save_path_prx2 = sprintf("dataset/uncertain_radar_pos/prx2_%d_rng3_sigma_%d.mat", nMC, sigma_pos_error);
    % save_path_prx3 = sprintf("dataset/uncertain_radar_pos/prx3_%d_rng3_sigma_%d.mat", nMC, sigma_pos_error);
    % save_path_rho = sprintf("dataset/uncertain_radar_pos/rho_%d_rng3_sigma_%d.mat", nMC, sigma_pos_error);
    % save_path_tar_pos = sprintf("dataset/uncertain_radar_pos/tarpos_%d_rng3_sigma_%d.mat", nMC, sigma_pos_error);
    % save_path_radar_pos = sprintf("dataset/uncertain_radar_pos/radar_%d_rng3_sigma_%d.mat", nMC, sigma_pos_error);
    % save(save_path_prx1, "prx1")
    % save(save_path_prx2, "prx2")
    % save(save_path_prx3, "prx3")
    % save(save_path_rho, "rho")
    % save(save_path_tar_pos, "target_pos")
    % save(save_path_radar_pos, "radar_pos")
end

%% Data Analysis
% ====== trajectory ====== %
% meanPosEst1 = mean(posEst_MC1, 3);
% meanPosEst2 = mean(posEst_MC2, 3);
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
% plot3(meanPosEst1(1,:),meanPosEst1(2,:),meanPosEst1(3,:), "-*")
% plot3(meanPosEst2(1,:),meanPosEst2(2,:),meanPosEst2(3,:), "-*")
% plot3(meanPosEst3(1,:),meanPosEst3(2,:),meanPosEst3(3,:), "-*")
% zlim([0.8*height, 1.2*height]) 
% % zlim([-100, 1.2*height]) 
% legend('radar','target', 'trajectory (policy1)',...
%     'trajectory (policy2)', 'trajectory (policy3)')
% title("Trajectory")
% saveas(h3,'figs/trajectory.fig')


% ===== Power Estimation ===== %
h10 = figure(10);
iPlot = 1;
plot((1:1:nTrack), pRx1(iPlot,:), "LineWidth", 5), hold on % true received power
plot((1:1:nTrack), pRxEst1(iPlot,:), "-o", "LineWidth", 2, "MarkerSize", 3) % estimation
legend("received power", "power estimation")
% % saveas(h10,'figs/power_estimation.fig')



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


