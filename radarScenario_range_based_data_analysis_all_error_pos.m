% implement radar scenario simulation
% the radar is tracking a self-defined trajectory in 3 dimensions
clear
close all
clc

%% Parameter Setting
rng(3)
% ====== radar characteristic ====== %
radarConfig;

% ====== method selection ====== %
adtest_enable = 1;
interdependency_enable = 0;

sigma_candidates = (0:1:20);
n_candidates = length(sigma_candidates);
P_fa_mi = zeros(1, n_candidates);
P_md_mi = zeros(1, n_candidates);
P_fa_dep = zeros(1, n_candidates);
P_md_dep = zeros(1, n_candidates);
for i_candidate = 1:n_candidates
%% Load data
% ====== load data ====== %
nDataset = 400;
sigma = sigma_candidates(i_candidate)

load_path_prx1 = sprintf("dataset/uncertain_radar_pos/prx1_%d_rng3_sigma_%d.mat", nDataset, sigma);
load_path_prx2 = sprintf("dataset/uncertain_radar_pos/prx2_%d_rng3_sigma_%d.mat", nDataset, sigma);
load_path_prx3 = sprintf("dataset/uncertain_radar_pos/prx3_%d_rng3_sigma_%d.mat", nDataset, sigma);
load_path_rho = sprintf("dataset/uncertain_radar_pos/rho_%d_rng3_sigma_%d.mat", nDataset, sigma);
load_path_tar_pos = sprintf("dataset/uncertain_radar_pos/tarpos_%d_rng3_sigma_%d.mat", nDataset, sigma);
load_path_radar_pos = sprintf("dataset/uncertain_radar_pos/radar_%d_rng3_sigma_%d.mat", nDataset, sigma);
prx1 = load(load_path_prx1);
prx1 = prx1.prx1;
prx2 = load(load_path_prx2);
prx2 = prx2.prx2;
prx3 = load(load_path_prx3);
prx3 = prx3.prx3;
rho123 = load(load_path_rho);
rho123 = rho123.rho;
radarPos123 = load(load_path_radar_pos);
radarPos123 = radarPos123.radar_pos;
% prx4 = load('dataset/policies_data/prx1_800_rng5.mat');
% prx4 = prx4.prx1;
% prx5 = load('dataset/policies_data/prx2_800_rng5.mat');
% prx5 = prx5.prx2;
% prx6 = load('dataset/policies_data/prx3_800_rng5.mat');
% prx6 = prx6.prx3;
% rho456 = load('dataset/policies_data/rho_800_rng5.mat');
% rho456 = rho456.rho;

count_fa_mi = zeros(1, nDataset);
count_md_mi = zeros(1, nDataset);
count_fa_dep = zeros(1, nDataset);
count_md_dep = zeros(1, nDataset);
for iDataset = 1:nDataset
    % ====== compute mutual information ======%  
    % load test data
    pr1_seq = prx1{iDataset};
    % pr2_seq = prx2{iDataset};
    pr3_seq = prx3{iDataset};
    rho123_seq = rho123{iDataset};
    % pr4_seq = prx4{iDataset};
    % pr5_seq = prx6{iDataset};
    % pr6_seq = prx6{iDataset};
    % rho456_seq = rho456{iDataset};
    
    % load surrogate data
    dist_surr = cell(nRadar, 1);
    surr_ind = randi([1, nDataset]);
    prx1_surr = prx1{surr_ind};
    rho_surr = rho123{surr_ind};
    
    % estimate marginal pdf
    pts1 = 0:0.005:1;
    pts2 = 0:0.005:1;
    % normalize power
    prx1EstNormalized = pr1_seq .* (rho123_seq.^2);
    prx1EstNormalized = prx1EstNormalized / max(prx1EstNormalized(:));
    % prx2EstNormalized = pr2_seq .* (rho123_seq.^2);
    % prx2EstNormalized = prx2EstNormalized / max(prx2EstNormalized(:));
    prx3EstNormalized = pr3_seq .* (rho123_seq.^2);
    prx3EstNormalized = prx3EstNormalized / max(prx3EstNormalized(:));
    % prx4EstNormalized = pr4_seq .* (rho456_seq.^2);
    % prx4EstNormalized = prx4EstNormalized / max(prx4EstNormalized(:));
    % prx5EstNormalized = pr5_seq .* (rho456_seq.^2);
    % prx5EstNormalized = prx5EstNormalized / max(prx5EstNormalized(:));
    % prx6EstNormalized = pr6_seq .* (rho456_seq.^2);
    % prx6EstNormalized = prx6EstNormalized / max(prx6EstNormalized(:));
    rho123Normalized = rho123_seq.^4;
    rho123Total = sum(rho123Normalized, 1);
    rho123Normalized = rho123Normalized ./ repmat(rho123Total, nRadar, 1);
    % rho456Normalized = rho456_seq.^4;
    % rho456Total = sum(rho456Normalized, 1);
    % rho456Normalized = rho456Normalized ./ repmat(rho456Total, nRadar, 1);
    prx1SurrNormalized = prx1_surr .* (rho_surr.^2);
    prx1SurrNormalized = prx1SurrNormalized / max(prx1SurrNormalized(:));
    
    %% Mutual inforamtion -- AD Test
    fprintf("Computing MI of dataset %d of %d \n",...
        iDataset, nDataset);
    mutualInfo1_seq = zeros(nRadar, nTrack-10);
    % mutualInfo2_seq = zeros(nRadar, nTrack-10);
    mutualInfo3_seq = zeros(nRadar, nTrack-10);
    % mutualInfo4_seq = zeros(nRadar, nTrack-10);
    % mutualInfo5_seq = zeros(nRadar, nTrack-10);
    % mutualInfo6_seq = zeros(nRadar, nTrack-10);
    mutualInfo_surr_seq = zeros(nRadar, nTrack-10);
    for iTrack = 1:nTrack - 10
        I_prx1_rho = zeros(nRadar,1);
        % I_prx2_rho = zeros(nRadar,1);
        I_prx3_rho = zeros(nRadar,1);
        % I_prx4_rho = zeros(nRadar,1);
        % I_prx5_rho = zeros(nRadar,1);
        % I_prx6_rho = zeros(nRadar,1);
        I_surr = zeros(nRadar,1);
        parfor iRadar = 1:nRadar
            I_prx1_rho(iRadar) = mutualInfo(prx1EstNormalized(iRadar,1:10 + iTrack),...
                rho123Normalized(iRadar,1:10 + iTrack), pts1, pts2);
            % I_prx2_rho(iRadar) = mutualInfo(prx2EstNormalized(iRadar,1:10 + iTrack),...
            %     rho123Normalized(iRadar,1:10 + iTrack), pts1, pts2);
            I_prx3_rho(iRadar) = mutualInfo(prx3EstNormalized(iRadar,1:10 + iTrack),...
                rho123Normalized(iRadar,1:10 + iTrack), pts1, pts2);
            % I_prx4_rho(iRadar) = mutualInfo(prx4EstNormalized(iRadar,1:10 + iTrack),...
            %     rho456Normalized(iRadar,1:10 + iTrack), pts1, pts2);
            % I_prx5_rho(iRadar) = mutualInfo(prx5EstNormalized(iRadar,1:10 + iTrack),...
            %     rho456Normalized(iRadar,1:10 + iTrack), pts1, pts2);
            % I_prx6_rho(iRadar) = mutualInfo(prx6EstNormalized(iRadar,1:10 + iTrack),...
            %     rho456Normalized(iRadar,1:10 + iTrack), pts1, pts2);
            I_surr(iRadar) = mutualInfo(prx1SurrNormalized(iRadar,1:10 + iTrack),...
                rho123Normalized(iRadar,1:10 + iTrack), pts1, pts2);
        end
        mutualInfo1_seq(:,iTrack) = I_prx1_rho;
        % mutualInfo2_seq(:,iTrack) = I_prx2_rho;
        mutualInfo3_seq(:,iTrack) = I_prx3_rho;
        % mutualInfo4_seq(:,iTrack) = I_prx4_rho;
        % mutualInfo5_seq(:,iTrack) = I_prx5_rho;
        % mutualInfo6_seq(:,iTrack) = I_prx6_rho;
        mutualInfo_surr_seq(:,iTrack) = I_surr;
    end
    
    % ====== use one of non-cognitive cases as surrogate data ======%  
    for iRadar = 1:nRadar
        mu = mean(mutualInfo_surr_seq(iRadar,:));
        sigma = std(mutualInfo_surr_seq(iRadar,:));
        dist_surr{iRadar} = makedist('Normal', 'mu',mu,'sigma',sigma);
    end
    
    % ====== use cognitive case as test data ======%  
    % p_fa2 = zeros(nRadar, nTrack-10);
    % hypo2 = zeros(nRadar, nTrack-10);
    p_value3 = zeros(nRadar, nTrack-10);
    hypo3_mi = zeros(nRadar, nTrack-10);
    % p_fa5 = zeros(nRadar, nTrack-10);
    % hypo5 = zeros(nRadar, nTrack-10);
    % p_fa6 = zeros(nRadar, nTrack-10);
    % hypo6 = zeros(nRadar, nTrack-10);
    for iTrack = 1:nTrack - 10
        for iRadar = 1:nRadar
            % [hypo2(iRadar, iTrack), p_fa2(iRadar, iTrack)] = ...
            %     adtest(I_prx2_rho(iRadar), 'Distribution', dist_surr{iRadar});
            [hypo3_mi(iRadar, iTrack), p_value3(iRadar, iTrack)] = ...
                adtest(mutualInfo3_seq(iRadar,iTrack), 'Distribution', dist_surr{iRadar});
            % [hypo5(iRadar, iTrack), p_fa5(iRadar, iTrack)] = ...
            %     adtest(I_prx5_rho(iRadar), 'Distribution', dist_surr{iRadar});
            % [hypo6(iRadar, iTrack), p_fa6(iRadar, iTrack)] = ...
            %     adtest(I_prx6_rho(iRadar), 'Distribution', dist_surr{iRadar});
        end
    end
    
    % ====== use non-cognitive case as test data ======%  
    p_value1 = zeros(nRadar, nTrack-10);
    hypo1_mi = zeros(nRadar, nTrack-10);
    
    for iTrack = 1:nTrack - 10
        for iRadar = 1:nRadar
            [hypo1_mi(iRadar, iTrack), p_value1(iRadar, iTrack)] = ...
                adtest(mutualInfo1_seq(iRadar,iTrack), 'Distribution', dist_surr{iRadar});
        end
    end
    
    %% Interdependency
    % ===== Initilization ===== %
    k = 10; % k closest neighbors
    S1_k_rho_given_power = zeros(1, nTrack - 10);
    S3_k_rho_given_power = zeros(1, nTrack - 10);
    thre = zeros(1, nTrack - 10);
    for iTrack = 1:nTrack - 10
        % ===== non-cognitive ===== %
        % search neighbiors for distance
        rho_point = rho123Normalized(:, 10 + iTrack).';
        [~, rho_dis] = knnsearch(rho123Normalized(:, 1:10 + iTrack).', rho_point, 'K', k);
        % rho_neighbors = rho123Normalized(:, k_ind).'; 
        R_k_rho = mean(rho_dis.^2);
        % search neighbiors for power
        power_point = prx1EstNormalized(:, 10 + iTrack).';
        [k_ind, ~] = knnsearch(prx1EstNormalized(:, 1:10 + iTrack).', power_point, 'K', k);
        rho_neighbors_given_power = rho123Normalized(:, k_ind).';
        R_k_rho_given_power = mean(vecnorm(rho_neighbors_given_power - repmat(rho_point, k, 1), 2, 2).^2);
        S1_k_rho_given_power(iTrack) = R_k_rho / R_k_rho_given_power;
    
        % ===== cognitive ===== %
        % search neighbiors for distance
        rho_point = rho123Normalized(:, 10 + iTrack).';
        [~, rho_dis] = knnsearch(rho123Normalized(:, 1:10 + iTrack).', rho_point, 'K', k);
        % rho_neighbors = rho123Normalized(:, k_ind).'; 
        R_k_rho = mean(rho_dis.^2);
        % search neighbiors for power
        power_point = prx3EstNormalized(:, 10 + iTrack).';
        [k_ind, ~] = knnsearch(prx3EstNormalized(:, 1:10 + iTrack).', power_point, 'K', k);
        rho_neighbors_given_power = rho123Normalized(:, k_ind).';
        R_k_rho_given_power = mean(vecnorm(rho_neighbors_given_power - repmat(rho_point, k, 1), 2, 2).^2);
        S3_k_rho_given_power(iTrack) = R_k_rho / R_k_rho_given_power;
    end
    % ===== Decision===== %
    hypo1_dep = zeros(1, nTrack-10);
    hypo3_dep = zeros(1, nTrack-10);
    % threshold
    for iTrack = 1:nTrack - 10
        thre(iTrack) = 0.99;%1.2*(k/(iTrack + 10))^(2/nRadar);
        if S1_k_rho_given_power(iTrack) >= thre(iTrack)
            hypo1_dep(iTrack) = 1;
        else
            hypo1_dep(iTrack) = 0;
        end
        if S3_k_rho_given_power(iTrack) >= thre(iTrack)
            hypo3_dep(iTrack) = 1;
        else
            hypo3_dep(iTrack) = 0;
        end
    end
    % ===== Statistics ===== %   
    count_fa_mi(iDataset) = sum(hypo1_mi == 1, "all");
    count_md_mi(iDataset) = sum(hypo3_mi == 0, "all");
    count_fa_dep(iDataset) = sum(hypo1_dep == 1, "all");
    count_md_dep(iDataset) = sum(hypo3_dep == 0, "all");
end
P_fa_mi(i_candidate) = sum(count_fa_mi) / (nRadar * nTrack * nDataset);
P_md_mi(i_candidate) = sum(count_md_mi) / (nRadar * nTrack * nDataset);
P_fa_dep(i_candidate) = sum(count_fa_dep) / (nTrack * nDataset);
P_md_dep(i_candidate) = sum(count_md_dep) / (nTrack * nDataset);
% fprintf("sigma = %d, Probability of false alarm   : MI %.4f, DEP %.4f \n",...
%     sigma, P_fa_mi(i_candidate), P_fa_dep(i_candidate))
% fprintf("sigma = %d, Probability of miss detection: MI %.4f, DEP %.4f \n",...
%     sigma, P_md_mi(i_candidate), P_md_dep(i_candidate))
end
% save_path_fa_mi = sprintf("result/P_fa_mi_%d_%d.mat", sigma_candidates(1), sigma_candidates(end));
% save_path_md_mi = sprintf("result/P_md_mi_%d_%d.mat", sigma_candidates(1), sigma_candidates(end));
% save_path_fa_dep = sprintf("result/P_fa_dep_%d_%d.mat", sigma_candidates(1), sigma_candidates(end));
% save_path_md_dep = sprintf("result/P_md_dep_%d_%d.mat", sigma_candidates(1), sigma_candidates(end));
% save(save_path_fa_mi, "P_fa_mi");
% save(save_path_md_mi, "P_md_mi");
% save(save_path_fa_dep, "P_fa_dep");
% save(save_path_md_dep, "P_md_dep");
%% plot
% ===== Mutual Information ===== %
% if adtest_enable == 1
% h11 = figure(11);
% x0_plot=500;
% y0_plot=400;
% width_plot = 1200;
% height_plot = 500;
% set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
% for iRadar = 1:nRadar
%     subplot(2, 4, iRadar)
%     yyaxis left
%     scatter(1:1:nTrack-10, hypo1_mi(iRadar,:))
%     ylim([0-0.01,1])
%     ylabel("Decision")
%     yyaxis right
%     plot(1:1:nTrack-10, p_value1(iRadar,:), 'LineWidth', 3)
%     ylim([0-0.01,1])
%     ylabel("$p$-value", "Interpreter","latex")
%     title(sprintf("Radar %d", iRadar))
% end
% sgtitle("Anderson-Darling for Policy 1")
% 
% h12 = figure(12);
% x0_plot=500;
% y0_plot=400;
% width_plot = 1200;
% height_plot = 500;
% set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
% for iRadar = 1:nRadar
%     subplot(2, 4, iRadar)
%     yyaxis left
%     scatter(1:1:nTrack-10, hypo2(iRadar,:))
%     ylim([0-0.01,1])
%     ylabel("Decision")
%     yyaxis right
%     plot(1:1:nTrack-10, p_value2(iRadar,:), 'LineWidth', 3)
%     ylim([0-0.01,1])
%     ylabel("$p$-value", "Interpreter","latex")
%     title(sprintf("Radar %d", iRadar))
% end
% sgtitle("Anderson-Darling for Policy 2")
% 
% h13 = figure(13);
% set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
% for iRadar = 1:nRadar
%     subplot(2, 4, iRadar)
%     yyaxis left
%     scatter(1:1:nTrack-10, hypo3_mi(iRadar,:))
%     ylim([0-0.01 1])
%     ylabel("Decision")
%     yyaxis right
%     plot(1:1:nTrack-10, p_value3(iRadar,:), 'LineWidth', 3)
%     ylim([0-0.01,1])
%     ylabel("$p$-value", "Interpreter","latex")
%     title(sprintf("Radar %d", iRadar))
% end
% sgtitle("Anderson-Darling for Policy 3")
% end
% 
% h14 = figure(14);
% plot(1:nTrack - 10, S1_k_rho_given_power);
% hold on
% plot(1:nTrack - 10, S3_k_rho_given_power);
% plot(1:nTrack - 10, thre)
% legend("policy 1", "policy 3", "threshold")
% 
% 
% x0_plot=500;
% y0_plot=400;
% width_plot = 500;
% height_plot =500;
% h15 = figure(15);
% set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
% scatter(1:1:nTrack-10, hypo1_dep)
% ylim([0-0.01 1])
% ylabel("Decision")
% sgtitle("Anderson-Darling for Policy 1")
% 
% h16 = figure(16);
% set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
% scatter(1:1:nTrack-10, hypo3_dep)
% ylim([0-0.01 1])
% ylabel("Decision")
% sgtitle("Anderson-Darling for Policy 3")

h17 = figure(17);
c1 = plot(sigma_candidates, P_fa_mi);
hold on 
c2 = plot(sigma_candidates, P_md_mi);
c3 = plot(sigma_candidates, P_fa_dep);
c4 = plot(sigma_candidates, P_md_dep);
legend(c1, "$P_{fa}$ MI", "Interpreter","latex")
legend(c2, "$P_{md}$ MI", "Interpreter","latex")
legend(c3, "$P_{fa}$ DEP", "Interpreter","latex")
legend(c4, "$P_{md}$ DEP", "Interpreter","latex")




