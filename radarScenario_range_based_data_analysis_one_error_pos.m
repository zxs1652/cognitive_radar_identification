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
%% Load data
% ====== load data ====== %
nDataset = 400;
sigma = 0;

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

% ====== compute mutual information ======%  
% load test data
iDataset = 1;
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
if adtest_enable == 1
mutualInfo1_seq = zeros(nRadar, nTrack-10);
% mutualInfo2_seq = zeros(nRadar, nTrack-10);
mutualInfo3_seq = zeros(nRadar, nTrack-10);
% mutualInfo4_seq = zeros(nRadar, nTrack-10);
% mutualInfo5_seq = zeros(nRadar, nTrack-10);
% mutualInfo6_seq = zeros(nRadar, nTrack-10);
mutualInfo_surr_seq = zeros(nRadar, nTrack-10);
for iTrack = 1:nTrack - 10
    fprintf("Computing MI of the %d step (dataset %d of %d) \n",...
        iTrack, iDataset, nDataset);
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
p_fa3_mi = zeros(nRadar, nTrack-10);
hypo3_mi = zeros(nRadar, nTrack-10);
% p_fa5 = zeros(nRadar, nTrack-10);
% hypo5 = zeros(nRadar, nTrack-10);
% p_fa6 = zeros(nRadar, nTrack-10);
% hypo6 = zeros(nRadar, nTrack-10);
for iTrack = 1:nTrack - 10
    for iRadar = 1:nRadar
        % [hypo2(iRadar, iTrack), p_fa2(iRadar, iTrack)] = ...
        %     adtest(I_prx2_rho(iRadar), 'Distribution', dist_surr{iRadar});
        [hypo3_mi(iRadar, iTrack), p_fa3_mi(iRadar, iTrack)] = ...
            adtest(mutualInfo3_seq(iRadar,iTrack), 'Distribution', dist_surr{iRadar});
        % [hypo5(iRadar, iTrack), p_fa5(iRadar, iTrack)] = ...
        %     adtest(I_prx5_rho(iRadar), 'Distribution', dist_surr{iRadar});
        % [hypo6(iRadar, iTrack), p_fa6(iRadar, iTrack)] = ...
        %     adtest(I_prx6_rho(iRadar), 'Distribution', dist_surr{iRadar});
    end
end

% ====== use non-cognitive case as test data ======%  
p_fa1_mi = zeros(nRadar, nTrack-10);
hypo1_mi = zeros(nRadar, nTrack-10);

for iTrack = 1:nTrack - 10
    for iRadar = 1:nRadar
        [hypo1_mi(iRadar, iTrack), p_fa1_mi(iRadar, iTrack)] = ...
            adtest(mutualInfo1_seq(iRadar,iTrack), 'Distribution', dist_surr{iRadar});
    end
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
    thre(iTrack) = 2*(k/(iTrack + 10))^(2/nRadar);
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
% count_fa_mi = 0;
% count_md_mi = 0;
% count_fa_dep = 0;
% count_md_dep = 0;

count_fa_mi = sum(hypo1_mi == 1, "all");
count_md_mi = sum(hypo3_mi == 0, "all");
count_fa_dep = sum(hypo1_dep == 1, "all");
count_md_dep = sum(hypo3_dep == 0, "all");
P_fa_mi = count_fa_mi / (nRadar * nTrack);
P_md_mi = count_md_mi / (nRadar * nTrack);
P_fa_dep = count_fa_dep / nTrack;
P_md_dep = count_md_dep / nTrack;
fprintf("Probability of false alarm   : MI %.4f, DEP %.4f \n", P_fa_mi, P_fa_dep)
fprintf("Probability of miss detection: MI %.4f, DEP %.4f \n", P_md_mi, P_md_dep)
%% plot
% ===== Mutual Information ===== %
% h11 = figure(11);
% iPlot = 1; 
% plot(1:1:nTrack-10, mutualInfo1_seq)
% legend('Radar 1', 'Radar 2', 'Radar 3', 'Radar 4',...
%     'Radar 5', 'Radar 6', 'Radar 7', 'Radar 8')
% h12 = figure(12);
% plot(1:1:nTrack-10, p_fa)
% legend('Radar 1', 'Radar 2', 'Radar 3', 'Radar 4',...
%     'Radar 5', 'Radar 6', 'Radar 7', 'Radar 8')
%

h11 = figure(11);
x0_plot=500;
y0_plot=400;
width_plot = 1000;
height_plot = 500;
set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
for iRadar = 1:nRadar
    subplot(2, 4, iRadar)
    yyaxis left
    scatter(1:1:nTrack-10, hypo1_mi(iRadar,:))
    ylim([0-0.01,1])
    ax = gca; % Get current axes
    ax.YTick = []; % Clear y-ticks
    ax.YTickLabel = [];
    if iRadar == 8
        ylabel("Decision", "FontSize",16)
    end
    yyaxis right
    plot(1:1:nTrack-10, p_fa1_mi(iRadar,:), 'LineWidth', 3)
    ylim([0-0.01,1])
    ax = gca; % Get current axes
    ax.YTick = []; % Clear y-ticks
    ax.YTickLabel = [];
    ax.FontSize = 12;
    if iRadar == 8
        ylabel("$p$-value", "Interpreter","latex", "FontSize",16)
        xlabel("Time Step")
    end
    title(sprintf("Radar #%d", iRadar), FontSize=12)
end
sgtitle("Anderson-Darling for Non-cognitive Radar$\,\,\,$     $\sigma_{\eta} = 0$", "fontSize", 16, "Interpreter","latex")
subplot(2,4,8)
yyaxis right
yticks([0, 0.5, 1])
yticklabels({'0','0.5','1'})
yyaxis left
yticks([0, 1])
yticklabels({'0','1'})
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
%     plot(1:1:nTrack-10, p_fa2(iRadar,:), 'LineWidth', 3)
%     ylim([0-0.01,1])
%     ylabel("$p$-value", "Interpreter","latex")
%     title(sprintf("Radar %d", iRadar))
% end
% sgtitle("Anderson-Darling for Policy 2")


h13 = figure(13);
set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
for iRadar = 1:nRadar
    subplot(2, 4, iRadar)
    yyaxis left
    scatter(1:1:nTrack-10, hypo3_mi(iRadar,:))
    ylim([0-0.01,1])
    ax = gca; % Get current axes
    ax.YTick = []; % Clear y-ticks
    ax.YTickLabel = [];
    if iRadar == 8
        ylabel("Decision", "FontSize",16)
    end
    yyaxis right
    plot(1:1:nTrack-10, p_fa3_mi(iRadar,:), 'LineWidth', 3)
    ylim([0-0.01,1])
    ax = gca; % Get current axes
    ax.YTick = []; % Clear y-ticks
    ax.YTickLabel = [];
    ax.FontSize = 12;
    if iRadar == 8
        ylabel("$p$-value", "Interpreter","latex", "FontSize",16)
        xlabel("Time Step")
    end
    title(sprintf("Radar #%d", iRadar), FontSize=12)
end
sgtitle("Anderson-Darling for Cognitive Radar$\,\,\,$     $\sigma_{\eta} = 0$", "fontSize", 16, "Interpreter","latex")
subplot(2,4,8)
yyaxis right
yticks([0, 0.5, 1])
yticklabels({'0','0.5','1'})
yyaxis left
yticks([0, 1])
yticklabels({'0','1'})

h14 = figure(14);
plot(1:nTrack - 10, S1_k_rho_given_power, LineWidth=1.5, Color="#0072bd");
hold on
plot(1:nTrack - 10, S3_k_rho_given_power, LineWidth=1.5, Color="#218f06");
plot(1:nTrack - 10, thre, LineWidth=1.5, LineStyle="--", Color="#d95319");
legend("Non-conitive Radar", "Cognitive Radar", "Threshold", fontsize=16)
xlabel("Time Step", FontSize=16)
ylabel("$S^k_l(X|Y)$", "Interpreter","latex", "FontSize",16)
title("Causality Inference $\sigma_{\eta} = 0$", "fontSize", 16, "Interpreter","latex")


x0_plot=500;
y0_plot=400;
width_plot = 500;
height_plot =500;
h15 = figure(15);
set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
scatter(1:1:nTrack-10, hypo1_dep)
ylim([0-0.01 1])
ylabel("Decision")
sgtitle("Anderson-Darling for Policy 1")

h16 = figure(16);
set(gcf,'position',[x0_plot,y0_plot,width_plot,height_plot])
scatter(1:1:nTrack-10, hypo3_dep)
ylim([0-0.01 1])
ylabel("Decision")
sgtitle("Anderson-Darling for Policy 3")



