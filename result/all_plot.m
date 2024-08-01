% P_fa_mi_0_20 = load("P_fa_mi_0_20.mat");
% P_fa_mi_0_20 = P_fa_mi_0_20.P_fa_mi;
% P_fa_mi_20_20 = load("P_fa_mi_20_20.mat");
% P_fa_mi_20_20 = P_fa_mi_20_20.P_fa_mi;
% P_fa_mi_0_20(end) = P_fa_mi_20_20;
% P_fa_mi_21_30 = load("P_fa_mi_21_30.mat");
% P_fa_mi_21_30 = P_fa_mi_21_30.P_fa_mi;
% P_fa_mi_31_40 = load("P_fa_mi_31_40.mat");
% P_fa_mi_31_40 = P_fa_mi_31_40.P_fa_mi;
% P_fa_mi_41_50 = load("P_fa_mi_41_50.mat");
% P_fa_mi_41_50 = P_fa_mi_41_50.P_fa_mi;
% P_fa_mi = [P_fa_mi_0_20, P_fa_mi_21_30, P_fa_mi_31_40, P_fa_mi_41_50];
% 
% P_md_mi_0_20 = load("P_md_mi_0_20.mat");
% P_md_mi_0_20 = P_md_mi_0_20.P_md_mi;
% P_md_mi_20_20 = load("P_md_mi_20_20.mat");
% P_md_mi_20_20 = P_md_mi_20_20.P_md_mi;
% P_md_mi_0_20(end) = P_md_mi_20_20;
% P_md_mi_21_30 = load("P_md_mi_21_30.mat");
% P_md_mi_21_30 = P_md_mi_21_30.P_md_mi;
% P_md_mi_31_40 = load("P_md_mi_31_40.mat");
% P_md_mi_31_40 = P_md_mi_31_40.P_md_mi;
% P_md_mi_41_50 = load("P_md_mi_41_50.mat");
% P_md_mi_41_50 = P_md_mi_41_50.P_md_mi;
% P_md_mi = [P_md_mi_0_20, P_md_mi_21_30, P_md_mi_31_40, P_md_mi_41_50];
% 
% P_fa_dep_0_20 = load("P_fa_dep_0_20.mat");
% P_fa_dep_0_20 = P_fa_dep_0_20.P_fa_dep;
% P_fa_dep_20_20 = load("P_fa_dep_20_20.mat");
% P_fa_dep_20_20 = P_fa_dep_20_20.P_fa_dep;
% P_fa_dep_0_20(end) = P_fa_dep_20_20;
% P_fa_dep_21_30 = load("P_fa_dep_21_30.mat");
% P_fa_dep_21_30 = P_fa_dep_21_30.P_fa_dep;
% P_fa_dep_31_40 = load("P_fa_dep_31_40.mat");
% P_fa_dep_31_40 = P_fa_dep_31_40.P_fa_dep;
% P_fa_dep_41_50 = load("P_fa_dep_41_50.mat");
% P_fa_dep_41_50 = P_fa_dep_41_50.P_fa_dep;
% P_fa_dep = [P_fa_dep_0_20, P_fa_dep_21_30, P_fa_dep_31_40, P_fa_dep_41_50];
% 
% P_md_dep_0_20 = load("P_md_dep_0_20.mat");
% P_md_dep_0_20 = P_md_dep_0_20.P_md_dep;
% P_md_dep_20_20 = load("P_md_dep_20_20.mat");
% P_md_dep_20_20 = P_md_dep_20_20.P_md_dep;
% P_md_dep_0_20(end) = P_md_dep_20_20;
% P_md_dep_21_30 = load("P_md_dep_21_30.mat");
% P_md_dep_21_30 = P_md_dep_21_30.P_md_dep;
% P_md_dep_31_40 = load("P_md_dep_31_40.mat");
% P_md_dep_31_40 = P_md_dep_31_40.P_md_dep;
% P_md_dep_41_50 = load("P_md_dep_41_50.mat");
% P_md_dep_41_50 = P_md_dep_41_50.P_md_dep;
% P_md_dep = [P_md_dep_0_20, P_md_dep_21_30, P_md_dep_31_40, P_md_dep_41_50];
%% plot
% P_fa_dep = P_fa_dep/2.5;
% sigma_candidates = (0:1:50);
% h17 = figure(17);
% c1 = plot(sigma_candidates, P_fa_mi,"LineWidth",2,"Color","#4DBEEE");
% hold on 
% c2 = plot(sigma_candidates, P_md_mi,"LineWidth",2,"Color","#0000FF");
% c3 = plot(sigma_candidates, P_fa_dep,"LineWidth",2,"Color","#69e80e");
% c4 = plot(sigma_candidates, P_md_dep,"LineWidth",2,"Color","#135702");
% legend("$P_{fa}$ MI", "$P_{md}$ MI", "$P_{fa}$ CI", "$P_{md}$ CI",...
%     "Interpreter","latex", "FontSize", 16)
% ylabel("Experimental Probability", "Interpreter","latex","FontSize", 16)
% xlabel("Uncertainty of Radar Position $\sigma_{\eta}$",...
%     "Interpreter","latex","FontSize", 16)

%%
x1 = load("individual_received_signal.mat");
x1 = x1.code0;
y1 = load("y_beam_1.mat");
y1 = y1.yBeam;
y2 = load("y_beam_2.mat");
y2 = y2.yBeam;
frx = 10e6;
nSampleRx = 1e4;
t = (1/frx:1/frx:nSampleRx/frx);
figure(18)
plot(t, real(code0(1,:)/1.21*100), "LineWidth",1.5);
xlabel('Time (s)', "Interpreter","latex","FontSize",16);
ylabel('Real($x_{r,m}$)', "Interpreter","latex","FontSize",16);
legend("Individual Received Signal",...
    "fontSize", 16)

figure(19) % plot beamformer output        
plot(t,real(y1(1 ,:)*100), "LineWidth",1.5)
hold on
plot(t,real(y2(1 ,:)/1.4), "LineWidth",1.5)
title('Output of Beamformer');
xlabel('Time (s)', "Interpreter","latex","FontSize",16);
ylabel('Real($y_{m}$)', "Interpreter","latex","FontSize",16);
legend("Single-Constraint Beamformer", "Multiple-Constraint Beamformer",...
    "fontSize", 16)

% figure(20)
% plot(t, real(code0(1,:)/1.21), "LineWidth",1.5);
% hold on
% plot(t,real(y1(1 ,:)*100), "LineWidth",1.5)
% hold on
% plot(t,real(y2(1 ,:)/1.4), "LineWidth",1.5)
% xlabel('Time (s)', "Interpreter","latex","FontSize",16);
% ylabel('Real Amplitude', "Interpreter","latex","FontSize",16);
% legend("Individual Received Signal", "Single-Constraint Beamformer", "Multiple-Constraint Beamformer",...
%     "fontSize", 16)
% 
% ylim([-1e3, 1e3])



