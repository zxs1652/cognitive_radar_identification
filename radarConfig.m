% This script includes all radar configurations

%%
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
hgt = 1000;
nTrack = 300;

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