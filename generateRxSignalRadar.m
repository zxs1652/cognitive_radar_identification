function [receivedSignal, georange, tau, alpha] = generateRxSignalRadar(p, radarInfo,tarInfo,s)
%This function genearates radar received signal according to 
%the given power allocation
%   radarInfo -- radar information
%       .nRadar -- number of radar
%       .radarPos -- radar positions
%   tarInfo -- target information
%       .tarPos -- target position
%       .tarVel -- target velocity
%   s -- trasmitted LFM waveform
    %% set parameter
    c = 299792458; % speed of light
    
    nRadar = size(radarInfo,1) - 5;
    radarPos = radarInfo(1:nRadar, 1:3);
    lambda = radarInfo(1:nRadar, 4);
    nSampleRx = radarInfo(nRadar + 1, 1);
    nSampleTx = radarInfo(nRadar + 2, 1);
    frx = radarInfo(nRadar + 3, 1);
    ftx = radarInfo(nRadar + 4, 1);
    fc = radarInfo(nRadar + 5, 1);

    tarPos = tarInfo(:,1:3);
    tarVel = tarInfo(:,4:6);

    Trx = 1/frx;

    %% set memory
    receivedSignal = zeros(nRadar, nSampleRx);
    %% generate signal
    georange = vecnorm(repmat(tarPos,nRadar,1) - radarPos, 2, 2);
    range = 2*georange;
    % delay
    tau = range./c;
    % alpha = (lambda./(4*pi*georange)).^4;
    ii = (1:1:nSampleRx);
    % phase
    phi = 2*pi*fc.*tau;
    % amplitude
    alpha = sqrt(p.* (1./(georange)).^4); 
    % generate received signal
    for iRadar = 1: nRadar
        sampleInd = 1 + mod(floor(tau(iRadar)/Trx + (ii-1)/(frx/ftx)), nSampleTx);
        code = s(iRadar,:);
        receivedSignal(iRadar,:) = alpha(iRadar) * code(sampleInd).*exp(-1i*phi(iRadar));
        % Noise
        % noise = sigmaW*(sqrt(1/2)*randn(1,nSampleRx) + 1i*sqrt(1/2)*randn(1,nSampleRx));
        % receivedSignal(iRadar,:) = alpha(iRadar)*receivedSignal(iRadar,:) + noise;
    end
    % snr = sum(abs(sqrt(alpha(iRadar))*receivedSignal(iRadar,:).^2))/sum(abs(noise.^2))
    % SNR = 10*log10(snr)
end

