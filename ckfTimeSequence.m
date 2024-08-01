function xEst = ckfTimeSequence(ckf, measurement, measurement_noise_std, T, radarPos)
%CKFTIMESEQUENCE use ckf to estimate state accordind to measurement
% in a time sequence
[nRadar, nTrack] = size(measurement);
xEst = zeros(6, nTrack);
for iTrack = 1:nTrack
    R = diag(measurement_noise_std(:, iTrack)).^2;
    ckf.MeasurementNoise = R;
    y = measurement(:,iTrack);
    [xPred, PPred] = predict(ckf, T);
    [xCorr, PCorr] = correct(ckf, y, nRadar, radarPos);

    xEst(:, iTrack) = xCorr;
end
end

