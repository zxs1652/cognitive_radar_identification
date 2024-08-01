function [cost, pSum] = objFunPCRB(p, tarPos, radarPos, sigmaW, fc, Jp)

nRadar = size(radarPos,1);
rho = vecnorm(repmat(tarPos,nRadar,1) - radarPos, 2, 2);
alpha = sqrt(p.* (1./(rho)).^4);
sigmaTau = sqrt(sigmaW^2./(8*pi^2*fc^2*alpha.^2));
R = diag(sigmaTau).^2;
dH = gradient_H(tarPos, radarPos);
Jd3 = dH.' * inv(R) * dH;
J = Jp + Jd3;
cost = trace(inv(J(1:3,1:3)));
pSum = sum(p);
end

function [H_gradient] = gradient_H(tarPos,radarPos)
    c = 299792458;
    nRadar = size(radarPos,1);
    diff = tarPos - radarPos;
    ranges = vecnorm(diff,2,2);
    H_gradient = [2/c * diff./repmat(ranges, 1, 3), zeros(nRadar, 3)];

end