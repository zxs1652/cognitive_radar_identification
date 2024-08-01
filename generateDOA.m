function [az, el] = generateDOA(radarPos, tarPos)
    deltaX = radarPos(:,1) - tarPos(1);
    deltaY = radarPos(:,2) - tarPos(2);
    deltaZ = radarPos(:,3) - tarPos(3);
    az = atan(deltaY./deltaX)/pi*180;
    for i = 1:length(deltaX)
        if deltaX(i)<0 && deltaY(i)>0
            az(i) = 180 + az(i);
        end
        if deltaX(i)<0 && deltaY(i)<0
            az(i) = az(i) - 180;
        end
    end
    el = atan(deltaZ./(sqrt(deltaX.^2 + deltaY.^2)))/pi*180;
end

