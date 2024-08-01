function point = randomPointInRectExclude(largerRect, smallerRect)
    % Define the boundaries of the larger rectangle
    xMinLarge = largerRect(1);
    yMinLarge = largerRect(2);
    widthLarge = largerRect(3);
    heightLarge = largerRect(4);
    xMaxLarge = xMinLarge + widthLarge;
    yMaxLarge = yMinLarge + heightLarge;

    % Define the boundaries of the smaller rectangle
    xMinSmall = smallerRect(1);
    yMinSmall = smallerRect(2);
    widthSmall = smallerRect(3);
    heightSmall = smallerRect(4);
    xMaxSmall = xMinSmall + widthSmall;
    yMaxSmall = yMinSmall + heightSmall;

    % Initialize point
    point = [0, 0];
    
    % Keep generating points until one is found outside the smaller rectangle
    while true
        % Generate a random point within the larger rectangle
        x = xMinLarge + rand * (xMaxLarge - xMinLarge);
        y = yMinLarge + rand * (yMaxLarge - yMinLarge);
        
        % Check if the point is outside the smaller rectangle
        if x < xMinSmall || x > xMaxSmall || y < yMinSmall || y > yMaxSmall
            point = [x, y];
            break;
        end
    end
end