function y = linearTransform(xStart, xEnd, yStart, yEnd, x)
%LINEARTRANSFORM transformiert Werte linear

    offset = xStart - yStart;
    scale = (yStart - yEnd) / (xStart - xEnd);
    
    y = (x - xStart - (xEnd - xStart)/2) * scale + (yEnd - yStart)/2 + yStart;
end

