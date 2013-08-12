function b = createVector(localF, geom, elements, f, xMin, xMax)
%CREATEVECTOR erzeugt die rechte Seite des Gleichungssystems
% 
% b = CREATEVECTOR(localF, geom, elements, func)
% b = CREATEVECTOR(localF, geom, elements, func, xMin, xMax)
    n = max(max(geom));
    elementCount = size(geom, 1);
    functionCount = length(localF);
    
    if not(exist('xMin', 'var'))
        xMin = 0;
    end
    
    if not(exist('xMax', 'var'))
        xMax = 1;
    end
    
    b = zeros(n, 1);
    
    for elementIndex = 1:elementCount
        left = elements(elementIndex);
        right = elements(elementIndex + 1);
            
        for functionIndex = 1:functionCount
            phiIndex = geom(elementIndex, functionIndex);
            
            func = @(x) localF{functionIndex}(linearTransform(left, right, xMin, xMax, x)) .* f(x);
            
            b(phiIndex) = b(phiIndex) + quadgk(func, left, right);
        end
    end
end

