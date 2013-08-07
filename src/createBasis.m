function Basis = createBasis(functions, elements, x, B, xMin, xMax)
%CREATEBASIS erstellt eine Basis-Matrix
% Basis = CREATEBASIS(functions, elements, x, B, xMin, xMax)

    if not(exist('xMin', 'var'))
        xMin = 0;
    end
    
    if not(exist('xMax', 'var'))
        xMax = 1;
    end

    n = max(max(B));
    elementCount = size(B,1);
    Basis = zeros(length(x), n);
    
    for elementIndex=1:elementCount
        row = B(elementIndex, :);
        
        for j=1:length(row)
            funcIndex = row(j);
            
            left = elements(elementIndex);
            right = elements(elementIndex+1);
                        
            interval = x >= left & x <= right;
            
            xInterval = linearTransform(left, right, xMin, xMax, x(interval));
            
            values = functions{j}(xInterval);
            
            Basis(interval,funcIndex) = values';
        end
    end
end

