function globalMatrix = createMatrix(localMatrix, geomMatrix, h)
% CREATEMATRIX erstellt globale Matrix
%
%   globalMatrix = CREATEMATRIX(localMatrix, geomMatrix)
%   globalMatrix = CREATEMATRIX(localMatrix, geomMatrix, h)
%
%
    if not(exist('h', 'var'))
        h = 1;
    end
    
    n = max(max(geomMatrix));
    m = size(geomMatrix, 1);
    
    if length(h) == 1
        hVector = h * ones(m,1);
    else
        hVector = h;
    end
    
    globalMatrix = zeros(n);
    
    for i = 1:n
        parfor j = 1:n
            [rows1, cols1] = find(geomMatrix == i);
            [rows2, cols2] = find(geomMatrix == j);
            
            rows = intersect(rows1, rows2);
            
            value = 0;

            for k = 1:length(rows)
                elementIndex = rows(k);

                iindex = find(geomMatrix(elementIndex,:) == i);
                jindex = find(geomMatrix(elementIndex,:) == j);

                value = value + hVector(elementIndex) * localMatrix(iindex, jindex);
            end

            globalMatrix(i, j) = value;
        end
    end
end

